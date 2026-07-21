#!/usr/bin/env python3
"""Extract CleaveLand4 target binding-site sequences and genome coordinates."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


RESULT_RE = re.compile(
    r"^5'\s+([ACGTUNacgtun-]+)\s+3'\s+Transcript:\s+([^:]+):(\d+)-(\d+)\s+Slice Site:(\d+)"
)


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.rstrip(";").split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def load_cds_map(gff_path: Path) -> dict[str, dict[str, object]]:
    cds_map: dict[str, dict[str, object]] = {}
    with gff_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            chrom, _source, _feature, start, end, _score, strand, _phase, attrs = fields
            parent = parse_attrs(attrs).get("Parent")
            if not parent:
                continue
            entry = cds_map.setdefault(
                parent,
                {"chrom": chrom, "strand": strand, "segments": []},
            )
            entry["segments"].append((int(start), int(end)))

    for entry in cds_map.values():
        segments = entry["segments"]
        if entry["strand"] == "+":
            segments.sort(key=lambda x: x[0])
        else:
            segments.sort(key=lambda x: x[0], reverse=True)
    return cds_map


def load_fasta(fasta_path: Path | None) -> dict[str, str]:
    if not fasta_path:
        return {}
    seqs: dict[str, list[str]] = {}
    current_id: str | None = None
    current_aliases: list[str] = []
    with fasta_path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
                current_id = header.split()[0]
                current_aliases = [current_id]
                for match in re.finditer(r"(?:^|\s)OriID=([^\s]+)", header):
                    current_aliases.append(match.group(1))
                seqs[current_id] = []
                for alias in current_aliases[1:]:
                    seqs[alias] = seqs[current_id]
            elif current_id is not None:
                seqs[current_id].append(line.upper())
    return {seq_id: "".join(parts) for seq_id, parts in seqs.items()}


def transcript_interval_to_genome(
    transcript_id: str,
    tx_start: int,
    tx_end: int,
    cds_map: dict[str, dict[str, object]],
) -> tuple[str, int, int, str, bool] | None:
    entry = cds_map.get(transcript_id)
    if not entry:
        return None

    chrom = entry["chrom"]
    strand = entry["strand"]
    coords: list[int] = []
    tx_pos = 1
    overlapped_segments = 0
    for seg_start, seg_end in entry["segments"]:
        seg_len = seg_end - seg_start + 1
        seg_tx_start = tx_pos
        seg_tx_end = tx_pos + seg_len - 1
        ov_start = max(tx_start, seg_tx_start)
        ov_end = min(tx_end, seg_tx_end)
        if ov_start <= ov_end:
            overlapped_segments += 1
            for pos in range(ov_start, ov_end + 1):
                offset = pos - seg_tx_start
                genome_pos = seg_start + offset if strand == "+" else seg_end - offset
                coords.append(genome_pos)
        tx_pos += seg_len

    expected_len = tx_end - tx_start + 1
    if len(coords) != expected_len:
        return None
    return chrom, min(coords), max(coords), strand, overlapped_segments > 1


def source_fields(result_file: Path, results_dir: Path) -> tuple[str, str, str]:
    rel = result_file.relative_to(results_dir.parent).as_posix()
    name = result_file.name
    if name.endswith(".res.txt"):
        stem = name[: -len(".res.txt")]
    else:
        stem = result_file.stem
    library, mirna = stem.split(".", 1)
    return library, mirna, rel


def load_id_map(path: Path | None) -> dict[str, str]:
    if not path:
        return {}
    mapping: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            old_id = row.get("old_id", "").strip()
            new_id = row.get("new_id", "").strip()
            if old_id:
                mapping[old_id] = new_id
    return mapping


def load_existing_detail(path: Path | None) -> tuple[list[dict[str, object]], Counter]:
    records: list[dict[str, object]] = []
    stats: Counter = Counter()
    if not path:
        return records, stats
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["binding_site_sequence_for_blast"] = (
                row["binding_site_sequence"].upper().replace("U", "T").replace("-", "")
            )
            row["transcript_start"] = int(row["transcript_start"])
            row["transcript_end"] = int(row["transcript_end"])
            row["slice_site"] = int(row["slice_site"])
            records.append(row)
    stats["existing_detail_records_loaded"] = len(records)
    return records, stats


def parse_results(
    results_dir: Path,
    id_map: dict[str, str] | None = None,
    skip_na_mapped_ids: bool = False,
) -> tuple[list[dict[str, object]], Counter]:
    records: list[dict[str, object]] = []
    stats: Counter = Counter()
    for result_file in sorted(results_dir.glob("*.res.txt")):
        stats["result_files_scanned"] += 1
        file_records = 0
        library_id, mirna_id, source_file = source_fields(result_file, results_dir)
        if id_map is not None and mirna_id in id_map:
            mapped_id = id_map[mirna_id]
            if skip_na_mapped_ids and mapped_id == "NA":
                stats["result_files_skipped_na_mapped_id"] += 1
                continue
            mirna_id = mapped_id
            stats["result_files_with_mapped_id"] += 1
        elif id_map is not None:
            stats["result_files_without_id_mapping"] += 1
        with result_file.open(errors="replace") as handle:
            for line in handle:
                match = RESULT_RE.match(line.strip())
                if not match:
                    continue
                binding_seq, transcript_id, tx_start, tx_end, slice_site = match.groups()
                target_gene = transcript_id.split(".m", 1)[0]
                records.append(
                    {
                        "library_id": library_id,
                        "miRNA_ID": mirna_id,
                        "target_gene": target_gene,
                        "binding_site_sequence": binding_seq.upper(),
                        "binding_site_sequence_for_blast": binding_seq.upper()
                        .replace("U", "T")
                        .replace("-", ""),
                        "transcript_id": transcript_id,
                        "transcript_start": int(tx_start),
                        "transcript_end": int(tx_end),
                        "slice_site": int(slice_site),
                        "source_file": source_file,
                    }
                )
                file_records += 1
        if file_records:
            stats["result_files_with_records"] += 1
            stats["raw_binding_site_records"] += file_records
        else:
            stats["result_files_without_records"] += 1
    return records, stats


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--transcriptome", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", default="results2025_239")
    parser.add_argument("--id-map", type=Path, help="Optional old_id/new_id two-column TSV for miRNA IDs.")
    parser.add_argument(
        "--skip-na-mapped-ids",
        action="store_true",
        help="Skip result files whose mapped miRNA ID is NA.",
    )
    parser.add_argument(
        "--existing-detail",
        type=Path,
        help="Optional existing binding_site_detail.tsv to include before summarizing.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    id_map = load_id_map(args.id_map)
    records, stats = parse_results(args.results_dir, id_map if args.id_map else None, args.skip_na_mapped_ids)
    existing_records, existing_stats = load_existing_detail(args.existing_detail)
    records = existing_records + records
    stats.update(existing_stats)
    if args.id_map:
        stats["id_mapping_entries_loaded"] = len(id_map)
    cds_map = load_cds_map(args.gff)
    transcriptome = load_fasta(args.transcriptome)

    detail_fields = [
        "library_id",
        "miRNA_ID",
        "target_gene",
        "binding_site_sequence",
        "transcript_id",
        "transcript_start",
        "transcript_end",
        "slice_site",
        "source_file",
    ]
    write_tsv(args.outdir / f"{args.prefix}_binding_site_detail.tsv", detail_fields, records)

    by_full_key: dict[tuple[str, str, str, str], dict[str, object]] = {}
    for row in records:
        key = (
            str(row["library_id"]),
            str(row["miRNA_ID"]),
            str(row["target_gene"]),
            str(row["binding_site_sequence"]),
        )
        by_full_key.setdefault(key, row)
    stats["deduplicated_by_library_mirna_target_site"] = len(by_full_key)

    unique: dict[tuple[str, str, str], dict[str, object]] = {}
    unique_support: dict[tuple[str, str, str], dict[str, set[str] | int | set[str]]] = {}
    for row in records:
        key = (
            str(row["miRNA_ID"]),
            str(row["target_gene"]),
            str(row["binding_site_sequence"]),
        )
        unique.setdefault(key, row)
        support = unique_support.setdefault(
            key,
            {"libraries": set(), "result_record_count": 0, "source_files": set()},
        )
        support["libraries"].add(str(row["library_id"]))
        support["source_files"].add(str(row["source_file"]))
        support["result_record_count"] = int(support["result_record_count"]) + 1

    unique_rows = []
    for key in sorted(unique):
        row = unique[key]
        support = unique_support[key]
        libraries = sorted(support["libraries"])
        unique_rows.append(
            {
                "miRNA_ID": row["miRNA_ID"],
                "target_gene": row["target_gene"],
                "binding_site_sequence": row["binding_site_sequence"],
                "library_count": len(libraries),
                "result_record_count": support["result_record_count"],
                "libraries": ";".join(libraries),
            }
        )

    unique_fields = [
        "miRNA_ID",
        "target_gene",
        "binding_site_sequence",
        "library_count",
        "result_record_count",
        "libraries",
    ]
    write_tsv(args.outdir / f"{args.prefix}_binding_site_unique.tsv", unique_fields, unique_rows)
    write_tsv(
        args.outdir / f"{args.prefix}_binding_site_unique_minimal.tsv",
        ["miRNA_ID", "target_gene", "binding_site_sequence"],
        unique_rows,
    )

    candidates = []
    missing_cds = 0
    bad_coords = 0
    for row in records:
        mapped = transcript_interval_to_genome(
            str(row["transcript_id"]),
            int(row["transcript_start"]),
            int(row["transcript_end"]),
            cds_map,
        )
        if mapped is None:
            if str(row["transcript_id"]) not in cds_map:
                missing_cds += 1
            else:
                bad_coords += 1
            continue
        chrom, genome_start, genome_end, strand, spans = mapped
        tx_seq = transcriptome.get(str(row["transcript_id"]), "")
        expected_seq = str(row["binding_site_sequence_for_blast"])
        if tx_seq:
            observed_seq = tx_seq[int(row["transcript_start"]) - 1 : int(row["transcript_end"])]
            sequence_matches = observed_seq == expected_seq
        else:
            sequence_matches = ""
        key = f"{row['miRNA_ID']}\x01{row['target_gene']}\x01{row['binding_site_sequence']}"
        candidates.append(
            {
                "key": key,
                "miRNA_ID": row["miRNA_ID"],
                "target_gene": row["target_gene"],
                "binding_site_sequence": row["binding_site_sequence"],
                "binding_site_sequence_for_blast": row["binding_site_sequence_for_blast"],
                "transcript_id": row["transcript_id"],
                "transcript_start": row["transcript_start"],
                "transcript_end": row["transcript_end"],
                "slice_site": row["slice_site"],
                "chromosome": chrom,
                "genome_start": genome_start,
                "genome_end": genome_end,
                "strand": strand,
                "spans_exon_junction": "TRUE" if spans else "FALSE",
                "sequence_matches_transcript": sequence_matches,
                "source_file": row["source_file"],
            }
        )

    candidate_fields = [
        "key",
        "miRNA_ID",
        "target_gene",
        "binding_site_sequence",
        "binding_site_sequence_for_blast",
        "transcript_id",
        "transcript_start",
        "transcript_end",
        "slice_site",
        "chromosome",
        "genome_start",
        "genome_end",
        "strand",
        "spans_exon_junction",
        "sequence_matches_transcript",
        "source_file",
    ]
    write_tsv(
        args.outdir / f"{args.prefix}_binding_site_genome_mapping_candidates.tsv",
        candidate_fields,
        candidates,
    )

    coord_support: dict[tuple[str, str, str, str, int, int], set[str]] = defaultdict(set)
    coords_by_key: dict[str, set[tuple[str, int, int]]] = defaultdict(set)
    for row in candidates:
        coord_key = (
            str(row["key"]),
            str(row["chromosome"]),
            int(row["genome_start"]),
            int(row["genome_end"]),
        )
        coord_support[(coord_key[0], coord_key[1], coord_key[2], coord_key[3], 0, 0)].add(
            str(row["source_file"])
        )
        coords_by_key[str(row["key"])].add(
            (str(row["chromosome"]), int(row["genome_start"]), int(row["genome_end"]))
        )

    ambiguous_rows = []
    for key, coord_set in sorted(coords_by_key.items()):
        if len(coord_set) <= 1:
            continue
        for chrom, start, end in sorted(coord_set):
            support = sum(
                1
                for row in candidates
                if row["key"] == key
                and row["chromosome"] == chrom
                and int(row["genome_start"]) == start
                and int(row["genome_end"]) == end
            )
            ambiguous_rows.append(
                {
                    "coord_key": f"{key}\x02{chrom}\x02{start}\x02{end}",
                    "supporting_result_records": support,
                    "key": key,
                    "chromosome": chrom,
                    "genome_start": start,
                    "genome_end": end,
                }
            )
    write_tsv(
        args.outdir / f"{args.prefix}_binding_site_genome_mapping_ambiguous.tsv",
        ["coord_key", "supporting_result_records", "key", "chromosome", "genome_start", "genome_end"],
        ambiguous_rows,
    )

    best_coord_rows = []
    candidate_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in candidates:
        candidate_groups[str(row["key"])].append(row)
    for key, group in sorted(candidate_groups.items()):
        coord_counts = Counter(
            (
                row["chromosome"],
                row["genome_start"],
                row["genome_end"],
                row["strand"],
            )
            for row in group
        )
        (chrom, start, end, strand), support = coord_counts.most_common(1)[0]
        first = group[0]
        best_coord_rows.append(
            {
                "miRNA_ID": first["miRNA_ID"],
                "target_gene": first["target_gene"],
                "binding_site_sequence": first["binding_site_sequence"],
                "chromosome": chrom,
                "genome_start": start,
                "genome_end": end,
                "strand": strand,
                "supporting_result_records": support,
                "candidate_coordinate_count": len(coord_counts),
            }
        )
    mapped_fields = [
        "miRNA_ID",
        "target_gene",
        "binding_site_sequence",
        "chromosome",
        "genome_start",
        "genome_end",
        "strand",
        "supporting_result_records",
        "candidate_coordinate_count",
    ]
    write_tsv(
        args.outdir / f"{args.prefix}_binding_site_unique_minimal.mapped.tsv",
        mapped_fields,
        best_coord_rows,
    )
    with (args.outdir / f"{args.prefix}_binding_site_unique_minimal.bed").open("w") as bed:
        for row in best_coord_rows:
            name = f"{row['miRNA_ID']}|{row['target_gene']}|{row['binding_site_sequence']}"
            bed.write(
                f"{row['chromosome']}\t{int(row['genome_start']) - 1}\t{row['genome_end']}\t"
                f"{name}\t{row['supporting_result_records']}\t{row['strand']}\n"
            )

    stats["unique_miRNA_target_site_records"] = len(unique_rows)
    stats["unique_miRNA_target_pairs"] = len({(r["miRNA_ID"], r["target_gene"]) for r in unique_rows})
    stats["unique_target_genes"] = len({r["target_gene"] for r in unique_rows})
    stats["records_with_genome_coordinates"] = len({r["key"] for r in candidates})
    stats["records_without_genome_coordinates"] = len(unique_rows) - stats["records_with_genome_coordinates"]
    stats["mapped_candidate_result_records"] = len(candidates)
    stats["records_with_multiple_candidate_coordinates"] = len({r["key"] for r in ambiguous_rows})
    if transcriptome:
        stats["candidate_records_not_matching_CDS_sequence"] = sum(
            1 for row in candidates if row["sequence_matches_transcript"] is False
        )
    stats["candidate_records_missing_CDS_transcript"] = missing_cds
    stats["candidate_records_bad_CDS_coordinates"] = bad_coords

    summary_rows = [{"metric": key, "value": value} for key, value in stats.items()]
    write_tsv(
        args.outdir / f"{args.prefix}_binding_site_extraction_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    write_tsv(
        args.outdir / f"{args.prefix}_binding_site_genome_mapping_summary.tsv",
        ["metric", "value"],
        [
            {"metric": "unique_miRNA_target_site_records", "value": stats["unique_miRNA_target_site_records"]},
            {"metric": "records_with_genome_coordinates", "value": stats["records_with_genome_coordinates"]},
            {"metric": "records_without_genome_coordinates", "value": stats["records_without_genome_coordinates"]},
            {"metric": "mapped_candidate_result_records", "value": stats["mapped_candidate_result_records"]},
            {
                "metric": "records_with_multiple_candidate_coordinates",
                "value": stats["records_with_multiple_candidate_coordinates"],
            },
            {
                "metric": "candidate_records_not_matching_CDS_sequence",
                "value": stats.get("candidate_records_not_matching_CDS_sequence", "NA"),
            },
            {"metric": "candidate_records_missing_CDS_transcript", "value": missing_cds},
            {"metric": "candidate_records_bad_CDS_coordinates", "value": bad_coords},
        ],
    )

    pairs_counter = Counter((row["miRNA_ID"], row["target_gene"]) for row in unique_rows)
    multi_pairs = [
        {"miRNA_ID": mirna, "target_gene": target, "binding_site_count": count}
        for (mirna, target), count in sorted(pairs_counter.items())
        if count > 1
    ]
    write_tsv(
        args.outdir / f"{args.prefix}_pairs_with_multiple_binding_sites.tsv",
        ["miRNA_ID", "target_gene", "binding_site_count"],
        multi_pairs,
    )

    id_map_note = ""
    if args.id_map:
        id_map_note = (
            f"\n老编号结果使用 `{args.id_map}` 转换为新编号；"
            f"映射为 NA 的 `.res.txt` 文件跳过 {stats['result_files_skipped_na_mapped_id']:,} 个。"
        )
    existing_note = ""
    if args.existing_detail:
        existing_note = (
            f"\n同时合并已有 detail 表 `{args.existing_detail}` 中的 "
            f"{stats['existing_detail_records_loaded']:,} 条记录。"
        )

    report = f"""# {args.prefix} mature miRNA 靶基因结合位点序列抽取结果

从 `{args.results_dir}` 中扫描 {stats['result_files_scanned']:,} 个 `.res.txt` 文件；其中 {stats['result_files_without_records']:,} 个没有 CleaveLand 命中记录，{stats['result_files_with_records']:,} 个包含结果。
{id_map_note}{existing_note}
从 CleaveLand Pretty 输出的 `5' ... 3' Transcript:` 行中解析出 {stats['raw_binding_site_records']:,} 条原始结合位点记录；按 `library_id + miRNA + target_gene + binding_site_sequence` 去重后为 {stats['deduplicated_by_library_mirna_target_site']:,} 条。
进一步按 `miRNA + target_gene + binding_site_sequence` 去冗余后，得到 {stats['unique_miRNA_target_site_records']:,} 条唯一记录，涉及 {stats['unique_miRNA_target_pairs']:,} 个 miRNA-target 对和 {stats['unique_target_genes']:,} 个靶基因。

主输出表：`{args.prefix}_binding_site_unique.tsv`，包含 `miRNA_ID`、`target_gene`、`binding_site_sequence`，并附带 `library_count`、`result_record_count` 和 `libraries` 便于后续筛选。
最小三列表：`{args.prefix}_binding_site_unique_minimal.tsv`。
基因组坐标表：`{args.prefix}_binding_site_genome_mapping_candidates.tsv`，坐标由 `{args.gff}` 中对应 CDS 坐标换算得到。
"""
    (args.outdir / f"{args.prefix}_binding_site_extraction_report.md").write_text(report)


if __name__ == "__main__":
    main()
