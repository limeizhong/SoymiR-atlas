#!/usr/bin/env python3
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import csv
import math
import re


SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
INPUT_DIR = MODULE_DIR / "input"
INTERMEDIATE_DIR = MODULE_DIR / "results" / "intermediate"
OVERLAP_RATIO = 0.8


def clean_chr(value: str) -> str:
    value = value.strip()
    value = re.sub(r"^chr", "", value, flags=re.I)
    value = re.sub(r"^Gma-Gm", "", value, flags=re.I)
    value = re.sub(r"^Gm", "", value, flags=re.I)
    if value.isdigit():
        return f"Chr{int(value):02d}"
    if value.startswith("Chr") and value[3:].isdigit():
        return f"Chr{int(value[3:]):02d}"
    return value


def clean_strand(value: str) -> str:
    value = value.strip().lower()
    if value in {"+", "plus", "1"}:
        return "+"
    if value in {"-", "minus", "-1"}:
        return "-"
    return value


def parse_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.inf if value.lower() == "inf" else 0.0


def interval_len(row: dict[str, str]) -> int:
    return int(row["Genome_end"]) - int(row["Genome_start"]) + 1


def overlap_bp(a: dict[str, str], b: dict[str, str]) -> int:
    if a["Chr"] != b["Chr"] or a["S_strand"] != b["S_strand"]:
        return 0
    return max(
        0,
        min(int(a["Genome_end"]), int(b["Genome_end"]))
        - max(int(a["Genome_start"]), int(b["Genome_start"]))
        + 1,
    )


def overlap_ratio(a: dict[str, str], b: dict[str, str]) -> float:
    ov = overlap_bp(a, b)
    if ov <= 0:
        return 0.0
    return ov / min(interval_len(a), interval_len(b))


def conflicts_with_used(candidate: dict[str, str], used: list[dict[str, str]]) -> bool:
    return any(overlap_ratio(candidate, u) >= OVERLAP_RATIO for u in used)


def hit_score(row: dict[str, str]) -> tuple[float, float, float, int]:
    return (
        parse_float(row["Bitscore"]),
        parse_float(row["Pident"]),
        parse_float(row["Qcovs"]),
        int(row["Aln_Length"]),
    )


def read_hits(path: Path) -> tuple[list[str], dict[str, list[dict[str, str]]]]:
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        hits = defaultdict(list)
        for row in reader:
            row = dict(row)
            s = int(row["Genome_start"])
            e = int(row["Genome_end"])
            if s > e:
                s, e = e, s
            row["Genome_start"] = str(s)
            row["Genome_end"] = str(e)
            hits[row["Query_ID"]].append(row)
    for q in hits:
        hits[q].sort(
            key=lambda r: (
                -parse_float(r["Bitscore"]),
                -parse_float(r["Pident"]),
                -parse_float(r["Qcovs"]),
                -int(r["Aln_Length"]),
                r["Chr"],
                int(r["Genome_start"]),
                int(r["Genome_end"]),
            )
        )
    return fieldnames, hits


def top_tie_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if not rows:
        return []
    best = hit_score(rows[0])
    return [r for r in rows if hit_score(r) == best]


def load_mirbase_original_positions(path: Path) -> dict[str, dict[str, str]]:
    mapping = {}
    with path.open() as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6:
                mapping[parts[0]] = {
                    "chr": clean_chr(parts[2]),
                    "strand": clean_strand(parts[5]),
                }
    return mapping


def load_pmiren_original_positions(path: Path) -> dict[str, dict[str, str]]:
    mapping = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["miRNA_locus_ID"]] = {
                "chr": clean_chr(row["Chromosome"]),
                "strand": clean_strand(row["Strand"]),
            }
    return mapping


def query_short(query_id: str) -> str:
    return query_id.split("|", 1)[0]


def query_chrom_from_header(query_id: str) -> str:
    parts = query_id.split("|")
    if len(parts) >= 3:
        m = re.match(r"([^:]+):", parts[2])
        if m:
            return clean_chr(m.group(1))
    return ""


def original_from_header(query_id: str) -> dict[str, str] | None:
    parts = query_id.split("|")
    if len(parts) >= 3:
        m = re.match(r"([^:]+):\d+-\d+\(([+-])\)", parts[2])
        if m:
            return {"chr": clean_chr(m.group(1)), "strand": clean_strand(m.group(2))}
    return None


def resolve_one(
    label: str,
    all_hits: Path,
    db_original_by_name: dict[str, dict[str, str]],
    output: Path,
    changed_output: Path,
    unresolved_output: Path,
):
    fieldnames, hits = read_hits(all_hits)
    selected: dict[str, dict[str, str]] = {}
    changed = []
    unresolved = []

    def db_original(q: str) -> dict[str, str]:
        return original_from_header(q) or db_original_by_name.get(query_short(q), {})

    def preferred_top_tie_candidates(q: str) -> list[dict[str, str]]:
        candidates = top_tie_candidates(hits[q])
        original = db_original(q)
        chrom = original.get("chr", "")
        strand = original.get("strand", "")
        if chrom and strand:
            compatible = [
                r for r in candidates
                if clean_chr(r["Chr"]) == chrom and clean_strand(r["S_strand"]) == strand
            ]
            if compatible:
                return compatible
        return candidates

    ordered_queries = sorted(hits, key=lambda q: (query_short(q), q))

    for q in ordered_queries:
        candidates = preferred_top_tie_candidates(q)
        original = db_original(q)
        chrom = original.get("chr", "")
        strand = original.get("strand", "")

        def candidate_key(row: dict[str, str]) -> tuple[int, str, int, int]:
            return (
                0 if chrom and clean_chr(row["Chr"]) == chrom and strand and clean_strand(row["S_strand"]) == strand else 1,
                row["Chr"],
                int(row["Genome_start"]),
                int(row["Genome_end"]),
            )

        chosen = sorted(candidates, key=candidate_key)[0]
        selected[q] = chosen
        original = hits[q][0]
        if (
            original["Chr"],
            original["Genome_start"],
            original["Genome_end"],
            original["S_strand"],
        ) != (
            chosen["Chr"],
            chosen["Genome_start"],
                chosen["Genome_end"],
                chosen["S_strand"],
        ):
            changed.append((q, chrom, strand, original, chosen, len(top_tie_candidates(hits[q]))))

    # Report residual conflicts that could not be removed by top-tier reassignment.
    selected_rows = list(selected.values())
    for i, a in enumerate(selected_rows):
        for b in selected_rows[i + 1 :]:
            if overlap_ratio(a, b) >= OVERLAP_RATIO:
                unresolved.append((a, b, overlap_bp(a, b), overlap_ratio(a, b)))

    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for q in sorted(selected):
            writer.writerow(selected[q])

    with changed_output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "Query_ID",
                "Database_chr",
                "Database_strand",
                "Top_tie_candidate_count",
                "Old_Chr",
                "Old_start",
                "Old_end",
                "Old_strand",
                "New_Chr",
                "New_start",
                "New_end",
                "New_strand",
            ]
        )
        for q, chrom, strand, old, new, n in changed:
            writer.writerow(
                [
                    q,
                    chrom,
                    strand,
                    n,
                    old["Chr"],
                    old["Genome_start"],
                    old["Genome_end"],
                    old["S_strand"],
                    new["Chr"],
                    new["Genome_start"],
                    new["Genome_end"],
                    new["S_strand"],
                ]
            )

    with unresolved_output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "Query_A",
                "Query_B",
                "Chr",
                "Strand",
                "A_start",
                "A_end",
                "B_start",
                "B_end",
                "Overlap_bp",
                "Overlap_over_min_len",
            ]
        )
        for a, b, ov, ratio in unresolved:
            writer.writerow(
                [
                    a["Query_ID"],
                    b["Query_ID"],
                    a["Chr"],
                    a["S_strand"],
                    a["Genome_start"],
                    a["Genome_end"],
                    b["Genome_start"],
                    b["Genome_end"],
                    ov,
                    f"{ratio:.4f}",
                ]
            )

    print(label)
    print(" queries:", len(selected))
    print(" changed_best_hit_selection:", len(changed))
    print(" residual_overlap80_pairs:", len(unresolved))
    print(" output:", output)
    print(" changed:", changed_output)
    print(" unresolved:", unresolved_output)


def main():
    resolve_one(
        "miRBase",
        INTERMEDIATE_DIR / "gma-hairpin_ZH13v2_genome_positions_all_hits.tsv",
        load_mirbase_original_positions(INPUT_DIR / "miRbase_gma_position.txt"),
        INTERMEDIATE_DIR / "gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv",
        INTERMEDIATE_DIR / "gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval_changed.tsv",
        INTERMEDIATE_DIR / "gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval_unresolved_overlap80.tsv",
    )
    resolve_one(
        "pmiREN",
        INTERMEDIATE_DIR / "pmiren_gmax_hairpin_ZH13v2_genome_positions_all_hits.tsv",
        load_pmiren_original_positions(INPUT_DIR / "pmiren_gma_position.txt"),
        INTERMEDIATE_DIR / "pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv",
        INTERMEDIATE_DIR / "pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval_changed.tsv",
        INTERMEDIATE_DIR / "pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval_unresolved_overlap80.tsv",
    )


if __name__ == "__main__":
    main()
