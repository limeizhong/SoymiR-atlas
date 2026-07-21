#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


OUT_FIELDS = [
    "Query_ID",
    "Chr",
    "Subject_ID",
    "Genome_start",
    "Genome_end",
    "S_strand",
    "Pident",
    "Qcovs",
    "Aln_Length",
    "Q_length",
    "Q_start",
    "Q_end",
    "S_start",
    "S_end",
    "Evalue",
    "Bitscore",
    "Subject_description",
]


def read_fasta_ids(path: Path) -> list[str]:
    ids = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def chr_from_subject(subject_id: str, subject_description: str) -> str:
    match = re.search(r"\bOriSeqID=(Chr\d+)", subject_description)
    if match:
        return match.group(1)
    match = re.search(r"GWHAAEV0*([0-9]+)\.1", subject_id)
    if match:
        return f"Chr{int(match.group(1)):02d}"
    return subject_id


def parse_blast(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for parts in reader:
            if not parts:
                continue
            if len(parts) not in {16, 17}:
                raise ValueError(f"Expected 16 or 17 BLAST columns, got {len(parts)} in: {path}")
            (
                qseqid,
                sseqid,
                pident,
                length,
                _mismatch,
                _gapopen,
                qstart,
                qend,
                sstart,
                send,
                evalue,
                bitscore,
                qlen,
                _slen,
                qcovs,
                sstrand,
                *rest,
            ) = parts[:17]
            stitle = rest[0] if rest else ""
            genome_start = min(int(sstart), int(send))
            genome_end = max(int(sstart), int(send))
            rows.append(
                {
                    "Query_ID": qseqid,
                    "Chr": chr_from_subject(sseqid, stitle),
                    "Subject_ID": sseqid,
                    "Genome_start": str(genome_start),
                    "Genome_end": str(genome_end),
                    "S_strand": "plus" if sstrand == "plus" else "minus",
                    "Pident": f"{float(pident):.3f}",
                    "Qcovs": str(int(round(float(qcovs)))),
                    "Aln_Length": str(int(length)),
                    "Q_length": str(int(qlen)),
                    "Q_start": str(int(qstart)),
                    "Q_end": str(int(qend)),
                    "S_start": str(int(sstart)),
                    "S_end": str(int(send)),
                    "Evalue": evalue,
                    "Bitscore": bitscore,
                    "Subject_description": stitle,
                }
            )
    rows.sort(
        key=lambda r: (
            r["Query_ID"],
            -float(r["Bitscore"]),
            -float(r["Pident"]),
            -float(r["Qcovs"]),
            -int(r["Aln_Length"]),
            r["Chr"],
            int(r["Genome_start"]),
            int(r["Genome_end"]),
        )
    )
    return rows


def best_score(row: dict[str, str]) -> tuple[float, float, float, int]:
    return (
        float(row["Bitscore"]),
        float(row["Pident"]),
        float(row["Qcovs"]),
        int(row["Aln_Length"]),
    )


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--all-hits", required=True, type=Path)
    parser.add_argument("--best-hit", required=True, type=Path)
    parser.add_argument("--full-length-pident100", required=True, type=Path)
    parser.add_argument("--no-hit", required=True, type=Path)
    args = parser.parse_args()

    query_ids = read_fasta_ids(args.query_fasta)
    rows = parse_blast(args.blast)

    by_query: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_query.setdefault(row["Query_ID"], []).append(row)

    best_rows = []
    for query_id in sorted(by_query):
        best_rows.append(sorted(by_query[query_id], key=lambda r: (-best_score(r)[0], -best_score(r)[1], -best_score(r)[2], -best_score(r)[3], r["Chr"], int(r["Genome_start"])))[0])

    full_rows = [
        row for row in rows
        if float(row["Pident"]) == 100.0
        and int(row["Qcovs"]) == 100
        and int(row["Aln_Length"]) == int(row["Q_length"])
        and int(row["Q_start"]) == 1
        and int(row["Q_end"]) == int(row["Q_length"])
    ]

    hit_ids = set(by_query)
    nohit_ids = [qid for qid in query_ids if qid not in hit_ids]

    write_table(args.all_hits, rows)
    write_table(args.best_hit, best_rows)
    write_table(args.full_length_pident100, full_rows)
    with args.no_hit.open("w") as f:
        for qid in nohit_ids:
            f.write(f"{qid}\n")

    print(f"All hits: {len(rows)} -> {args.all_hits}")
    print(f"Best-hit queries: {len(best_rows)} -> {args.best_hit}")
    print(f"Full-length pident100 hits: {len(full_rows)} -> {args.full_length_pident100}")
    print(f"No-hit queries: {len(nohit_ids)} -> {args.no_hit}")


if __name__ == "__main__":
    main()
