#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import csv


SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
INPUT_DIR = MODULE_DIR / "input"
INTERMEDIATE_DIR = MODULE_DIR / "results" / "intermediate"

INPUT_FILE = INPUT_DIR / "2814_precusor_miRNAs.txt"
OUTPUT_FILE = INTERMEDIATE_DIR / "2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv"
SUMMARY_FILE = INTERMEDIATE_DIR / "2814_precusor_miRNAs_mature_position_variant_clusters_overlap80_summary.txt"
OVERLAP_THRESHOLD = 0.8


def interval_overlap_ratio(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if overlap <= 0:
        return 0.0
    return overlap / min(a_end - a_start + 1, b_end - b_start + 1)


def chrom_key(chrom: str):
    return (0, int(chrom)) if chrom.isdigit() else (1, chrom)


def row_sort_key(row: dict[str, str]):
    return (
        chrom_key(row["Chr"]),
        row["Strand"],
        int(row["M_start"]),
        int(row["M_end"]),
        int(row["H_start"]),
        int(row["H_end"]),
        row["Seq-ID"],
    )


def connected_components(indices: list[int], rows: list[dict[str, str]]) -> list[list[int]]:
    parent = {i: i for i in indices}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for pos, i in enumerate(indices):
        a = rows[i]
        for j in indices[pos + 1:]:
            b = rows[j]
            ratio = interval_overlap_ratio(
                int(a["M_start"]),
                int(a["M_end"]),
                int(b["M_start"]),
                int(b["M_end"]),
            )
            if ratio >= OVERLAP_THRESHOLD:
                union(i, j)

    groups = defaultdict(list)
    for i in indices:
        groups[find(i)].append(i)
    components = [sorted(members, key=lambda x: row_sort_key(rows[x])) for members in groups.values()]
    components.sort(key=lambda members: row_sort_key(rows[members[0]]))
    return components


def main():
    INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

    with INPUT_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    by_chr_strand = defaultdict(list)
    for i, row in enumerate(rows):
        by_chr_strand[(row["Chr"], row["Strand"])].append(i)

    clusters = []
    for key in sorted(by_chr_strand, key=lambda x: (chrom_key(x[0]), x[1])):
        for component in connected_components(by_chr_strand[key], rows):
            if len(component) >= 2:
                clusters.append(component)

    clusters.sort(key=lambda members: row_sort_key(rows[members[0]]))

    fieldnames = [
        "Variant_Cluster",
        "Cluster_Size",
        "Seq-ID",
        "Sequences_Mature",
        "Chr",
        "M_start",
        "M_end",
        "H_start",
        "H_end",
        "Strand",
    ]
    with OUTPUT_FILE.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for n, members in enumerate(clusters, 1):
            cluster_id = f"V{n}"
            for i in members:
                row = rows[i]
                writer.writerow({
                    "Variant_Cluster": cluster_id,
                    "Cluster_Size": len(members),
                    "Seq-ID": row["Seq-ID"],
                    "Sequences_Mature": row["Sequences_Mature"],
                    "Chr": row["Chr"],
                    "M_start": row["M_start"],
                    "M_end": row["M_end"],
                    "H_start": row["H_start"],
                    "H_end": row["H_end"],
                    "Strand": row["Strand"],
                })

    size_counts = Counter(len(members) for members in clusters)
    records_in_clusters = sum(len(members) for members in clusters)
    with SUMMARY_FILE.open("w") as handle:
        handle.write("2814_precusor_miRNAs.txt mature-position variant cluster summary\n")
        handle.write(f"Input: {INPUT_FILE.name}\n")
        handle.write(f"Total_records: {len(rows)}\n")
        handle.write(f"Overlap_threshold_of_shorter_mature_interval: {OVERLAP_THRESHOLD}\n")
        handle.write(f"Variant_clusters_size>=2: {len(clusters)}\n")
        handle.write(f"Records_in_variant_clusters_size>=2: {records_in_clusters}\n")
        handle.write(f"Singleton_records_not_output: {len(rows) - records_in_clusters}\n\n")
        handle.write("Cluster_Size\tClusters\n")
        for size in sorted(size_counts):
            handle.write(f"{size}\t{size_counts[size]}\n")

    print(OUTPUT_FILE)
    print(SUMMARY_FILE)


if __name__ == "__main__":
    main()
