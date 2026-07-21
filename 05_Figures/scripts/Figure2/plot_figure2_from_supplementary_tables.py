#!/usr/bin/env python3
"""Regenerate Figure 2 from the normalized miRNA annotation categories."""

from pathlib import Path
import subprocess

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parents[1]
R_SCRIPT = MODULE_DIR / "scripts/misc/plot_miRNA_annotation_source_status_100pct_bar.R"
SIMPLIFIED_ANNOTATION = MODULE_DIR.parent / "02_miRNA_annotation/results/simplified_annotation.xlsx"
COUNT_INPUT = MODULE_DIR / "input/plotting_data/misc/2814_annotation_hierarchical_counts_source_annotation_category.tsv"

CATEGORY_ORDER = [
    "reference_matched",
    "reference_locus_variant",
    "known_family_new_locus",
    "known_family_new_locus_variant",
    "novel_family_new_locus",
    "novel_family_new_locus_variant",
    "unannotated_opposite_arm_product",
    "unannotated_opposite_arm_variant",
]
SOURCE_ORDER = ["miRbase", "pmiren", "soymir"]


def refresh_source_annotation_category_counts() -> None:
    """Refresh Figure 2 plotting data from the latest simplified annotation."""
    category_rank = {category: idx for idx, category in enumerate(CATEGORY_ORDER)}
    source_rank = {source: idx for idx, source in enumerate(SOURCE_ORDER)}

    df = pd.read_excel(SIMPLIFIED_ANNOTATION, dtype=str).fillna("")
    counts = df.groupby(["Source", "annotation_category"]).size().reset_index(name="Count")
    counts = counts.sort_values(
        by=["Source", "annotation_category"],
        key=lambda col: col.map(source_rank if col.name == "Source" else category_rank),
    )

    COUNT_INPUT.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(COUNT_INPUT, sep="\t", index=False)


if __name__ == "__main__":
    refresh_source_annotation_category_counts()
    subprocess.run(["Rscript", str(R_SCRIPT)], check=True)
