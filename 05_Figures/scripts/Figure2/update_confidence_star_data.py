#!/usr/bin/env python3
"""
Update Panel B input data (confidence_star_distribution_input.xlsx) from the
latest miRNA annotation workflow results.

Data flow:
  S7 data (副本S7(1).xlsx)
  ├── annotation_category  ←  updated from workflow TSV (this script)
  ├── Confidence (3-7★)   ←  carried over unchanged from S7
  └── other columns        ←  carried over unchanged from S7
       ↓
  confidence_star_distribution_input.xlsx  →  Figure 2 Panel B

Matching strategy (coordinate-first):
  1. Match by Seq-ID + hairpin coordinate (Chr_H_start_H_end) — the most
     stable record identifier across annotation re-runs.
  2. Fall back to Annotation name matching for any remaining unmatched records.

Only the annotation_category column is updated. Confidence stars and all other
columns from the S7 table are preserved as-is. When the annotation pipeline
regenerates 副本S7(1).xlsx with updated Confidence values, copy that file to
confidence_star_distribution_input.xlsx first, then run this script to sync
annotation_category.

Input:
  - 02_miRNA_annotation/results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv
  - 05_Figures/input/plotting_data/Figure2/confidence_star_distribution_input.xlsx

Output:
  - 05_Figures/input/plotting_data/Figure2/confidence_star_distribution_input.xlsx (updated in-place)

Usage:
  cd SoymiR-atlas
  python3 05_Figures/scripts/Figure2/update_confidence_star_data.py
"""

from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parents[1]  # 05_Figures/
ATLAS_ROOT = MODULE_DIR.parent      # SoymiR-atlas/

WORKFLOW_TSV = ATLAS_ROOT / "02_miRNA_annotation/results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv"
CONFIDENCE_XLSX = MODULE_DIR / "input/plotting_data/Figure2/confidence_star_distribution_input.xlsx"


def main() -> None:
    if not WORKFLOW_TSV.exists():
        raise FileNotFoundError(f"Workflow TSV not found: {WORKFLOW_TSV}")
    if not CONFIDENCE_XLSX.exists():
        raise FileNotFoundError(f"Confidence input not found: {CONFIDENCE_XLSX}")

    wf = pd.read_csv(WORKFLOW_TSV, sep="\t")
    conf = pd.read_excel(CONFIDENCE_XLSX)

    # Build coordinate key: Seq-ID + hairpin position
    wf["hairpin_key"] = (
        wf["Chr"].astype(str)
        + "_"
        + wf["H_start"].astype(str)
        + "_"
        + wf["H_end"].astype(str)
    )
    wf["coord_key"] = wf["Seq-ID"].astype(str) + "|" + wf["hairpin_key"]
    conf["coord_key"] = (
        conf["Seq-ID"].astype(str) + "|" + conf["Precusor_Coordinate"].astype(str)
    )

    # Step 1: coordinate matching (most stable)
    wf_coord_map = wf.groupby("coord_key")["annotation_category"].first().to_dict()
    conf["new_cat"] = conf["coord_key"].map(wf_coord_map)
    coord_matched = conf["new_cat"].notna().sum()

    # Step 2: Annotation name fallback
    wf_ann_map = dict(
        zip(wf["Annotation"].astype(str), wf["annotation_category"])
    )
    missing_mask = conf["new_cat"].isna()
    ann_filled = 0
    for idx in conf[missing_mask].index:
        ann = str(conf.at[idx, "Annotation"])
        if ann in wf_ann_map:
            conf.at[idx, "new_cat"] = wf_ann_map[ann]
            ann_filled += 1

    still_missing = missing_mask.sum() - ann_filled
    changed = (conf["annotation_category"] != conf["new_cat"]).sum()

    # Update and clean up
    conf["annotation_category"] = conf["new_cat"]
    conf = conf.drop(columns=["new_cat", "coord_key"])

    # Verify
    updated = conf["annotation_category"].value_counts()
    expected = wf["annotation_category"].value_counts()
    compare = (
        pd.DataFrame({"updated": updated, "expected": expected})
        .fillna(0)
        .astype(int)
    )
    compare["diff"] = compare["updated"] - compare["expected"]
    all_match = (compare["diff"] == 0).all()

    # Save
    conf.to_excel(CONFIDENCE_XLSX, index=False)

    # Report
    print(f"Coordinate matched : {coord_matched}/{len(conf)}")
    print(f"Annotation fallback: {ann_filled}")
    print(f"Still unmatched    : {still_missing}")
    print(f"Category changed   : {changed}/{len(conf)}")
    print()
    if all_match:
        print("✓ All annotation_category counts match the workflow summary.")
    else:
        print("⚠ Mismatches remain:")
        print(compare[compare["diff"] != 0].to_string())
    print(f"\nSaved: {CONFIDENCE_XLSX}")


if __name__ == "__main__":
    main()
