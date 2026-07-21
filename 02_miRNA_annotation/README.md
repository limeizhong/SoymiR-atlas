# 02_miRNA_annotation

This module reproduces the comparative annotation of the 2,814 predicted soybean miRNA precursor-location records.

The current workflow uses database precursor overlap on the Zhonghuang 13 genome to define known precursor anchors.

## Inputs

```text
input/1588_mature_miRNAs.fasta
input/2814_precusor_miRNAs.txt
input/miRbase-mature.fa
input/pmiren_all_mature_clean.fa
input/gma-hairpin.fa
input/pmiren_gmax_hairpin.fa
input/miRbase_gma_position.txt
input/pmiren_gma_position.txt
```

## Run Order

```bash
bash scripts/run_mature_miRNA_search_and_clustering.sh
python3 scripts/build_position_variant_clusters.py
GENOME=/path/to/GWHAAEV00000000.1.genome.fasta bash scripts/run_hairpin_ZH13_blast_unique_anchor.sh
python3 scripts/annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 scripts/normalize_annotation_categories.py
```

`run_hairpin_ZH13_blast_unique_anchor.sh` uses BLAST E-value `1e-10`. For each database precursor, one global best genome hit is selected. If several top-score hits tie, the hit matching the original database chromosome and strand is preferred.

## Intermediate Evidence

```text
results/intermediate/1588_mature_miRNAs_vs_miRbase_mature_blastn_e1e-4.tsv
results/intermediate/1588_mature_miRNAs_vs_pmiren_mature_blastn_e1e-4.tsv
results/intermediate/1588_mature_miRNAs_cdhit_est_c0.8.clstr
results/intermediate/2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv
results/intermediate/gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
```

## Final Outputs

```text
results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow.tsv
results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv
results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_summary.txt
results/intermediate/2814_annotation_category_rename_count_comparison.tsv
results/intermediate/2814_annotation_category_normalization_issues.tsv
results/simplified_annotation.xlsx
results/简化结果.xlsx
../05_Figures/input/source_tables/2026_merged_supplementary_tables.xlsx
```

## Current Summary

```text
Input records: 2814

annotation_category counts:
reference_matched: 515
reference_locus_variant: 550
known_family_new_locus: 313
known_family_new_locus_variant: 79
novel_family_new_locus: 1040
novel_family_new_locus_variant: 80
unannotated_opposite_arm_product: 179
unannotated_opposite_arm_variant: 58

locus_class counts:
reference: 1280
known-family-new: 408
novel-family-new: 1126

variant_status counts:
anchor: 2047
variant: 767

arm_status counts:
anchor: 2577
opposite: 237

Source counts:
miRbase: 1177
pmiren: 511
soymir: 1126
```

The input filename keeps the project spelling `precusor`.
