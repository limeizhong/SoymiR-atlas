# SoymiR-atlas

This repository archives the analysis workflow and key outputs for soybean miRNA identification, annotation, and downstream analyses.

The repository follows a module-based structure:

```text
01_miRNA_identification/
02_miRNA_annotation/
03_miRNA_targets/
04_miRNA_expression/
05_Figures/
```

The current completed modules are:

```text
01_miRNA_identification/  (raw-data example workflow)
02_miRNA_annotation/
03_miRNA_targets/
04_miRNA_expression/
05_Figures/
```

## Completed module

### 01_miRNA_identification

This module provides a migrated raw-data example workflow for identifying
soybean miRNAs from a representative sRNA-seq library with miRDP2. It documents
the upstream discovery step and archives a small example mature-miRNA output.

Run example:

```bash
bash 01_miRNA_identification/scripts/raw_data_example/02_mirna_download_qc.sh SRR27718796
bash 01_miRNA_identification/scripts/raw_data_example/03_mirna_mirdeep2.sh SRR27718796
bash 01_miRNA_identification/scripts/raw_data_example/04_mirna_mature.sh SRR27718796
```

### 02_miRNA_annotation

This module classifies and annotates 2814 non-redundant predicted soybean miRNA precursor-location records corresponding to 1588 mature miRNA sequences. It uses miRBase, pmiREN, position-variant clusters, and CD-HIT mature-sequence clusters to assign reported miRNAs, known member variants, known-family new members, unreported mature arms, and no-hit soymir candidates.

Run order:

```bash
bash 02_miRNA_annotation/scripts/run_mature_miRNA_search_and_clustering.sh
python3 02_miRNA_annotation/scripts/build_position_variant_clusters.py
GENOME=/path/to/GWHAAEV00000000.1.genome.fasta bash 02_miRNA_annotation/scripts/run_hairpin_ZH13_blast_unique_anchor.sh
python3 02_miRNA_annotation/scripts/annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 02_miRNA_annotation/scripts/update_simplified_result_and_supplementary_s2_s4.py
```

Final outputs are stored in:

```text
02_miRNA_annotation/results/
```

### 03_miRNA_targets

This module archives the final degradome-supported target prediction tables
used downstream, including 8,293 non-redundant miRNA-target pairs. It also
contains a migrated raw-data example workflow for processing representative
degradome libraries with CleaveLand4.

Final target outputs are stored in:

```text
03_miRNA_targets/results/results2025_239_binding_site_sequences/
```

### 04_miRNA_expression

This module rebuilds Supplementary Table S3 from the cleaned miRNA expression
table. For each library in which a miRNA is detected, read abundance is
normalized as reads per million total clean reads. Tissue-level expression is
calculated as the average RPM across detected libraries in the same tissue
category, and S3 reports `log2(mean RPM + 1)` values for 1,588 mature miRNAs
across 10 tissue categories.

Run:

```bash
python3 04_miRNA_expression/scripts/build_s3_tissue_expression_profiles.py
```

Final outputs are stored in:

```text
04_miRNA_expression/results/
```

### 05_Figures

This module archives the plotting scripts, plotting inputs, and final files for
manuscript Figures 1-4. Figures 1 and 2 can be regenerated directly from the
final supplementary workbook.

Run:

```bash
python3 05_Figures/scripts/Figure1/plot_figure1_from_supplementary_tables.py
python3 05_Figures/scripts/Figure2/plot_figure2_from_supplementary_tables.py
```

Final outputs are stored in:

```text
05_Figures/results/final_figures/
```

## Documentation

```text
docs/miRNA_annotation_workflow_summary.md
docs/miRNA_annotation_workflow_summary_zh.md
docs/materials_methods_miRNA_annotation_zh.md
```

## Current annotation summary

```text
Input precursor-location records: 2814
miRNA loci: 1943
miRNA families: 911
miRBase annotations: 1198
pmiREN annotations: 518
soymir no-hit candidates: 1098
conserved records: 851
specific records: 1963
```
