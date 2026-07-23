# Reproducibility

## Environment setup

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate soymir-atlas
```

Or install Python packages only:

```bash
pip install -r requirements.txt
```

The core analysis requires:
- Python 3.10+ with `openpyxl`, `pandas`, `numpy`, `scipy`, `matplotlib`, `Pillow`
- BLAST+ (`makeblastdb`, `blastn`) and CD-HIT-EST
- R 4.3+ with `ggplot2`, `dplyr`, `readr`, `scales` (for Figure 2)

The optional raw-data example workflow additionally requires:
- `bowtie`, `fastp`, `seqkit`, `samtools`, `sra-tools`, `viennarna`, `perl`
- miRDP2 v1.1.4 and CleaveLand4 v4.5 (manual installation; see `environment.yml`)

See `environment/software_versions.txt` and `environment/raw_data_example_software_versions.txt`
for recorded versions used during development.

## Core modules

The completed `02_miRNA_annotation` module can be reproduced with:

```bash
bash 02_miRNA_annotation/scripts/run_mature_miRNA_search_and_clustering.sh
python3 02_miRNA_annotation/scripts/build_position_variant_clusters.py
GENOME=/path/to/GWHAAEV00000000.1.genome.fasta bash 02_miRNA_annotation/scripts/run_hairpin_ZH13_blast_unique_anchor.sh
python3 02_miRNA_annotation/scripts/annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 02_miRNA_annotation/scripts/normalize_annotation_categories.py
```

The completed `04_miRNA_expression` module can be reproduced with:

```bash
python3 04_miRNA_expression/scripts/build_s3_tissue_expression_profiles.py
```

## Figures

Figures 1 and 2 can be regenerated from the supplementary workbook:

```bash
python3 05_Figures/scripts/Figure1/plot_figure1_from_supplementary_tables.py
python3 05_Figures/scripts/Figure2/compose_figure2_annotation_and_confidence.py
```

Figures 3 and 4 require additional plotting data under `05_Figures/input/plotting_data/`:

```bash
python3 05_Figures/scripts/Figure3/analyze_global_isomiR_expression_library_overlap.py
python3 05_Figures/scripts/Figure3/update_isomiR_composition_detection_rich_families_with_global_panel_e.py
python3 05_Figures/scripts/Figure4/update_figure4_from_supplementary_tables.py
```
