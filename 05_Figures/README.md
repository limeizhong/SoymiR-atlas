# 05_Figures

This module archives the plotting scripts, plotting inputs, and final figure
files for the SoymiR manuscript figures.

## Directory contents

```text
scripts/
input/
results/
```

## Scripts

```text
scripts/Figure1/plot_figure1_from_supplementary_tables.py
scripts/Figure2/compose_figure2_annotation_and_confidence.py
scripts/Figure2/plot_figure2_from_supplementary_tables.py
scripts/Figure2/plot_status_star_distribution.py
scripts/Figure3/analyze_global_isomiR_expression_library_overlap.py
scripts/Figure3/update_isomiR_composition_detection_rich_families_with_global_panel_e.py
scripts/Figure4/update_figure4_from_supplementary_tables.py
scripts/misc/plot_miRNA_annotation_source_status_100pct_bar.R
```

## Inputs

```text
input/source_tables/2026_merged_supplementary_tables.xlsx
input/plotting_data/Figure2/
input/plotting_data/Figure3/
input/plotting_data/Figure4/
input/plotting_data/misc/
```

`Figure 1` is regenerated directly from the final supplementary workbook. For
the manuscript comparison panel, miRBase and PmiREN retain the database-native
counts used in the original figure, while SoymiR values are updated from the
current supplementary tables.

## Run

From the repository root:

```bash
python3 05_Figures/scripts/Figure1/plot_figure1_from_supplementary_tables.py
python3 05_Figures/scripts/Figure2/compose_figure2_annotation_and_confidence.py
```

Figures 1 and 2 are regenerated from the final supplementary workbook and the
latest normalized annotation table in `../02_miRNA_annotation/results/`.
Figure 2 is a two-panel composite:

```text
Panel A: scripts/Figure2/plot_figure2_from_supplementary_tables.py
Panel B: scripts/Figure2/plot_status_star_distribution.py
Composite: scripts/Figure2/compose_figure2_annotation_and_confidence.py
```

Panel B uses
`input/plotting_data/Figure2/confidence_star_distribution_input.xlsx`.

Figure 3 is regenerated from `input/source_tables/2026_merged_supplementary_tables.xlsx`
plus the Figure 3 plotting inputs, including
`input/plotting_data/Figure3/68282.1588.expression.rawdata.txt`,
`input/plotting_data/Figure3/gma_mature_wo_U.fa`, and
`input/plotting_data/Figure3/global_expression_library_overlap/`.

```bash
python3 05_Figures/scripts/Figure3/analyze_global_isomiR_expression_library_overlap.py
python3 05_Figures/scripts/Figure3/update_isomiR_composition_detection_rich_families_with_global_panel_e.py
```

This command writes both the intermediate working figure and the manuscript
final figure:

```text
results/intermediate_figures/Figure3/Figure_isomiR_composition_detection_rich_families_global_panel_E_600dpi.png
results/final_figures/Figure_3.png
```

Figure 4 is regenerated from the same final supplementary workbook and writes
both the intermediate working figure and final manuscript figure:

```bash
python3 05_Figures/scripts/Figure4/update_figure4_from_supplementary_tables.py
```

```text
results/intermediate_figures/Figure4/Figure_miRNA_target_support_count_overlap_ABCD_600dpi.png
results/final_figures/Figure_4.png
```

Only scripts and input/output files required to reproduce the current four
final manuscript figures are retained.

## Final outputs

```text
results/final_figures/Figure_1.png
results/final_figures/Figure_2.png
results/final_figures/Figure_3.png
results/final_figures/Figure_4.png
```
