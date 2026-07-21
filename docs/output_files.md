# Output files

The completed annotation module produces:

```text
02_miRNA_annotation/results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow.tsv
02_miRNA_annotation/results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv
02_miRNA_annotation/results/2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_summary.txt
02_miRNA_annotation/results/simplified_annotation.xlsx
```

The full TSV includes evidence columns. The short TSV and Excel file are intended for reporting and downstream summary analyses.

The completed expression module produces:

```text
04_miRNA_expression/results/supplementary_table_S3_tissue_expression_profiles.tsv
04_miRNA_expression/results/2026_merged_supplementary_tables_S3_updated.xlsx
```

The S3 TSV reports `log2(mean RPM + 1)` tissue expression values for 1,588
mature miRNAs across 10 major tissue categories. The mean RPM is calculated
from libraries in which the miRNA is detected; tissue categories with no
detected libraries for that miRNA are reported as 0.
