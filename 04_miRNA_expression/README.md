# 04_miRNA_expression

This module rebuilds Supplementary Table S3 from the final cleaned miRNA
expression table.

## Input

```text
input/68282.1588.expression.rawdata.txt
```

Columns:

```text
sample  mature_seq  Tissue  reads  Seq-ID  Reported_Status  Total_Reads
```

The input table contains 68,282 expression records for 1,588 non-redundant
mature miRNA sequences. Tissue names are normalized to the 10 tissue categories
used in S3.

The workbook template is:

```text
input/2026_merged_supplementary_tables_template.xlsx
```

## Method

For each library in which a miRNA is detected, miRNA abundance is normalized as
reads per million total clean reads:

```text
RPM = read_count / Total_Reads * 1,000,000
```

The read count is parsed from the `reads` field, for example
`read00023764_x471` contributes 471 reads. Tissue-level expression is calculated
as the average RPM across detected libraries in the same tissue category.
Undetected libraries are not included in the average; if a miRNA is not detected
in any library from a tissue category, that tissue value is set to 0.

S3 values are:

```text
Xi = log2(mean_RPM + 1)
```

Values are rounded to two decimal places.

## Run

From the repository root:

```bash
python3 04_miRNA_expression/scripts/build_s3_tissue_expression_profiles.py
```

## Output

```text
results/supplementary_table_S3_tissue_expression_profiles.tsv
results/2026_merged_supplementary_tables_S3_updated.xlsx
```

The output S3 table contains 1,588 mature miRNA IDs across 10 tissue categories:

```text
cotyledon, flower, leaves, mixed, nodules, pod, root, seed, shoot, stem
```
