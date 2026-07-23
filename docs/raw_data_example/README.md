# Raw-data example workflow

This directory documents the optional raw-data example workflow migrated from
the working analysis folder.

The example covers two upstream components:

1. sRNA-seq FASTQ to candidate mature miRNAs with miRDP2.
2. Degradome FASTQ to CleaveLand4 target predictions.

Main scripts:

```text
workflow/raw_data_example/
workflow/raw_data_example/
```

Example outputs:

```text
01_miRNA_identification/results/example/SRR27718796/
03_miRNA_targets/results/example_degradome/
```

Large raw sequencing files, reference FASTA files, Bowtie indexes, miRDP2, and
CleaveLand4 are not committed. Configure local paths in
`workflow/raw_data_example/config.sh`.

Additional notes:

```text
workflow_overview.md
input_files.md
output_files.md
```
