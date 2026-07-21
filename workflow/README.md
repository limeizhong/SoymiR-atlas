# workflow

This directory stores workflow-level entry points.

## Raw-data example workflow

`raw_data_example/` contains the migrated example workflow for reproducing the
two raw-data analysis components that precede the curated atlas tables:

1. sRNA-seq to miRNA discovery with miRDP2.
2. Degradome FASTQ to miRNA target prediction with CleaveLand4.

Run order:

```bash
bash workflow/raw_data_example/01_prepare_reference.sh
bash workflow/raw_data_example/02_mirna_download_qc.sh SRR27718796
bash workflow/raw_data_example/03_mirna_mirdeep2.sh SRR27718796
bash workflow/raw_data_example/04_mirna_mature.sh SRR27718796
bash workflow/raw_data_example/05_deg_download_qc.sh SRR23932132 SRR29504124
bash workflow/raw_data_example/06_deg_mode1.sh SRR23932132
bash workflow/raw_data_example/06_deg_mode1.sh SRR29504124
bash workflow/raw_data_example/07_deg_mode4.sh SRR23932132 SRR29504124
```

Or use:

```bash
bash workflow/raw_data_example/run_all.sh
```

Large FASTQ, SRA, genome, CDS, index, and third-party software files are not
committed. Configure their paths in `workflow/raw_data_example/config.sh`.
