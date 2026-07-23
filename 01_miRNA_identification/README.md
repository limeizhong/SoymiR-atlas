# 01_miRNA_identification

This module contains an optional raw-data example workflow for identifying
soybean miRNAs from a representative sRNA-seq library with miRDP2.

The final atlas annotation used in the manuscript is stored in
`02_miRNA_annotation/`; this module documents how a representative miRNA
discovery step can be reproduced from raw sequencing data.

## Example workflow

Scripts are archived under:

```text
../workflow/raw_data_example/
```

Main steps:

```bash
bash ../workflow/raw_data_example/02_mirna_download_qc.sh SRR27718796
bash ../workflow/raw_data_example/03_mirna_mirdeep2.sh SRR27718796
bash ../workflow/raw_data_example/04_mirna_mature.sh SRR27718796
```

The scripts expect reference files and third-party tools to be available as
described in `../workflow/raw_data_example/config.sh`.

## Example result

```text
results/example/SRR27718796/SRR27718796_miRNA_query.fa
```

This example file contains the 15 non-redundant mature miRNA sequences extracted
from the miRDP2 high-confidence prediction output in the migrated example
workflow.

## Inputs not committed

Large primary files are not stored in the repository:

- raw sRNA-seq FASTQ/SRA files
- soybean genome FASTA and Bowtie indexes
- locally installed miRDP2 software

Place those files under the paths defined in `scripts/raw_data_example/config.sh`
or override the paths with environment variables before running the scripts.
