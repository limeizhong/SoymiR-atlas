# Input files

The completed annotation module uses the following inputs:

```text
02_miRNA_annotation/input/1588_mature_miRNAs.fasta
02_miRNA_annotation/input/2814_precusor_miRNAs.txt
02_miRNA_annotation/input/miRbase-mature.fa
02_miRNA_annotation/input/pmiren-mature.fa
02_miRNA_annotation/input/miRbase_gma_position.txt
02_miRNA_annotation/input/pmiren_gma_position.txt
```

BLAST and CD-HIT intermediate evidence files are stored under:

```text
02_miRNA_annotation/results/intermediate/
```

The position-variant cluster file is generated from `input/2814_precusor_miRNAs.txt` by:

```bash
python3 02_miRNA_annotation/scripts/build_position_variant_clusters.py
```
