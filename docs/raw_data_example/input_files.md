# 输入文件说明

## miRNA-seq

| 文件 | 来源 | 格式 | 说明 |
|---|---|---|---|
| SRR27718796.fastq.gz | NCBI SRA | FASTQ | PE150 |
| SRR27718796.clean.fastq.gz | fastp | FASTQ | 19-24nt 过滤后 |

## Degradome

| 文件 | 来源 | 格式 | 说明 |
|---|---|---|---|
| SRR23932132_1.fastq.gz | NCBI SRA | FASTQ | SE47 |
| SRR23932132.clean.fastq.gz | fastp | FASTQ | 无长度截断 |

## 参考数据

| 文件 | 格式 | 说明 |
|---|---|---|
| data/reference/genome/GWHAAEV...genome.fasta | FASTA | 大豆基因组 (58 chr) |
| data/reference/cds/GWHAAEV...CDS.simple.fasta | FASTA | 清洗后 CDS (CleaveLand4用) |
