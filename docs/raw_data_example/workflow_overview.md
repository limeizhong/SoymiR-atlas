# 分析流程概览

```
raw FASTQ
  │
  ├── miRNA-seq (SRR27718796)
  │     │ 02_mirna_download_qc.sh
  │     ├── fastp (19-24nt) → clean FASTQ
  │     │ 03_mirna_mirdeep2.sh
  │     ├── seqkit fq2fa → FA
  │     ├── filter.fasta2input.format.py → 去重+19-24nt 筛选
  │     └── miRDP2 (--rpm 5) → miRNA 预测
  │     │ 04_mirna_mature.sh
  │     └── 从 filter_P_prediction 提取 → 15条 unique miRNA
  │
  └── Degradome (SRR23932132, SRR29504124)
        │ 05_deg_download_qc.sh  
        ├── fastp (不截断, 保留47bp) → clean FASTQ
        │ 06_deg_mode1.sh
        ├── FASTQ→FASTA → CleaveLand4 Mode1 → 靶标 + TPlot
        │       ├── dd.txt → 04_results/
        │       └── BAM → 04_results/{SRR}_sorted.bam
        │ 07_deg_mode4.sh
        └── dd.txt + GSTAr.txt → CleaveLand4 Mode4
```

## 模块说明

- **01_miRNA_identification**: sRNA-seq → miRDP2 → miRNA 预测
- **03_miRNA_targets**: 降解组 → CleaveLand4 → 靶标预测
- **workflow/**: 流水线脚本
- **workflow/raw_data_example/config.sh**: 配置参数
- **environment/**: Conda 环境配置
- **data/reference/**: 参考基因组 + CDS
- **metadata/**: 样本信息
