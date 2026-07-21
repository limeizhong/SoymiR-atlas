# 输出文件说明

## 01_miRNA_identification

| 文件 | 说明 |
|---|---|
| results/SRR27718796_miRNA_query.fa | 15条 unique 成熟体 miRNA |
| results/mirdp2/*/filter_filter_P_prediction | 26条高置信度预测 |
| results/mirdp2/*/filter_predictions | 2747条全部预测 |

## 03_miRNA_targets — Mode1

| 文件 | 说明 |
|---|---|
| 04_results/{SRR}_mode1.txt | 靶标预测结果 + GSTAr/RNAplex 统计 |
| 04_results/{SRR}_sorted.bam | 降解组比对 BAM |
| 04_results/TPlots/ | TPlot 图 (显著靶标) |

## 03_miRNA_targets — Mode4

| 文件 | 说明 |
|---|---|
| 04_results/{SRR}_mode4.txt | 降解组密度 (按类别0-4排列) |
