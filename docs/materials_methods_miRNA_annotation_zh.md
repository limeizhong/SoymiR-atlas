# 材料与方法：miRNA 分类注释

## miRNA 数据集与参考数据库

本研究共获得 1,588 条非冗余成熟 miRNA 序列及其在大豆基因组上的 2,814 个候选前体定位记录。成熟 miRNA 序列保存于 `1588_mature_miRNAs.fasta`，候选前体及成熟序列定位信息保存于 `2814_precusor_miRNAs.txt`。参考数据库包括 miRBase mature miRNA、pmiREN mature miRNA、miRBase soybean hairpin precursor 和 pmiREN soybean hairpin precursor。

## mature miRNA 比对与序列聚类

以 1,588 条成熟 miRNA 序列为 query，分别对 miRBase 和 pmiREN mature 数据库进行 blastn 搜索，E-value 阈值为 `1e-4`。大豆 subject 命中用于已报道成员、已知成员变体和已知家族新成员判定；非大豆 subject 命中用于 family 证据和保守性判断。

使用 CD-HIT-EST 对 1,588 条成熟 miRNA 序列进行聚类，相似性阈值为 0.8。CD-HIT 聚类主要用于数据库证据耗尽后的 no-hit 序列命名和家族归并，不覆盖 miRBase 或 pmiREN 已支持的注释。

## 位置变体簇

根据 `2814_precusor_miRNAs.txt` 中成熟序列的预测位置构建位置变体簇。同一染色体、同一链、成熟区间重叠比例不低于 80% 的记录归为同一变体簇。位置变体簇仅用于同一 mature 区域或同一 mature arm 内的 `.v` 变体传播；同一前体上成熟区间不重叠的两组序列优先作为 5p/3p arm 关系检查。

## 数据库 precursor 的 ZH13 锚定

由于本研究统一使用中黄 13 基因组，miRBase 和 pmiREN 的 soybean hairpin precursor 序列先分别比对到 ZH13.v2 genome，BLAST E-value 为 `1e-10`。每个数据库 precursor 只保留一个 genome hit：先选择全局最优 hit；若多个最优 hit 并列，优先选择与数据库原始记录染色体和链方向一致的 hit；若仍并列，则按染色体和坐标顺序选择一个。

后续注释中，候选前体必须与数据库 precursor 的 ZH13 best-hit 区间同染色体、同链且发生前体区间重叠，才可视为同一 precursor anchor。

## 注释优先级与类别

注释按 miRBase、pmiREN、SoymiR no-hit 的顺序进行。miRBase 依次使用大豆 full-length 100% 命中、大豆 100% 非 full-length 命中和大豆非 100% 相似命中；pmiREN 依次使用大豆 full-length 100% 命中、大豆 100% 非 full-length 命中和其它物种 mature family 证据。高优先级注释不被低优先级证据覆盖。

`reported` 表示 mature 序列与大豆数据库 mature full-length 100% 一致，且前体与对应数据库 precursor 的 ZH13 锚定区间重叠。`known_member_variant` 表示已知数据库成员的 isomiR 或位置变体。`unreported_mature_arm` 表示同一预测前体或已注释 precursor anchor 上未被数据库 mature 记录覆盖的另一条成熟臂，其位置变体记为 `unreported_mature_arm_variant`。若 family 在对应数据库的大豆记录中存在但无可挂接的已知 precursor 成员，则记为 `known_family_new_member`；其位置变体记为 `known_family_new_member_variant`。miRBase 和 pmiREN 均不能解释的 no-hit 候选按 SoymiR 规则命名为 `new_family_member` 或 `new_family_member_variant`。

SoymiR no-hit 主编号采用 `gma_miRNC0001`、`gma_miRNC0002` 等格式，并按最终 SoymiR 家族成员数降序分配。同一家族内多个独立成员用 `a/b/c...` 区分；同一成熟区域的位置变体按 read abundance 降序编号为 `.v1/.v2...`。

## 报道状态与保守性

最终表中仅 `Status=reported` 的记录记为 `Reported_Status=reported`，其它变体、新 arm、新成员和 SoymiR no-hit 候选均记为 `unreported`。若 query 在 miRBase 或 pmiREN mature BLAST 结果中存在非大豆 subject 命中，则记为 `Conservation=conserved`；否则记为 `specific`。非大豆 subject 的去重物种数保存于 `Conserved_Species_Count`。

## 输出

注释流程由 `annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py` 实现，主要输出为：

```text
2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow.tsv
2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv
2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_summary.txt
simplified_annotation.xlsx
```

当前 2,814 条候选前体定位记录中，miRBase、pmiREN 和 SoymiR 来源分别为 1,198、518 和 1,098 条。状态统计为：`reported` 514 条，`known_member_variant` 549 条，`known_family_new_member` 310 条，`known_family_new_member_variant` 73 条，`unreported_mature_arm` 209 条，`unreported_mature_arm_variant` 66 条，`new_family_member` 1,015 条，`new_family_member_variant` 78 条。
