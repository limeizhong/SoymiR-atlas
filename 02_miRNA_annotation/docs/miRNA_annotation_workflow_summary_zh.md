# miRNA 注释流程

本文档用于复现 `2814_precusor_miRNAs.txt` 中 2814 条预测 miRNA 前体定位记录的分类注释。流程以 mature 序列相似性、数据库 hairpin 在 ZH13.v2 基因组上的前体区间重叠、位置变体簇和 CD-HIT 序列簇为主要证据，依次判断参考库匹配位点、参考位点变体、未注释 opposite-arm 产物、已知家族新位点和新家族新位点。

## 1. 输入文件

核心输入：

```text
2814_precusor_miRNAs.txt
1588_mature_miRNAs.fasta
miRbase-mature.fa
pmiren_all_mature_clean.fa
miRbase_gma_position.txt
gma-hairpin.fa
pmiren_gmax_hairpin.fa
2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv
1588_mature_miRNAs_cdhit_est_c0.8.clstr
../../01_reference_genome/GWHAAEV00000000.1.genome.fasta
```

成熟序列 BLAST 结果作为已生成输入使用：

```text
1588_mature_miRNAs_vs_miRbase_mature_blastn_e1e-4.tsv
1588_mature_miRNAs_vs_pmiren_mature_blastn_e1e-4.tsv
```

其中 mature BLAST 使用全物种 miRBase 和 pmiREN mature 数据库，因此既用于大豆数据库注释，也用于非大豆命中支持的保守性判断。pmiREN 前体锚定直接使用当前的 `pmiren_gmax_hairpin.fa`，不再从 `pmiren-hairpin.fa` 重新抽取覆盖。

## 2. 数据库前体在 ZH13 上的锚定

miRBase 和 pmiREN soybean hairpin 序列分别 BLAST 到 ZH13.v2 基因组：

```bash
bash run_hairpin_ZH13_blast_unique_anchor.sh
```

该脚本执行以下步骤：

```text
1. 为 ZH13.v2 genome 建立 BLAST 数据库。
2. 用 gma-hairpin.fa 和 pmiren_gmax_hairpin.fa 分别进行 blastn 搜索，E-value = 1e-10。
3. 将 BLAST 结果整理为 all-hit、best-hit、full-length 100% hit 和 no-hit 文件。
4. 调用 resolve_hairpin_best_hits_unique_intervals.py，为每个数据库 precursor 确定一个最终 ZH13 best-hit 区间。
```

下游注释使用：

```text
gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
```

best-hit 选择规则：

```text
1. 每个数据库 precursor 只选择一个全局最高分 genome hit。
2. 若多个最高分 hit 并列，优先选择与数据库原始染色体和链方向一致的 hit。
3. 若仍并列，按染色体、起点和终点排序，选择排序最靠前的一条。
4. 不使用较低分 hit 替代全局最高分 hit。
```

## 3. 数据库 precursor anchor 判定

预测前体记录与数据库 precursor 同时满足以下条件时，才可建立已知成员 anchor：

```text
1. 染色体相同
2. 链方向相同
3. 预测前体区间 [H_start, H_end] 与数据库 hairpin 的 ZH13 区间 [Genome_start, Genome_end] 至少重叠 1 bp
```

若多个同 family 数据库 precursor 均可作为 anchor，优先选择与预测前体重叠 bp 数最大的 precursor；若仍相同，再按数据库坐标和名称排序。

最终输出前再次执行数据库 anchor 一致性检查：`reference_matched` 和 `reference_locus_variant` 的最终注释名称若对应一个大豆数据库 precursor，则该预测前体必须与该数据库 precursor 的 ZH13 区间至少重叠 1 bp。若不重叠，该记录不能保留为该参考成员或其 `.v` 变体，而应撤回为同 family 的新位点/新位点变体。

## 4. 注释优先级

注释按证据优先级逐层进行。低优先级结果不覆盖高优先级结果；miRBase 先于 pmiREN，pmiREN 先于本研究 no-hit 候选。

```text
1. miRBase 大豆 full-length 100% mature 命中
2. miRBase 大豆 pident = 100 但非 full-length mature 命中
3. miRBase 大豆非 100% mature 相似命中
4. pmiREN 大豆 full-length 100% mature 命中
5. pmiREN 大豆 pident = 100 但非 full-length mature 命中
6. pmiREN 大豆非 100% mature 相似命中
7. SoymiR no-hit 候选命名
```

miRBase 和 pmiREN 的其它物种 mature 命中不作为 family 注释层级使用；它们仅用于 `Conservation` 和 `Conserved_Species_Count`。

每类数据库证据完成后，立即执行共有检查：

```text
1. 位置变体传播：同一位置变体簇内 mature 区间高度重叠的未注释记录，随当前 anchor 注释为 .v 变体。
2. 同前体 arm 查询：同一预测前体或高度重叠预测前体内，若已有明确 anchor，另一组 mature 区域只有同时满足“mature 区域不重叠、双方 mature 区域均落入对方预测前体范围内、相对前体位置分处两端”并有 family、arm 或 CD-HIT 证据支持时，才注释为未报道 opposite-arm 产物。
```

SoymiR no-hit 命名和 CD-HIT family 回收后，再执行一次共有检查。若同一 mature 序列在任一预测位置已被注释为明确已知 family，则其它完全相同 mature 序列的 no-hit 位点继承该 family；没有对应数据库 precursor anchor 时，注释为 `known_family_new_locus`，而不是 `novel_family_new_locus`。随后再次检查同一或高度重叠预测前体上的两端 mature 区域，避免 CD-HIT 回收之后遗漏 5p/3p opposite-arm 关系。

位置变体簇仅来自 `2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv`，不使用数据库物理坐标判断 mature 区域是否重叠。

局部 tandem miRNA 的顺序映射只用于多个预测前体和多个数据库 precursor 在同一区域内都存在实际 ZH13 区间重叠的情况；该步骤不能把不重叠的远距离预测前体强行映射到某个数据库 precursor。

## 5. 注释分类和派生属性

最终输出同时保留 `annotation_category_original` 和规范化后的 `annotation_category`。`annotation_category_original` 记录原始流程内部类别，便于追溯；正式统计和作图使用 `annotation_category`。

`reference_matched`：mature 序列与大豆数据库 mature 全长 100% 一致，且预测前体与对应数据库 precursor 在 ZH13 上满足同染色体、同链和前体区间重叠。

`reference_locus_variant`：query mature 与参考 mature 属于同一 mature 区域或同一 arm，可解释为参考位点的 isomiR 或成熟序列位置变体；同时预测前体可由同 family、同染色体、同链、前体区间重叠的数据库 precursor 解释。

`unannotated_opposite_arm_product`：同一预测前体或高度重叠预测前体中已有非 opposite-arm anchor，另一组 mature 与 anchor mature 区域不重叠，双方 mature 区域均落在对方预测前体范围内，且相对位置分处前体两端，并可由 family、arm 或 CD-HIT 证据解释为该 precursor 的另一条 5p/3p arm。

`known_family_new_locus`：query 有已知 family 证据，且该 family 在当前注释来源的大豆数据库中存在，但没有合适的参考成员 anchor。

`novel_family_new_locus`：query 无法被 miRBase 或 pmiREN 的大豆已知 family 可靠解释，最终作为本研究新发现 family 成员保留。

对应变体类别为：

```text
known_family_new_locus_variant
novel_family_new_locus_variant
unannotated_opposite_arm_variant
```

`.v` 只表示同一 mature 区域或同一 arm 的 isomiR/位置变体，按 `Seq-ID` 中 `_R数字_` 的 reads 数降序编号。

三个可独立统计的派生字段为：

```text
locus_class: reference / known-family-new / novel-family-new
variant_status: anchor / variant
arm_status: anchor / opposite
```

前六类的 `locus_class` 由类别直接确定：`reference_*` 为 `reference`，`known_family_new_*` 为 `known-family-new`，`novel_family_new_*` 为 `novel-family-new`。`*_variant` 的 `variant_status = variant`，其它为 `anchor`；前六类的 `arm_status = anchor`。

`unannotated_opposite_arm_product` 和 `unannotated_opposite_arm_variant` 的 `arm_status = opposite`，`variant_status` 分别为 `anchor` 和 `variant`；其 `locus_class` 优先继承同一 precursor locus 上对应 anchor mature product 的 `locus_class`。若未能回溯到 anchor，但记录来源为 miRBase 或 pmiREN，则按数据库来源赋为 `reference`；仍无法确定时才写为 `NA` 并输出到 `2814_annotation_category_normalization_issues.tsv` 供人工检查。

## 6. 命名规则

通用命名规则：

```text
1. 注释名称使用 mature miRNA 名称，不使用 precursor 名称。
2. miRBase 阶段只避开大豆 miRBase 已占用成员名。
3. pmiREN 阶段只避开大豆 pmiREN 已占用成员名。
4. miRBase 与 pmiREN 的同名编号不强行一一对应；跨来源同名冲突只顺延自动生成的新成员名，reported 数据库名称保持不变。
5. 成员名比较时忽略大小写差异，避免 Gma-MIRxxx 与 Gma-miRxxx 被当作两个独立编号。
6. 输出表中的 Family 列统一为 mature-family 形式并忽略大小写差异：Gma/gma 统一写为 gma，MIR/MIRN 统一写为 miR/miRN；例如 Gma-MIR12410、gma-MIR12410 和 Gma-miR12410 均写为 gma-miR12410。
7. 同一位置变体簇内，reads 最高者保留基础名，其余按 reads 降序写 .v1、.v2 ...
8. 同一 predicted precursor 或高度重叠 predicted precursors 上若有两个符合 5p/3p 两端位置的不重叠 mature arms，二者属于同一 locus；已有 anchor 保留原名，另一端写为 -5p 或 -3p。
```

未报道 mature arm 的补充限制：

```text
1. 仅凭“同一前体或前体高度重叠且 mature 区间不重叠”不足以判定为 mature arm；必须同时满足前体 overlap >= 80%、双方 mature 区域均落入对方预测前体范围内、两组 mature 相对前体位置分处两端，并有 family、arm 或 CD-HIT 支持。
2. 若同一预测前体或高度重叠预测前体内存在多个互不重叠 mature 区域，只保留与 anchor 构成典型 5p/3p 两端关系且证据最充分的区域作为 `unannotated_opposite_arm_product`；其它区域不自动作为 arm 注释。
3. 同一 inferred arm 内，只有属于同一位置变体簇或同一 CD-HIT cluster 的记录，才作为 .v variant。
4. 若记录不在同一位置变体簇，也不在同一 CD-HIT cluster，则不强行合并为同一 arm variant，而拆分为独立 known_family_new_locus。
5. 最终输出前检查同一 source、同一 anchor、同一 arm 名称下是否存在重复基础名；若可由位置变体簇或 CD-HIT cluster 支持，则统一为 base + .v，否则拆分为独立成员名。
```

SoymiR no-hit 命名规则：

```text
1. 剩余 no-hit 记录按 CD-HIT cluster 归入 SoymiR family。
2. family 主编号按家族成员数降序分配：gma_miRNC0001、gma_miRNC0002 ...
3. 同一 CD-HIT family 内的独立基因组位点用 a/b/c 区分。
4. 同一位置变体位点内，reads 最高者保留基础名，其余写 .v1、.v2 ...
5. 同一 CD-HIT family 内两个 no-hit 位点只有在同染色体、同链、预测前体区间 overlap >= 80%、mature 区间不重叠、双方 mature 区域均落入对方预测前体范围内，且相对前体位置分处两端时，才合并为同一 locus 的两个 arms。
6. 若 no-hit 位点所在 CD-HIT cluster 只包含一个明确 miRBase/pmiREN family，则继承该 family；该 family 存在于对应大豆数据库时写 known_family_new_locus，否则写 novel_family_new_locus。
```

## 7. 报道状态与保守性

`Reported_Status = reported` 等价于 `annotation_category = reference_matched`；其它 `annotation_category` 均为 `unreported`。

`Conservation`：

```text
conserved：miRBase 或 pmiREN mature BLAST 中存在非大豆 subject 命中
specific：无非大豆 subject 命中
```

`Conserved_Species_Count` 为非大豆命中物种数；specific 记为 0。

## 8. 输出文件

核心注释输出包括：

```text
完整注释表：保留 evidence、matched mature、matched precursor 和 Matched_precursor_overlap_bp
简化注释表：保留补充表和作图所需核心字段
注释统计表：记录 source、annotation_category、reported status 和 conservation 频数
```

当前归档表：

```text
简化结果.xlsx
2026_合并_supplementary_tables.xlsx
2814_annotation_hierarchical_counts_source_reported_status.tsv
2814_annotation_category_rename_count_comparison.tsv
2814_annotation_category_normalization_issues.tsv
2814_annotation_category_normalization_summary.txt
```

简化注释表核心字段：

```text
Seq-ID
Sequences_Mature
Chr
M_start
M_end
H_start
H_end
Strand
annotation_category_original
annotation_category
locus_class
variant_status
arm_status
Annotation
miRNA_Locus
Reported_Status
Conservation
Conserved_Species_Count
Source
Family
Position_Variant_Cluster
CDHIT_Cluster_ID
```

## 9. 当前结果统计

```text
Input records: 2814
Family count after case-insensitive mature-family normalization: 928
miRBase ZH13 hairpin anchors loaded: 682
pmiREN ZH13 hairpin anchors loaded: 1347
```

Source counts:

```text
miRbase: 1177
pmiren: 511
soymir: 1126
```

annotation_category counts:

```text
reference_matched: 515
reference_locus_variant: 550
known_family_new_locus: 313
known_family_new_locus_variant: 79
unannotated_opposite_arm_product: 179
unannotated_opposite_arm_variant: 58
novel_family_new_locus: 1040
novel_family_new_locus_variant: 80
```

Derived attribute counts:

```text
locus_class:
reference: 1280
known-family-new: 408
novel-family-new: 1126

variant_status:
anchor: 2047
variant: 767

arm_status:
anchor: 2577
opposite: 237
```

Reported and conservation counts:

```text
reported: 515
unreported: 2299
conserved: 851
specific: 1963
```

## 10. 复现命令

```bash
cd SoymiR-atlas/02_miRNA_annotation
bash scripts/run_mature_miRNA_search_and_clustering.sh
python3 scripts/build_position_variant_clusters.py
GENOME=/path/to/GWHAAEV00000000.1.genome.fasta bash scripts/run_hairpin_ZH13_blast_unique_anchor.sh
python3 scripts/annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 scripts/normalize_annotation_categories.py
```
