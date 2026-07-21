# results2025_239 mature miRNA 靶基因结合位点序列抽取结果

从 `results2025` 中扫描 13,634 个 `.res.txt` 文件；其中 0 个没有 CleaveLand 命中记录，12,911 个包含结果。

老编号结果使用 `2026-miRNA/combine_data/id_transversion/1406_id_with_new.txt` 转换为新编号；映射为 NA 的 `.res.txt` 文件跳过 723 个。
同时合并已有 detail 表 `/private/tmp/results2025_239_binding_site_detail.before_merge.tsv` 中的 1,717 条记录。
从 CleaveLand Pretty 输出的 `5' ... 3' Transcript:` 行中解析出 77,016 条原始结合位点记录；按 `library_id + miRNA + target_gene + binding_site_sequence` 去重后为 42,547 条。
进一步按 `miRNA + target_gene + binding_site_sequence` 去冗余后，得到 8,301 条唯一记录，涉及 8,293 个 miRNA-target 对和 5,496 个靶基因。

主输出表：`results2025_239_binding_site_unique.tsv`，包含 `miRNA_ID`、`target_gene`、`binding_site_sequence`，并附带 `library_count`、`result_record_count` 和 `libraries` 便于后续筛选。
最小三列表：`results2025_239_binding_site_unique_minimal.tsv`。
基因组坐标表：`results2025_239_binding_site_genome_mapping_candidates.tsv`，坐标由 `2026-miRNA/01_reference_genome/GWHAAEV00000000.1.gff` 中对应 CDS 坐标换算得到。
