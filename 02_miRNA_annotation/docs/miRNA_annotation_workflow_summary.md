# miRNA Annotation Workflow: ZH13 Precursor-Overlap Unique-Anchor Version

This document describes the workflow used to annotate the 2814 predicted soybean miRNA precursor-location records in `2814_precusor_miRNAs.txt`. This version uses precursor overlap on the Zhonghuang 13 genome as the genomic anchor criterion. It no longer uses the former same-chromosome, same-strand, precursor-midpoint distance <= 10 Mb rule. The final output uses normalized `annotation_category` labels and keeps the original internal category in `annotation_category_original`.

## 1. Inputs

```text
2814_precusor_miRNAs.txt
1588_mature_miRNAs.fasta
miRbase-mature.fa
miRbase_gma_position.txt
gma-hairpin.fa
pmiren_all_mature_clean.fa
pmiren_gma_position.txt
pmiren_gmax_hairpin.fa
2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv
1588_mature_miRNAs_cdhit_est_c0.8.clstr
../../data/reference/genome/GWHAAEV00000000.1.genome.fasta
```

The mature BLAST results are used as precomputed inputs:

```text
1588_mature_miRNAs_vs_miRbase_mature_blastn_e1e-4.tsv
1588_mature_miRNAs_vs_pmiren_mature_blastn_e1e-4.tsv
```

## 2. Database Precursor Anchoring on ZH13

miRBase and pmiREN soybean hairpin sequences are mapped to the ZH13.v2 genome with:

```text
run_hairpin_ZH13_blast_unique_anchor.sh
```

The script runs `blastn` with `E-value = 1e-10`, formats all hits, selects best hits, reports full-length 100% hits and no-hit records, and then resolves a final unique best-hit interval for each database precursor.

The downstream annotation uses:

```text
gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
```

Best-hit selection:

```text
1. Select one global best genome hit for each database precursor.
2. If multiple equivalent best hits exist, prefer the hit matching the original database chromosome and strand.
3. If ties remain, sort by chromosome, start, and end, and retain the first hit.
4. Lower-scoring hits are not used to replace a global best hit.
```

## 3. Precursor-Overlap Criterion

A predicted precursor can be anchored to a database precursor only when all conditions are met:

```text
1. Same chromosome
2. Same strand
3. Predicted precursor interval [H_start, H_end] overlaps the database hairpin ZH13 interval [Genome_start, Genome_end] by at least 1 bp
```

If multiple same-family database precursors are compatible, the one with the largest overlap in bp is preferred; remaining ties are resolved by database coordinate and name.

Before final output, database-anchor consistency is checked again. If the final name of a `reference_matched` or `reference_locus_variant` record corresponds to a soybean database precursor, the predicted precursor must overlap that database precursor on ZH13 by at least 1 bp. Otherwise, the record cannot remain assigned to that reference member or its `.v` variant and is withdrawn to a same-family new-locus/new-locus variant.

## 4. Annotation Priority

Lower-priority evidence does not overwrite higher-priority annotations. miRBase is applied before pmiREN, and pmiREN before SoymiR no-hit naming.

```text
1. miRBase soybean full-length 100% mature hit
2. miRBase soybean pident = 100 but non-full-length mature hit
3. miRBase soybean non-100% mature similarity hit
4. pmiREN soybean full-length 100% mature hit
5. pmiREN soybean pident = 100 but non-full-length mature hit
6. pmiREN soybean non-100% mature similarity hit
7. SoymiR no-hit candidate naming
```

After each evidence class, two shared checks are performed:

```text
1. Position-variant propagation: unannotated records in the same position-variant cluster with highly overlapping mature intervals inherit the current anchor as .v variants. If a later step promotes the cluster anchor from a SoymiR novel-family assignment to a miRBase/pmiREN known-family assignment, remaining low-priority variants in the same position cluster inherit the promoted source, family, and name. This synchronized propagation is not applied when the propagated database member would lack precursor overlap with the corresponding database precursor.
2. Same-precursor arm search: a second mature region is annotated as an unannotated opposite-arm product only when it has a clear anchor on the same or highly overlapping predicted precursor, the two mature regions do not overlap, each mature region falls within the other predicted precursor interval, the two mature regions occupy opposite precursor ends, and family, arm, or CD-HIT evidence supports the relationship.
```

After SoymiR no-hit naming and CD-HIT-based known-family rescue, the shared checks are repeated. If an identical mature sequence has already been assigned to a clear known family at any predicted location, other no-hit locations with the same mature sequence inherit that family. Without a matching database precursor anchor, these records are annotated as `known_family_new_locus` rather than `novel_family_new_locus`. Same-precursor opposite-arm detection is then repeated so that 5p/3p relationships created by late CD-HIT rescue are not missed.

Position-variant clusters are taken only from `2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv`; database physical coordinates are not used to judge mature-region overlap.

Tandem-order mapping is allowed only when each paired predicted precursor and database precursor in the local tandem block has a real ZH13 interval overlap. This step cannot force a distant non-overlapping predicted precursor onto a database precursor.

## 5. Annotation Categories and Derived Attributes

The final tables keep `annotation_category_original` for traceability and use the normalized `annotation_category` for statistics and figures.

`reference_matched`: the mature sequence is a full-length 100% match to a soybean database mature miRNA, and the predicted precursor overlaps the corresponding database precursor on the same ZH13 chromosome and strand.

`reference_locus_variant`: the query mature sequence belongs to the same mature region or arm as a reference member and can be explained by an overlapping same-family database precursor. Assignment to a specific reference mature-arm `.v` variant requires not only precursor overlap but also high mature-sequence similarity to that reference arm, or membership in the same mature-region interval. A same-family mature hit on an overlapping precursor is not sufficient when the query sequence is clearly dissimilar to the named reference arm.

The same restriction is applied during later local tandem-order mapping and unused database-precursor anchor rescue. These rescue steps may create `reference_locus_variant` records only when mature sequences are similar; otherwise the records are withdrawn to same-family `known_family_new_locus` or its variant category.

`unannotated_opposite_arm_product`: a non-overlapping mature product from the opposite arm of the same or highly overlapping predicted precursor where an anchor mature product already exists; both mature regions must fall within the paired predicted precursor intervals and occupy opposite precursor ends.

`known_family_new_locus`: known-family evidence is present and the family exists in the current soybean reference database, but no suitable reference-locus anchor is available.

`novel_family_new_locus`: the family is absent from miRBase and pmiREN soybean annotations, or the record remains a study-specific no-hit candidate.

The corresponding variant categories are `known_family_new_locus_variant`, `novel_family_new_locus_variant`, and `unannotated_opposite_arm_variant`. The `.v` suffix denotes an isomiR or mature-position variant within the same mature region or arm and is numbered by descending `_R` read count in `Seq-ID`.

Three derived attributes are added for independent summaries:

```text
locus_class: reference / known-family-new / novel-family-new
variant_status: anchor / variant
arm_status: anchor / opposite
```

For the six non-opposite-arm categories, `locus_class` is determined directly from the category: `reference_*` is `reference`, `known_family_new_*` is `known-family-new`, and `novel_family_new_*` is `novel-family-new`. Variant categories have `variant_status = variant`; all other non-opposite categories have `variant_status = anchor` and `arm_status = anchor`.

For `unannotated_opposite_arm_product` and `unannotated_opposite_arm_variant`, `arm_status = opposite` and `variant_status` is `anchor` or `variant`, respectively. Their `locus_class` is preferentially inherited from the corresponding anchor mature product on the same precursor locus. If no anchor can be traced but the record source is miRBase or pmiREN, the record is assigned `locus_class = reference`. Only records still lacking a source-supported assignment are written as `NA` and exported to `2814_annotation_category_normalization_issues.tsv`.

## 6. Naming Rules

```text
1. Use mature miRNA names rather than precursor names.
2. During the miRBase stage, avoid soybean miRBase member names already occupied.
3. During the pmiREN stage, avoid soybean pmiREN member names already occupied.
4. miRBase and pmiREN member suffixes are not forced to correspond one-to-one.
5. Family labels in output tables are normalized to case-insensitive mature-family form: Gma/gma is written as gma, and MIR/MIRN is written as miR/miRN. For example, Gma-MIR12410, gma-MIR12410, and Gma-miR12410 are all reported as gma-miR12410 in the Family column.
6. Within a position-variant cluster, the highest-read member keeps the base name; the others are named .v1, .v2, ...
7. If two non-overlapping mature regions occur on the same or highly overlapping predicted precursors and satisfy the opposite-end arm rule, they are treated as one locus; the anchored arm keeps the base name and the other arm is named -5p or -3p.
```

SoymiR no-hit naming:

```text
1. Remaining no-hit records are grouped by CD-HIT cluster.
2. SoymiR family numbers are assigned by descending family size: gma_miRNC0001, gma_miRNC0002, ...
3. Independent genomic loci within the same CD-HIT family receive a/b/c suffixes.
4. Mature-position variants are named .v1, .v2, ... by descending read count.
5. Two no-hit loci in the same CD-HIT family are merged as two arms of one locus only when they have the same chromosome and strand, precursor overlap >= 80%, non-overlapping mature intervals, mutual mature-in-precursor containment, and opposite-end mature positions.
6. If a no-hit locus belongs to a CD-HIT cluster containing one unambiguous miRBase/pmiREN family, it inherits that family. The category is `known_family_new_locus` only when the family exists in the corresponding soybean database; otherwise it remains `novel_family_new_locus`.
```

## 7. Reported Status and Conservation

`Reported_Status` is `reported` only when `annotation_category = reference_matched`; all other records are `unreported`.

`Conservation` is `conserved` when a non-soybean subject is detected in miRBase or pmiREN mature BLAST results, otherwise `specific`. `Conserved_Species_Count` is the number of non-soybean species with mature BLAST hits. The mature BLAST inputs are full multi-species databases, so this field reflects non-soybean miRBase/pmiREN hits.

## 8. Outputs

```text
2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow.tsv
2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv
2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_summary.txt
```

The simplified output contains:

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

## 9. Current Result Summary

```text
Input records: 2814
Family count after case-insensitive mature-family normalization: 925
miRBase ZH13 hairpin anchors loaded: 682
pmiREN ZH13 hairpin anchors loaded: 1347
```

Source counts:

```text
miRbase: 1178
pmiren: 521
soymir: 1115
```

annotation_category counts:

```text
reference_matched: 515
reference_locus_variant: 473
known_family_new_locus: 371
known_family_new_locus_variant: 109
unannotated_opposite_arm_product: 179
unannotated_opposite_arm_variant: 58
novel_family_new_locus: 1036
novel_family_new_locus_variant: 73
```

Derived attribute counts:

```text
locus_class:
reference: 1198
known-family-new: 501
novel-family-new: 1115

variant_status:
anchor: 2101
variant: 713

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

## 10. Reproducibility

```bash
cd SoymiR-atlas/02_miRNA_annotation
bash scripts/run_mature_miRNA_search_and_clustering.sh
python3 scripts/build_position_variant_clusters.py
GENOME=/path/to/GWHAAEV00000000.1.genome.fasta bash scripts/run_hairpin_ZH13_blast_unique_anchor.sh
python3 scripts/annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 scripts/normalize_annotation_categories.py
```
