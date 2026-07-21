# miRNA Annotation Workflow: ZH13 Precursor-Overlap Unique-Anchor Version

This document describes the workflow used to annotate the 2814 predicted soybean miRNA precursor-location records in `2814_precusor_miRNAs.txt`. This version uses precursor overlap on the Zhonghuang 13 genome as the genomic anchor criterion. It no longer uses the former same-chromosome, same-strand, precursor-midpoint distance <= 10 Mb rule.

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
../../01_reference_genome/GWHAAEV00000000.1.genome.fasta
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

Before final output, database-anchor consistency is checked again. If the final name of a `reported` or `known_member_variant` record corresponds to a soybean database precursor, the predicted precursor must overlap that database precursor on ZH13 by at least 1 bp. Otherwise, the record cannot remain assigned to that known member or its `.v` variant and is withdrawn to a same-family new member/new-member variant.

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
1. Position-variant propagation: unannotated records in the same position-variant cluster with highly overlapping mature intervals inherit the current anchor as .v variants.
2. Same-precursor arm search: non-overlapping mature regions on the same or highly overlapping predicted precursor are preferentially annotated as unreported mature arms.
```

Position-variant clusters are taken only from `2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv`; database physical coordinates are not used to judge mature-region overlap.

Tandem-order mapping is allowed only when each paired predicted precursor and database precursor in the local tandem block has a real ZH13 interval overlap. This step cannot force a distant non-overlapping predicted precursor onto a database precursor.

## 5. Status Definitions

`reported`: the mature sequence is a full-length 100% match to a soybean database mature miRNA, and the predicted precursor overlaps the corresponding database precursor on the same ZH13 chromosome and strand.

`known_member_variant`: the query mature sequence belongs to the same mature region or arm as a known database member and can be explained by an overlapping same-family database precursor.

`unreported_mature_arm`: another non-overlapping mature arm from the same or highly overlapping predicted precursor where a non-arm anchor already exists.

`known_family_new_member`: known-family evidence is present and the family exists in the current soybean reference database, but no suitable known-member anchor is available.

`new_family_member`: the family is absent from the current soybean reference database, or the record remains a study-specific no-hit candidate.

Position variants are written as `known_family_new_member_variant`, `new_family_member_variant`, or `unreported_mature_arm_variant`. The `.v` suffix denotes an isomiR or mature-position variant within the same mature region or arm and is numbered by descending `_R` read count in `Seq-ID`.

## 6. Naming Rules

```text
1. Use mature miRNA names rather than precursor names.
2. During the miRBase stage, avoid soybean miRBase member names already occupied.
3. During the pmiREN stage, avoid soybean pmiREN member names already occupied.
4. miRBase and pmiREN member suffixes are not forced to correspond one-to-one.
5. Family labels in output tables are normalized to case-insensitive mature-family form: Gma/gma is written as gma, and MIR/MIRN is written as miR/miRN. For example, Gma-MIR12410, gma-MIR12410, and Gma-miR12410 are all reported as gma-miR12410 in the Family column.
6. Within a position-variant cluster, the highest-read member keeps the base name; the others are named .v1, .v2, ...
7. If two non-overlapping mature arms occur on the same predicted precursor, they are treated as the same locus; the anchored arm keeps the base name and the other arm is named -5p or -3p.
```

SoymiR no-hit naming:

```text
1. Remaining no-hit records are grouped by CD-HIT cluster.
2. SoymiR family numbers are assigned by descending family size: gma_miRNC0001, gma_miRNC0002, ...
3. Independent genomic loci within the same CD-HIT family receive a/b/c suffixes.
4. Mature-position variants are named .v1, .v2, ... by descending read count.
5. Two no-hit loci in the same CD-HIT family are merged as two arms of one locus when they have the same chromosome and strand, precursor overlap >= 80%, and non-overlapping mature intervals.
6. If a no-hit locus belongs to a CD-HIT cluster containing one unambiguous miRBase/pmiREN family, it inherits that family. The status is `known_family_new_member` only when the family exists in the corresponding soybean database; otherwise it remains `new_family_member`.
```

## 7. Reported Status and Conservation

`Reported_Status` is `reported` only when `Status = reported`; all other records are `unreported`.

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
Status
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
Family count after case-insensitive mature-family normalization: 911
miRBase ZH13 hairpin anchors loaded: 682
pmiREN ZH13 hairpin anchors loaded: 1347
```

Source counts:

```text
miRbase: 1198
pmiren: 518
soymir: 1098
```

Status counts:

```text
reported: 514
known_member_variant: 549
known_family_new_member: 310
known_family_new_member_variant: 73
unreported_mature_arm: 209
unreported_mature_arm_variant: 66
new_family_member: 1015
new_family_member_variant: 78
```

Reported and conservation counts:

```text
reported: 514
unreported: 2300
conserved: 851
specific: 1963
```

## 10. Reproducibility

```bash
cd 2026-miRNA/combine_data/current_database
bash run_hairpin_ZH13_blast_unique_anchor.sh
python3 annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 update_simplified_result_and_supplementary_s2_s4.py
```
