# File format description

## Annotation short table

Main columns in `2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv`:

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

`miRNA_Locus` removes `.v` suffixes and `-5p/-3p` arm suffixes from `Annotation`.

`Reported_Status` is `reported` only for records with `Status=reported`; all other records are `unreported`.

`Conservation` is `conserved` when the query has at least one non-soybean mature-miRNA BLAST hit in miRBase or pmiREN; otherwise it is `specific`.
