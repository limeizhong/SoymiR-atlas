# Reproducibility

The completed `02_miRNA_annotation` module can be reproduced with:

```bash
bash 02_miRNA_annotation/scripts/run_mature_miRNA_search_and_clustering.sh
python3 02_miRNA_annotation/scripts/build_position_variant_clusters.py
GENOME=/path/to/GWHAAEV00000000.1.genome.fasta bash 02_miRNA_annotation/scripts/run_hairpin_ZH13_blast_unique_anchor.sh
python3 02_miRNA_annotation/scripts/annotate_2814_precusor_miRNAs_precursor_overlap_unique_anchor_workflow.py
python3 02_miRNA_annotation/scripts/update_simplified_result_and_supplementary_s2_s4.py
```

The BLAST and CD-HIT steps require BLAST+ and CD-HIT-EST. `run_hairpin_ZH13_blast_unique_anchor.sh` uses `blastn` from `PATH`, or `conda run -n bioinfo blastn` when available. The Excel table builders require Python with `openpyxl`.

The completed `04_miRNA_expression` module can be reproduced with:

```bash
python3 04_miRNA_expression/scripts/build_s3_tissue_expression_profiles.py
```

This step requires Python with `openpyxl`.
