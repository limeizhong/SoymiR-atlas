# 03_miRNA_targets

This module archives target prediction results for the 239 newly added mature miRNAs.

## Workflow

1. CleaveLand4 was run on the server using `run_CleaveLand4_239_targets_fixed.sh`.
2. The fixed result directory was:
   `results2025_239_fixed`
3. Binding-site sequences and transcript/genome coordinates were extracted with:
   `extract_cleaveland_binding_sites.py`
4. Output tables were written to:
   `results/results2025_239_binding_site_sequences/`

The extraction keeps the original `results2025_239` file prefix for compatibility with downstream plotting scripts, but the contents are based on `results2025_239_fixed`.

## Main Outputs

- `results2025_239_binding_site_unique.tsv`: non-redundant miRNA-target-site records with library support.
- `results2025_239_binding_site_unique_minimal.tsv`: minimal three-column table.
- `results2025_239_binding_site_genome_mapping_candidates.tsv`: binding-site genome coordinate candidates converted from CDS coordinates.
- `results2025_239_binding_site_extraction_report.md`: extraction summary.

Current summary after merging the original 1,406-miRNA CleaveLand results and the 239 newly added miRNA results:

- old-numbered `.res.txt` files scanned from `results2025`: 13,634
- files skipped because the new ID was `NA`: 723
- old-numbered files with CleaveLand records after ID conversion: 12,911
- existing 239-miRNA binding-site records merged: 1,717
- total parsed binding-site result records after merging: 78,733
- unique miRNA-target-site records: 8,301
- unique miRNA-target pairs: 8,293
- unique target genes: 5,496
- genome-coordinate candidate records matching the CDS sequence: 78,733/78,733

## Optional raw-data example workflow

The migrated example scripts from the raw degradome workflow are stored under:

```text
../../workflow/raw_data_example/
```

They show how representative degradome libraries can be processed from raw
FASTQ through CleaveLand4:

```bash
bash ../../workflow/raw_data_example/05_deg_download_qc.sh SRR23932132 SRR29504124
bash ../../workflow/raw_data_example/06_deg_mode1.sh SRR23932132
bash ../../workflow/raw_data_example/06_deg_mode1.sh SRR29504124
bash ../../workflow/raw_data_example/07_deg_mode4.sh SRR23932132 SRR29504124
```

These scripts require local installations of CleaveLand4, GSTAr, Bowtie,
fastp, seqkit, and the soybean CDS reference files. They are provided as a
small raw-data example and are separate from the final merged target tables
listed above.

Example outputs from two degradome libraries are archived in:

```text
results/example_degradome/SRR23932132/04_results/
results/example_degradome/SRR29504124/04_results/
```
