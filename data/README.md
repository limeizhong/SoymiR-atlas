# data

Shared example inputs or external data links can be placed here.

Large primary data files should preferably be documented in `external_links/` rather than committed directly, unless they are required small reference files for reproducibility.

The raw-data example workflow expects reference files under:

```text
data/reference/genome/GWHAAEV00000000.1.genome.fasta
data/reference/cds/GWHAAEV00000000.1.CDS.simple.fasta
```

These large reference files are not committed by default. Users can either place
them at the paths above or override `GENOME`, `GENOME_INDEX`, and `CDS` in
`workflow/raw_data_example/config.sh`.
