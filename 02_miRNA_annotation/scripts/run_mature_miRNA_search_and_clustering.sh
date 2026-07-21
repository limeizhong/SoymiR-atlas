#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="$MODULE_DIR/input"
INTERMEDIATE_DIR="$MODULE_DIR/results/intermediate"
mkdir -p "$INTERMEDIATE_DIR"

QUERY="$INPUT_DIR/1588_mature_miRNAs.fasta"
MIRBASE_DB_FASTA="$INPUT_DIR/miRbase-mature.fa"
PMIREN_DB_FASTA="$INPUT_DIR/pmiren_all_mature_clean.fa"

MIRBASE_DB_PREFIX="$INTERMEDIATE_DIR/miRbase-mature.blastdb"
PMIREN_DB_PREFIX="$INTERMEDIATE_DIR/pmiren_all_mature_clean.blastdb"

MIRBASE_BLAST_OUT="$INTERMEDIATE_DIR/1588_mature_miRNAs_vs_miRbase_mature_blastn_e1e-4.tsv"
PMIREN_BLAST_OUT="$INTERMEDIATE_DIR/1588_mature_miRNAs_vs_pmiren_mature_blastn_e1e-4.tsv"

CDHIT_OUT_PREFIX="$INTERMEDIATE_DIR/1588_mature_miRNAs_cdhit_est_c0.8"

BLAST_OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq stitle"

require_file() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "Missing required input file: $path" >&2
    exit 1
  fi
}

require_command() {
  local command_name="$1"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "Required command not found in PATH: $command_name" >&2
    exit 1
  fi
}

require_file "$QUERY"
require_file "$MIRBASE_DB_FASTA"
require_file "$PMIREN_DB_FASTA"

require_command makeblastdb
require_command blastn
require_command cd-hit-est

echo "[1/5] Building miRBase all-species mature BLAST database"
makeblastdb \
  -in "$MIRBASE_DB_FASTA" \
  -dbtype nucl \
  -parse_seqids \
  -out "$MIRBASE_DB_PREFIX"

echo "[2/5] Searching mature miRNAs against miRBase"
blastn \
  -task blastn-short \
  -query "$QUERY" \
  -db "$MIRBASE_DB_PREFIX" \
  -evalue 1e-4 \
  -outfmt "$BLAST_OUTFMT" \
  -out "$MIRBASE_BLAST_OUT"

echo "[3/5] Building pmiREN all-species mature BLAST database"
makeblastdb \
  -in "$PMIREN_DB_FASTA" \
  -dbtype nucl \
  -parse_seqids \
  -out "$PMIREN_DB_PREFIX"

echo "[4/5] Searching mature miRNAs against pmiREN"
blastn \
  -task blastn-short \
  -query "$QUERY" \
  -db "$PMIREN_DB_PREFIX" \
  -evalue 1e-4 \
  -outfmt "$BLAST_OUTFMT" \
  -out "$PMIREN_BLAST_OUT"

echo "[5/5] Clustering mature miRNAs with CD-HIT-EST at c=0.8"
cd-hit-est \
  -i "$QUERY" \
  -o "$CDHIT_OUT_PREFIX" \
  -c 0.8 \
  -n 4 \
  -d 0 \
  -T 0 \
  -M 0

echo "Done."
echo "Generated:"
echo "  $MIRBASE_BLAST_OUT"
echo "  $PMIREN_BLAST_OUT"
echo "  ${CDHIT_OUT_PREFIX}"
echo "  ${CDHIT_OUT_PREFIX}.clstr"
