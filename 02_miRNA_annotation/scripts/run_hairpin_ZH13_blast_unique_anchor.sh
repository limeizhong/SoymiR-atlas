#!/usr/bin/env bash
set -euo pipefail

# Map miRBase and pmiREN soybean hairpin precursors to the ZH13.v2 genome,
# format genome hits, and resolve one preferred ZH13 interval per precursor.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_DIR="$(cd "${MODULE_DIR}/.." && pwd)"
INPUT_DIR="${MODULE_DIR}/input"
INTERMEDIATE_DIR="${MODULE_DIR}/results/intermediate"
GENOME="${PROJECT_DIR}/01_reference_genome/GWHAAEV00000000.1.genome.fasta"
DB_PREFIX="${INTERMEDIATE_DIR}/blastdb/ZH13v2_genome"
EVALUE="${EVALUE:-1e-10}"
THREADS="${THREADS:-8}"

MIRBASE_HAIRPIN="${INPUT_DIR}/gma-hairpin.fa"
PMIREN_HAIRPIN="${INPUT_DIR}/pmiren_gmax_hairpin.fa"

BLAST_OUTFMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sstrand stitle'

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $1" >&2
    exit 1
  fi
}

usage() {
  cat <<'USAGE'
Usage:
  bash run_hairpin_ZH13_blast_unique_anchor.sh

Environment variables:
  THREADS=8       number of BLAST threads
  EVALUE=1e-10    BLAST e-value cutoff

Inputs:
  input/gma-hairpin.fa
  input/pmiren_gmax_hairpin.fa
  ../01_reference_genome/GWHAAEV00000000.1.genome.fasta
  input/miRbase_gma_position.txt
  input/pmiren_gma_position.txt

Outputs:
  results/intermediate/gma-hairpin_vs_ZH13v2_genome_blastn_e1e-10.tsv
  results/intermediate/gma-hairpin_ZH13v2_genome_positions_all_hits.tsv
  results/intermediate/gma-hairpin_ZH13v2_genome_positions_best_hit.tsv
  results/intermediate/gma-hairpin_ZH13v2_genome_positions_full_length_pident100.tsv
  results/intermediate/gma-hairpin_ZH13v2_genome_no_blast_hit.txt
  results/intermediate/gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
  results/intermediate/gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval_changed.tsv
  results/intermediate/gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval_unresolved_overlap80.tsv

  results/intermediate/pmiren_gmax_hairpin_vs_ZH13v2_genome_blastn_e1e-10.tsv
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_all_hits.tsv
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit.tsv
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_full_length_pident100.tsv
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_no_blast_hit.txt
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval_changed.tsv
  results/intermediate/pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval_unresolved_overlap80.tsv
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

require_command python3

mkdir -p "${INTERMEDIATE_DIR}/blastdb"

echo "Using pmiREN soybean hairpin file without overwriting: ${PMIREN_HAIRPIN}"

if [[ ! -s "${GENOME}" ]]; then
  echo "ERROR: genome fasta not found: ${GENOME}" >&2
  exit 1
fi

if [[ ! -s "${DB_PREFIX}.nsq" ]]; then
  require_command makeblastdb
  makeblastdb -in "${GENOME}" -dbtype nucl -parse_seqids -out "${DB_PREFIX}"
fi

run_one() {
  local query_fasta="$1"
  local prefix="$2"
  local raw="${INTERMEDIATE_DIR}/${prefix}_vs_ZH13v2_genome_blastn_e1e-10.tsv"
  local legacy_raw="${raw}"
  local all_hits="${INTERMEDIATE_DIR}/${prefix}_ZH13v2_genome_positions_all_hits.tsv"
  local best="${INTERMEDIATE_DIR}/${prefix}_ZH13v2_genome_positions_best_hit.tsv"
  local full="${INTERMEDIATE_DIR}/${prefix}_ZH13v2_genome_positions_full_length_pident100.tsv"
  local nohit="${INTERMEDIATE_DIR}/${prefix}_ZH13v2_genome_no_blast_hit.txt"

  if [[ ! -s "${query_fasta}" ]]; then
    echo "ERROR: query fasta not found: ${query_fasta}" >&2
    exit 1
  fi

  if command -v blastn >/dev/null 2>&1; then
    blastn \
      -task blastn \
      -query "${query_fasta}" \
      -db "${DB_PREFIX}" \
      -evalue "${EVALUE}" \
      -num_threads "${THREADS}" \
      -outfmt "${BLAST_OUTFMT}" \
      -out "${raw}"
  elif command -v conda >/dev/null 2>&1 && conda run -n bioinfo blastn -version >/dev/null 2>&1; then
    conda run -n bioinfo blastn \
      -task blastn \
      -query "${query_fasta}" \
      -db "${DB_PREFIX}" \
      -evalue "${EVALUE}" \
      -num_threads "${THREADS}" \
      -outfmt "${BLAST_OUTFMT}" \
      -out "${raw}"
  elif [[ -s "${raw}" ]]; then
    echo "WARNING: blastn not found; reusing existing ${raw##*/}" >&2
  elif [[ -s "${legacy_raw}" ]]; then
    echo "WARNING: blastn not found; filtering existing ${legacy_raw##*/} to EVALUE <= ${EVALUE}" >&2
    python3 "${SCRIPT_DIR}/filter_blast_by_evalue.py" \
      --input "${legacy_raw}" \
      --output "${raw}" \
      --max-evalue "${EVALUE}"
  else
    echo "ERROR: blastn not found and fallback raw file is missing: ${legacy_raw}" >&2
    exit 1
  fi

  python3 "${SCRIPT_DIR}/format_hairpin_ZH13_blast_hits.py" \
    --query-fasta "${query_fasta}" \
    --blast "${raw}" \
    --all-hits "${all_hits}" \
    --best-hit "${best}" \
    --full-length-pident100 "${full}" \
    --no-hit "${nohit}"
}

run_one "${MIRBASE_HAIRPIN}" "gma-hairpin"
run_one "${PMIREN_HAIRPIN}" "pmiren_gmax_hairpin"

python3 "${SCRIPT_DIR}/resolve_hairpin_best_hits_unique_intervals.py"
