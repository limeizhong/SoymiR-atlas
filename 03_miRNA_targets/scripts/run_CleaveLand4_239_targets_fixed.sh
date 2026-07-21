#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$BASE_DIR"

CONDA="/home/zhonglm092001/zhonglm092001/miniconda3/bin/conda"
CLEAVELAND_ENV="cleaveland_fixed"
DEGRADOME_DIR_ORIG="fa_dd"
DEGRADOME_DIR="fa_dd_abs_transcriptome"
GSTAR_DIR="fas_file_239_not_in_1406_tabular"
RESULTS_DIR="results2025_239_fixed"
PLOTS_DIR="plots2025_239_fixed"
CLEAVELAND="../CleaveLand4.pl"
TRANSCRIPTS="ZH13.v2.CDS.fasta"
TRANSCRIPTS_ABS="$BASE_DIR/$TRANSCRIPTS"
COMMANDS="run_239_fixed.all.sh"
SUMMARY="run_CleaveLand4_239_fixed.summary.txt"
JOBS="${JOBS:-1}"
MODE="${1:-run}"

if [[ "$MODE" != "run" && "$MODE" != "commands-only" ]]; then
  echo "Usage: $0 [run|commands-only]" >&2
  exit 1
fi

for required in "$CLEAVELAND" "$TRANSCRIPTS"; do
  if [[ ! -s "$required" ]]; then
    echo "ERROR: required file not found or empty: $required" >&2
    exit 1
  fi
done

if [[ ! -x "$CONDA" ]]; then
  echo "ERROR: conda not found or not executable: $CONDA" >&2
  exit 1
fi

"$CONDA" run -n "$CLEAVELAND_ENV" perl -e 'print "ok\n"' >/dev/null 2>&1

mkdir -p "$DEGRADOME_DIR" "$RESULTS_DIR" "$PLOTS_DIR"

while IFS= read -r dd; do
  out="$DEGRADOME_DIR/$(basename "$dd")"
  if [[ -s "$out" ]] && grep -q "^# Transcriptome:$TRANSCRIPTS_ABS$" "$out"; then
    continue
  fi
  awk -v transcriptome="$TRANSCRIPTS_ABS" '
    /^# Transcriptome:/ {
      print "# Transcriptome:" transcriptome
      next
    }
    { print }
  ' "$dd" > "$out"
done < <(find "$DEGRADOME_DIR_ORIG" -maxdepth 1 -type f -name "*.fa_dd.txt" | sort)

degradome_count=$(find "$DEGRADOME_DIR" -maxdepth 1 -type f -name "*.fa_dd.txt" | wc -l | awk '{print $1}')
gstar_count=$(find "$GSTAR_DIR" -maxdepth 1 -type f -name "*.fa_GSTAr.txt" ! -name "*.tmp.*" | wc -l | awk '{print $1}')
tabular_count=$(grep -l "^# Output Format: Tabular" "$GSTAR_DIR"/*.fa_GSTAr.txt | wc -l | awk '{print $1}')

if [[ "$degradome_count" -eq 0 ]]; then
  echo "ERROR: no fixed degradome files found in $DEGRADOME_DIR" >&2
  exit 1
fi

if [[ "$gstar_count" -ne 239 || "$tabular_count" -ne 239 ]]; then
  echo "ERROR: expected 239 tabular GSTAr files; found result=$gstar_count tabular=$tabular_count" >&2
  exit 1
fi

{
  echo "#!/usr/bin/env bash"
  echo "set -euo pipefail"
  echo "export LC_ALL=C"
  echo 'cd "$(cd "$(dirname "$0")" && pwd)"'

  while IFS= read -r dd; do
    sample=$(basename "$dd" .fa_dd.txt)
    mkdir -p "$PLOTS_DIR/$sample"
    while IFS= read -r gstar; do
      mirna=$(basename "$gstar" .fa_GSTAr.txt)
      out="$RESULTS_DIR/${sample}.${mirna}.res.txt"
      plot_dir="$PLOTS_DIR/$sample"
      log="${out}.log"
      printf 'mkdir -p %q && tmp=%q.tmp.$$ && log=%q && if [[ -s %q ]]; then exit 0; fi && if %q run -n %q perl %q -d %q -g %q -p 0.05 -c 1 -o %q > "$tmp" 2> "$log"; then mv "$tmp" %q; else rm -f "$tmp"; exit 1; fi\n' \
        "$plot_dir" "$out" "$log" "$out" "$CONDA" "$CLEAVELAND_ENV" "$CLEAVELAND" "$dd" "$gstar" "$plot_dir" "$out"
    done < <(find "$GSTAR_DIR" -maxdepth 1 -type f -name "*.fa_GSTAr.txt" ! -name "*.tmp.*" | sort)
  done < <(find "$DEGRADOME_DIR" -maxdepth 1 -type f -name "*.fa_dd.txt" | sort)
} > "$COMMANDS"
chmod +x "$COMMANDS"

command_count=$(grep -c "CleaveLand4.pl" "$COMMANDS")

{
  echo "Original_degradome_dir: $DEGRADOME_DIR_ORIG"
  echo "Fixed_degradome_dir: $DEGRADOME_DIR"
  echo "Transcript_FASTA: $TRANSCRIPTS_ABS"
  echo "Degradome_count: $degradome_count"
  echo "GSTAr_dir: $GSTAR_DIR"
  echo "GSTAr_count: $gstar_count"
  echo "GSTAr_tabular_count: $tabular_count"
  echo "Expected_CleaveLand4_jobs: $((degradome_count * gstar_count))"
  echo "Command_count: $command_count"
  echo "Results_dir: $RESULTS_DIR"
  echo "Plots_dir: $PLOTS_DIR"
  echo "CleaveLand4_env: $CLEAVELAND_ENV"
  echo "Command_list: $COMMANDS"
  echo "Parallel_jobs: $JOBS"
  echo "Mode: $MODE"
} > "$SUMMARY"

cat "$SUMMARY"

if [[ "$MODE" == "commands-only" ]]; then
  exit 0
fi

grep "CleaveLand4.pl" "$COMMANDS" |
  xargs -I {} -P "$JOBS" bash -c '{}'
