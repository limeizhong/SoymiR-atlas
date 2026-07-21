#!/bin/bash
# 03_mirna_mirdeep2.sh
# 格式转换 + miRDeep-P2 预测
# 用法：bash workflow/03_mirna_mirdeep2.sh SRR27718796

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

SRR="$1"
[ -z "$SRR" ] && echo "用法: $0 <SRR号>" && exit 1

DATA_DIR="${PROJECT_DIR}/01_miRNA_identification/${SRR}"
FASTP_DIR="${DATA_DIR}/fastp"
RESULT_DIR="${DATA_DIR}/results"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "${RESULT_DIR}"

INPUT="${FASTP_DIR}/${SRR}.clean.fastq.gz"
FA_GZ="${RESULT_DIR}/${SRR}.fa.gz"
FILTER_FA="${RESULT_DIR}/${SRR}.filter.fa"
OUTDIR="${RESULT_DIR}/mirdp2"

[ -f "${INPUT}" ] || { echo "[错误] 输入不存在: ${INPUT}"; exit 1; }

echo "=== 1/3 fastq → fa.gz ==="
"${SEQKIT}" fq2fa "${INPUT}" -o "${FA_GZ}"

echo "=== 2/3 去重+过滤（19-24nt） ==="
python3 "${SCRIPT_DIR}/filter.fasta2input.format.py" "${FA_GZ}"

echo "=== 3/3 miRDeep-P2 预测 ==="

mkdir -p "${OUTDIR}"
bash "${MIRDEP2}" \
  -g "${GENOME}" \
  -x "${GENOME_INDEX}" \
  -f \
  -i "${FILTER_FA}" \
  --rpm 5 \
  -o "${OUTDIR}" \
  -p 1

echo "[完成] ${SRR} → ${OUTDIR}"
