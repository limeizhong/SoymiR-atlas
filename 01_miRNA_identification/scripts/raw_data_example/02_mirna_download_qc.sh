#!/bin/bash
# 02_mirna_download_qc.sh
# miRNA-seq 下载 + fastp 质控
# 需先 conda activate mirna_env，或通过 run_all.sh 调用
# 用法：bash workflow/02_mirna_download_qc.sh SRR27718796

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

SRR="$1"
[ -z "$SRR" ] && echo "用法: $0 <SRR号>" && exit 1

DATA_DIR="${PROJECT_DIR}/01_miRNA_identification/${SRR}"
RAW_DIR="${DATA_DIR}/raw"
FASTP_DIR="${DATA_DIR}/fastp"
mkdir -p "${RAW_DIR}" "${FASTP_DIR}"

echo "=== 1/2 下载 ${SRR} ==="
"${PREFETCH}" "${SRR}" -O "${RAW_DIR}"
SRA_FILE="${RAW_DIR}/${SRR}/${SRR}.sra"
if [ -f "${SRA_FILE}" ]; then
    "${FASTERQ}" "${SRA_FILE}" -O "${RAW_DIR}"
    # PE 输出 SRR_1.fastq + SRR_2.fastq，SE 输出 SRR.fastq
    for f in "${RAW_DIR}/${SRR}_1.fastq" "${RAW_DIR}/${SRR}_2.fastq" "${RAW_DIR}/${SRR}.fastq"; do
        [ -f "$f" ] && gzip -f "$f"
    done 2>/dev/null
    rm -rf "${RAW_DIR}/${SRR}"
fi

echo "=== 2/2 fastp 质控 ==="
# 检测 PE 或 SE 的 fastq 文件名
RAW_GZ="${RAW_DIR}/${SRR}_1.fastq.gz"
[ ! -f "${RAW_GZ}" ] && RAW_GZ="${RAW_DIR}/${SRR}.fastq.gz"
"${FASTP}" \
  -i "${RAW_GZ}" \
  -o "${FASTP_DIR}/${SRR}.clean.fastq.gz" \
  -h "${FASTP_DIR}/${SRR}.fastp.html" \
  -j "${FASTP_DIR}/${SRR}.fastp.json" \
  -l 19 \
  --length_limit 24 \
  --detect_adapter_for_pe \
  --thread 8

echo "[完成] ${SRR} 质控完成"
