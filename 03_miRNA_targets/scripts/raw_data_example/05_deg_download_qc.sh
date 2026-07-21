#!/bin/bash
# 05_deg_download_qc.sh
# 降解组 SRA 下载 + fastp 质控
# 需先 conda activate mirna_env，或通过 run_all.sh 调用
# 用法：bash workflow/05_deg_download_qc.sh SRR23932132 [SRR29504124 ...]

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

[ $# -eq 0 ] && echo "用法: $0 <SRR> [SRR ...]" && exit 1

for SRR in "$@"; do
    DATA_DIR="${PROJECT_DIR}/03_miRNA_targets/${SRR}"
    RAW_DIR="${DATA_DIR}/01_raw"
    FASTP_DIR="${DATA_DIR}/02_fastp"
    mkdir -p "${RAW_DIR}" "${FASTP_DIR}"

    echo "=== 1/2 下载 ${SRR} ==="
    "${PREFETCH}" "${SRR}" -O "${RAW_DIR}"
    SRA_FILE="${RAW_DIR}/${SRR}/${SRR}.sra"
    if [ -f "${SRA_FILE}" ]; then
        "${FASTERQ}" "${SRA_FILE}" -O "${RAW_DIR}"
        # fasterq-dump 输出可能是 SRR.fastq 或 SRR_1.fastq
        for f in "${RAW_DIR}/${SRR}.fastq" "${RAW_DIR}/${SRR}_1.fastq"; do
            [ -f "$f" ] && gzip -f "$f"
        done 2>/dev/null
        rm -rf "${RAW_DIR}/${SRR}"
    fi

    # 检测 fastq 文件名（_1.fastq.gz 或 .fastq.gz）
    RAW_GZ="${RAW_DIR}/${SRR}_1.fastq.gz"
    [ ! -f "${RAW_GZ}" ] && RAW_GZ="${RAW_DIR}/${SRR}.fastq.gz"
    [ ! -f "${RAW_GZ}" ] && { echo "[错误] 未找到原始 FASTQ"; exit 1; }

    echo "=== 2/2 fastp 质控 ==="
    "${FASTP}" \
      -i "${RAW_GZ}" \
      -o "${FASTP_DIR}/${SRR}.clean.fastq.gz" \
      -h "${FASTP_DIR}/${SRR}_fastp.html" \
      -j "${FASTP_DIR}/${SRR}_fastp.json" \
      -l 19 \
      --detect_adapter_for_pe \
      --thread 8

    echo "[完成] ${SRR}"
done
echo "=== 全部完成 ==="
