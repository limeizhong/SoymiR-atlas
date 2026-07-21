#!/bin/bash
# 06_deg_mode1.sh
# CleaveLand4 Mode 1：全自动降解组靶基因预测
# 用法：bash workflow/06_deg_mode1.sh SRR23932132 [SRR29504124 ...]

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

[ $# -eq 0 ] && echo "用法: $0 <SRR> [SRR ...]" && exit 1

export PATH="${GSTAR_DIR}:${PATH}"

# R 脚本软链接（CleaveLand4 从 CWD 查找）
ln -sf "${PROJECT_DIR}/software/CleaveLand4-4.5/CleaveLand4_Tplotter.R" \
       "${PROJECT_DIR}/03_miRNA_targets/" 2>/dev/null || true

for SRR in "$@"; do
    SRR_DIR="${PROJECT_DIR}/03_miRNA_targets/${SRR}"
    FASTP_DIR="${SRR_DIR}/02_fastp"
    CLEAN_FASTA="${FASTP_DIR}/${SRR}.clean.fasta"
    CLEAN_FASTQ="${FASTP_DIR}/${SRR}.clean.fastq"
    CLEAN_GZ="${FASTP_DIR}/${SRR}.clean.fastq.gz"

    # 转换为 FASTA（CleaveLand4 需要 FASTA 格式）
    if [ ! -f "${CLEAN_FASTA}" ]; then
        if [ -f "${CLEAN_FASTQ}" ]; then
            "${SEQKIT}" fq2fa "${CLEAN_FASTQ}" -o "${CLEAN_FASTA}"
        elif [ -f "${CLEAN_GZ}" ]; then
            zcat "${CLEAN_GZ}" | "${SEQKIT}" fq2fa -o "${CLEAN_FASTA}"
        else
            echo "[错误] 输入不存在: ${CLEAN_FASTQ} / ${CLEAN_GZ}"; exit 1
        fi
    fi

    OUTDIR="${SRR_DIR}/04_results"
    mkdir -p "${OUTDIR}"
    TPlot_DIR="${OUTDIR}/TPlots"
    RESULT_FILE="${OUTDIR}/${SRR}_mode1.txt"

    echo "=== Mode 1: ${SRR} ==="
    perl "${CLEAVELAND4}" \
      -e "${CLEAN_FASTA}" \
      -u "${MIRNA_QUERY}" \
      -n "${CDS}" \
      -p 0.05 \
      -c 1 \
      -o "${TPlot_DIR}" \
      > "${RESULT_FILE}" 2>&1

    echo "[完成] ${SRR} → ${RESULT_FILE}"

    # 清理：挪 dd.txt 和 BAM 到 04_results/
    DD_FILE="${FASTP_DIR}/${SRR}.clean.fasta_dd.txt"
    if [ -f "${DD_FILE}" ]; then
        mv "${DD_FILE}" "${OUTDIR}/"
        echo "[info] dd.txt → ${OUTDIR}/"
    fi
    BAM_FILE="$(dirname "${CDS}")/$(basename "${CDS}")_sorted.bam"
    if [ -f "${BAM_FILE}" ]; then
        mv "${BAM_FILE}" "${OUTDIR}/${SRR}_sorted.bam"
        echo "[info] BAM → ${OUTDIR}/${SRR}_sorted.bam"
    fi
done
echo "=== 全部完成 ==="
