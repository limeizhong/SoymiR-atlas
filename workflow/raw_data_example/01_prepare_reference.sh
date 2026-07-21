#!/bin/bash
# 01_prepare_reference.sh
# 构建参考基因组和 CDS 的 bowtie 索引
# 用法：bash workflow/01_prepare_reference.sh

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

echo "=== 1/2 构建基因组索引 ==="
"${BOWTIE_BUILD}" -f "${GENOME}" "${GENOME_INDEX}"

echo "=== 2/2 构建 CDS 索引 ==="
CDS_DIR="$(dirname "${CDS}")"
"${BOWTIE_BUILD}" -f "${CDS}" "${CDS_DIR}/GWHAAEV00000000.1.CDS.simple.fasta"

echo "[完成]"
echo "  基因组索引: $(dirname ${GENOME_INDEX})/"
echo "  CDS 索引: ${CDS_DIR}/"
