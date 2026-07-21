#!/bin/bash
# run_all.sh
# 全自动运行：miRNA 鉴定 + 降解组靶标预测
# 用法：bash workflow/run_all.sh

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
else
    echo "[警告] conda 不在 PATH 中；继续运行并假定所需软件已在 PATH 中。"
fi

echo "=========================================="
echo " Soybean miRNA 鉴定 + 降解组靶标预测"
echo "=========================================="

# === 1. 参考索引 ===
echo ""
echo ">>> [1/6] 构建参考索引"
bash "$(dirname "$0")/01_prepare_reference.sh"

# === 2. miRNA-seq ===
echo ""
echo ">>> [2/6] miRNA-seq 下载 + 质控"
bash "$(dirname "$0")/02_mirna_download_qc.sh" SRR27718796

echo ""
echo ">>> [3/6] miRDeep-P2 预测"
bash "$(dirname "$0")/03_mirna_mirdeep2.sh" SRR27718796

echo ""
echo ">>> [4/6] 提取成熟体"
bash "$(dirname "$0")/04_mirna_mature.sh" SRR27718796

# === 3. 降解组 ===
echo ""
echo ">>> [5/6] 降解组下载 + 质控"
bash "$(dirname "$0")/05_deg_download_qc.sh" SRR23932132 SRR29504124

echo ""
echo ">>> [6/6] CleaveLand4 靶标预测（串行）"
bash "$(dirname "$0")/06_deg_mode1.sh" SRR23932132
bash "$(dirname "$0")/06_deg_mode1.sh" SRR29504124
bash "$(dirname "$0")/07_deg_mode4.sh" SRR23932132
bash "$(dirname "$0")/07_deg_mode4.sh" SRR29504124

echo ""
echo "=========================================="
echo " 全部完成！"
echo "=========================================="
echo " miRNA: 01_miRNA_identification/SRR27718796/results/SRR27718796_miRNA_query.fa"
echo " Mode1: 03_miRNA_targets/SRR*/04_results/*_mode1.txt"
echo " Mode4: 03_miRNA_targets/SRR*/04_results/*_mode4.txt"
echo " TPlot: 03_miRNA_targets/SRR*/04_results/TPlots/"
