#!/bin/bash
# config.sh - 项目路径配置
# 所有脚本 source 此文件获取公共路径

# 项目根目录。该配置可被 workflow/raw_data_example/、01_miRNA_identification/
# 或 03_miRNA_targets/ 下的示例脚本 source。
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "${SCRIPT_DIR}/../../02_miRNA_annotation" ]; then
    PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
elif [ -d "${SCRIPT_DIR}/../../../02_miRNA_annotation" ]; then
    PROJECT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
else
    echo "[错误] 无法自动定位 SoymiR-atlas 根目录，请检查脚本位置。" >&2
    exit 1
fi

# Conda 环境（所有工具统一装在此环境）
CONDA_ENV="mirna_env"

# 工具命令（依赖 conda 环境 PATH）
PREFETCH="prefetch"
FASTERQ="fasterq-dump"
FASTP="fastp"
SEQKIT="seqkit"
BOWTIE_BUILD="bowtie-build"

# 软件路径。miRDP2 和 CleaveLand4 通常不随仓库发布，需要用户按 README
# 下载或安装到 software/ 下，或在这里改为本机绝对路径。
MIRDEP2="${MIRDEP2:-${PROJECT_DIR}/software/miRDP2-v1.1.4/miRDP2-v1.1.4_pipeline.bash}"
CLEAVELAND4="${CLEAVELAND4:-${PROJECT_DIR}/software/CleaveLand4-4.5/CleaveLand4.pl}"
GSTAR_DIR="${GSTAR_DIR:-${PROJECT_DIR}/software/CleaveLand4-4.5/GSTAr_v1-0}"

# 参考数据
GENOME="${GENOME:-${PROJECT_DIR}/data/reference/genome/GWHAAEV00000000.1.genome.fasta}"
GENOME_INDEX="${GENOME_INDEX:-${PROJECT_DIR}/data/reference/genome/GWHAAEV00000000.1}"
CDS="${CDS:-${PROJECT_DIR}/data/reference/cds/GWHAAEV00000000.1.CDS.simple.fasta}"

# miRNA 输入文件（来自 filter_P_prediction）
MIRNA_QUERY="${MIRNA_QUERY:-${PROJECT_DIR}/01_miRNA_identification/SRR27718796/results/SRR27718796_miRNA_query.fa}"
