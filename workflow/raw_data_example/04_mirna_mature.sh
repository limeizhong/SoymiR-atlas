#!/bin/bash
# 04_mirna_mature.sh
# 从 miRDP2 结果提取成熟体并去冗余（生成 miRNA_query.fa）
# 用法：bash workflow/04_mirna_mature.sh SRR27718796

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

SRR="$1"
[ -z "$SRR" ] && echo "用法: $0 <SRR号>" && exit 1

RESULT_DIR="${PROJECT_DIR}/01_miRNA_identification/${SRR}/results"
MIRDP_DIR="${RESULT_DIR}/mirdp2"
PREDICTION="${MIRDP_DIR}/${SRR}.filter/${SRR}.filter_filter_P_prediction"
OUTPUT="${RESULT_DIR}/${SRR}_miRNA_query.fa"

[ -f "${PREDICTION}" ] || { echo "[错误] 预测结果不存在: ${PREDICTION}"; exit 1; }

echo "=== 提取成熟体序列 ==="
python3 - "${PREDICTION}" "${OUTPUT}" <<'PYEOF'
import sys
from collections import OrderedDict

infile = sys.argv[1]
outfile = sys.argv[2]

seen = OrderedDict()
with open(infile) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 7:
            qid = parts[2]    # read6129_x151
            seq = parts[6]    # mature sequence
            if seq not in seen:
                seen[seq] = qid

with open(outfile, 'w') as f:
    for seq, qid in seen.items():
        f.write(f'>{qid}\n{seq}\n')

print(f'[完成] {len(seen)} 条 unique 序列 → {outfile}')
PYEOF
