#!/bin/bash
# 07_deg_mode4.sh
# CleaveLand4 Mode 4：用已有 dd.txt + GSTAr.txt 分析降解组密度
# 用法：bash workflow/07_deg_mode4.sh SRR23932132 [SRR29504124 ...]

set -e
source "$(cd "$(dirname "$0")" && pwd)/config.sh"

[ $# -eq 0 ] && echo "用法: $0 <SRR> [SRR ...]" && exit 1

export PATH="${GSTAR_DIR}:${PATH}"

for SRR in "$@"; do
    SRR_DIR="${PROJECT_DIR}/03_miRNA_targets/${SRR}"
    OUTDIR="${SRR_DIR}/04_results"
    BOWTIE_DIR="${SRR_DIR}/03_bowtie"
    TPlot_DIR="${OUTDIR}/TPlots"
    mkdir -p "${TPlot_DIR}"

    DD_RAW="${OUTDIR}/${SRR}.clean.fasta_dd.txt"
    GSTAR_FILE="${BOWTIE_DIR}/${SRR}_GSTAr.txt"

    # 如果 04_results/ 没有 dd.txt，尝试 02_fastp/（旧版）
    [ -f "${DD_RAW}" ] || DD_RAW="${SRR_DIR}/02_fastp/${SRR}.clean.fasta_dd.txt"
    [ -f "${DD_RAW}" ] || DD_RAW="${BOWTIE_DIR}/${SRR}_dd.txt"
    [ -f "${DD_RAW}" ] || { echo "[错误] dd.txt 未找到"; exit 1; }
    [ -f "${GSTAR_FILE}" ] || { echo "[错误] GSTAr.txt 不存在: ${GSTAR_FILE}"; exit 1; }

    # 从 GSTAr.txt 读转录组路径，补到 dd.txt header 中
    DD_INPUT="${OUTDIR}/${SRR}_mode4_dd.txt"
    python3 - "${DD_RAW}" "${GSTAR_FILE}" "${DD_INPUT}" <<'PY'
import sys, os, re
dd_raw, gstar, dd_out = sys.argv[1], sys.argv[2], sys.argv[3]

tx_path = ""
with open(gstar) as f:
    for line in f:
        if line.startswith("# Transcripts:"):
            tx_path = line.split(":", 1)[1].strip()
            break

with open(dd_raw) as f:
    content = f.read()

if re.search(r'^# Transcriptome:', content, re.MULTILINE):
    content = re.sub(r'^# Transcriptome:.*', f'# Transcriptome:{tx_path}', content, flags=re.MULTILINE)
else:
    lines = content.split("\n")
    header_end = 0
    for i, line in enumerate(lines):
        if line.startswith("#"):
            header_end = i + 1
        else:
            break
    lines.insert(header_end, f"# Transcriptome:{tx_path}")
    content = "\n".join(lines)

with open(dd_out, "w") as f:
    f.write(content)
PY

    echo "=== Mode 4: ${SRR} ==="
    perl "${CLEAVELAND4}" \
      -d "${DD_INPUT}" \
      -g "${GSTAR_FILE}" \
      -p 0.05 \
      -c 1 \
      -o "${TPlot_DIR}" \
      > "${OUTDIR}/${SRR}_mode4.txt" 2> "${OUTDIR}/${SRR}_mode4.log"

    rm -f "${DD_INPUT}"
    echo "[完成] ${SRR}_mode4.txt"
done
echo "=== 全部完成 ==="
