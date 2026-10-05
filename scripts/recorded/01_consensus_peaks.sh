#!/usr/bin/env bash
# Recorded Week 2 command body, recovered from the saved July 2026 workflow.
# Institutional BASE path replaced with ATAC_RECORDED_BASE for public release.
# This retains the historical unchanged BED-to-SAF start; see ../README.md.
# This is a historical reference, separate from the corrected main scripts.
set -euo pipefail

: "${ATAC_RECORDED_BASE:?Set ATAC_RECORDED_BASE to your local teaching-data project directory}"
BASE="$ATAC_RECORDED_BASE"
PEAK_DIR="$BASE/P07_Regulatory_T_vs_Naive_Teffs"
OUT_DIR="$BASE/out_dir"

SAMPLES=(
    "SRR7650753"
    "SRR7650755"
    "SRR7650836"
    "SRR7650838"
    "SRR7650870"
    "SRR7650874"
    "SRR7650796"
    "SRR7650798"
)

[[ ! -e "$OUT_DIR/consensus_peaks.bed" && ! -e "$OUT_DIR/consensus_peaks.saf" ]] || { printf "%s\n" "Consensus outputs already exist." >&2; exit 1; }
mkdir -p "$OUT_DIR"

# 1. Extract, sort and merge overlapping regions
for sample in "${SAMPLES[@]}"; do
    awk '{print $1"\t"$2"\t"$3}' \
        "${PEAK_DIR}/${sample}_REP1.mLb.clN_peaks.broadPeak"
done | sort -k1,1 -k2,2n |
    bedtools merge -i - \
    > "${OUT_DIR}/consensus_peaks.bed"

# 2. Convert the BED file to SAF format for featureCounts
awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     {print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}' \
    "${OUT_DIR}/consensus_peaks.bed" \
    > "${OUT_DIR}/consensus_peaks.saf"

