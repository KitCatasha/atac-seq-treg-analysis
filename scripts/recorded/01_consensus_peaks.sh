#!/usr/bin/env bash
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

for sample in "${SAMPLES[@]}"; do
    awk '{print $1"\t"$2"\t"$3}' \
        "${PEAK_DIR}/${sample}_REP1.mLb.clN_peaks.broadPeak"
done | sort -k1,1 -k2,2n |
    bedtools merge -i - \
    > "${OUT_DIR}/consensus_peaks.bed"

awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     {print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}' \
    "${OUT_DIR}/consensus_peaks.bed" \
    > "${OUT_DIR}/consensus_peaks.saf"

