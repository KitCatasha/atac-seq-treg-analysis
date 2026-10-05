#!/usr/bin/env bash
# Recorded Week 2 command body, recovered from the saved July 2026 workflow.
# Institutional BASE path replaced with ATAC_RECORDED_BASE for public release.
# This retains the historical unchanged BED-to-SAF start; see ../README.md.
# This is a historical reference, separate from the corrected main scripts.
set -euo pipefail

: "${ATAC_RECORDED_BASE:?Set ATAC_RECORDED_BASE to your local teaching-data project directory}"
BASE="$ATAC_RECORDED_BASE"
BAM_DIR="$BASE/P07_Regulatory_T_vs_Naive_Teffs"
OUT_DIR="$BASE/out_dir"
SAF="$OUT_DIR/consensus_peaks.saf"

BAMS=(
    "$BAM_DIR/SRR7650753_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650755_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650836_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650838_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650870_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650874_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650796_REP1.mLb.clN.sorted.bam"
    "$BAM_DIR/SRR7650798_REP1.mLb.clN.sorted.bam"
)

[[ ! -e "$OUT_DIR/counts_matrix.txt" && ! -e "$OUT_DIR/counts_clean.txt" ]] || { printf "%s\n" "Count outputs already exist." >&2; exit 1; }
mkdir -p "$OUT_DIR"

# 1. Check the SAF file
if [[ ! -f "$SAF" ]]; then
    echo "ERROR: SAF file not found: $SAF"
    exit 1
fi

# 2. Check the BAM files and their indices
for bam in "${BAMS[@]}"; do
    if [[ ! -f "$bam" ]]; then
        echo "ERROR: BAM file not found: $bam"
        exit 1
    fi

    if [[ ! -f "${bam}.bai" ]] && [[ ! -f "${bam%.bam}.bai" ]]; then
        echo "ERROR: BAM index not found for: $bam"
        echo "Run: samtools index \"$bam\""
        exit 1
    fi
done

# 3. Count paired-end fragments overlapping the consensus peaks
featureCounts \
    -F SAF \
    -a "$SAF" \
    -o "${OUT_DIR}/counts_matrix.txt" \
    -p \
    --countReadPairs \
    -B \
    -C \
    -T 2 \
    "${BAMS[@]}"

# 4. Retain PeakID and all sample-count columns
#    Extract the SRR accessions automatically from the BAM filenames
awk '
BEGIN {
    OFS = "\t"
}

NR == 1 && /^#/ {
    next
}

NR == 2 {
    printf "PeakID"

    for (i = 7; i <= NF; i++) {
        sample = $i
        sub(/^.*\//, "", sample)
        sub(/_REP1\.mLb\.clN\.sorted\.bam$/, "", sample)
        printf OFS sample
    }

    print ""
    next
}

{
    printf "%s", $1

    for (i = 7; i <= NF; i++) {
        printf OFS "%s", $i
    }

    print ""
}
' "${OUT_DIR}/counts_matrix.txt" \
> "${OUT_DIR}/counts_clean.txt"

