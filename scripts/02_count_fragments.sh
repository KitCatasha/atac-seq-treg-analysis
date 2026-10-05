#!/usr/bin/env bash
# Public-release adaptation of the saved July 2026 Week 2 counting workflow,
# cross-checked with the SysBio Internship Schedule (SS2026, Week 2, pp. 5-6).
# The recovered historical command body is in recorded/02_count_fragments.sh.
# Usage: bash scripts/02_count_fragments.sh /path/to/cleaned_BAM_directory
set -euo pipefail
[[ $# -eq 1 ]] || { printf '%s\n' 'Usage: bash scripts/02_count_fragments.sh BAM_DIR' >&2; exit 1; }
bam_dir=$(realpath -e -- "$1")
source "$(dirname -- "${BASH_SOURCE[0]}")/project_paths.sh"
command -v featureCounts >/dev/null || die 'Install Subread/featureCounts before counting fragments.'
[[ -d "$bam_dir" ]] || die 'BAM_DIR must be a directory.'
[[ -s "$run_dir/input/consensus_peaks.saf" ]] || die 'Run 01_consensus_peaks.sh first, using the same ATAC_RUN_DIR.'
[[ ! -e "$run_dir/input/counts_clean.txt" && ! -e "$run_dir/qc/counts_matrix.txt" ]] || die 'Count outputs already exist; choose a new ATAC_RUN_DIR.'
stage=$(mktemp -d "$run_dir/.counting.XXXXXX")
trap 'rm -rf -- "$stage"' EXIT
bams=()
for sample in "${samples[@]}"; do
    files=("$bam_dir/$sample"*.bam)
    [[ ${#files[@]} -eq 1 && -s "${files[0]}" ]] || die "Need exactly one nonempty cleaned $sample*.bam file."
    # As in the course example, accept either standard BAM-index filename.
    [[ -s "${files[0]}.bai" || -s "${files[0]%.bam}.bai" ]] \
        || die "BAM index missing for $sample; run samtools index first."
    bams+=("${files[0]}")
    printf '%s\t%s\n' "$sample" "${files[0]}" >> "$stage/bam_order.tsv"
done
# Paired fragments, both ends mapped, chimeric pairs excluded, two threads.
# -B requires both ends to be mapped; it does not require both to overlap a peak.
featureCounts -F SAF -p --countReadPairs -B -C -T 2 \
    -a "$run_dir/input/consensus_peaks.saf" \
    -o "$stage/counts_matrix.txt" "${bams[@]}" \
    > "$stage/featureCounts.log" 2>&1
[[ -s "$stage/counts_matrix.txt.summary" ]] || die 'featureCounts assignment summary is missing.'
# Label columns in the same order as the actual BAM arguments, not saved matrix order.
awk -F '\t' 'BEGIN {OFS="\t"}
    FNR==NR {ids[++n]=$1; paths[n]=$2; next}
    /^#/ {next}
    !header {
        if (NF!=6+n || $1!="Geneid") exit 1;
        printf "PeakID";
        for (i=1;i<=n;i++) {
            if ($(i+6)!=paths[i]) exit 1;
            printf "\t%s", ids[i]
        }
        printf "\n"; header=1; next
    }
    {
        if (NF!=6+n) exit 1;
        printf "%s", $1;
        for (i=7;i<=NF;i++) {
            if ($i !~ /^[0-9]+$/) exit 1;
            printf "\t%s", $i
        }
        printf "\n"; rows++
    }
    END {if (!header || !rows) exit 1}' \
    "$stage/bam_order.tsv" "$stage/counts_matrix.txt" > "$stage/counts_clean.txt" \
    || die 'Unexpected featureCounts columns or noninteger counts.'
mkdir -p -- "$run_dir/qc"
cp -- "$stage/counts_clean.txt" "$run_dir/input/"
cp -- "$stage/counts_matrix.txt" "$stage/counts_matrix.txt.summary" \
    "$stage/featureCounts.log" "$stage/bam_order.tsv" "$run_dir/qc/"
printf 'Rebuilt count matrix: %s/input/counts_clean.txt\n' "$run_dir"
