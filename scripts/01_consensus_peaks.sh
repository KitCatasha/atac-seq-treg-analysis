#!/usr/bin/env bash
# Public-release adaptation of the saved July 2026 Week 2 consensus workflow,
# cross-checked with the SysBio Internship Schedule (SS2026, Week 2, pp. 5-6).
# The recovered historical command body is in recorded/01_consensus_peaks.sh.
# Usage: bash scripts/01_consensus_peaks.sh /path/to/broadPeak_directory
set -euo pipefail
export LC_ALL=C
[[ $# -eq 1 ]] || { printf '%s\n' 'Usage: bash scripts/01_consensus_peaks.sh PEAK_DIR' >&2; exit 1; }
peak_dir=$(realpath -e -- "$1")
source "$(dirname -- "${BASH_SOURCE[0]}")/project_paths.sh"
command -v bedtools >/dev/null || die 'Install bedtools before rebuilding consensus peaks.'
[[ -d "$peak_dir" ]] || die 'PEAK_DIR must be a directory.'
for name in consensus_peaks.bed consensus_peaks.saf samples.tsv; do
    [[ ! -e "$run_dir/input/$name" ]] || die "Output already exists: $run_dir/input/$name"
done
stage=$(mktemp -d "$run_dir/.consensus.XXXXXX")
trap 'rm -rf -- "$stage"' EXIT
: > "$stage/combined.bed"
for sample in "${samples[@]}"; do
    files=("$peak_dir/$sample"*.broadPeak)
    [[ ${#files[@]} -eq 1 && -s "${files[0]}" ]] || die "Need exactly one nonempty $sample*.broadPeak file."
    awk 'BEGIN {OFS="\t"}
         NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 {exit 1}
         {print $1, $2, $3}' "${files[0]}" >> "$stage/combined.bed" || die "Invalid peak coordinates for $sample."
done
# Keep the supplied chromosome names; BED coordinates stay zero-based.
sort -k1,1 -k2,2n -k3,3n "$stage/combined.bed" > "$stage/sorted.bed"
bedtools merge -i "$stage/sorted.bed" > "$stage/merged.bed"
[[ -s "$stage/merged.bed" ]] || die 'No consensus intervals were produced.'
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "peak_" NR, 0, "."}' \
    "$stage/merged.bed" > "$stage/consensus_peaks.bed"
# SAF uses one-based, inclusive coordinates: BED start + 1, unchanged end.
# This corrects the schedule's unchanged BED start. The saved original counts
# used the unchanged start; rebuilding therefore need not match them exactly.
awk 'BEGIN {OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"}
     {print $4, $1, $2+1, $3, "."}' \
    "$stage/consensus_peaks.bed" > "$stage/consensus_peaks.saf"
mkdir -p -- "$run_dir/input"
cp -- "$stage/consensus_peaks.bed" "$stage/consensus_peaks.saf" input/samples.tsv "$run_dir/input/"
printf 'Rebuilt consensus inputs: %s/input\n' "$run_dir"
