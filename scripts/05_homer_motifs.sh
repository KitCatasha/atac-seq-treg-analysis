#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 0 ]] || { printf '%s\n' 'Usage: bash scripts/05_homer_motifs.sh' >&2; exit 1; }
source "$(dirname -- "${BASH_SOURCE[0]}")/project_paths.sh"
command -v findMotifsGenome.pl >/dev/null || die 'Install HOMER and its hg38 genome first.'
threads=${HOMER_THREADS:-2}
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die 'HOMER_THREADS must be a positive integer.'
sets=(treg naive treg_promoter naive_promoter treg_nonpromoter naive_nonpromoter)
backgrounds=(background background background_promoter background_promoter background_nonpromoter background_nonpromoter)
for i in "${!sets[@]}"; do
    target="$run_dir/homer_input/${sets[$i]}_DARs.bed"
    background="$run_dir/homer_input/${backgrounds[$i]}_peaks.bed"
    [[ -s "$target" && -s "$background" ]] || die "Missing or empty target/background BED for ${sets[$i]}."
    [[ ! -e "$run_dir/homer_results/${sets[$i]}_200" ]] || die "HOMER output already exists for ${sets[$i]}."
done
mkdir -p -- "$run_dir/homer_results"
for i in "${!sets[@]}"; do
    output="$run_dir/homer_results/${sets[$i]}_200"
    findMotifsGenome.pl "$run_dir/homer_input/${sets[$i]}_DARs.bed" hg38 "$output" \
        -size 200 -bg "$run_dir/homer_input/${backgrounds[$i]}_peaks.bed" -p "$threads" \
        > "$run_dir/homer_results/${sets[$i]}_200.log" 2>&1
    [[ -s "$output/knownResults.txt" ]] || die "Known-motif output missing for ${sets[$i]}; check its log."
done
printf 'HOMER outputs: %s/homer_results\n' "$run_dir"

