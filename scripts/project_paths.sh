#!/usr/bin/env bash
# Shared paths for the public-release Bash scripts. Requires Bash 4+.

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
cd -- "$repo_root"
run_dir=$(realpath -m -- "${ATAC_RUN_DIR:-runs/local}")
case "$run_dir" in
    "$repo_root"/runs/*) ;;
    *) die 'ATAC_RUN_DIR must be a subfolder of this repository’s runs/ directory.' ;;
esac
mapfile -t samples < <(awk -F '\t' 'NR > 1 {sub(/\r$/, ""); print $1}' input/samples.tsv)
[[ ${#samples[@]} -eq 8 ]] || die 'Expected eight samples in input/samples.tsv.'
for sample in "${samples[@]}"; do
    [[ "$sample" =~ ^SRR[0-9]+$ ]] || die 'Invalid sample accession.'
done
[[ $(printf '%s\n' "${samples[@]}" | sort -u | wc -l) -eq 8 ]] || die 'Duplicate sample accessions.'
mkdir -p -- "$run_dir"
shopt -s nullglob
