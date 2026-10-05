#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Adapted from the saved Week 1 workflow of 26 July 2026.
# Processing commands below retain the recorded parameters; paths and guards
# are adapted for this public repository. This is a chr22 teaching exercise,
# separate from the eight-sample whole-genome comparison.
# Blacklist exclusion is absent from the recovered Week 1 command document.
# It is not reconstructed here; this script is not the complete execution log.
# Usage: bash scripts/00_week1_chr22_training.sh R1.fastq.gz R2.fastq.gz BT2_CHR22_PREFIX
set -euo pipefail
[[ $# -eq 3 ]] || {
    printf '%s\n' 'Usage: bash scripts/00_week1_chr22_training.sh R1.fastq.gz R2.fastq.gz BT2_CHR22_PREFIX' >&2
    exit 1
}
R1=$(realpath -e -- "$1")
R2=$(realpath -e -- "$2")
BT2_INDEX=$(realpath -m -- "$3")
export ATAC_RUN_DIR="${ATAC_RUN_DIR:-runs/week1}"
source "$(dirname -- "${BASH_SOURCE[0]}")/project_paths.sh"
for command_name in fastqc cutadapt bowtie2 samtools picard macs2 multiqc gzip; do
    command -v "$command_name" >/dev/null || die "Required program missing: $command_name"
done
[[ -s "$R1" && -s "$R2" ]] || die 'Both raw FASTQ files must be nonempty.'
index_files=("$BT2_INDEX"*.bt2 "$BT2_INDEX"*.bt2l)
[[ ${#index_files[@]} -eq 6 ]] || die 'Provide one complete six-file Bowtie2 chromosome-22 index.'
for index_file in "${index_files[@]}"; do
    [[ -s "$index_file" ]] || die "Empty index file: $index_file"
done
PROJECT="$run_dir/week1_chr22"
[[ ! -e "$PROJECT" ]] || die 'Week 1 output directory already exists; choose a new ATAC_RUN_DIR.'
SAMPLE=SRR7650763
TRIMMED="$PROJECT/trimmed"
ALIGN="$PROJECT/alignment"
FILTERED="$PROJECT/filtered_bam"
PEAKS="$PROJECT/peaks"
QC="$PROJECT/qc"
MULTIQC_INPUTS="$PROJECT/multiqc_inputs"
MULTIQC_OUTPUT="$PROJECT/multiqc_output"
TRIM_R1="$TRIMMED/${SAMPLE}_1_chr22.trimmed.fastq.gz"
TRIM_R2="$TRIMMED/${SAMPLE}_2_chr22.trimmed.fastq.gz"
SAM="$ALIGN/alignment.sam"
DEDUP_BAM="$FILTERED/${SAMPLE}.filtered.sorted.dedup.bam"
mkdir -p "$TRIMMED" "$ALIGN" "$FILTERED" "$PEAKS"
mkdir -p "$QC/fastqc_raw" "$QC/fastqc_trimmed"
mkdir -p "$MULTIQC_INPUTS" "$MULTIQC_OUTPUT"
gzip -t "$R1" "$R2"

# Recorded Week 1 workflow, section 5.
fastqc \
  -t 2 \
  -o "$QC/fastqc_raw" \
  "$R1" "$R2"

# Recorded Week 1 workflow, section 6.
cutadapt \
  -j 2 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -O 5 \
  -m 20 \
  -o "$TRIM_R1" \
  -p "$TRIM_R2" \
  "$R1" "$R2" \
  > "$TRIMMED/${SAMPLE}.cutadapt.log" 2>&1

# Recorded Week 1 workflow, section 7.
fastqc \
  -t 2 \
  -o "$QC/fastqc_trimmed" \
  "$TRIM_R1" "$TRIM_R2"

# Recorded Week 1 workflow, section 8.
{ time bowtie2 \
  -p 2 \
  -X 2000 \
  -x "$BT2_INDEX" \
  -1 "$TRIM_R1" \
  -2 "$TRIM_R2" \
  -S "$SAM"; } \
  2> "$ALIGN/bowtie2.log"

# Recorded Week 1 workflow, section 9.
samtools view \
  -@ 2 \
  -b \
  -q 30 \
  -F 1804 \
  -f 2 \
  -o "$FILTERED/${SAMPLE}.filtered.bam" \
  "$SAM"

# Recorded Week 1 workflow, section 10.
samtools sort \
  -@ 2 \
  -o "$FILTERED/${SAMPLE}.filtered.sorted.bam" \
  "$FILTERED/${SAMPLE}.filtered.bam"

# Recorded Week 1 workflow, section 11.
picard AddOrReplaceReadGroups \
  I="$FILTERED/${SAMPLE}.filtered.sorted.bam" \
  O="$FILTERED/${SAMPLE}.filtered.sorted.rg.bam" \
  RGID="$SAMPLE" \
  RGLB="ATACseq" \
  RGPL="ILLUMINA" \
  RGPU="unit1" \
  RGSM="$SAMPLE" \
  SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=SILENT

# Recorded Week 1 workflow, section 12.
picard MarkDuplicates \
  I="$FILTERED/${SAMPLE}.filtered.sorted.rg.bam" \
  O="$DEDUP_BAM" \
  M="$FILTERED/${SAMPLE}.duplicate_metrics.txt" \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=SILENT

# Recorded Week 1 workflow, section 13.
samtools index -@ 2 "$DEDUP_BAM"

# Recorded Week 1 workflow, section 14.
samtools flagstat \
  -@ 2 \
  "$DEDUP_BAM" \
  > "$QC/${SAMPLE}.final.flagstat.txt"

# Recorded Week 1 workflow, section 14.
samtools stats \
  -@ 2 \
  "$DEDUP_BAM" \
  > "$QC/${SAMPLE}.final.samtools_stats.txt"

# Recorded Week 1 workflow, section 14.
samtools idxstats \
  "$DEDUP_BAM" \
  > "$QC/${SAMPLE}.final.idxstats.txt"

# Recorded Week 1 workflow, section 14.
picard CollectInsertSizeMetrics \
  I="$DEDUP_BAM" \
  O="$QC/${SAMPLE}.insert_size_metrics.txt" \
  H="$QC/${SAMPLE}.insert_size_histogram.pdf" \
  M=0.5 \
  VALIDATION_STRINGENCY=SILENT

# Recorded Week 1 workflow, section 15.
macs2 callpeak \
  -t "$DEDUP_BAM" \
  -f BAMPE \
  -g 50818468 \
  -n "${SAMPLE}_ATAC_chr22" \
  --outdir "$PEAKS" \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  > "$PEAKS/macs2.log" 2>&1

# Recorded Week 1 workflow, section 16.
mkdir -p "$MULTIQC_INPUTS/fastqc"
mkdir -p "$MULTIQC_INPUTS/cutadapt"
mkdir -p "$MULTIQC_INPUTS/bowtie2"
mkdir -p "$MULTIQC_INPUTS/picard"
mkdir -p "$MULTIQC_INPUTS/samtools"
mkdir -p "$MULTIQC_INPUTS/macs2"

# Recorded Week 1 workflow, section 16.
cp "$QC/fastqc_raw/"* "$MULTIQC_INPUTS/fastqc/"
cp "$QC/fastqc_trimmed/"* "$MULTIQC_INPUTS/fastqc/"
cp "$TRIMMED/${SAMPLE}.cutadapt.log" "$MULTIQC_INPUTS/cutadapt/"
cp "$ALIGN/bowtie2.log" "$MULTIQC_INPUTS/bowtie2/"
cp "$FILTERED/${SAMPLE}.duplicate_metrics.txt" "$MULTIQC_INPUTS/picard/"
cp "$QC/${SAMPLE}.insert_size_metrics.txt" "$MULTIQC_INPUTS/picard/"
cp "$QC/${SAMPLE}.final.flagstat.txt" "$MULTIQC_INPUTS/samtools/"
cp "$QC/${SAMPLE}.final.samtools_stats.txt" "$MULTIQC_INPUTS/samtools/"
cp "$QC/${SAMPLE}.final.idxstats.txt" "$MULTIQC_INPUTS/samtools/"
cp "$PEAKS/macs2.log" "$MULTIQC_INPUTS/macs2/"
cp "$PEAKS/${SAMPLE}_ATAC_chr22_peaks.narrowPeak" "$MULTIQC_INPUTS/macs2/"

# Recorded Week 1 workflow, section 16.
printf "Sample\tPeaks\n%s\t%s\n" \
  "${SAMPLE}.final" \
  "$(wc -l < "$PEAKS/${SAMPLE}_ATAC_chr22_peaks.narrowPeak")" \
  > "$MULTIQC_INPUTS/macs2/macs2_peaks_mqc.tsv"

# Recorded Week 1 workflow, section 17.
multiqc \
  "$PROJECT/multiqc_inputs" \
  -o "$PROJECT/multiqc_output" \
  --force

printf 'Week 1 teaching outputs: %s\n' "$PROJECT"

