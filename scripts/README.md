# Running the analysis

Run these scripts from the repository root. They use the eight-sample,
paired-donor design listed in `input/samples.tsv`.

## R analysis from the supplied inputs

Use R with the packages loaded by the scripts:
SummarizedExperiment, DESeq2, apeglm, vsn, ggplot2, pheatmap, GenomicRanges,
GenomicFeatures, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db,
clusterProfiler and enrichplot. Bioconductor packages can be installed through
[BiocManager](https://www.bioconductor.org/install/).
The recorded analysis used R 4.3.1 and Bioconductor 3.18.

In a fresh terminal:

```bash
Rscript scripts/differential_accessibility.R
Rscript scripts/peak_annotation_GO_and_HOMER_input.R
Rscript scripts/HOMER_summary_and_plots.R
```

The first two commands use `input/`. The last summarises HOMER outputs.
New tables and plots go under `runs/local/`; R scripts can replace files there.

For a different output folder, set `ATAC_RUN_DIR` to a subfolder of `runs/`.
`ATAC_INPUT_DIR` optionally selects a different input folder for R.
To return to defaults:

```bash
unset ATAC_RUN_DIR ATAC_INPUT_DIR
```

## Consensus peaks, counts and motif enrichment

Requirements: Bash 4+, GNU coreutils, bedtools, Subread/featureCounts,
and HOMER with its hg38 genome package. Use WSL/Linux.
The saved counting output records featureCounts 2.1.1.
Input BAMs must be cleaned, duplicate-removed paired-end hg38 alignments.

Provide exactly one nonempty broadPeak file and one BAM per accession in
`input/samples.tsv`, with filenames beginning with that accession. Each BAM
must have a `.bam.bai` or `.bai` index. Peak and BAM chromosome names must agree.
These large upstream files and reference genomes are not included.

```bash
export ATAC_RUN_DIR=runs/rebuilt
bash scripts/01_consensus_peaks.sh /path/to/broadPeak_directory
bash scripts/02_count_fragments.sh /path/to/cleaned_BAM_directory
export ATAC_INPUT_DIR=runs/rebuilt/input
Rscript scripts/differential_accessibility.R
Rscript scripts/peak_annotation_GO_and_HOMER_input.R
bash scripts/05_homer_motifs.sh
Rscript scripts/HOMER_summary_and_plots.R
```

Replace the placeholder paths with your input directories. Rebuilt inputs,
raw counting output, logs, tables and plots stay under `runs/rebuilt/`.
Bash output guards reject existing outputs. Choose a new run folder to repeat
those stages. `runs/` is ignored by Git.

The motif runner uses fixed 200-bp windows, separately testing Treg and naive
DARs against non-differential accessible peaks. Promoter and non-promoter
analyses use the corresponding region-matched backgrounds.
`HOMER_THREADS` sets the thread count; the default is two.

Software, genome, annotation and motif-database versions can affect new results.
The Bash scripts use configurable paths and input checks. The HOMER runner
was reconstructed from the analysis methods; the original shell script was
not available. These Bash scripts have not been rerun on the upstream data.

## Coordinates and saved results

BED uses zero-based starts; SAF uses one-based, inclusive coordinates.
`01_consensus_peaks.sh` adds one to each BED start and keeps the end unchanged.
See the [bedtools coordinate description](https://bedtools.readthedocs.io/en/latest/content/general-usage.html)
and [Subread user's guide, section 6.2.2](https://subread.sourceforge.net/SubreadUsersGuide.pdf).

For example, BED `chr1 100 200` becomes SAF start `101`, end `200`.
The saved counts used unchanged BED starts. Recounting with corrected
coordinates may change counts and DARs; this effect has not been measured.

`recorded/` keeps the Week 2 commands for comparison, including the original
SAF conversion and sample order. Set `ATAC_RECORDED_BASE` to their data folder.
Use the main scripts above when rebuilding counts.

## BATF/LEF1 gene follow-up

The saved motif-hit tables identify BATF-positive Treg peaks and LEF1-positive
naive peaks and their nearest genes. Motif enrichment does not recreate these scans.

With rebuilt inputs, the final R script makes enrichment plots and skips
this additional gene follow-up unless newly scanned tables are present as:

```text
runs/rebuilt/motif_peak_analysis/treg_BATF_motif_nearby_genes.txt
runs/rebuilt/motif_peak_analysis/naive_LEF1_motif_nearby_genes.txt
```

Use the same columns as the supplied examples, including `NearestGene`.
Nearest-gene associations do not establish direct regulatory targets.

## Separate Week 1 chromosome-22 teaching exercise

The Week 1 script runs FastQC, Cutadapt, Bowtie2, SAMtools, Picard, MACS2
and MultiQC on teaching sample SRR7650763. It is separate from the
eight-sample whole-genome comparison.

Activate an environment containing those tools, then supply the two original
chromosome-22 FASTQ files and a complete chromosome-22 Bowtie2 index:

```bash
export ATAC_RUN_DIR=runs/teaching
bash scripts/00_week1_chr22_training.sh \
  /path/to/SRR7650763_1_chr22.raw.fastq.gz \
  /path/to/SRR7650763_2_chr22.raw.fastq.gz \
  /path/to/GRCh38_chr22
```

Outputs go under `runs/teaching/week1_chr22/`. Without `ATAC_RUN_DIR`, the
default is `runs/week1/week1_chr22/`. Existing teaching outputs are protected.
This teaching script does not include blacklist filtering or whole-genome
preprocessing.
