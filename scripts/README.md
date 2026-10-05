# Running the analysis

These scripts support the Treg/naive CD4+ T-cell case study in the main README.
They are open source and can be adapted to other experiments by updating
sample metadata, input matching, count filters and statistical contrasts.
The supplied runners expect this project's eight-sample paired design.

Browsing or publishing the repository does not require running the analysis.
Use the main scripts below for a new run; `recorded/` retains the earlier
Week 2 command bodies for comparison with the saved results.

## R analysis from the supplied inputs

Run commands from the repository root, the folder containing the main README.
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

The first two commands use the supplied count matrix, consensus BED and
metadata in `input/`. The last command summarises saved HOMER outputs; it does
not run HOMER. Newly generated tables and plots go under `runs/local/`.
The R scripts may replace outputs within that selected run directory.
Saved files in `input/` and `results/` remain separate.

For a different output folder, set `ATAC_RUN_DIR` to a subfolder of `runs/`.
`ATAC_INPUT_DIR` optionally selects a different input folder for R.
Both settings apply to the current terminal session. To return to defaults:

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
`HOMER_THREADS` optionally changes the public runner's two-thread default.

Genome, annotation, motif-database and software versions can affect new results.
The public Bash runners have been checked for syntax and wrapper behaviour
but have not been scientifically rerun against the original upstream data.

## Coordinate correction and historical commands

The recovered Week 2 script and the course example pass BED starts unchanged
into SAF. The original featureCounts table records these same start positions.
The main public script adds one to each zero-based BED start for one-based,
inclusive SAF coordinates, while keeping the end unchanged.
See the [bedtools coordinate description](https://bedtools.readthedocs.io/en/latest/content/general-usage.html)
and [Subread user's guide, section 6.2.2](https://subread.sourceforge.net/SubreadUsersGuide.pdf).

For example, BED `chr1 100 200` becomes SAF start `101`, end `200` for the same
100-base interval. Retaining start `100` in SAF adds one base to its left edge.
Saved counts and results are the original analysis snapshot. Corrected counts
may differ, and their effect on downstream DARs has not been measured.

`recorded/01_consensus_peaks.sh` and `recorded/02_count_fragments.sh` retain
the July Week 2 processing command bodies, including the historical SAF
conversion and sample order. Only the institutional base path is replaced by
`ATAC_RECORDED_BASE`, and existing-output checks are added. These files are
historical reference code, rather than the recommended corrected entry points.

## BATF/LEF1 gene follow-up

The supplied motif-hit tables identify BATF-positive Treg peaks and LEF1-positive
naive peaks and their nearest genes. They belong to the saved analysis.
The motif runner performs enrichment but does not recreate these motif scans.

With rebuilt inputs, the final R script makes enrichment plots and skips
this additional gene follow-up unless newly scanned tables are present as:

```text
runs/rebuilt/motif_peak_analysis/treg_BATF_motif_nearby_genes.txt
runs/rebuilt/motif_peak_analysis/naive_LEF1_motif_nearby_genes.txt
```

Use the same columns as the supplied examples, including `NearestGene`.
Nearest-gene associations do not establish direct regulatory targets.

## Separate Week 1 chromosome-22 teaching exercise

The preprocessing script adapts the commands in the saved workflow of
26 July 2026. It includes FastQC, Cutadapt, Bowtie2, SAMtools, Picard, MACS2
and MultiQC. It uses teaching sample SRR7650763 and is independent of the
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
The public version changes file routing and uses the saved standard MultiQC
command without an unavailable institution-specific configuration.
Blacklist exclusion was not present in the recovered command document and
is not reconstructed in this script. This is a preserved teaching workflow,
not a complete historical execution log or a whole-genome preprocessing recipe.
