# Public-release file routing; analytical operations retained from the supplied script.
source("scripts/project_paths.R")

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)


consensus_bed <- read.delim(
  input_path("consensus_peaks.bed"),
  header = FALSE,
  col.names = c("Chr", "Start", "End", "PeakID", "Score", "Strand"),
  stringsAsFactors = FALSE
)

missing_chr <- !grepl("^chr", consensus_bed$Chr)
consensus_bed$Chr[missing_chr] <- paste0("chr", consensus_bed$Chr[missing_chr])
consensus_bed$Chr[consensus_bed$Chr == "chrMT"] <- "chrM"

all_results <- read.delim(read_result_path("all_DA_results.txt"), header = TRUE, check.names = FALSE)
significant_results <- read.delim(read_result_path("significant_DARs.txt"), header = TRUE, check.names = FALSE)
treg_results <- read.delim(read_result_path("regulatory_accessible_DARs.txt"), header = TRUE, check.names = FALSE)
naive_results <- read.delim(read_result_path("naive_accessible_DARs.txt"), header = TRUE, check.names = FALSE)

cat("Consensus peaks:", nrow(consensus_bed), "\n")
cat("All tested peaks:", nrow(all_results), "\n")
cat("Significant DARs:", nrow(significant_results), "\n")
cat("Treg-accessible DARs:", nrow(treg_results), "\n")
cat("Naive-accessible DARs:", nrow(naive_results), "\n")

add_coordinates <- function(result_table, bed_table) {
  matched_rows <- match(result_table$PeakID, bed_table$PeakID)
  stopifnot(!anyNA(matched_rows))
  output <- cbind(result_table, bed_table[matched_rows, c("Chr", "Start", "End")])
  rownames(output) <- NULL
  output
}

all_results_coords <- add_coordinates(all_results, consensus_bed)
significant_results_coords <- add_coordinates(significant_results, consensus_bed)
treg_results_coords <- add_coordinates(treg_results, consensus_bed)
naive_results_coords <- add_coordinates(naive_results, consensus_bed)

write.table(
  significant_results_coords,
  result_path("significant_DARs_with_coordinates.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  treg_results_coords,
  result_path("treg_DARs_with_coordinates.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  naive_results_coords,
  result_path("naive_DARs_with_coordinates.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

make_granges <- function(peak_table) {
  makeGRangesFromDataFrame(
    peak_table,
    seqnames.field = "Chr",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    starts.in.df.are.0based = TRUE
  )
}

all_gr <- make_granges(all_results_coords)
treg_gr <- make_granges(treg_results_coords)
naive_gr <- make_granges(naive_results_coords)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

annotate_peaks <- function(peaks) {
  annotatePeak(
    peaks,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
}

all_annotation_df <- as.data.frame(annotate_peaks(all_gr))
treg_annotation_df <- as.data.frame(annotate_peaks(treg_gr))
naive_annotation_df <- as.data.frame(annotate_peaks(naive_gr))

write.table(all_annotation_df, result_path("all_peaks_annotated.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(treg_annotation_df, result_path("treg_DARs_annotated.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(naive_annotation_df, result_path("naive_DARs_annotated.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

simplify_annotation <- function(annotation) {
  category <- rep("Other or unannotated", length(annotation))
  category[grepl("Promoter", annotation)] <- "Promoter"
  category[grepl("UTR", annotation)] <- "UTR"
  category[grepl("Exon", annotation)] <- "Exon"
  category[grepl("Intron", annotation)] <- "Intron"
  category[grepl("Downstream", annotation)] <- "Downstream"
  category[grepl("Intergenic", annotation)] <- "Distal intergenic"
  category
}

treg_annotation_df$Category <- simplify_annotation(treg_annotation_df$annotation)
naive_annotation_df$Category <- simplify_annotation(naive_annotation_df$annotation)

annotation_comparison <- rbind(
  data.frame(Cell_type = "Treg-accessible", Category = treg_annotation_df$Category),
  data.frame(Cell_type = "Naive-accessible", Category = naive_annotation_df$Category)
)

annotation_comparison$Cell_type <- factor(
  annotation_comparison$Cell_type,
  levels = c("Treg-accessible", "Naive-accessible")
)

annotation_comparison$Category <- factor(
  annotation_comparison$Category,
  levels = c(
    "Promoter",
    "UTR",
    "Exon",
    "Intron",
    "Downstream",
    "Distal intergenic",
    "Other or unannotated"
  )
)

annotation_summary <- as.data.frame(
  prop.table(table(annotation_comparison$Cell_type, annotation_comparison$Category), margin = 1) * 100
)
colnames(annotation_summary) <- c("Cell_type", "Category", "Percentage")

write.table(annotation_summary, result_path("genomic_annotation_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

p_annotation <- ggplot(annotation_summary, aes(x = Category, y = Percentage, fill = Cell_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f", Percentage)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3
  ) +
  scale_fill_manual(values = c("Treg-accessible" = "#D55E00", "Naive-accessible" = "#0072B2")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Genomic Distribution of Differentially Accessible Regions",
    subtitle = "Regulatory T cells versus naive conventional CD4 T cells",
    x = "Genomic annotation",
    y = "Percentage of DARs",
    fill = "Accessibility direction"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

ggsave(
  plot_path("07_genomic_annotation_distribution.png"),
  plot = p_annotation,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

treg_candidates <- c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "ENTPD1")
naive_candidates <- c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "MAL")

candidate_columns <- c(
  "PeakID",
  "SYMBOL",
  "GENENAME",
  "annotation",
  "distanceToTSS",
  "log2FoldChange",
  "padj"
)

treg_candidate_results <- treg_annotation_df[
  treg_annotation_df$SYMBOL %in% treg_candidates,
  candidate_columns
]
treg_candidate_results <- treg_candidate_results[
  order(treg_candidate_results$SYMBOL, treg_candidate_results$padj),
]

naive_candidate_results <- naive_annotation_df[
  naive_annotation_df$SYMBOL %in% naive_candidates,
  candidate_columns
]
naive_candidate_results <- naive_candidate_results[
  order(naive_candidate_results$SYMBOL, naive_candidate_results$padj),
]

write.table(
  treg_candidate_results,
  result_path("treg_candidate_gene_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  naive_candidate_results,
  result_path("naive_candidate_gene_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

rank_genes <- function(annotation_table) {
  gene_peaks <- annotation_table[
    !is.na(annotation_table$SYMBOL) &
      annotation_table$SYMBOL != "" &
      !is.na(annotation_table$padj),
  ]

  peaks_by_gene <- split(gene_peaks, gene_peaks$SYMBOL)
  gene_summary <- lapply(peaks_by_gene, function(peaks) {
    best_peak <- peaks[which.min(peaks$padj), ]
    data.frame(
      SYMBOL = best_peak$SYMBOL,
      GENENAME = best_peak$GENENAME,
      Number_of_DARs = nrow(peaks),
      Promoter_DARs = sum(grepl("Promoter", peaks$annotation)),
      Best_peak = best_peak$PeakID,
      Best_annotation = best_peak$annotation,
      Distance_to_TSS = best_peak$distanceToTSS,
      Log2_fold_change = best_peak$log2FoldChange,
      Minimum_padj = best_peak$padj
    )
  })

  gene_summary <- do.call(rbind, gene_summary)
  gene_summary <- gene_summary[
    order(
      gene_summary$Minimum_padj,
      -abs(gene_summary$Log2_fold_change),
      -gene_summary$Number_of_DARs
    ),
  ]
  rownames(gene_summary) <- NULL
  gene_summary
}

treg_gene_ranking <- rank_genes(treg_annotation_df)
naive_gene_ranking <- rank_genes(naive_annotation_df)

write.table(treg_gene_ranking, result_path("treg_gene_ranking.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(naive_gene_ranking, result_path("naive_gene_ranking.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

background_gene_ids <- unique(na.omit(as.character(all_annotation_df$geneId)))
treg_gene_ids <- unique(na.omit(as.character(treg_annotation_df$geneId)))
naive_gene_ids <- unique(na.omit(as.character(naive_annotation_df$geneId)))

background_gene_ids <- background_gene_ids[background_gene_ids != ""]
treg_gene_ids <- intersect(treg_gene_ids[treg_gene_ids != ""], background_gene_ids)
naive_gene_ids <- intersect(naive_gene_ids[naive_gene_ids != ""], background_gene_ids)

go_treg <- enrichGO(
  gene = treg_gene_ids,
  universe = background_gene_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_naive <- enrichGO(
  gene = naive_gene_ids,
  universe = background_gene_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_treg_df <- as.data.frame(go_treg)
go_naive_df <- as.data.frame(go_naive)

write.table(go_treg_df, result_path("treg_GO_biological_process.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(go_naive_df, result_path("naive_GO_biological_process.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

p_go_treg <- dotplot(go_treg, showCategory = 15, font.size = 10, label_format = 50) +
  labs(
    title = "GO Enrichment of Treg-Accessible DAR-Associated Genes",
    subtitle = "Biological Process",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p_go_naive <- dotplot(go_naive, showCategory = 15, font.size = 10, label_format = 50) +
  labs(
    title = "GO Enrichment of Naive-Accessible DAR-Associated Genes",
    subtitle = "Biological Process",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

ggsave(plot_path("08_Treg_GO_enrichment.png"), plot = p_go_treg, width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(plot_path("09_Naive_GO_enrichment.png"), plot = p_go_naive, width = 10, height = 7, units = "in", dpi = 300, bg = "white")

standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
treg_homer_bed <- treg_results_coords[treg_results_coords$Chr %in% standard_chromosomes, c("Chr", "Start", "End", "PeakID")]
naive_homer_bed <- naive_results_coords[naive_results_coords$Chr %in% standard_chromosomes, c("Chr", "Start", "End", "PeakID")]

background_results_coords <- all_results_coords[
  !all_results_coords$PeakID %in% significant_results_coords$PeakID,
]
background_homer_bed <- background_results_coords[
  background_results_coords$Chr %in% standard_chromosomes,
  c("Chr", "Start", "End", "PeakID")
]

peak_width_summary <- data.frame(
  Peak_set = c("Treg DARs", "Naive DARs", "Background"),
  Number_of_peaks = c(nrow(treg_homer_bed), nrow(naive_homer_bed), nrow(background_homer_bed)),
  Median_width = c(
    median(treg_homer_bed$End - treg_homer_bed$Start),
    median(naive_homer_bed$End - naive_homer_bed$Start),
    median(background_homer_bed$End - background_homer_bed$Start)
  ),
  Mean_width = c(
    mean(treg_homer_bed$End - treg_homer_bed$Start),
    mean(naive_homer_bed$End - naive_homer_bed$Start),
    mean(background_homer_bed$End - background_homer_bed$Start)
  ),
  Minimum_width = c(
    min(treg_homer_bed$End - treg_homer_bed$Start),
    min(naive_homer_bed$End - naive_homer_bed$Start),
    min(background_homer_bed$End - background_homer_bed$Start)
  ),
  Maximum_width = c(
    max(treg_homer_bed$End - treg_homer_bed$Start),
    max(naive_homer_bed$End - naive_homer_bed$Start),
    max(background_homer_bed$End - background_homer_bed$Start)
  )
)

write.table(peak_width_summary, result_path("HOMER_peak_width_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(treg_homer_bed, homer_input_path("treg_DARs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(naive_homer_bed, homer_input_path("naive_DARs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(background_homer_bed, homer_input_path("background_peaks.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

is_promoter <- function(annotation) {
  !is.na(annotation) & grepl("^Promoter", annotation)
}

treg_promoter_ids <- treg_annotation_df$PeakID[is_promoter(treg_annotation_df$annotation)]
naive_promoter_ids <- naive_annotation_df$PeakID[is_promoter(naive_annotation_df$annotation)]
background_promoter_ids <- all_annotation_df$PeakID[
  all_annotation_df$PeakID %in% background_homer_bed$PeakID & is_promoter(all_annotation_df$annotation)
]

treg_promoter_bed <- treg_homer_bed[treg_homer_bed$PeakID %in% treg_promoter_ids, ]
treg_nonpromoter_bed <- treg_homer_bed[!treg_homer_bed$PeakID %in% treg_promoter_ids, ]
naive_promoter_bed <- naive_homer_bed[naive_homer_bed$PeakID %in% naive_promoter_ids, ]
naive_nonpromoter_bed <- naive_homer_bed[!naive_homer_bed$PeakID %in% naive_promoter_ids, ]
background_promoter_bed <- background_homer_bed[background_homer_bed$PeakID %in% background_promoter_ids, ]
background_nonpromoter_bed <- background_homer_bed[!background_homer_bed$PeakID %in% background_promoter_ids, ]

region_counts <- data.frame(
  Peak_set = c(
    "Treg promoter",
    "Treg non-promoter",
    "Naive promoter",
    "Naive non-promoter",
    "Background promoter",
    "Background non-promoter"
  ),
  Number_of_peaks = c(
    nrow(treg_promoter_bed),
    nrow(treg_nonpromoter_bed),
    nrow(naive_promoter_bed),
    nrow(naive_nonpromoter_bed),
    nrow(background_promoter_bed),
    nrow(background_nonpromoter_bed)
  )
)

write.table(region_counts, result_path("HOMER_region_counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(treg_promoter_bed, homer_input_path("treg_promoter_DARs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(treg_nonpromoter_bed, homer_input_path("treg_nonpromoter_DARs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(naive_promoter_bed, homer_input_path("naive_promoter_DARs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(naive_nonpromoter_bed, homer_input_path("naive_nonpromoter_DARs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(background_promoter_bed, homer_input_path("background_promoter_peaks.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(background_nonpromoter_bed, homer_input_path("background_nonpromoter_peaks.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

candidate_symbols <- c("FOXP3", "IL2RA", "CTLA4", "TCF7", "LEF1", "CCR7")
candidate_entrez <- mapIds(
  org.Hs.eg.db,
  keys = candidate_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
candidate_entrez <- unname(na.omit(candidate_entrez))

transcripts_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
candidate_transcript_list <- transcripts_by_gene[candidate_entrez]
candidate_transcripts <- unlist(candidate_transcript_list, use.names = FALSE)
candidate_transcripts$geneId <- rep(
  names(candidate_transcript_list),
  S4Vectors::elementNROWS(candidate_transcript_list)
)
candidate_transcripts$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = candidate_transcripts$geneId,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

candidate_promoter_gr <- promoters(candidate_transcripts, upstream = 3000, downstream = 1000)

collect_promoter_hits <- function(peaks, promoters_gr, accessibility) {
  hits <- findOverlaps(peaks, promoters_gr, ignore.strand = TRUE)
  output <- as.data.frame(peaks)[queryHits(hits), , drop = FALSE]
  rownames(output) <- NULL
  output$CandidateGene <- promoters_gr$SYMBOL[subjectHits(hits)]
  output <- output[!duplicated(output[c("PeakID", "CandidateGene")]), ]
  output$Accessibility <- accessibility
  output
}

treg_candidate_promoter_hits <- collect_promoter_hits(
  treg_gr,
  candidate_promoter_gr,
  "Treg-accessible"
)
naive_candidate_promoter_hits <- collect_promoter_hits(
  naive_gr,
  candidate_promoter_gr,
  "Naive-accessible"
)

candidate_promoter_results <- rbind(
  treg_candidate_promoter_hits,
  naive_candidate_promoter_hits
)
candidate_promoter_results <- candidate_promoter_results[
  ,
  c(
    "PeakID",
    "CandidateGene",
    "seqnames",
    "start",
    "end",
    "log2FoldChange",
    "padj",
    "Accessibility"
  )
]
rownames(candidate_promoter_results) <- NULL

write.table(
  candidate_promoter_results,
  result_path("candidate_gene_promoter_DARs.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(go_treg, result_path("go_treg.rds"))
saveRDS(go_naive, result_path("go_naive.rds"))
capture.output(sessionInfo(), file = result_path("week3_R_session_info.txt"))
