# Public-release file routing; analytical operations retained from the supplied script.
source("scripts/project_paths.R")

library(SummarizedExperiment)
library(DESeq2)
library(apeglm)
library(vsn)
library(ggplot2)
library(pheatmap)


count_table <- read.delim(input_path("counts_clean.txt"), header = TRUE, check.names = FALSE)
metadata <- read.delim(input_path("samples.tsv"), header = TRUE, stringsAsFactors = FALSE)

stopifnot(!anyDuplicated(count_table$PeakID), !anyDuplicated(metadata$ID))
stopifnot(setequal(metadata$ID, colnames(count_table)[-1]))
stopifnot(all(is.finite(as.matrix(count_table[, -1, drop = FALSE]))))
stopifnot(all(as.matrix(count_table[, -1, drop = FALSE]) >= 0))
stopifnot(all(as.matrix(count_table[, -1, drop = FALSE]) == floor(as.matrix(count_table[, -1, drop = FALSE]))))

count_matrix <- as.matrix(count_table[, -1, drop = FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- count_table$PeakID

rownames(metadata) <- metadata$ID
metadata <- metadata[colnames(count_matrix), , drop = FALSE]
stopifnot(identical(rownames(metadata), colnames(count_matrix)))

metadata$donor <- factor(metadata$donor)
metadata$cell_type <- factor(metadata$cell_type, levels = c("Naive_Teffs", "Regulatory_T"))

se <- SummarizedExperiment(assays = list(counts = count_matrix), colData = metadata)
keep <- rowSums(assay(se) >= 10) >= 4
se_filtered <- se[keep, ]

dds <- DESeqDataSet(se_filtered, design = ~ donor + cell_type)
dds <- DESeq(dds)

normalized_log_counts <- log2(counts(dds, normalized = TRUE) + 1)
mean_sd_before <- meanSdPlot(normalized_log_counts, plot = FALSE)

p_mean_sd_before <- mean_sd_before$gg +
  labs(
    title = "Mean-SD Relationship Before VST",
    subtitle = "Log2-transformed normalized fragment counts",
    x = "Rank of mean signal",
    y = "Standard deviation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5)
  )

ggsave(
  plot_path("01_meanSD_before_VST.png"),
  plot = p_mean_sd_before,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

vsd <- vst(dds, blind = FALSE)
mean_sd_after <- meanSdPlot(assay(vsd), plot = FALSE)

p_mean_sd_after <- mean_sd_after$gg +
  labs(
    title = "Mean-SD Relationship After VST",
    subtitle = "Variance-stabilized fragment counts",
    x = "Rank of mean signal",
    y = "Standard deviation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5)
  )

ggsave(
  plot_path("02_meanSD_after_VST.png"),
  plot = p_mean_sd_after,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

pca_data <- plotPCA(vsd, intgroup = c("cell_type", "donor"), returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_type, shape = donor)) +
  geom_point(size = 4) +
  geom_text(aes(label = name), vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(
    values = c("Naive_Teffs" = "#2CA58D", "Regulatory_T" = "#E66101"),
    labels = c("Naive conventional CD4 T", "Regulatory T")
  ) +
  labs(
    title = "PCA of VST-Transformed ATAC-seq Counts",
    subtitle = "Colour represents cell type; shape represents donor",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    color = "Cell type",
    shape = "Donor"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5)
  )

ggsave(
  plot_path("03_PCA_VST.png"),
  plot = p_pca,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

sample_distances <- dist(t(assay(vsd)))
sample_distance_matrix <- as.matrix(sample_distances)
annotation_col <- as.data.frame(colData(vsd)[, c("donor", "cell_type")])

pheatmap(
  sample_distance_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  clustering_distance_rows = sample_distances,
  clustering_distance_cols = sample_distances,
  main = "Sample-to-Sample Distances After VST",
  filename = plot_path("04_sample_distance_heatmap.png"),
  width = 9,
  height = 8
)

coefficient <- "cell_type_Regulatory_T_vs_Naive_Teffs"
stopifnot(coefficient %in% resultsNames(dds))

res_shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm")
all_results <- data.frame(PeakID = rownames(res_shrunk), as.data.frame(res_shrunk), check.names = FALSE)
all_results <- all_results[order(all_results$padj, na.last = TRUE), ]

significant_results <- all_results[
  !is.na(all_results$padj) &
    all_results$padj < 0.05 &
    abs(all_results$log2FoldChange) >= 1,
]

treg_results <- significant_results[significant_results$log2FoldChange >= 1, ]
naive_results <- significant_results[significant_results$log2FoldChange <= -1, ]

write.table(all_results, result_path("all_DA_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(significant_results, result_path("significant_DARs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(treg_results, result_path("regulatory_accessible_DARs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(naive_results, result_path("naive_accessible_DARs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

analysis_summary <- data.frame(
  Category = c(
    "Consensus peaks",
    "Peaks tested by DESeq2",
    "Significant DARs",
    "Treg-accessible DARs",
    "Naive-accessible DARs"
  ),
  Number_of_peaks = c(
    nrow(count_matrix),
    nrow(dds),
    nrow(significant_results),
    nrow(treg_results),
    nrow(naive_results)
  )
)

write.table(analysis_summary, result_path("differential_accessibility_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

png(plot_path("05_MA_plot_DA.png"), width = 1800, height = 1400, res = 200)
DESeq2::plotMA(
  res_shrunk,
  alpha = 0.05,
  ylim = c(-5, 5),
  main = "Differential Accessibility: Regulatory T vs Naive T"
)
dev.off()

volcano_df <- all_results[
  !is.na(all_results$padj) & !is.na(all_results$log2FoldChange),
]

smallest_nonzero_padj <- min(volcano_df$padj[volcano_df$padj > 0], na.rm = TRUE)
volcano_df$padj_for_plot <- ifelse(
  volcano_df$padj == 0,
  smallest_nonzero_padj / 10,
  volcano_df$padj
)

volcano_df$minus_log10_padj <- -log10(volcano_df$padj_for_plot)
volcano_df$Direction <- "Not significant"
volcano_df$Direction[
  volcano_df$padj < 0.05 & volcano_df$log2FoldChange >= 1
] <- "More accessible in Regulatory T"
volcano_df$Direction[
  volcano_df$padj < 0.05 & volcano_df$log2FoldChange <= -1
] <- "More accessible in Naive T"

volcano_df$Direction <- factor(
  volcano_df$Direction,
  levels = c(
    "Not significant",
    "More accessible in Regulatory T",
    "More accessible in Naive T"
  )
)

p_volcano <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = minus_log10_padj, color = Direction)
) +
  geom_point(size = 1.3, alpha = 0.65) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey75",
      "More accessible in Regulatory T" = "#E66101",
      "More accessible in Naive T" = "#0072B2"
    )
  ) +
  labs(
    title = "Differential Chromatin Accessibility",
    subtitle = "Regulatory T cells versus naive conventional CD4 T cells",
    x = "Shrunken log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = "Peak classification",
    caption = "Significant: adjusted p-value < 0.05 and |log2 fold change| >= 1"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey35", hjust = 0.5),
    plot.caption = element_text(size = 9, color = "grey40", hjust = 1),
    plot.title.position = "plot",
    legend.position = "right"
  )

ggsave(
  plot_path("06_volcano_plot_adjusted_pvalue.png"),
  plot = p_volcano,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

capture.output(sessionInfo(), file = result_path("R_session_info.txt"))
