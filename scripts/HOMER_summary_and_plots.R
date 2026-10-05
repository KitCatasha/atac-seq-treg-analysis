# Public-release file routing; analytical operations retained from the supplied script.
source("scripts/project_paths.R")

library(ggplot2)


extract_motif <- function(result_file, motif_pattern, motif_label, dar_set, context = "All DARs") {
  motif_results <- read.delim(result_file, header = TRUE, check.names = FALSE)
  target_column <- grep("^% of Target Sequences with Motif", colnames(motif_results), value = TRUE)
  background_column <- grep("^% of Background Sequences with Motif", colnames(motif_results), value = TRUE)
  motif_row <- motif_results[grepl(motif_pattern, motif_results[["Motif Name"]]), , drop = FALSE]

  stopifnot(length(target_column) == 1, length(background_column) == 1, nrow(motif_row) >= 1)
  motif_row <- motif_row[1, , drop = FALSE]

  target_percentage <- as.numeric(sub("%", "", motif_row[[target_column]]))
  background_percentage <- as.numeric(sub("%", "", motif_row[[background_column]]))

  data.frame(
    Motif = motif_label,
    DAR_set = dar_set,
    Context = context,
    Target_percentage = target_percentage,
    Background_percentage = background_percentage,
    Fold_enrichment = target_percentage / background_percentage,
    P_value = motif_row[["P-value"]]
  )
}

primary_motifs <- rbind(
  extract_motif(homer_result_path("treg_200/knownResults.txt"), "^JunB\\(bZIP\\)", "JunB", "Treg-accessible"),
  extract_motif(homer_result_path("treg_200/knownResults.txt"), "^Fra2\\(bZIP\\)", "Fra2", "Treg-accessible"),
  extract_motif(homer_result_path("treg_200/knownResults.txt"), "^BATF\\(bZIP\\)", "BATF", "Treg-accessible"),
  extract_motif(homer_result_path("treg_200/knownResults.txt"), "^Atf3\\(bZIP\\)", "ATF3", "Treg-accessible"),
  extract_motif(homer_result_path("naive_200/knownResults.txt"), "^LEF1\\(HMG\\)", "LEF1", "Naive-accessible"),
  extract_motif(homer_result_path("naive_200/knownResults.txt"), "^RUNX1\\(Runt\\)", "RUNX1", "Naive-accessible"),
  extract_motif(homer_result_path("naive_200/knownResults.txt"), "^Tcf3\\(HMG\\)", "TCF3", "Naive-accessible"),
  extract_motif(homer_result_path("naive_200/knownResults.txt"), "^Tcf7\\(HMG\\)", "TCF7", "Naive-accessible")
)

primary_motifs$Motif <- factor(
  primary_motifs$Motif,
  levels = rev(c("JunB", "Fra2", "BATF", "ATF3", "LEF1", "RUNX1", "TCF3", "TCF7"))
)

write.table(
  primary_motifs,
  result_path("HOMER_primary_motif_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_primary_motifs <- ggplot(primary_motifs, aes(x = Motif, y = Fold_enrichment, fill = DAR_set)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Fold_enrichment)), hjust = -0.15, size = 3.5, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Treg-accessible" = "#7B2CBF", "Naive-accessible" = "#2E9F75")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Distinct TF Motifs Separate Treg and Naive DARs",
    subtitle = "Known motifs analysed using fixed 200-bp windows",
    x = NULL,
    y = "Fold enrichment: target / background",
    fill = "DAR set"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  plot_path("primary_motifs.png"),
  plot = p_primary_motifs,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

region_motifs <- rbind(
  extract_motif(
    homer_result_path("treg_nonpromoter_200/knownResults.txt"),
    "^AP-1\\(bZIP\\)",
    "AP-1\n(Treg)",
    "Treg-accessible",
    "Non-promoter"
  ),
  extract_motif(
    homer_result_path("treg_promoter_200/knownResults.txt"),
    "^AP-1\\(bZIP\\)",
    "AP-1\n(Treg)",
    "Treg-accessible",
    "Promoter"
  ),
  extract_motif(
    homer_result_path("treg_nonpromoter_200/knownResults.txt"),
    "^BATF\\(bZIP\\)",
    "BATF\n(Treg)",
    "Treg-accessible",
    "Non-promoter"
  ),
  extract_motif(
    homer_result_path("treg_promoter_200/knownResults.txt"),
    "^BATF\\(bZIP\\)",
    "BATF\n(Treg)",
    "Treg-accessible",
    "Promoter"
  ),
  extract_motif(
    homer_result_path("naive_nonpromoter_200/knownResults.txt"),
    "^LEF1\\(HMG\\)",
    "LEF1\n(Naive)",
    "Naive-accessible",
    "Non-promoter"
  ),
  extract_motif(
    homer_result_path("naive_promoter_200/knownResults.txt"),
    "^LEF1\\(HMG\\)",
    "LEF1\n(Naive)",
    "Naive-accessible",
    "Promoter"
  ),
  extract_motif(
    homer_result_path("naive_nonpromoter_200/knownResults.txt"),
    "^Tcf7\\(HMG\\)",
    "TCF7\n(Naive)",
    "Naive-accessible",
    "Non-promoter"
  ),
  extract_motif(
    homer_result_path("naive_promoter_200/knownResults.txt"),
    "^Tcf7\\(HMG\\)",
    "TCF7\n(Naive)",
    "Naive-accessible",
    "Promoter"
  ),
  extract_motif(
    homer_result_path("naive_nonpromoter_200/knownResults.txt"),
    "^RUNX1\\(Runt\\)",
    "RUNX1\n(Naive)",
    "Naive-accessible",
    "Non-promoter"
  ),
  extract_motif(
    homer_result_path("naive_promoter_200/knownResults.txt"),
    "^RUNX1\\(Runt\\)",
    "RUNX1\n(Naive)",
    "Naive-accessible",
    "Promoter"
  )
)

region_motifs$Motif <- factor(
  region_motifs$Motif,
  levels = c("AP-1\n(Treg)", "BATF\n(Treg)", "LEF1\n(Naive)", "TCF7\n(Naive)", "RUNX1\n(Naive)")
)
region_motifs$Context <- factor(region_motifs$Context, levels = c("Non-promoter", "Promoter"))

write.table(
  region_motifs,
  result_path("HOMER_region_motif_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_region_motifs <- ggplot(region_motifs, aes(x = Motif, y = Fold_enrichment, fill = Context)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.68) +
  geom_text(
    aes(label = sprintf("%.2f", Fold_enrichment)),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    size = 3.4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("Non-promoter" = "#B39DDB", "Promoter" = "#5A189A")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Motif Enrichment Persists Outside Promoters",
    subtitle = "Known motifs using fixed 200-bp windows and region-matched backgrounds",
    x = NULL,
    y = "Fold enrichment",
    fill = "Genomic context"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  plot_path("region_motifs.png"),
  plot = p_region_motifs,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# The supplied BATF/LEF1 motif-hit tables describe the saved analysis only.
# Rebuilt inputs require new motif scans before this optional gene follow-up.
if (!use_saved_input && !all(file.exists(file.path(run_dir, "motif_peak_analysis", c(
  "treg_BATF_motif_nearby_genes.txt", "naive_LEF1_motif_nearby_genes.txt"
))))) {
  message("Motif plots complete. Skipping BATF/LEF1 gene follow-up: new motif-hit tables are required for rebuilt inputs.")
  capture.output(sessionInfo(), file = result_path("HOMER_followup_R_session_info.txt"))
  quit(save = "no", status = 0)
}

batf_peaks <- read.delim(
  motif_peak_path("treg_BATF_motif_nearby_genes.txt"),
  header = TRUE,
  check.names = FALSE
)
lef1_peaks <- read.delim(
  motif_peak_path("naive_LEF1_motif_nearby_genes.txt"),
  header = TRUE,
  check.names = FALSE
)

count_nearby_genes <- function(peak_table) {
  gene_counts <- as.data.frame(sort(table(peak_table$NearestGene), decreasing = TRUE))
  colnames(gene_counts) <- c("Gene", "NumberOfPeaks")
  gene_counts[!is.na(gene_counts$Gene) & gene_counts$Gene != "", ]
}

batf_gene_counts <- count_nearby_genes(batf_peaks)
lef1_gene_counts <- count_nearby_genes(lef1_peaks)
candidate_symbols <- c("FOXP3", "IL2RA", "CTLA4", "TCF7", "LEF1", "CCR7")
batf_candidate_peaks <- batf_peaks[batf_peaks$NearestGene %in% candidate_symbols, ]
lef1_candidate_peaks <- lef1_peaks[lef1_peaks$NearestGene %in% candidate_symbols, ]

write.table(
  batf_gene_counts,
  result_path("BATF_positive_Treg_nearby_gene_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  lef1_gene_counts,
  result_path("LEF1_positive_naive_nearby_gene_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  batf_candidate_peaks,
  result_path("BATF_Treg_candidate_gene_peaks.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  lef1_candidate_peaks,
  result_path("LEF1_naive_candidate_gene_peaks.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

motif_peak_summary <- data.frame(
  Motif = c("BATF in Treg DARs", "LEF1 in naive DARs"),
  Motif_positive_peaks = c(nrow(batf_peaks), nrow(lef1_peaks)),
  Unique_nearest_genes = c(nrow(batf_gene_counts), nrow(lef1_gene_counts))
)

write.table(motif_peak_summary, result_path("HOMER_motif_peak_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = result_path("HOMER_followup_R_session_info.txt"))
