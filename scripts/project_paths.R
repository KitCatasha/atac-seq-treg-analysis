# Run the analysis scripts from the repository root.
if (!file.exists("input/samples.tsv")) {
  stop("Run Rscript from the repository root (the folder containing README.md).")
}
run_dir <- Sys.getenv("ATAC_RUN_DIR", unset = "runs/local")
if (dir.exists(run_dir)) {
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
} else {
  parent <- dirname(run_dir)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(run_dir))
}
repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
allowed_root <- paste0(file.path(repo_root, "runs"), "/")
if (!startsWith(paste0(run_dir, "/"), allowed_root)) {
  stop("ATAC_RUN_DIR must point to a subfolder of this repository's runs/ directory.")
}
for (subdir in c("results", "plots", "homer_input", "homer_results")) {
  dir.create(file.path(run_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}
analysis_input <- Sys.getenv("ATAC_INPUT_DIR", unset = "input")
use_saved_input <- identical(
  normalizePath(analysis_input, winslash = "/", mustWork = TRUE),
  normalizePath("input", winslash = "/", mustWork = TRUE)
)
input_path <- function(name) file.path(analysis_input, name)
result_path <- function(name) file.path(run_dir, "results", name)
plot_path <- function(name) file.path(run_dir, "plots", name)
homer_input_path <- function(name) file.path(run_dir, "homer_input", name)
read_result_path <- function(name) {
  rerun <- result_path(name)
  saved <- file.path("results", "tables", name)
  selected <- if (file.exists(rerun) || !use_saved_input) rerun else saved
  if (!file.exists(selected)) stop("Required result is missing: ", selected)
  message("Reading analysis table: ", selected)
  selected
}
homer_result_path <- function(name) {
  rerun <- file.path(run_dir, "homer_results", name)
  saved <- file.path("results", "motifs", name)
  selected <- if (file.exists(rerun) || !use_saved_input) rerun else saved
  if (!file.exists(selected)) stop("Required HOMER result is missing: ", selected)
  message("Reading HOMER result: ", selected)
  selected
}
motif_peak_path <- function(name) {
  selected <- if (use_saved_input) file.path("results", "motif_peak_analysis", name) else
    file.path(run_dir, "motif_peak_analysis", name)
  if (!file.exists(selected)) stop("Required motif-peak table is missing: ", selected)
  selected
}
message("New outputs will be written under ", run_dir)
