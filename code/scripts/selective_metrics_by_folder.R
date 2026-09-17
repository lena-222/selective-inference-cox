#!/usr/bin/env Rscript
# scripts/selective_metrics_by_folder.R
# Discover runs under root_dir and compute (A) CI error metrics and (B) selection summaries.
# Both are merged into one per-run CSV row set.

options(encoding = "UTF-8")

ul <- Sys.getenv("R_LIBS_USER")
if (!nzchar(ul)) {
  ul <- file.path(path.expand("~"), "Rlibs")
  Sys.setenv(R_LIBS_USER = ul)
}
if (!dir.exists(ul)) dir.create(ul, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(ul, .libPaths()))

ensure_packages <- function(pkgs, lib = .libPaths()[1]) {
  ip <- rownames(installed.packages(lib.loc = .libPaths()))
  miss <- setdiff(pkgs, ip)
  if (length(miss)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    message("Installing packages: ", paste(miss, collapse = ", "))
    install.packages(miss, lib = lib, dependencies = TRUE,
                     Ncpus = 1, INSTALL_opts = c("--no-lock"))
  }
}

ensure_packages(c("data.table", "dplyr"), lib = .libPaths()[1])

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("--file=", args)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", args[m[1]]))))
  normalizePath(getwd())
}
ROOT <- get_script_dir()

src <- function(f) source(file.path(ROOT, f), chdir = TRUE)
src("helpers/meta_io.R")
src("helpers/ci_metrics.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

DEFAULT_ROOT_DIR <- file.path(ROOT, "results", "_processed")
args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1L && nzchar(args[1])) args[1] else DEFAULT_ROOT_DIR

OUT_DIR <- file.path(ROOT, "selective_metrics_by_folder")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

folder_name <- basename(normalizePath(root_dir))
OUT_FILE <- file.path(OUT_DIR, paste0("selective_metrics_", folder_name, ".csv"))
DIAG_LOG <- file.path(OUT_DIR, "selective_diagnostics_by_folder.log")

log_diag <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", paste(..., collapse = " "))
  message(msg)
  cat(msg, file = DIAG_LOG, sep = "\n", append = TRUE)
}

log_diag("Start selective metrics for root_dir =", root_dir)

if (file.exists(OUT_FILE)) {
  log_diag("Removing existing output file:", OUT_FILE)
  file.remove(OUT_FILE)
}

meta_files <- list.files(root_dir, pattern = "meta\\.txt$", recursive = TRUE, full.names = TRUE)
meta_files <- meta_files[!grepl("/[^/]*meta[^/]*/", meta_files, ignore.case = TRUE)]

log_diag("Found meta files:", length(meta_files))
if (!length(meta_files)) quit(status = 0)

all_runs <- list()

for (mpath in meta_files) {
  cat("==> meta:", mpath, "\n")
  
  meta_row <- parse_meta_file(mpath)
  if (is.null(meta_row)) next
  
  run_dir <- meta_row$run_dir[1]
  rds_files <- list.files(run_dir, pattern = "_outs\\.rds$", full.names = TRUE)
  if (!length(rds_files)) next
  
  outs_path <- rds_files[1]
  outs_raw <- try(readRDS(outs_path), silent = TRUE)
  if (inherits(outs_raw, "try-error") || is.null(outs_raw)) next
  
  outs_list <- Filter(Negate(is.null), outs_raw)
  if (!length(outs_list)) next
  
  n_sims_total  <- length(outs_raw)
  n_sims_ok     <- length(outs_list)
  n_sims_failed <- n_sims_total - n_sims_ok
  
  methods_vec <- attr(meta_row, "methods_vec") %||%
    if ("methods" %in% names(meta_row)) strsplit(as.character(meta_row$methods[1]), "\\s*,\\s*")[[1]] else NULL
  methods_vec <- methods_vec[nzchar(methods_vec)]
  if (!length(methods_vec)) methods_vec <- names(.res_name_map)
  
  alpha_ci <- meta_row$alpha[1] %||% 0.10
  
  ci_err <- try(
    compute_ci_error_metrics(outs_list, methods = methods_vec, alpha_for_ci = alpha_ci),
    silent = TRUE
  )
  if (inherits(ci_err, "try-error") || is.null(ci_err) || !nrow(ci_err)) next
  
  sel_sum <- try(
    compute_selection_summary_metrics(outs_list, methods = methods_vec, alpha_for_ci = alpha_ci),
    silent = TRUE
  )
  if (inherits(sel_sum, "try-error") || is.null(sel_sum) || !nrow(sel_sum)) {
    sel_sum <- data.frame()
  }
  
  # Merge both metric blocks on (method, term, j, beta_true, is_signal) when possible
  metrics <- dplyr::left_join(
    ci_err,
    sel_sum,
    by = intersect(names(ci_err), names(sel_sum))
  )
  
  run_id <- sub("_outs\\.rds$", "", basename(outs_path))
  
  keep_meta_cols <- intersect(
    c("timestamp","n","p","n_sim","dist","beta_type",
      "lambda_choice","lambda_fix","target_censoring","rho",
      "correlated","bin_cov","adaptive","sim_design"),
    names(meta_row)
  )
  meta_small <- meta_row[rep(1, nrow(metrics)), keep_meta_cols, drop = FALSE]
  
  metrics_run <- cbind(
    data.frame(
      run_id        = run_id,
      run_dir       = run_dir,
      meta_path     = mpath,
      n_sims_total  = n_sims_total,
      n_sims_ok     = n_sims_ok,
      n_sims_failed = n_sims_failed,
      stringsAsFactors = FALSE
    ),
    meta_small,
    metrics
  )
  
  all_runs[[length(all_runs) + 1L]] <- metrics_run
}

if (!length(all_runs)) quit(status = 0)

final_dt <- dplyr::bind_rows(all_runs)
data.table::fwrite(final_dt, OUT_FILE)

log_diag("Done. Wrote:", OUT_FILE, "rows:", nrow(final_dt))
cat("Output:", OUT_FILE, "\n")
