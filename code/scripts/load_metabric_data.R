# scripts/load_metabric_data.R
#
# Download and cache clinical data for the METABRIC breast cancer cohort
# from cBioPortal. The data are stored locally to ensure reproducibility
# and to avoid repeated downloads.
#
# Primary data source:
#   https://www.cbioportal.org/study/summary?id=brca_metabric
#
# The function returns a cleaned data.frame suitable for downstream
# survival analyses.

options(encoding = "UTF-8")

# ---------------------------------------------------------------
# Package setup
# ---------------------------------------------------------------
ensure_packages <- function(pkgs, lib = .libPaths()[1]) {
  ip <- rownames(installed.packages(lib.loc = .libPaths()))
  miss <- setdiff(pkgs, ip)
  if (length(miss)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages(miss, lib = lib, dependencies = TRUE)
  }
}

ensure_packages(c("cgdsr", "dplyr", "janitor"))

suppressPackageStartupMessages({
  library(cgdsr)
  library(dplyr)
  library(janitor)
})

# ---------------------------------------------------------------
# Paths
# ---------------------------------------------------------------
DATA_DIR  <- file.path("data_raw", "cbioportal")
CACHE_RDS <- file.path(DATA_DIR, "metabric_clinical.rds")

if (!dir.exists(DATA_DIR)) {
  dir.create(DATA_DIR, recursive = TRUE)
}

# ---------------------------------------------------------------
# Main loader
# ---------------------------------------------------------------
load_metabric_data <- function(force = FALSE) {
  
  # -----------------------------------------------------------
  # Use cached version if available
  # -----------------------------------------------------------
  if (file.exists(CACHE_RDS) && !isTRUE(force)) {
    message("[METABRIC] Using cached data: ", CACHE_RDS)
    return(readRDS(CACHE_RDS))
  }
  
  message("[METABRIC] Downloading clinical data from cBioPortal …")
  message("[METABRIC] This may take several minutes.")
  
  # -----------------------------------------------------------
  # Connect to cBioPortal
  # -----------------------------------------------------------
  mycgds <- CGDS("https://www.cbioportal.org/")
  
  study_id <- "brca_metabric"
  case_list_id <- "brca_metabric_all"
  
  # -----------------------------------------------------------
  # Fetch clinical data
  # -----------------------------------------------------------
  clin_raw <- try(
    getClinicalData(mycgds, case_list_id),
    silent = TRUE
  )
  
  if (inherits(clin_raw, "try-error") || is.null(clin_raw)) {
    stop(
      "Failed to download METABRIC data from cBioPortal.\n",
      "If this persists, use the manual download option described in the README."
    )
  }
  
  # -----------------------------------------------------------
  # Basic cleaning / normalization
  # -----------------------------------------------------------
  clin <- clin_raw %>%
    as.data.frame() %>%
    janitor::clean_names()
  
  # Keep patient identifier explicitly
  if (!"patient_id" %in% names(clin)) {
    clin$patient_id <- rownames(clin)
  }
  
  # -----------------------------------------------------------
  # Cache result
  # -----------------------------------------------------------
  saveRDS(clin, CACHE_RDS)
  message("[METABRIC] Data cached at: ", CACHE_RDS)
  
  clin
}

# ---------------------------------------------------------------
# Example usage (interactive)
# ---------------------------------------------------------------
# metabric_df <- load_metabric_data()
# metabric_df <- load_metabric_data(force = TRUE)
