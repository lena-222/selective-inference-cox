# scripts/00_setup_packages.R
# Installs required CRAN + Bioconductor packages .

install_if_missing <- function(pkgs) {
  pkgs <- unique(pkgs)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install.packages(missing, dependencies = TRUE)
  }
  invisible(missing)
}

# CRAN packages
cran_pkgs <- c(
  "dplyr", "janitor", "tidyr", "survival", "glmnet",
  "ggplot2", "scales", "patchwork", "future", "doFuture",
  "foreach", "MASS", "selectiveInference", "flexsurv", "doRNG",
  "data.table", "rstudioapi", "knitr"
)
install_if_missing(cran_pkgs)

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c("cBioPortalData")
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

message("Package setup complete.")
