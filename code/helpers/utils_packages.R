options(encoding = "UTF-8")

# Package setup used by run scripts.
# This version:
# - prefers user library (R_LIBS_USER)
# - can install missing packages if desired

setup_user_library <- function() {
  ul <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(ul)) {
    ul <- file.path(path.expand("~"), "Rlibs")
    Sys.setenv(R_LIBS_USER = ul)
  }
  if (!dir.exists(ul)) dir.create(ul, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(ul, .libPaths()))
  invisible(ul)
}

ensure_packages <- function(pkgs, lib = .libPaths()[1], install_if_missing = TRUE) {
  ip <- rownames(installed.packages(lib.loc = .libPaths()))
  miss <- setdiff(pkgs, ip)
  
  if (length(miss) && isTRUE(install_if_missing)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    message("Installing missing packages to: ", lib)
    install.packages(
      miss,
      lib = lib,
      dependencies = TRUE,
      Ncpus = 1,
      INSTALL_opts = c("--no-lock")
    )
  }
  
  invisible(miss)
}
