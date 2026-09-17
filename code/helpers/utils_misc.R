options(encoding = "UTF-8")

# helpers/utils_misc.R


# Small utilities used across run scripts and analysis cores.

# Return directory of the currently running script in various environments:
# - Rscript --file=...
# - RStudio (active document)
# - knitr
# Fallback: getwd()
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("--file=", args)
  if (length(m)) {
    return(dirname(normalizePath(sub("^--file=", "", args[m[1]]))))
  }

  if (nzchar(Sys.getenv("RSTUDIO", ""))) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }

  if (isTRUE(getOption("knitr.in.progress"))) {
    p <- tryCatch(knitr::current_input(), error = function(e) NULL)
    if (!is.null(p)) return(dirname(normalizePath(p)))
  }

  normalizePath(getwd())
}

# Parse CLI arguments of the form --key=value.
# Uses type.convert() so numbers and TRUE/FALSE become proper types.
get_arg <- function(key, default, args = commandArgs(trailingOnly = TRUE)) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  val <- sub(paste0("^--", key, "="), "", hit[1])
  type.convert(val, as.is = TRUE)
}

# Safe division with NA on invalid inputs.
safe_div <- function(a, b) {
  ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)
}

# Small helper to select from a named vector by index set j.
pick_named <- function(named_vec, jset) {
  if (!length(named_vec)) return(rep(NA_real_, length(jset)))
  out <- named_vec[as.character(jset)]
  out[is.na(out)] <- NA_real_
  unname(out)
}

# Aggregate binary outcomes success_vec by group ids in id_vec:
# returns k = sum(success), n = number of finite entries.
agg_bin <- function(success_vec, id_vec) {
  if (!length(success_vec)) return(list(k = numeric(0), n = numeric(0)))
  k <- tapply(success_vec, id_vec, sum, na.rm = TRUE)
  n <- tapply(success_vec, id_vec, function(x) sum(is.finite(x)))
  if (length(k) == 0L) { k <- numeric(0); n <- numeric(0) }
  list(k = k, n = n)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

make_surv_formula <- function(Xs) {
  if (is.null(Xs) || ncol(Xs) == 0) return(NULL)
  rhs <- paste0(colnames(Xs), collapse = " + ")
  as.formula(paste("Surv(time, status) ~", rhs))
}
# Small helper to check whether exact_posi output seems usable.
.extract_kkt_flag <- function(exact_posi_obj) {
  if (is.null(exact_posi_obj) || inherits(exact_posi_obj, "try-error")) return(0L)
  kk <- try(exact_posi_obj$kkt_ok, silent = TRUE)
  if (!inherits(kk, "try-error") && length(kk) == 1) return(as.integer(isTRUE(kk)))
  pv <- try(exact_posi_obj$pv, silent = TRUE)
  if (!inherits(pv, "try-error") && is.numeric(pv) && length(pv) > 0) return(as.integer(all(is.finite(pv))))
  1L
}

# Create an empty result table aligned to all_coef.
.init_res_df <- function(all_coef) {
  out <- data.frame(
    var     = all_coef,
    bhat    = NA_real_,
    se      = NA_real_,
    lower   = NA_real_,
    upper   = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )
  rownames(out) <- all_coef
  out
}