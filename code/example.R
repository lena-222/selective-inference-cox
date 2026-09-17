# Source local functions
source("analyze_all_meta.R")

required <- c("dplyr","janitor","survival","glmnet","ggplot2","scales","patchwork","cBioPortalData")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nRun: source('scripts/00_setup_packages.R')")
}

# Load packages
library(dplyr)
library(janitor)
library(survival)
library(glmnet)
library(ggplot2)
library(scales)

# Create output folder
dir.create("example", showWarnings = FALSE, recursive = TRUE)

# RNG setup (parallel-safe)
RNGkind("L'Ecuyer-CMRG")

# Load and clean METABRIC
data_fp <- Sys.getenv(
  "METABRIC_CSV",
  unset = file.path("data_raw", "cbioportal", "Breast Cancer METABRIC.csv")
)
if (!file.exists(data_fp)) {
  stop(
    "METABRIC CSV not found: ", data_fp, "\n",
    "Set METABRIC_CSV or place the downloaded clinical CSV under data_raw/cbioportal/."
  )
}
metabric_raw <- read.csv(data_fp, stringsAsFactors = TRUE, check.names = FALSE) |>
  janitor::clean_names()

metabric_os <- metabric_raw |>
  dplyr::filter(!is.na(overall_survival_months)) |>
  dplyr::mutate(
    os_time  = overall_survival_months,
    os_event = ifelse(grepl("^Die", patients_vital_status), 1L, 0L),
    tumor_stage_f = factor(tumor_stage, levels = sort(unique(tumor_stage))),
    er_status_clean = dplyr::na_if(er_status, ""),
    er_status_f = dplyr::case_when(
      er_status_clean == "Positive" ~ "ER_pos",
      er_status_clean == "Negative" ~ "ER_neg",
      TRUE ~ NA_character_
    ) |> factor(),
    her2_status_f = dplyr::case_when(
      her2_status %in% c("Positive", "Pos") ~ "HER2_pos",
      her2_status %in% c("Negative", "Neg") ~ "HER2_neg",
      TRUE ~ NA_character_
    ) |> factor()
  ) |>
  dplyr::filter(os_time > 0)

# Define covariates used for the design matrix
metabric_covariates <- c(
  "age_at_diagnosis",
  "tumor_stage_f",
  "er_status_f",
  "her2_status_f",
  "pr_status",
  "chemotherapy",
  "hormone_therapy",
  "radio_therapy",
  "neoplasm_histologic_grade",
  "tumor_size",
  "lymph_nodes_examined_positive",
  "nottingham_prognostic_index"
)

X_df <- metabric_os |>
  dplyr::select(dplyr::all_of(metabric_covariates))

idx_complete <- complete.cases(X_df, metabric_os$os_time, metabric_os$os_event)

metabric_obj <- list(
  time   = metabric_os$os_time[idx_complete],
  status = metabric_os$os_event[idx_complete],
  X      = X_df[idx_complete, , drop = FALSE]
)

# Build a dummy-coded design matrix without intercept
get_design_matrix <- function(X_df) {
  stats::model.matrix(~ . , data = X_df)[, -1, drop = FALSE]
}

# Choose lambda by AIC/BIC along glmnet path using Cox partial log-likelihood with offset
choose_lambda_ic_cox <- function(time, status, X,
                                 alpha = 1,
                                 standardize = FALSE,
                                 ic = c("aic","bic"),
                                 ties = "efron",
                                 max_warn = 5) {
  ic <- match.arg(ic)
  y <- survival::Surv(time, status)
  
  fit <- glmnet::glmnet(
    x = X, y = y,
    family = "cox",
    alpha = alpha,
    standardize = standardize
  )
  
  eta_mat <- stats::predict(fit, newx = X, type = "link")
  if (is.vector(eta_mat)) eta_mat <- matrix(eta_mat, ncol = 1)
  
  n <- length(time)
  L <- ncol(eta_mat)
  df <- fit$df
  if (length(df) != L) df <- rep(df[1], L)
  
  loglik <- rep(NA_real_, L)
  warn_ct <- 0L
  
  for (j in seq_len(L)) {
    eta <- as.numeric(eta_mat[, j])
    if (any(!is.finite(eta))) next
    if (sd(eta) == 0) next
    
    m <- try(survival::coxph(y ~ stats::offset(eta), ties = ties), silent = TRUE)
    if (inherits(m, "try-error")) {
      warn_ct <- warn_ct + 1L
      if (warn_ct <= max_warn) message("AIC/BIC: coxph failed at lambda index ", j, " (skipping).")
      next
    }
    loglik[j] <- m$loglik[2]
  }
  
  ok <- is.finite(loglik) & is.finite(df)
  if (!any(ok)) stop("AIC/BIC: no valid log-likelihood values were computed.")
  
  crit <- rep(Inf, L)
  if (ic == "aic") crit[ok] <- -2 * loglik[ok] + 2 * df[ok]
  if (ic == "bic") crit[ok] <- -2 * loglik[ok] + log(n) * df[ok]
  
  j_best <- which.min(crit)
  list(lambda = fit$lambda[j_best], fit = fit, ic = ic, crit = crit)
}

# Select support using glmnet on the full dataset under a lambda rule
glmnet_support <- function(time, status, X_df,
                           lambda_choice = c("min","1se","aic","bic"),
                           nfolds = 10,
                           seed = NULL) {
  lambda_choice <- match.arg(lambda_choice)
  if (!is.null(seed)) set.seed(seed)
  
  X <- get_design_matrix(X_df)
  y <- survival::Surv(time, status)
  
  if (lambda_choice %in% c("min","1se")) {
    cvfit <- glmnet::cv.glmnet(
      x = X, y = y,
      family = "cox",
      alpha = 1,
      nfolds = nfolds,
      standardize = FALSE
    )
    s <- if (lambda_choice == "min") "lambda.min" else "lambda.1se"
    b <- as.matrix(stats::coef(cvfit, s = s))
    sel <- rownames(b)[as.numeric(b) != 0]
    return(list(sel = sel, lambda = cvfit[[s]]))
  }
  
  icfit <- choose_lambda_ic_cox(
    time = time, status = status, X = X,
    alpha = 1, standardize = FALSE,
    ic = lambda_choice
  )
  b <- as.matrix(stats::coef(icfit$fit, s = icfit$lambda))
  sel <- rownames(b)[as.numeric(b) != 0]
  list(sel = sel, lambda = icfit$lambda)
}

# Select support using glmnet on a training split
split_support <- function(time, status, X_df,
                          splitratio = 0.5,
                          lambda_choice = c("min","1se","aic","bic"),
                          nfolds = 10,
                          seed = NULL) {
  lambda_choice <- match.arg(lambda_choice)
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(time)
  idx_train <- sample.int(n, size = floor(splitratio * n), replace = FALSE)
  
  X <- get_design_matrix(X_df)
  y <- survival::Surv(time, status)
  
  Xtr <- X[idx_train, , drop = FALSE]
  ytr <- y[idx_train]
  time_tr <- time[idx_train]
  status_tr <- status[idx_train]
  
  if (lambda_choice %in% c("min","1se")) {
    cvfit <- glmnet::cv.glmnet(
      x = Xtr, y = ytr,
      family = "cox",
      alpha = 1,
      nfolds = nfolds,
      standardize = FALSE
    )
    s <- if (lambda_choice == "min") "lambda.min" else "lambda.1se"
    b <- as.matrix(stats::coef(cvfit, s = s))
    sel <- rownames(b)[as.numeric(b) != 0]
    return(list(sel = sel, idx_train = idx_train, lambda = cvfit[[s]]))
  }
  
  icfit <- choose_lambda_ic_cox(
    time = time_tr, status = status_tr, X = Xtr,
    alpha = 1, standardize = FALSE,
    ic = lambda_choice
  )
  b <- as.matrix(stats::coef(icfit$fit, s = icfit$lambda))
  sel <- rownames(b)[as.numeric(b) != 0]
  list(sel = sel, idx_train = idx_train, lambda = icfit$lambda)
}
