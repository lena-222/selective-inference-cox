## run_analyze_toy.R
options(encoding = "UTF-8")

## Load helper modules (project utilities, simulation, IO, lambda rules, etc.)
source("helpers/utils_misc.R")
source("helpers/utils_packages.R")
source("helpers/ci_metrics.R")
source("helpers/simulation_generic.R")
source("helpers/meta_io.R")
source("helpers/lambda_rules.R")
source("helpers/funs.design.R")
source("helpers/funs.debiased_cox.R")
source("analyze_toy_sim.R")

## Resolve project root (based on script location / RStudio / knitr / wd)
ROOT <- get_script_dir()

## Install/setup packages in user library (useful on clusters / restricted systems)
setup_user_library()

ensure_packages(c(
  "survival","glmnet","future","doFuture","foreach","MASS","selectiveInference",
  "rstudioapi","knitr"
), install_if_missing = TRUE)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(future)
  library(selectiveInference)
  library(doFuture)
  library(foreach)
  library(MASS)
})

RNGkind("Mersenne-Twister")

## Command-line arguments: all parameters can be overridden via --key=value
args <- commandArgs(trailingOnly = TRUE)

## Simulation parameters
beta_type    <- get_arg("beta_type", "realistic")   # allones | highcontrast | realistic | sparse
n            <- get_arg("n", 70)
p            <- get_arg("p", 10)
n_sim        <- get_arg("n_sim", 10)
dist_name    <- get_arg("dist", "weibull")
correlated   <- get_arg("correlated", FALSE)
rho          <- get_arg("rho", 0.0)
bin_cov      <- get_arg("bin_cov", FALSE)
target_cens  <- get_arg("target_censoring", 0.10)

## Inference / tuning parameters
alpha_fit     <- get_arg("alpha", 0.10)     # nominal alpha, CI level = 1 - alpha
lambda_choice <- get_arg("lambda_choice", "fix")
nfolds        <- get_arg("nfolds", 5)
alpha_glm     <- get_arg("alpha_glm", 1)
gamma         <- get_arg("gamma", 1)
eps           <- get_arg("eps", 1e-6)
seed0         <- get_arg("seed", 1000)
splitratio    <- get_arg("splitratio", 0.5)
support_tol   <- get_arg("support_tol", 0)
robust_cox    <- get_arg("robust_cox", FALSE)
ties          <- get_arg("ties", "efron")

## Method flags (toggle components inside analyze_surv_dataset)
full_flag      <- get_arg("full", TRUE)
oracle_flag    <- get_arg("oracle", TRUE)
split_flag     <- get_arg("split", TRUE)
cox_flag       <- get_arg("cox", TRUE)
cox0_flag      <- get_arg("cox0", TRUE)
exact_posi_flag       <- get_arg("exact_posi", TRUE)
debiased_flag  <- get_arg("debiased", TRUE)
## Debiased Cox tuning defaults (used inside analyze_surv_dataset)
debiased_lambda     <- get_arg("debiased_lambda", "lambda.min")  # "lambda.min" or "lambda.1se"
debiased_node_lambda<- get_arg("debiased_node_lambda", NA)       # numeric or NA
adaptive       <- get_arg("adaptive", FALSE)
ic_logn        <- get_arg("ic_logn", "events")
methods_arg    <- get_arg("methods", "auto")

## Helper: read lambda_fix from CSV based on simulation settings
## dist should match the CSV column values (e.g., weibull | exponential)
get_lambda_fix <- function(
    file  = "config/lambda_values_used_in_plots_n100000.csv",
    beta_type = "realistic",
    dist  = "weibull",
    p     = 10,
    rho   = 0.0,
    target_censoring = 0.10,
    use   = c("mean", "median", "min")
) {
  use <- match.arg(use)
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  ## Robust numeric parsing (handles commas and spaces)
  numize <- function(x) as.numeric(gsub(",", ".", gsub("\\s", "", x)))
  for (nm in c("rho","target_censoring","p","mean_lambda","median_lambda","sd_lambda")) {
    if (nm %in% names(df) && is.character(df[[nm]])) {
      df[[nm]] <- if (nm == "p") as.integer(df[[nm]]) else numize(df[[nm]])
    }
  }
  
  sub <- df[
    df$beta_type == beta_type &
      df$dist == dist &
      df$p == p &
      df$rho == rho &
      df$target_censoring == target_censoring,
    , drop = FALSE
  ]
  
  if (nrow(sub) == 0) stop(sprintf(
    "No CSV row found for beta_type=%s, dist=%s, p=%d, rho=%.3f, target_censoring=%.2f.",
    beta_type, dist, p, rho, target_censoring
  ))
  
  mean_val   <- mean(as.numeric(sub$mean_lambda),   na.rm = TRUE)
  median_val <- mean(as.numeric(sub$median_lambda), na.rm = TRUE)
  
  lambda <- switch(
    use,
    mean   = mean_val,
    median = median_val,
    min    = min(mean_val, median_val, na.rm = TRUE)
  )
  
  lambda <- as.numeric(lambda)[1]
  if (!is.finite(lambda) || lambda <= 0) stop("lambda_fix from CSV is not positive/finite.")
  lambda
}

if (identical(tolower(lambda_choice), "fix")) {
  if (!is.finite(lambda_fix)) {
    f_csv <- file.path(ROOT, "config/lambda_values_used_in_plots_n100000.csv")
    stopifnot(file.exists(f_csv))
    
    lambda_fix <- get_lambda_fix(
      file = f_csv,
      beta_type = beta_type,
      dist =  dist_name,
      p = p,
      rho =  rho,
      target_censoring = target_cens,
      use = "mean"#,
      #match_p = !grepl("^metabric_", sim_design)
    )
  }
}

cat(sprintf(
  "RUN: beta_type=%s n=%d p=%d dist=%s rho=%.2f target_cens=%.2f lambda_choice=%s lambda_fix=%s\n",
  beta_type, n, p, dist_name, rho, target_cens, lambda_choice,
  if (is.na(lambda_fix)) "NA" else sprintf("%.10f", lambda_fix)
))

## Beta generator for toy simulations
make_beta <- function(beta_type, p) {
  stopifnot(p >= 1)
  
  base_vals <- switch(
    beta_type,
    "allones"      = rep(1, p),
    "highcontrast" = c(0.3, 1.0, 0.3, 1.0, rep(0, p - 4)),
    "realistic"    = c(0.8, 0.7, 0.5, 0.8, rep(0, p - 4)),
    "sparse"       = { out <- rep(0, p); out[1:min(2, p)] <- 1; out },
    c(0.8, 0.7, 0.5, 0.8, rep(0, p - 4))
  )
  
  beta <- base_vals[1:p]
  names(beta) <- paste0("X", seq_len(p))
  beta
}

## Ground-truth coefficients (used for oracle + coverage/power metrics)
beta_true <- make_beta(beta_type, p)
if (length(beta_true) != p) {
  stop(sprintf("beta_true length (%d) != p (%d). Adjust beta_true or p.", length(beta_true), p))
}
names(beta_true) <- paste0("X", seq_along(beta_true))

## Map method identifiers to result-table names inside analyze_surv_dataset output
.res_name_map <- list(
  full     = "res_full",
  oracle   = "res_oracle",
  refit    = "res_refit",
  refit0   = "res_refit0",
  split    = "res_split",
  exact_posi      = "res_exact_posi",
  debiased = "res_debiased"
)

## Run simulation and analysis
set.seed(seed0)

if (identical(methods_arg, "auto")) {
  methods <- c("split", "full", "oracle", "refit", "refit0", "exact_posi", "debiased")
} else {
  methods <- strsplit(methods_arg, ",")[[1]]
  methods <- trimws(methods)
}

plan(multisession, workers = parallelly::availableCores())

registerDoFuture()

outs <- foreach(i = 1:n_sim, .packages = c("survival", "glmnet")) %dopar% {
  S <- simulate_fun(
    n = n, beta = beta_true,
    dist = dist_name,
    correlated = correlated, rho = rho,
    bin_cov = bin_cov, target_censoring = target_cens
  )
  
  res <- tryCatch(
    analyze_surv_dataset(
      obj = S, name = paste0("sim", i),
      beta_true = beta_true,
      full   = full_flag,
      oracle = oracle_flag,
      split  = split_flag,
      cox    = cox_flag,
      cox0   = cox0_flag,
      exact_posi    = exact_posi_flag,
      debiased        = debiased_flag,
      debiased_lambda = debiased_lambda,
      debiased_nodewise_lambda =
        if (is.na(debiased_node_lambda)) NULL else as.numeric(debiased_node_lambda),
      adaptive  = adaptive,
      nfolds    = nfolds,
      alpha_glm = alpha_glm,
      lambda_choice = lambda_choice,
      lambda_fix    = lambda_fix,
      ties      = ties,
      ic_logn   = ic_logn,
      robust_cox = robust_cox,
      splitratio = splitratio,
      support_tol = support_tol,
      gamma = gamma,
      eps   = eps,
      seed  = seed0 + i,
      alpha = alpha_fit
    ),
    error = function(e) {
      message("Error in sim ", i, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(res)) res$beta_true <- beta_true
  res
}

outs <- Filter(Negate(is.null), outs)
if (!length(outs)) {
  warning(sprintf(
    "All analyze_surv_dataset() runs returned NULL; skipping config beta_type=%s, n=%d, p=%d, dist=%s, lambda_choice=%s, target_cens=%.2f, rho=%.2f.",
    beta_type, n, p, dist_name, lambda_choice, target_cens, rho
  ))
  q("no", status = 0)
}

perf_beta <- perf_per_beta_split(outs, methods = methods, alpha_for_ci = alpha_fit)

## Save results
datetag <- format(Sys.Date(), "%d-%m-%Y")
outdir  <- file.path(ROOT, "results", "new", datetag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

timetag <- format(Sys.time(), "%H-%M")

tag_base <- sprintf(
  "new-%s_n%dp%d_sim%d_dist-%s_alpha-%.2f_lambda-%s",
  beta_type, n, p, n_sim, dist_name, alpha_fit, lambda_choice
)
tag <- paste0(tag_base, "_", timetag)

file_perf <- file.path(outdir, paste0(tag, "_perf_beta.csv"))
file_outs <- file.path(outdir, paste0(tag, "_outs.rds"))
file_meta <- file.path(outdir, paste0(tag, "_meta.txt"))

write.csv(perf_beta, file_perf, row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(outs, file_outs)

meta <- list(
  timestamp = as.character(Sys.time()),
  n = n, p = p, n_sim = n_sim,
  dist = dist_name,
  correlated = correlated, rho = rho,
  bin_cov = bin_cov,
  target_censoring = target_cens,
  adaptive = adaptive,
  alpha = alpha_fit,
  lambda_choice = lambda_choice,
  beta_type = beta_type,
  splitratio = splitratio,
  exact_posi = exact_posi_flag,
  debiased = debiased_flag,
  debiased_lambda = debiased_lambda,
  debiased_nodewise_lambda = debiased_node_lambda,
  methods = paste(methods, collapse = ", "),
  seed0 = seed0,
  ties = ties,
  lambda_fix = if (is.na(lambda_fix)) NA_real_ else lambda_fix
)

capture.output(str(meta), file = file_meta)
capture.output(sessionInfo(), file = file_meta, append = TRUE)

cat("Saved:\n  - ", file_perf, "\n  - ", file_outs, "\n  - ", file_meta, "\n", sep = "")

