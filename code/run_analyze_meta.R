## run_analyze_meta.R
## End-to-end: simulieren -> analysieren -> Metriken -> speichern
options(encoding = "UTF-8")

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("--file=", args)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", args[m[1]]))))
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
ROOT <- get_script_dir()

## Pakete
ul <- Sys.getenv("R_LIBS_USER")
if (!nzchar(ul)){
  ul <- file.path(path.expand("~"), "Rlibs")
  Sys.setenv(R_LIBS_USER=ul)
}
if(!dir.exists(ul)) dir.create(ul,recursive = TRUE, showWarnings = FALSE)
.libPaths(c(ul,.libPaths()))

ensure_packages <- function(pkgs, lib = .libPaths()[1]) {
  ip <- rownames(installed.packages(lib.loc = .libPaths()))
  miss <- setdiff(pkgs, ip)
  if (length(miss)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    message("Installing to user library: ", lib)
    install.packages(
      miss, lib = lib, dependencies = TRUE,
      Ncpus = 1, INSTALL_opts = c("--no-lock")
    )
  }
}

ensure_packages(
  c(
    "survival","glmnet","future","doFuture","foreach","MASS","selectiveInference",
    "flexsurv","dplyr","janitor"
  ),
  lib = .libPaths()[1]
)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(future)
  library(selectiveInference)
  library(doFuture)
  library(foreach)
  library(MASS)
  library(flexsurv)
  library(dplyr)
  library(janitor)
})

RNGkind("Mersenne-Twister")

## Projekt-Helfer laden
src <- function(f) source(file.path(ROOT, f), chdir = TRUE)
#src("funs.inf2.R")
#src("funs.fixedCox2.R")
#src("funs.fixed.R")
#src("funs.common.R")
#src("funs.randomizedCox.R")
src("analyze_all_meta.R")     # enthält analyze_surv_dataset(), perf_per_beta_split(), ...
# Build the METABRIC Weibull-PH reference model from the public clinical data.
# Set the environment variable METABRIC_CSV or place the CSV at
# data_raw/cbioportal/Breast Cancer METABRIC.csv (see README).
src("helpers/funs.metabric_weibph.R") # METABRIC-Weibull-PH + sim_metabric_weibph(), get_metabric_beta_true()

## CLI-Args mit Defaults
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  val <- sub(paste0("^--", key, "="), "", hit[1])
  type.convert(val, as.is = TRUE)
}

beta_type    <- get_arg("beta_type",     "realistic")
n <- get_arg("n", 250)
n <- as.integer(n)
if (is.na(n) || n <= 0) {
  stop(sprintf("--n muss positive ganze Zahl sein, bekommen: %s", n))
}

p            <- get_arg("p",              10)
n_sim        <- get_arg("n_sim",          10)
dist_name    <- get_arg("dist",           "weibull")
correlated   <- get_arg("correlated",     FALSE)
rho          <- get_arg("rho",            0.0)
bin_cov      <- get_arg("bin_cov",        FALSE)
target_cens  <- get_arg("target_censoring", 0.10)
alpha_fit    <- get_arg("alpha",          0.10)
lambda_choice<- get_arg("lambda_choice",  "min")
nfolds       <- get_arg("nfolds",         5)
alpha_glm    <- get_arg("alpha_glm",      1)
gamma        <- get_arg("gamma",          1)
eps          <- get_arg("eps",            1e-6)
seed0        <- get_arg("seed",           1000)
splitratio   <- get_arg("splitratio",     0.5)
support_tol  <- get_arg("support_tol",    0)
robust_cox   <- get_arg("robust_cox",     FALSE)
ties         <- get_arg("ties",           "efron")
lambda_fix <- NA_real_

## Design-Schalter: "generic" vs "metabric_weibph" vs "metabric_weibph_truth"
sim_design   <- get_arg("sim_design",     "generic")

## Modell-Flags
## Modell-Flags
full_flag    <- get_arg("full",           TRUE)
oracle_flag  <- get_arg("oracle",         TRUE)
split_flag   <- get_arg("split",          TRUE)
cox_flag     <- get_arg("cox",            TRUE)
cox0_flag    <- get_arg("cox0",           TRUE)
fli_flag     <- get_arg("fli",            TRUE)
adaptive     <- get_arg("adaptive",       FALSE)
ic_logn      <- get_arg("ic_logn",        "events")
methods_arg  <- get_arg("methods",        "auto")

## Debiased-Optionen
debiased_flag        <- get_arg("debiased",        TRUE)

# --- Helper: lambda aus CSV ----------------------------------------------
get_lambda_fix <- function(
    file  = "lambda_values_used_in_plots_n100000.csv",
    beta_type = "realistic",
    dist  = "weibull",
    p     = 10,
    rho   = 0.0,
    target_censoring = 0.10,
    use   = c("mean", "median", "min"),
    match_p = TRUE              # <--- NEU: steuert, ob nach p gefiltert wird
) {
  use <- match.arg(use)
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  numize <- function(x) as.numeric(gsub(",", ".", gsub("\\s", "", x)))
  for (nm in c("rho","target_censoring","p","mean_lambda","median_lambda","sd_lambda")) {
    if (nm %in% names(df) && is.character(df[[nm]])) {
      df[[nm]] <- if (nm == "p") as.integer(df[[nm]]) else numize(df[[nm]])
    }
  }
  
  ## Basis-Filter ohne p
  rows <- df$beta_type == beta_type &
    df$dist      == dist &
    df$rho       == rho &
    df$target_censoring == target_censoring
  
  ## Optional p dazunehmen
  if (isTRUE(match_p)) {
    rows <- rows & df$p == p
  }
  
  sub <- df[rows, , drop = FALSE]
  
  if (nrow(sub) == 0) {
    stop(sprintf(
      "Keine Zeile für beta_type=%s, dist=%s, rho=%.3f, target_censoring=%.2f%s gefunden.",
      beta_type, dist, rho, target_censoring,
      if (match_p) sprintf(", p=%d", p) else " (ohne p-Filter)"
    ))
  }
  
  mean_val   <- mean(as.numeric(sub$mean_lambda),   na.rm = TRUE)
  median_val <- mean(as.numeric(sub$median_lambda), na.rm = TRUE)
  lambda <- switch(use,
                   mean   = mean_val,
                   median = median_val,
                   min    = min(mean_val, median_val, na.rm = TRUE))
  lambda <- as.numeric(lambda)[1]
  if (!is.finite(lambda) || lambda <= 0) stop("lambda_fix aus CSV ist nicht positiv/endlich.")
  lambda
}


## Design-Schalter: "generic" vs "metabric_weibph" vs "metabric_weibph_truth"
sim_design   <- get_arg("sim_design",     "generic")

## --- lambda_choice behandeln: nur für generic-Design fix aus CSV holen ----
if (identical(sim_design, "generic") && identical(tolower(lambda_choice), "fix")) {
  f_csv <- file.path(ROOT, "config", "lambda_values_used_in_plots_n100000.csv")
  cat("ROOT =", ROOT, "\n")
  cat(sprintf("CSV-Pfad: %s  exists=%s\n", f_csv, file.exists(f_csv)))
  stopifnot(file.exists(f_csv))
  
  lambda_fix <- get_lambda_fix(
    file = f_csv,
    beta_type = beta_type,
    dist = dist_name,
    p = p, 
    rho = rho,
    target_censoring = target_cens,
    use = "mean",
    match_p = TRUE    # hier darf p gefiltert werden, weil generic
  )
  
  stopifnot(
    is.numeric(lambda_fix),
    length(lambda_fix) == 1L,
    is.finite(lambda_fix),
    lambda_fix > 0
  )
}

## Für METABRIC-Design: lambda_choice='fix' nicht verwenden
if (sim_design %in% c("metabric_weibph", "metabric_weibph_truth") &&
    identical(tolower(lambda_choice), "fix")) {
  warning("lambda_choice='fix' wird für sim_design=metabric_* ignoriert; setze lambda_choice='negahban'.")
  lambda_choice <- "negahban"
  lambda_fix    <- NA_real_
}


cat(sprintf(
  "RUN: beta_type=%s n=%d p=%d dist=%s rho=%.2f target_cens=%.2f lambda_choice=%s lambda_fix=%s sim_design=%s\n",
  beta_type, n, p, dist_name, rho, target_cens, lambda_choice,
  if (is.na(lambda_fix)) "NA" else sprintf("%.10f", lambda_fix),
  sim_design
))

## Beta-Generator (nur für generic-Design)
make_beta <- function(beta_type, p) {
  stopifnot(p >= 1)
  base_vals <- switch(
    beta_type,
    "allones"     = rep(1, p),
    "highcontrast"= c(0.3, 1.0, 0.3, 1.0, rep(0, p-4)),
    "realistic"   = c(0.8, 0.7, 0.5, 0.8, rep(0, p-4)),
    "sparse"      = { out <- rep(0, p); out[1:min(2,p)] <- 1; out },
    { c(0.8, 0.7, 0.5, 0.8, rep(0, p-4)) }
  )
  beta <- base_vals[1:p]
  names(beta) <- paste0("X", seq_len(p))
  beta
}

## Wahr-Parameter je nach Design
if (identical(sim_design, "generic")) {
  
  beta_true <- make_beta(beta_type, p)
  if (length(beta_true) != p) {
    stop(sprintf("beta_true length (%d) != p (%d). Passe beta_true oder p an.", length(beta_true), p))
  }
  names(beta_true) <- paste0("X", seq_along(beta_true))
  
} else if (identical(sim_design, "metabric_weibph")) {
  
  ## reine METABRIC-Simulation ohne Coverage/Power
  beta_true <- NULL
  oracle_flag <- FALSE
  
} else if (identical(sim_design, "metabric_weibph_truth")) {
  
  ## METABRIC-Simulation + "wahre" Betas = log-HRs aus m_os_weibph
  beta_true <- get_metabric_beta_true()   # Namen: age10, tumor_stage_f1, ...
  oracle_flag <- TRUE
  
} else {
  stop(sprintf("Unbekannter sim_design-Wert: %s", sim_design))
}

# -------- Mapping etc. ----------------------------------------------------
.res_name_map <- list(
  full                   = "res_full",
  oracle                 = "res_oracle",
  refit                  = "res_refit",
  refit0                 = "res_refit0",
  split                  = "res_split",
  fli                    = "res_fli",
  debiased               = "res_debiased"
)

safe_div <- function(a,b) ifelse(is.finite(a) & is.finite(b) & b > 0, a/b, NA_real_)
pick_named <- function(named_vec, jset){
  if (!length(named_vec)) return(rep(NA_real_, length(jset)))
  out <- named_vec[as.character(jset)]
  out[is.na(out)] <- NA_real_
  unname(out)
}
agg_bin <- function(success_vec, id_vec){
  if (!length(success_vec)) return(list(k = numeric(0), n = numeric(0)))
  k <- tapply(success_vec, id_vec, sum, na.rm = TRUE)
  n <- tapply(success_vec, id_vec, function(x) sum(is.finite(x)))
  if (length(k) == 0L) { k <- numeric(0); n <- numeric(0) }
  list(k = k, n = n)
}

## get_ci_df_present(), normalize_ci_names(), perf_per_beta_split()
## -> kommen aus analyze_all_meta.R

## Simulation & Analyse
set.seed(seed0)

if (identical(methods_arg, "auto")) {
  methods <- c(
    "split", "full", "oracle",
    "refit","refit0","fli",
    "debiased"
  )
} else {
  methods <- strsplit(methods_arg, ",")[[1]]
  methods <- trimws(methods)
}

#plan(multisession, workers = 16)
#registerDoFuture()

### Wrapper, der je nach sim_design simuliert
simulate_one_dataset <- function(n, beta_true, dist_name,
                                 correlated, rho, bin_cov, target_cens) {
  if (sim_design %in% c("metabric_weibph", "metabric_weibph_truth")) {
    sim_metabric_weibph(n = n)
  } else {
    simulate_fun(
      n = n, beta = beta_true,
      dist = dist_name,
      correlated = correlated, rho = rho,
      bin_cov = bin_cov, target_censoring = target_cens
    )
  }
}
plan(sequential)
registerDoSEQ()
#plan(multisession, workers = 16)
#registerDoFuture()
outs <- foreach(
  i = 1:n_sim,
  .packages = c("survival","glmnet","flexsurv","dplyr","janitor")
) %dopar% {
  
  tryCatch({
    ## 1) Simulation
    S <- simulate_one_dataset(
      n         = n,
      beta_true = beta_true,
      dist_name = dist_name,
      correlated = correlated,
      rho        = rho,
      bin_cov    = bin_cov,
      target_cens = target_cens
    )
    
    ## 2) Analyse
    res <- analyze_surv_dataset(
      obj = S, name = paste0("sim", i),
      beta_true = beta_true,          # bei METABRIC_truth = log-HR-Vektor, sonst NULL
      full   = full_flag,
      oracle = oracle_flag,
      split  = split_flag,
      cox    = cox_flag,
      cox0   = cox0_flag,
      fli    = fli_flag,
      debiased        = debiased_flag,
      debiased_lambda = debiased_lambda,
      debiased_nodewise_lambda =
        if (is.na(debiased_node_lambda)) NULL else as.numeric(debiased_node_lambda),
      adaptive                = adaptive,
      nfolds                  = nfolds,
      alpha_glm               = alpha_glm,
      lambda_choice           = lambda_choice,
      ties                    = ties,
      ic_logn                 = ic_logn,
      robust_cox              = robust_cox,
      splitratio              = splitratio,
      support_tol             = support_tol,
      gamma                   = gamma,
      eps                     = eps,
      seed                    = seed0 + i,
      alpha                   = alpha_fit,
      lambda_fix              = lambda_fix
    )
    
    # nach res <- analyze_surv_dataset(...)
    
    if (!is.null(res) && identical(sim_design, "metabric_weibph_truth")) {
      bt_un <- get_metabric_beta_true()   # unskaliert (Original-X)
      s     <- res$s_vec                  # Skalierung aus analyze_surv_dataset()
      
      # auf gemeinsame Namen bringen
      common <- intersect(names(bt_un), names(s))
      bt_un  <- bt_un[common]
      s      <- s[common]
      
      # auf skaliertes X umrechnen (X_scaled = X / s)
      bt_sc <- bt_un * s
      
      res$beta_true <- bt_sc
    }
    
    res
  }, error = function(e) {
    message(sprintf("Error in sim %d: %s", i, conditionMessage(e)))
    NULL
  })
}

outs <- Filter(Negate(is.null), outs)
if (!length(outs)) {
  warning(sprintf(
    "All analyze_surv_dataset() runs returned NULL; skipping config beta_type=%s, n=%d, p=%d, dist=%s, lambda_choice=%s, target_cens=%.2f, rho=%.2f, sim_design=%s.",
    beta_type, n, p, dist_name, lambda_choice, target_cens, rho, sim_design
  ))
  q("no", status = 0)
}#

### Perf nur, wenn beta_true vorhanden
has_truth <- length(outs) > 0 && !is.null(outs[[1]]$beta_true)

if (has_truth) {
  perf_beta <- perf_per_beta_split(outs, methods = methods, alpha_for_ci = alpha_fit)
} else {
  perf_beta <- data.frame()
}


head(perf_beta)
table(perf_beta$method)
#subset(perf_beta, var == "age10")


## -------- Speichern ------------------------------------------------------
datetag <- format(Sys.Date(), "%Y-%m-%d")
timetag <- format(Sys.time(), "%H-%M")

outdir  <- file.path(ROOT, "results", "meta")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tag_base <- sprintf(
  "met2-%s_n%dcen%02d_sim%d_dist-%s_alpha-%.2f_lambda-%s_design-%s_%s",
  beta_type, n, round(100 * target_cens), n_sim, alpha_fit, alpha_fit,
  lambda_choice, sim_design, datetag
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
  adaptive=adaptive,
  alpha = alpha_fit, 
  lambda_choice = lambda_choice, 
  beta_type = beta_type,
  splitratio = splitratio,
  fli = fli_flag, 
  debiased = debiased_flag,
  #debiased_lambda = debiased_lambda, 
  #debiased_nodewise_lambda = debiased_node_lambda,
  methods = paste(methods, collapse = ", "),
  seed0 = seed0, 
  ties = ties, 
  robust_cox = robust_cox,
  sim_design = sim_design,
  lambda_fix = if (is.na(lambda_fix)) NA_real_ else lambda_fix
)
capture.output(str(meta), file = file_meta)
capture.output(sessionInfo(), file = file_meta, append = TRUE)

cat("Saved:\n  - ", file_perf, "\n  - ", file_outs, "\n  - ", file_meta, "\n", sep = "")

