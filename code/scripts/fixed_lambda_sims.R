#!/usr/bin/env Rscript
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
    message("Installing to user library: ", lib)
    install.packages(
      miss, lib = lib, dependencies = TRUE,
      Ncpus = 1, INSTALL_opts = c("--no-lock")
    )
  }
}

ensure_packages(c("survival", "glmnet", "future", "doFuture", "foreach", "MASS", "doRNG"))

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(future)
  library(doFuture)
  library(foreach)
  library(MASS)
  library(doRNG)
})

RNGkind("L'Ecuyer-CMRG")

src <- function(f) source(file.path(ROOT, f), chdir = TRUE)
#src("funs.inf.R")  # only needed if selective inference package is not downloadable
#src("funs.fixedCox.R")
#src("funs.fixed.R")
#src("funs.common.R")
src("analyze_toy_sim.R")
src("helpers/funs.metabric_weibph.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  val <- sub(paste0("^--", key, "="), "", hit[1])
  type.convert(val, as.is = TRUE)
}

sim_design    <- get_arg("sim_design", "generic")

beta_type     <- get_arg("beta_type", "realistic")
n             <- get_arg("n", 100)
p             <- get_arg("p", 10)
n_sim         <- get_arg("n_sim", 200)
dist_name     <- get_arg("dist", "weibull")
correlated    <- get_arg("correlated", FALSE)
rho           <- get_arg("rho", 0.0)
bin_cov       <- get_arg("bin_cov", FALSE)
target_cens   <- get_arg("target_censoring", 0.10)

lambda_choice <- get_arg("lambda_choice", "min")
nfolds        <- get_arg("nfolds", 5)
alpha_glm     <- get_arg("alpha_glm", 1)
seed0         <- get_arg("seed", 1000)

out_dir       <- get_arg("out_dir", file.path(ROOT, "out"))
job_tag       <- get_arg("tag", "run")

nlambda           <- get_arg("nlambda", 100L)
lambda_min_ratio  <- get_arg("lambda_min_ratio", 1e-4)
glmnet_thresh     <- get_arg("glmnet_thresh", 1e-7)
glmnet_maxit      <- get_arg("glmnet_maxit", 1e+5)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_beta <- function(beta_type, p) {
  stopifnot(p >= 1)
  
  bt <- switch(
    beta_type,
    "allones" = rep(1, p),
    "highcontrast" = {
      out <- rep(0, p)
      idx <- seq_len(min(4, p))
      out[idx] <- c(0.3, 1.0, 0.3, 1.0)[idx]
      out
    },
    "realistic" = {
      out <- rep(0, p)
      idx <- seq_len(min(4, p))
      out[idx] <- c(0.8, 0.7, 0.5, 0.8)[idx]
      out
    },
    "sparse" = {
      out <- rep(0, p)
      out[1] <- 1
      if (p >= 2) out[2] <- 1
      out
    },
    {
      out <- rep(0, p)
      idx <- seq_len(min(4, p))
      out[idx] <- c(0.8, 0.7, 0.5, 0.8)[idx]
      out
    }
  )
  
  names(bt) <- paste0("X", seq_len(p))
  bt
}

beta_true <- make_beta(beta_type, p)
if (length(beta_true) != p) stop(sprintf("beta_true length (%d) != p (%d).", length(beta_true), p))
names(beta_true) <- paste0("X", seq_along(beta_true))

to_mm <- function(x) {
  all_numeric <- function(z) {
    if (is.matrix(z)) return(is.numeric(z))
    if (is.data.frame(z)) return(all(vapply(z, is.numeric, logical(1))))
    is.numeric(z)
  }
  
  if (all_numeric(x)) {
    M <- if (is.matrix(x)) x else as.matrix(x)
  } else {
    df <- if (is.matrix(x)) as.data.frame(x, check.names = TRUE) else as.data.frame(x, check.names = TRUE)
    if (is.null(names(df))) names(df) <- paste0("V", seq_len(ncol(df)))
    names(df) <- make.names(names(df), unique = TRUE)
    df[] <- lapply(df, function(z) if (is.character(z)) factor(z) else z)
    M <- model.matrix(~ . - 1, data = df, na.action = stats::na.pass)
  }
  
  storage.mode(M) <- "double"
  attr(M, "scaled:center") <- NULL
  attr(M, "scaled:scale")  <- NULL
  M
}

l2_normalize <- function(X, return_scale = TRUE) {
  n  <- nrow(X)
  s_vec <- sqrt(colSums(X^2)) / sqrt(n)
  s_vec[!is.finite(s_vec) | s_vec == 0] <- 1
  Xn <- sweep(X, 2, s_vec, "/")
  if (return_scale) return(list(X = Xn, scale = s_vec))
  Xn
}

make_unique <- function(v, eps = 1e-8) {
  dup <- duplicated(v) | duplicated(v, fromLast = TRUE)
  if (!any(dup)) return(v)
  
  v2 <- v
  for (group in split(which(dup), v[dup])) {
    k <- length(group)
    off <- seq_len(k) - (k + 1) / 2
    v2[group] <- v2[group] + off * eps * max(1, diff(range(v)))
  }
  v2
}

simulate_fun <- function(
    n, beta, rho = 1, dist = "weibull", k_shape = 2, lambda = 1,
    correlated = FALSE, bin_cov = FALSE, target_censoring = NULL,
    jitter_eps_T = 1e-8, jitter_eps_C = 5e-9
) {
  p <- length(beta)
  
  if (correlated) {
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  }
  
  if (bin_cov) {
    binomial_cols <- sample.int(p, size = floor(p / 3), replace = FALSE)
    X[, binomial_cols] <- rbinom(n * length(binomial_cols), size = 1, prob = 0.5)
  }
  
  X <- scale(X)
  attr(X, "scaled:center") <- NULL
  attr(X, "scaled:scale")  <- NULL
  
  lin <- drop(X %*% beta)
  
  if (dist == "exponential") {
    T <- -log(runif(n)) / exp(lin)
  } else if (dist == "weibull") {
    T <- (-(log(runif(n)) / (lambda * exp(lin))))^(1 / k_shape)
  } else if (dist == "lognormal") {
    T <- exp(rnorm(n, mean = lin, sd = 1))
  } else if (dist == "loglogistic") {
    Uu <- runif(n)
    T <- (Uu / (1 - Uu))^(1 / k_shape) * exp(lin)
  } else if (dist == "gompertz") {
    gamma <- 0.1
    T <- (1 / gamma) * log(1 - (log(runif(n)) * gamma / exp(lin)))
    T[T <= 0] <- min(T[T > 0]) * 0.5
  } else {
    stop("unknown distribution")
  }
  
  Uc <- runif(n)
  
  if (is.null(target_censoring)) {
    rateC <- 1 / runif(n, 1, 3)
    C <- rexp(n, rate = rateC)
  } else {
    r_lo <- 1e-6
    r_hi <- 1e+2
    for (it in 1:40) {
      r_mid <- sqrt(r_lo * r_hi)
      C_try <- -log(Uc) / r_mid
      cens_rate <- mean(C_try < T)
      if (cens_rate > target_censoring) r_hi <- r_mid else r_lo <- r_mid
    }
    rateC <- r_mid
    C <- -log(Uc) / rateC
  }
  
  T <- make_unique(T, eps = jitter_eps_T)
  C <- make_unique(C, eps = jitter_eps_C)
  
  Y <- pmin(T, C)
  status <- as.integer(T <= C)
  censoring_rate <- mean(status == 0)
  
  list(
    time = Y,
    status = status,
    X = X,
    T = T,
    C = C,
    censoring_rate = censoring_rate,
    target_censoring = target_censoring,
    rateC = if (is.null(target_censoring)) NA_real_ else rateC
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

slurm_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
slurm_cpus    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(slurm_cpus) || slurm_cpus < 1) slurm_cpus <- 1L

doFuture::registerDoFuture()
future::plan(future::multisession, workers = slurm_cpus)

sim_ids   <- seq_len(n_sim)
doRNGseed <- seed0 + slurm_task_id

simulate_one_dataset <- function(n, beta_true, dist_name,
                                 correlated, rho, bin_cov, target_cens) {
  if (identical(sim_design, "metabric_weibph")) {
    sim_metabric_weibph(n = n)
  } else {
    simulate_fun(
      n = n, beta = beta_true,
      dist = dist_name,
      correlated = correlated, rho = rho,
      bin_cov = bin_cov,
      target_censoring = target_cens
    )
  }
}

res_list <- foreach(
  i = sim_ids,
  .packages = c("glmnet", "survival"),
  .options.RNG = doRNGseed
) %dorng% {
  obj <- simulate_one_dataset(
    n = n,
    beta_true = beta_true,
    dist_name = dist_name,
    correlated = correlated,
    rho = rho,
    bin_cov = bin_cov,
    target_cens = target_cens
  )
  
  X_mm <- to_mm(obj$X)
  tmp  <- l2_normalize(X_mm, return_scale = TRUE)
  X    <- tmp$X
  
  y <- survival::Surv(obj$time, obj$status)
  
  cvfit <- glmnet::cv.glmnet(
    X, y, family = "cox",
    nfolds = nfolds,
    standardize = FALSE,
    alpha = alpha_glm,
    nlambda = nlambda,
    lambda.min.ratio = lambda_min_ratio,
    thresh = glmnet_thresh,
    maxit = glmnet_maxit
  )
  
  lambda_pick <- switch(
    lambda_choice,
    "1se" = cvfit$lambda.1se,
    "min" = cvfit$lambda.min,
    cvfit$lambda.min
  )
  
  list(
    lambda     = lambda_pick,
    cvm_min    = min(cvfit$cvm),
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    cens       = obj$censoring_rate %||% NA_real_
  )
}

lams    <- vapply(res_list, function(z) z$lambda,     numeric(1))
cvm     <- vapply(res_list, function(z) z$cvm_min,    numeric(1))
lam_min <- vapply(res_list, function(z) z$lambda_min, numeric(1))
lam_1se <- vapply(res_list, function(z) z$lambda_1se, numeric(1))
cens    <- vapply(res_list, function(z) z$cens,       numeric(1))

summary_df <- data.frame(
  sim_design = sim_design,
  lambda_choice = lambda_choice,
  beta_type = beta_type,
  n = n,
  p = p,
  dist = dist_name,
  correlated = correlated,
  rho = rho,
  bin_cov = bin_cov,
  target_censoring = target_cens,
  alpha_glm = alpha_glm,
  nfolds = nfolds,
  n_sim = n_sim,
  seed0 = seed0,
  array_id = slurm_task_id,
  mean_lambda   = mean(lams, na.rm = TRUE),
  sd_lambda     = sd(lams, na.rm = TRUE),
  median_lambda = median(lams, na.rm = TRUE),
  q05_lambda    = unname(quantile(lams, 0.05, na.rm = TRUE)),
  q95_lambda    = unname(quantile(lams, 0.95, na.rm = TRUE)),
  mean_cvm       = mean(cvm, na.rm = TRUE),
  mean_censoring = mean(cens, na.rm = TRUE),
  stringsAsFactors = FALSE
)

ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
prefix <- sprintf(
  "bestlambda_%s_%s_design-%s_beta%s_n%d_p%d_seed%d_arr%d",
  ts, job_tag, sim_design, beta_type, n, p, seed0, slurm_task_id
)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res_list, file = file.path(out_dir, paste0(prefix, "_raw.rds")))
write.csv(summary_df, file = file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)

cat("Done\n")
print(summary_df)
