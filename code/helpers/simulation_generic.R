# simulation_generic.R

options(encoding = "UTF-8")

# Add tiny jitter to duplicate values to avoid ties.
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

safe_y_cox <- function(time, status) {
  time <- as.numeric(time)
  status <- as.integer(status > 0)
  
  if (any(time <= 0, na.rm = TRUE)) {
    tpos <- suppressWarnings(min(time[time > 0], na.rm = TRUE))
    if (!is.finite(tpos)) tpos <- 1
    eps <- min(1e-8, tpos * 1e-6)
    time[time <= 0] <- eps
  }
  
  list(time = time, status = status)
}


# Generic simulator (kept as-is, but cleaned).
simulate_fun <- function(
    n, beta, rho = 0.0,
    dist = "weibull", k_shape = 2, lambda = 1,
    correlated = FALSE, bin_cov = FALSE,
    target_censoring = NULL,
    jitter_eps_T = 1e-8, jitter_eps_C = 5e-9
) {
  p <- length(beta)
  
  if (correlated) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS required for correlated X (mvrnorm).")
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  }
  
  if (bin_cov) {
    bin_cols <- sample.int(p, size = floor(p / 3), replace = FALSE)
    X[, bin_cols] <- rbinom(n * length(bin_cols), size = 1, prob = 0.5)
  }
  
  X <- scale(X)
  attr(X, "scaled:center") <- NULL
  attr(X, "scaled:scale")  <- NULL
  
  lin <- as.numeric(X %*% beta)
  
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
    g <- 0.1
    T <- (1 / g) * log(1 - (log(runif(n)) * g / exp(lin)))
    T[T <= 0] <- min(T[T > 0]) * 0.5
  } else {
    stop("Unknown dist.")
  }
  
  Uc <- runif(n)
  if (is.null(target_censoring)) {
    rateC <- 1 / runif(n, 1, 3)
    C <- rexp(n, rate = rateC)
    rateC_used <- NA_real_
  } else {
    r_lo <- 1e-6; r_hi <- 1e2
    r_mid <- NA_real_
    for (it in 1:40) {
      r_mid <- sqrt(r_lo * r_hi)
      C_try <- -log(Uc) / r_mid
      cens_rate <- mean(C_try < T)
      if (cens_rate > target_censoring) r_hi <- r_mid else r_lo <- r_mid
    }
    rateC_used <- r_mid
    C <- -log(Uc) / rateC_used
  }
  
  T <- make_unique(T, eps = jitter_eps_T)
  C <- make_unique(C, eps = jitter_eps_C)
  
  Y <- pmin(T, C)
  status <- as.integer(T <= C)
  
  list(
    time = Y,
    status = status,
    X = X,
    T = T,
    C = C,
    censoring_rate = mean(status == 0),
    target_censoring = target_censoring,
    rateC = rateC_used
  )
}
