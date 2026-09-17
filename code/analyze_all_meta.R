
# analyze_all_meta.R
to_mm <- function(x) {
  # quick check: all columns numeric?
  all_numeric <- function(z) {
    if (is.matrix(z)) return(is.numeric(z))
    if (is.data.frame(z)) return(all(vapply(z, is.numeric, logical(1))))
    is.numeric(z)
  }
  
  if (all_numeric(x)) {
    M <- if (is.matrix(x)) x else as.matrix(x)   # keep numeric as-is
  } else {
    df <- if (is.matrix(x)) as.data.frame(x, check.names = TRUE) else as.data.frame(x, check.names = TRUE)
    if (is.null(names(df))) names(df) <- paste0("V", seq_len(ncol(df)))
    names(df) <- make.names(names(df), unique = TRUE)
    df[] <- lapply(df, function(z) if (is.character(z)) factor(z) else z)
    M <- model.matrix(~ . - 1, data = df, na.action = stats::na.pass)  # one-hot, no intercept
  }
  
  storage.mode(M) <- "double"         # ensure double
  #M[!is.finite(M)] <- 0               # replace NA/Inf
  attr(M, "scaled:center") <- NULL    # drop scale attrs
  attr(M, "scaled:scale")  <- NULL
  M
}


# L2-normalize and RETURN the scale s_vec = ||x_j||2 / sqrt(n) (stabilizes numerics)
l2_normalize <- function(X, return_scale = TRUE) {
  n  <- nrow(X)
  s_vec <- sqrt(colSums(X^2)) / sqrt(n)
  s_vec[!is.finite(s_vec) | s_vec == 0] <- 1
  Xn <- sweep(X, 2, s_vec, "/")
  if (return_scale) return(list(X = Xn, scale = s_vec))
  Xn
}
#  Negahban-Lambda per Monte-Carlo
negahban_like_lambda_cox <- function(X, time, status, B = 1000,
                                     standardize = TRUE, seed = 1) {
  stopifnot(is.matrix(X))
  if (standardize) X <- scale(X)
  n_events <- sum(status)
  
  # Null-Cox (nur Grundhazard), Martingale-Residuals ~ Score-Anteil je Person
  fit0 <- coxph(Surv(time, status) ~ 1)
  mres <- residuals(fit0, type = "martingale")   # Länge n
  
  set.seed(seed)
  # Multiplier-Bootstrap (Gaussian)
  G <- matrix(rnorm(length(mres) * B), nrow = length(mres), ncol = B)
  # gewichtete Residuen
  WG <- mres * G                                 # n x B
  XT_WG <- t(X) %*% WG                           # p x B
  norms <- apply(XT_WG, 2, function(v) max(abs(v)))
  E_norm <- mean(norms)
  
  lam <- 2 * E_norm / n_events                   # glmnet-kompatible Skalierung
  #list(lambda = as.numeric(lam)#, E_norm = E_norm, norms = norms, n_events = n_events )
  return(as.numeric(lam))
}

make_unique <- function(v, eps = 1e-8){
  dup <- duplicated(v) | duplicated(v, fromLast = TRUE)
  if (!any(dup)) return(v)
  v2 <- v
  for (group in split(which(dup), v[dup])) {
    k <- length(group)
    off <- seq_len(k) - (k + 1)/2
    v2[group] <- v2[group] + off * eps * max(1, diff(range(v)))
  }
  v2
}

# 3) Ensure valid y for Cox (strictly positive times, 0/1 status)
safe_y_cox <- function(time, status) {
  if (any(time <= 0, na.rm = TRUE)) {
    tpos <- suppressWarnings(min(y$time[time > 0], na.rm = TRUE)); if (!is.finite(tpos)) tpos <- 1
    eps <- min(1e-8, tpos * 1e-6); time <- ifelse(time <= 0, eps, time)
  }
  y$time <- time
  y$status <- as.integer(status)  # 0/1
  y
}

## get_ci_df_present, normalize_ci_names, perf_per_beta_split etc. bleiben unverändert
## (ich lasse sie aus Platzgründen hier weg, du hast sie ja oben schon)

## ------------------------------------------------------------
## Helper: Namen bereinigen / normalisieren
## ------------------------------------------------------------
## ------------------------------------------------------------
## Helper: Namen bereinigen / normalisieren
## ------------------------------------------------------------
normalize_ci_names <- function(x) {
  if (is.null(x)) return(x)
  # Backticks weg: `age10` -> age10
  x <- sub("^`(.*)`$", "\\1", x)
  # führende/trailende Leerzeichen
  x <- trimws(x)
  x
}

## ------------------------------------------------------------
## CI-Sammlung über alle Simulationen
## nimmt outs (Liste von analyze_surv_dataset-Listen)
## und res_slot, z.B. "res_fli", "res_wb_std", ...
## ------------------------------------------------------------
get_ci_df_present <- function(outs, res_slot) {
  if (!length(outs)) return(data.frame())
  
  out_list <- list()
  idx <- 1L
  
  for (i in seq_along(outs)) {
    obj <- outs[[i]]
    if (!is.list(obj)) next
    if (!res_slot %in% names(obj)) next
    
    df <- obj[[res_slot]]
    if (is.null(df) || !is.data.frame(df)) next
    if (!all(c("lower", "upper") %in% names(df))) next
    
    ## Variablennamen: bevorzugt Spalte "var", sonst rownames
    var_names <- NULL
    if ("var" %in% names(df)) {
      var_names <- as.character(df$var)
    } else if (!is.null(rownames(df))) {
      var_names <- rownames(df)
    } else {
      var_names <- paste0("V", seq_len(nrow(df)))
    }
    
    ## Namen normalisieren (Backticks etc. entfernen)
    var_names <- normalize_ci_names(var_names)
    
    ## beta_true aus dem Objekt holen; kann NULL sein (METABRIC ohne Truth)
    bt <- obj$beta_true
    if (is.null(bt)) {
      beta_true_vec <- rep(NA_real_, length(var_names))
    } else {
      bt2 <- bt
      if (!is.null(names(bt2))) {
        names(bt2) <- normalize_ci_names(names(bt2))
      }
      beta_true_vec <- bt2[var_names]
      ## für nicht gematchte Namen -> NA
      beta_true_vec[is.na(beta_true_vec)] <- NA_real_
    }
    
    tmp <- data.frame(
      sim       = i,
      var       = var_names,
      lower     = as.numeric(df$lower),
      upper     = as.numeric(df$upper),
      beta_true = as.numeric(beta_true_vec),
      stringsAsFactors = FALSE
    )
    
    out_list[[idx]] <- tmp
    idx <- idx + 1L
  }
  
  out_list <- Filter(Negate(is.null), out_list)
  if (!length(out_list)) return(data.frame())
  do.call(rbind, out_list)
}

## ------------------------------------------------------------
## Performance-Kennzahlen pro Beta / Methode
##  - coverage         = P(truth in CI)
##  - rejection_rate   = P(CI schließt 0 aus)
##  - selective_power  = P(reject | selected, beta_true != 0)
##  - selective_type1  = P(reject | selected, beta_true = 0)
##  - mean_ci_width    = mittlere CI-Breite
## ------------------------------------------------------------
perf_per_beta_split <- function(outs, methods, alpha_for_ci = 0.10) {
  stopifnot(length(outs) > 0)
  methods <- intersect(methods, names(.res_name_map))
  if (!length(methods)) return(data.frame())
  
  res_all <- list()
  idx <- 1L
  
  for (mth in methods) {
    res_slot <- .res_name_map[[mth]]   # z.B. "res_fli", "res_wb_std", ...
    
    ## CIs über alle Replikate einsammeln
    ci_df <- get_ci_df_present(outs, res_slot)
    if (is.null(ci_df) || !nrow(ci_df)) next
    
    needed <- c("var", "lower", "upper", "beta_true")
    if (!all(needed %in% names(ci_df))) next
    
    vars <- sort(unique(ci_df$var))
    
    for (v in vars) {
      sub <- ci_df[ci_df$var == v, , drop = FALSE]
      if (!nrow(sub)) next
      
      truth_vec <- sub$beta_true
      truth <- unique(truth_vec[is.finite(truth_vec)])
      truth <- if (length(truth)) truth[1L] else NA_real_
      
      lower <- sub$lower
      upper <- sub$upper
      ci_width <- upper - lower
      
      ## Coverage (nur, falls truth bekannt)
      in_ci <- if (is.finite(truth)) {
        (lower <= truth & upper >= truth)
      } else {
        rep(NA, length(lower))
      }
      
      ## Ablehnung H0: beta = 0 (0 außerhalb CI)
      reject <- (lower > 0 | upper < 0)
      
      ## „Selektiert“ = CI existiert (beide Grenzen endlich)
      selected <- is.finite(lower) & is.finite(upper)
      
      ## Unbedingte Kennzahlen
      coverage <- if (all(is.na(in_ci))) NA_real_ else mean(in_ci, na.rm = TRUE)
      rejection_rate <- if (all(is.na(reject))) NA_real_ else mean(reject, na.rm = TRUE)
      mean_ci_width <- if (all(is.na(ci_width))) NA_real_ else mean(ci_width[is.finite(ci_width)], na.rm = TRUE)
      
      ## Selektive Kennzahlen
      idx_sel <- which(selected & is.finite(reject))
      n_sel   <- length(idx_sel)
      
      selective_power <- NA_real_
      selective_type1 <- NA_real_
      
      if (n_sel > 0L && is.finite(truth)) {
        if (truth != 0) {
          ## Nicht-Null-Beta: selektive Power
          selective_power <- mean(reject[idx_sel], na.rm = TRUE)
        } else {
          ## Null-Beta: selektiver Typ-I-Fehler
          selective_type1 <- mean(reject[idx_sel], na.rm = TRUE)
        }
      }
      
      ## Gruppe (optional): signal / null / unknown
      group <- if (!is.finite(truth)) {
        "unknown"
      } else if (truth == 0) {
        "null"
      } else {
        "signal"
      }
      
      res_all[[idx]] <- data.frame(
        method          = mth,
        var             = v,
        group           = group,
        beta_true       = truth,
        coverage        = coverage,
        rejection_rate  = rejection_rate,
        selective_power = selective_power,
        selective_type1 = selective_type1,
        mean_ci_width   = mean_ci_width,
        n_rep           = sum(is.finite(lower) & is.finite(upper)),
        n_sel           = n_sel,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  if (!length(res_all)) return(data.frame())
  do.call(rbind, res_all)
}


# nach der Simulation:
#T <- make_unique(T, eps = 1e-8)
#C <- make_unique(C, eps = 5e-9)   # anderer Schritt, reduziert T==C
#Y <- pmin(T, C)
#status <- as.integer(T <= C)

# helper (außerhalb oder ganz oben in analyze_surv_dataset definieren)
.extract_kkt_flag <- function(fli_obj) {
  if (is.null(fli_obj) || inherits(fli_obj, "try-error")) return(0L)
  kk <- try(fli_obj$kkt_ok, silent = TRUE)
  if (!inherits(kk, "try-error") && length(kk) == 1) return(as.integer(isTRUE(kk)))
  pv <- try(fli_obj$pv, silent = TRUE)
  if (!inherits(pv, "try-error") && is.numeric(pv) && length(pv) > 0) {
    return(as.integer(all(is.finite(pv))))
  }
  1L
}

##  Beispiel-Simulator: Cox-Daten mit AR(1)-Kovarianz 
simulate_fun_ar1_cox <- function(n, beta, p = length(beta),
                                 rho = 0.3, h0 = 0.1, censor_rate = 0.35) {
  stopifnot(length(beta) == p)
  # AR(1)-Sigma
  idx <- 1:p
  Sigma <- outer(idx, idx, function(i,j) rho^abs(i - j))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)
  eta <- as.vector(X %*% beta)
  
  # Ereigniszeiten ~ Exp(h0 * exp(eta))
  T_event <- rexp(n, rate = h0 * exp(eta))
  # Zensierung ~ Exp(h0 * censor_rate)
  C <- rexp(n, rate = h0 * censor_rate)
  
  time  <- pmin(T_event, C)
  status <- as.integer(T_event <= C)
  list(X = X, time = time, status = status)
}

simulate_fun1 <- function(n, beta, rho=1, dist = "weibull", k_shape = 2, lambda = 1, correlated = FALSE, bin_cov = FALSE) {
  p <- length(beta)
  
  # Kovariaten-Matrix mit oder ohne Korrelation
  if (correlated) {
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  }
  
  # Falls bin_cov TRUE ist, ersetze 1/3 der Kovariaten durch binomial verteilte Variablen
  if (bin_cov) {
    binomial_cols <- sample(1:p, size = floor(p / 3), replace = FALSE)
    X[, binomial_cols] <- rbinom(n * length(binomial_cols), size = 1, prob = 0.5)
  }
  X <- scale(X)
  # Ãœberlebenszeiten basierend auf gewÃ¤hlter Verteilung
  if (dist == "exponential") {
    T <- -log(runif(n)) / exp(X %*% beta)
  } else if (dist == "weibull") {
    T <- (-(log(runif(n)) / (lambda * exp(X %*% beta))))^(1 / k_shape)
    #arg_T1 <- beta_0*A + colSums(rep(c(1,0,1,0,0,1),2)*log(2) * t(Z))
    #beta_0 <- 2
    #A <- rbinom(n, 1, 0.5)
    #arg_T1 <- beta_0*A + colSums(beta * t(X))
    #T <-  rweibull(n, shape = 2, scale = 10/sqrt(exp(arg_T2)))
    #T <- rweibull(n, shape = 2, scale = 10/sqrt(exp(arg_T1)))
    
  } else if (dist == "lognormal") {
    T <- exp(rnorm(n, mean = X %*% beta, sd = 1))
  } else if (dist == "loglogistic") {
    T <- (runif(n) / (1 - runif(n)))^(1 / k_shape) * exp(X %*% beta)
  } else if (dist == "gompertz") {
    gamma <- 0.1  # Gompertz-Shape-Parameter
    T <- (1 / gamma) * log(1 - (log(runif(n)) * gamma / exp(X %*% beta)))
  } else {
    stop("Unbekannte Verteilung")
  }
  
  # Zensierungszeiten
  U <- runif(n, 1, 3)
  C <- rexp(n, rate = 1/U)
  
  # no ties
  #T <- make_unique(T, eps = 1e-8)
  #C <- make_unique(C, eps = 5e-9)   # anderer Schritt, reduziert T==C
  
  # Beobachtete Zeiten und Statusvariablen
  Y <- pmin(T, C)
  status <- as.integer(T <= C)
  
  return(list(time = Y, status = status, X=X))
}
simulate_one_dataset <- function(n, beta_true, dist_name,
                                 correlated, rho, bin_cov, target_cens) {
  if (identical(sim_design, "metabric_weibph")) {
    # param. METABRIC-Simulation, beta_true wird hier nicht verwendet
    sim_metabric_weibph(n = n)
  } else {
    # bisherige, generische Simulation
    simulate_fun(
      n = n, beta = beta_true,
      dist = dist_name,
      correlated = correlated, rho = rho,
      bin_cov = bin_cov, target_censoring = target_cens
    )
  }
}

simulate_fun <- function(n, beta, rho = 1,
                         dist = "weibull", k_shape = 2, lambda = 1,
                         correlated = FALSE, bin_cov = FALSE,
                         target_censoring = NULL,     # z.B. 0.30; wenn NULL: keine Steuerung
                         jitter_eps_T = 1e-8, jitter_eps_C = 5e-9) {
  p <- length(beta)
  
  # Kovariaten
  if (correlated) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Bitte MASS installieren (mvrnorm).")
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  }
  if (bin_cov) {
    binomial_cols <- sample.int(p, size = floor(p / 3), replace = FALSE)
    X[, binomial_cols] <- rbinom(n * length(binomial_cols), size = 1, prob = 0.5)
  }
  #X <- as.matrix(X)
  #storage.mode(X) <- "double"
  #beta <- as.numeric(beta)
  
  # Standardisieren
  X <- scale(X)
  # nur die scale-Attribute löschen
  attr(X, "scaled:center") <- NULL
  attr(X, "scaled:scale")  <- NULL
  
  # Linearprädiktor – KEIN drop(); as.numeric ist sicher
  lin <- as.numeric(X %*% beta)
  # Überlebenszeiten T
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
  } else stop("Unbekannte Verteilung")
  
  # censoring C
  # Falls target_censoring gesetzt ist, rate per Bisektion suchen (C = -log(Uc)/rateC).
  # Vorab Uniforms ziehen, damit Monotonie in rateC gewährleistet ist.
  Uc <- runif(n)
  if (is.null(target_censoring)) {
    rateC <- 1 / runif(n, 1, 3)   # heterogenes rates
    C <- rexp(n, rate = rateC)
  } else {
    # Bisection auf gemeinsamer Exponentialrate rC in [r_lo, r_hi]
    r_lo <- 1e-6; r_hi <- 1e+2
    for (it in 1:40) {
      r_mid <- sqrt(r_lo * r_hi)              # geometric mean (more stable)
      C_try <- -log(Uc) / r_mid
      cens_rate <- mean(C_try < T)            # amount of censoring = P(C < T)
      if (cens_rate > target_censoring) {
        # zu viel zensiert -> C wird kleiner -> rate zu groß -> senken
        r_hi <- r_mid
      } else {
        r_lo <- r_mid
      }
    }
    rateC <- r_mid
    C <- -log(Uc) / rateC
  }
  
  # no ties
  T <- make_unique(T, eps = jitter_eps_T)
  C <- make_unique(C, eps = jitter_eps_C)
  #
  #observed time/status
  Y <- pmin(T, C)
  status <- as.integer(T <= C)
  
  # calculate censoring rate 
  censoring_rate <- mean(status == 0)
  
  return(list(
    time = Y,
    status = status,
    X = X,
    T = T,
    C = C,
    censoring_rate = censoring_rate,
    target_censoring = target_censoring,
    rateC = if (is.null(target_censoring)) NA_real_ else rateC
  ))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
# sichere Formel-Erstellung
make_surv_formula <- function(Xs) {
  if (is.null(Xs) || ncol(Xs) == 0) return(NULL)
  rhs <- paste0(colnames(Xs), collapse = " + ")
  as.formula(paste("Surv(time, status) ~", rhs))
}

#' Fit penalized Cox (glmnet) with optional adaptive weights
#' Supports lambda selection via min, 1se, aic, bic, negahban.
#' If adaptive = TRUE, performs a second glmnet fit using adaptive penalty factors.
fit_cox_glmnet <- function(
    x_full, time, status, all_coef,
    alpha = 1,
    nfolds = 5,
    lambda_choice = c("min","1se","aic","bic","negahban"),
    adaptive = FALSE,
    gamma = 1,
    eps = 1e-6,
    seed = 123
) {
  lambda_choice <- match.arg(lambda_choice)
  set.seed(seed)
  
  n <- nrow(x_full)
  ok_glm <- FALSE
  err_glm <- NULL
  lambda_pick <- NA_real_
  lambda_pick_ada <- NULL
  beta <- setNames(rep(0, length(all_coef)), all_coef)
  
  res <- tryCatch({
    
    cvfit <- glmnet::cv.glmnet(
      x = x_full,
      y = survival::Surv(time, status),
      family = "cox",
      alpha = alpha,
      nfolds = nfolds,
      standardize = FALSE
    )
    
    glm_fit <- cvfit$glmnet.fit
    
    n_obs <- length(status)
    dev <- glm_fit$dev.ratio
    loglik <- glm_fit$nulldev * dev / (-2)
    df <- glm_fit$df
    
    # Lambda selection
    lambda_pick <- switch(lambda_choice,
                          "min" = cvfit$lambda.min,
                          "1se" = cvfit$lambda.1se,
                          "aic" = {
                            aic <- -2 * loglik + 2 * df
                            glm_fit$lambda[which.min(aic)]
                          },
                          "bic" = {
                            bic <- -2 * loglik + log(n_obs) * df
                            glm_fit$lambda[which.min(bic)]
                          },
                          "negahban" = {
                            p <- length(all_coef)
                            lam_seq <- glm_fit$lambda
                            lam_target <- sqrt(log(p) / n_obs)
                            lam_seq[which.min(abs(lam_seq - lam_target))]
                          }
    )
    
    # Coefficients
    beta <- drop(as.matrix(coef(glm_fit, s = lambda_pick)))
    names(beta) <- all_coef
    ok_glm <- TRUE
    
    if (isTRUE(adaptive)) {
      penalty.factor <- 1 / (abs(beta) + eps)^gamma
      
      cvfit_ada <- glmnet::cv.glmnet(
        x = x_full,
        y = survival::Surv(time, status),
        family = "cox",
        alpha = alpha,
        nfolds = nfolds,
        standardize = FALSE,
        penalty.factor = penalty.factor
      )
      
      # adaptive lambda selection
      lambda_pick_ada <- switch(lambda_choice,
                                "min" = cvfit_ada$lambda.min,
                                "1se" = cvfit_ada$lambda.1se,
                                "aic" = {
                                  dev_path <- deviance(cvfit_ada$glmnet.fit, newx = x_full,
                                                       y = survival::Surv(time, status))
                                  df_path <- cvfit_ada$glmnet.fit$df
                                  aic <- dev_path + 2 * df_path
                                  cvfit_ada$glmnet.fit$lambda[which.min(aic)]
                                },
                                "bic" = {
                                  dev_path <- deviance(cvfit_ada$glmnet.fit, newx = x_full,
                                                       y = survival::Surv(time, status))
                                  df_path <- cvfit_ada$glmnet.fit$df
                                  bic <- dev_path + log(sum(status)) * df_path
                                  cvfit_ada$glmnet.fit$lambda[which.min(bic)]
                                },
                                "negahban" = {
                                  p <- length(all_coef)
                                  lam_seq <- cvfit_ada$glmnet.fit$lambda
                                  lam_target <- sqrt(log(p) / n_obs)
                                  lam_seq[which.min(abs(lam_seq - lam_target))]
                                }
      )
      
      beta <- drop(as.matrix(coef(cvfit_ada$glmnet.fit, s = lambda_pick_ada)))
      names(beta) <- all_coef
    }
    
    list(cvfit = cvfit, lambda_pick_ada = lambda_pick_ada, beta = beta)
  }, error = function(e) {
    err_glm <- conditionMessage(e)
    NULL
  })
  
  list(
    beta = beta,
    lambda_pick = lambda_pick,
    lambda_pick_ada = lambda_pick_ada,
    ok_glm = ok_glm,
    err_glm = err_glm,
    used_cols = all_coef,
    cvfit = if (!is.null(res)) res$cvfit else NULL
  )
}

# Kleiner Helfer: backticks entfernen; 1..p -> X1..Xp
norm_names <- function(nms) {
  nms <- sub("^`(.*)`$", "\\1", nms)
  if (all(grepl("^[0-9]+$", nms))) nms <- paste0("X", as.integer(nms))
  nms
}

.riskset_sums <- function(time, status, X, eta) {
  w <- exp(eta)
  S0 <- rev(cumsum(rev(w)))
  WX <- w * X
  S1 <- apply(WX, 2, function(col) rev(cumsum(rev(col))))
  list(S0 = S0, S1 = S1)
}

.compute_event_contribs <- function(time, status, X, eta) {
  ord <- order(time, decreasing = FALSE)
  time <- time[ord]; status <- status[ord]; X <- X[ord, , drop = FALSE]
  eta <- as.numeric(eta)[ord]
  rs <- .riskset_sums(time, status, X, eta)
  xbar <- sweep(rs$S1, 1, rs$S0, "/")
  ev_local <- which(status == 1L)
  z_events <- X[ev_local, , drop = FALSE] - xbar[ev_local, , drop = FALSE]
  ev_idx_global <- ord[ev_local]
  list(z_events = z_events, ev_idx_global = ev_idx_global)
}

.nodewise_precision <- function(Z, lambda = NULL, standardize = TRUE, cv_folds = 5) {
  p <- ncol(Z)
  Theta <- matrix(0, p, p)
  colnames(Theta) <- rownames(Theta) <- colnames(Z)
  for (j in seq_len(p)) {
    Zj <- Z[, j]; Zm <- Z[, -j, drop = FALSE]
    if (ncol(Zm) == 0L) { Theta[j, j] <- 1 / var(Zj); next }
    if (is.null(lambda)) {
      cvfit <- glmnet::cv.glmnet(Zm, Zj, family = "gaussian",
                                 intercept = TRUE, standardize = standardize)
      fit <- glmnet::glmnet(Zm, Zj, family = "gaussian",
                            intercept = TRUE, standardize = standardize,
                            lambda = cvfit$lambda.min)
    } else {
      fit <- glmnet::glmnet(Zm, Zj, family = "gaussian",
                            intercept = TRUE, standardize = standardize,
                            lambda = lambda)
    }
    g <- as.numeric(coef(fit)); r <- Zj - (g[1] + Zm %*% g[-1])
    tau2 <- mean(r^2)
    Theta[j, j] <- 1 / tau2
    Theta[-j, j] <- -g[-1] / tau2
  }
  0.5 * (Theta + t(Theta))
}

debiased_cox_lasso <- function(X, time, status,
                               lambda_cox = c("lambda.min","lambda.1se"),
                               nodewise_lambda = NULL,
                               nfolds = 5) {
  X <- as.matrix(X)
  p <- ncol(X)
  # glmnet Cox
  cvfit <- glmnet::cv.glmnet(X, survival::Surv(time, status), family = "cox",
                             nfolds = nfolds, standardize = FALSE, intercept = FALSE)
  lam_use <- if (match.arg(lambda_cox) == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min
  fit <- glmnet::glmnet(X, survival::Surv(time, status), family = "cox",
                        lambda = lam_use, standardize = FALSE, intercept = FALSE)
  bhat <- as.numeric(coef(fit)); names(bhat) <- colnames(X)
  eta  <- as.numeric(X %*% bhat)
  
  ev <- .compute_event_contribs(time, status, X, eta)
  Z <- ev$z_events; m <- nrow(Z)
  U <- colSums(Z); Sigma_z <- crossprod(Z) / m
  colnames(Z) <- names(bhat)
  Theta <- .nodewise_precision(Z, lambda = nodewise_lambda, standardize = TRUE, cv_folds = nfolds)
  
  # debias & se
  beta_tilde <- bhat + as.numeric(Theta %*% (U / m))
  V <- (Theta %*% Sigma_z %*% t(Theta)) / m
  se <- sqrt(pmax(diag(V), 0))
  
  list(beta_debiased = setNames(beta_tilde, colnames(X)),
       se = setNames(se, colnames(X)),
       Theta = Theta, Z = Z, lambda_used = lam_use, events = m)
}
debiased_cox_from_beta <- function(X, time, status, beta_hat,
                                   nodewise_lambda = NULL, nfolds = 5) {
  X <- as.matrix(X)
  stopifnot(length(status) == nrow(X), length(time) == nrow(X))
  # Namenssicherheit
  if (is.null(names(beta_hat))) names(beta_hat) <- colnames(X)
  beta_hat <- as.numeric(beta_hat); names(beta_hat) <- colnames(X)
  
  eta  <- as.numeric(X %*% beta_hat)
  ev   <- .compute_event_contribs(time, status, X, eta)
  Z    <- ev$z_events
  colnames(Z) <- colnames(X)
  
  m        <- nrow(Z)
  U        <- colSums(Z)
  Sigma_z  <- crossprod(Z) / m
  Theta    <- .nodewise_precision(Z, lambda = nodewise_lambda,
                                  standardize = TRUE, cv_folds = nfolds)
  
  beta_tilde <- beta_hat + as.numeric(Theta %*% (U / m))
  V          <- (Theta %*% Sigma_z %*% t(Theta)) / m
  se         <- sqrt(pmax(diag(V), 0))
  
  list(beta_debiased = setNames(beta_tilde, colnames(X)),
       se            = setNames(se,          colnames(X)),
       Theta = Theta, Z = Z, events = m)
}

.draw_weights_event <- function(m, type = c("rademacher","normal","mammen",
                                            "weird_poisson","weird_exponential")) {
  type <- match.arg(type)
  if (type == "rademacher") sample(c(-1,1), m, TRUE)
  else if (type == "normal") rnorm(m)
  else if (type == "mammen") {
    p1 <- (sqrt(5)+1)/(2*sqrt(5)); a <- (1 - sqrt(5))/2; b <- (1 + sqrt(5))/2
    ifelse(runif(m) < p1, a, b)
  } else if (type == "weird_poisson") rpois(m, 1) - 1
  else rexp(m, 1) - 1
}

.draw_weights_subject <- function(n, type = c("normal","rademacher")) {
  type <- match.arg(type); if (type == "normal") rnorm(n) else sample(c(-1,1), n, TRUE)
}

wildbootstrap_debiased_cox <- function(fit_db, X, time, status, B = 500,
                                       scheme = c("event_rademacher","event_normal","event_mammen",
                                                  "weird_poisson","weird_exponential",
                                                  "subject_normal","subject_rademacher")) {
  scheme <- match.arg(scheme)
  Z <- fit_db$Z; Theta <- fit_db$Theta; m <- nrow(Z); p <- ncol(Z); n <- nrow(X)
  beta_hat <- (as.numeric(coef(glmnet::glmnet(X, survival::Surv(time, status),
                                              family="cox",
                                              lambda = fit_db$lambda_used,
                                              standardize = FALSE, intercept = FALSE))))
  # Score draws
  U_boot <- matrix(0, p, B)
  if (startsWith(scheme, "subject_")) {
    g_subj <- replicate(B, .draw_weights_subject(n, sub("^subject_", "", scheme)))
    g <- g_subj[fit_db$ev_idx_global, , drop = FALSE]  # m x B
    U_boot <- t(Z) %*% g                                # p x B
  } else {
    g <- replicate(B, .draw_weights_event(m, sub("^(event_|weird_)", "", scheme)))
    U_boot <- t(Z) %*% g                                # p x B
  }
  Delta <- Theta %*% (U_boot / m)                       # p x B
  beta_star <- matrix(fit_db$beta_debiased, p, B) + Delta
  rownames(beta_star) <- colnames(X)
  beta_star
}


# Kompakte Zusammenfassung (Bias/SE) aus Boot-Koeffizienten
.summarize_lasso_boot <- function(boot_coefs, beta_hat_named) {
  stopifnot(is.matrix(boot_coefs))
  bh <- as.numeric(beta_hat_named)
  names(bh) <- names(beta_hat_named)
  boot_bias <- colMeans(boot_coefs, na.rm = TRUE) - bh
  boot_se   <- apply(boot_coefs, 2, function(z) sd(z, na.rm = TRUE))
  data.frame(variable = names(bh),
             beta_hat = bh,
             bias_boot = as.numeric(boot_bias),
             se_boot   = as.numeric(boot_se),
             row.names = NULL, check.names = FALSE)
}

.summarize_lasso_boot_all <- function(boot_coefs, beta_hat_named, alpha = 0.10) {
  stopifnot(is.matrix(boot_coefs))
  bh <- as.numeric(beta_hat_named); names(bh) <- names(beta_hat_named)
  z  <- qnorm(1 - alpha/2)
  se_boot   <- apply(boot_coefs, 2, sd, na.rm = TRUE)
  mu_boot   <- colMeans(boot_coefs, na.rm = TRUE)
  bias_boot <- mu_boot - bh
  # Guard gegen se==0
  se_eff <- ifelse(is.finite(se_boot) & se_boot > 0, se_boot, NA_real_)
  # Perzentil (abs) und studentisiert
  Tstar   <- sweep(boot_coefs, 2, bh, "-")                    # B x p (Differenzen)
  q_abs   <- apply(abs(Tstar), 2, quantile, 1 - alpha/2, na.rm = TRUE, type = 1)
  Tstd    <- sweep(Tstar, 2, se_eff, "/")
  q_std   <- apply(abs(Tstd),  2, quantile, 1 - alpha/2, na.rm = TRUE, type = 1)
  
  make_df <- function(lo, up) {
    data.frame(
      var = names(bh),
      bhat = unname(bh),
      se_boot = unname(se_boot),
      bias_boot = unname(bias_boot),
      lower = unname(lo),
      upper = unname(up),
      check.names = FALSE
    )
  }
  list(
    df_abs  = make_df(bh - q_abs,          bh + q_abs),
    df_std  = make_df(bh - q_std*se_eff,   bh + q_std*se_eff),
    df_wald = make_df(bh - z*se_eff,       bh + z*se_eff)
  )
}

# Residual-Bootstrap (Martingale + offset) – verwendet X/time/status aus analyze()-Kontext
.cox_lasso_boot_residual_inline <- function(X, time, status,
                                            lambda_opt, B = 200, a_n = 0.01,
                                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Basisfit für Residuen und linearen Prädiktor (gleiche Daten!)
  df <- data.frame(time = time, status = status, X, check.names = FALSE)
  f  <- as.formula(paste0("Surv(time, status) ~ ",
                          paste(colnames(X), collapse = " + ")))
  fit_cox <- survival::coxph(f, data = df, ties = "efron")
  lp  <- as.numeric(predict(fit_cox, type = "lp"))
  M   <- as.numeric(residuals(fit_cox, type = "martingale"))
  
  p <- ncol(X)
  boot_coefs <- matrix(NA_real_, nrow = B, ncol = p,
                       dimnames = list(NULL, colnames(X)))
  for (b in seq_len(B)) {
    resampled <- sample(M - mean(M), replace = TRUE)
    offset_lp <- lp + resampled
    bf <- glmnet::glmnet(x = X,
                         y = survival::Surv(time, status),
                         family = "cox",
                         lambda = lambda_opt,
                         standardize = FALSE,
                         offset = offset_lp)
    boot_coefs[b, ] <- as.numeric(coef(bf))
  }
  boot_coefs
}

# Pairs-Bootstrap (Case-Resampling); optional pro Bootstrap Lambda re-fitten
.cox_lasso_boot_pairs_inline <- function(X, time, status,
                                         lambda_base,
                                         B = 200,
                                         refit_lambda_each_boot = FALSE,
                                         alpha_glm = 1, nfolds = 5,
                                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  boot_coefs <- matrix(NA_real_, nrow = B, ncol = p,
                       dimnames = list(NULL, colnames(X)))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    xb  <- X[idx, , drop = FALSE]
    yb  <- survival::Surv(time[idx], status[idx])
    
    lb <- lambda_base
    if (refit_lambda_each_boot) {
      cvb <- glmnet::cv.glmnet(xb, yb, family = "cox",
                               nfolds = nfolds, standardize = FALSE, alpha = alpha_glm)
      lb  <- cvb$lambda.min
    }
    bf <- glmnet::glmnet(x = xb, y = yb, family = "cox",
                         lambda = lb, standardize = FALSE, alpha = alpha_glm)
    boot_coefs[b, ] <- as.numeric(coef(bf))
  }
  boot_coefs
}

.estimate_H0_offset <- function(time, status, eta) {
  fit_off <- survival::coxph(
    survival::Surv(time, status) ~ 1 + offset(eta),
    ties = "breslow"
  )
  bh <- survival::basehaz(fit_off, centered = FALSE)
  list(time = bh$time, H0 = bh$hazard)
}

.H0_inv_step <- function(H0_time, H0_vals, y) {
  maxH <- max(H0_vals)
  if (y > maxH) return(Inf)
  idx <- which(H0_vals >= y)[1]
  H0_time[idx]
}

.estimate_censor_sampler <- function(time, status) {
  statusC <- 1L - status
  fitC <- survival::survfit(survival::Surv(time, statusC) ~ 1)
  times <- c(0, fitC$time)
  Gvals <- c(1, fitC$surv)         # G(t) = P(C > t)
  Fc <- 1 - Gvals                  # Verteilungsfunktion der Zensurzeiten
  
  function(n) {
    u <- runif(n)
    out <- numeric(n)
    for (i in seq_len(n)) {
      idx <- which(Fc >= u[i])[1]
      if (is.na(idx)) {
        # u liegt über letzter Stufe -> sehr große Zensurzeit
        # einfache Extrapolation
        dt <- diff(times)
        dt_pos <- dt[dt > 0]
        step <- if (length(dt_pos)) median(dt_pos) else 1
        out[i] <- max(times) + rexp(1, rate = 1 / (step + 1e-6))
      } else {
        out[i] <- times[idx]
      }
    }
    out
  }
}

.simulate_T_given_H0_eta <- function(eta, H0_time, H0_vals) {
  n <- length(eta)
  u <- runif(n)
  y <- -log(u) / exp(eta)
  vapply(y, function(yi) .H0_inv_step(H0_time, H0_vals, yi), numeric(1))
}

.cox_lasso_boot_parametric_inline <- function(
    X, time, status,
    lambda_base,
    B = 200,
    alpha_glm = 1,
    nfolds = 5,
    seed = NULL,
    use_cv_lambda = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  y_surv <- survival::Surv(time, status)
  n <- length(time); p <- ncol(X)
  
  ## 1) Basis-Lasso-Fit
  if (use_cv_lambda) {
    cvfit <- glmnet::cv.glmnet(
      X, y_surv, family = "cox",
      alpha = alpha_glm, nfolds = nfolds,
      standardize = FALSE
    )
    lambda_use <- cvfit$lambda.min
  } else {
    lambda_use <- lambda_base
  }
  fit <- glmnet::glmnet(
    X, y_surv, family = "cox",
    alpha = alpha_glm,
    lambda = lambda_use,
    standardize = FALSE
  )
  beta_hat <- drop(as.matrix(coef(fit)))
  names(beta_hat) <- rownames(coef(fit))
  
  eta_hat <- as.numeric(X %*% beta_hat)
  
  ## 2) H0 und Censurmodell schätzen
  H0 <- .estimate_H0_offset(time, status, eta_hat)
  censor_sampler <- .estimate_censor_sampler(time, status)
  
  ## 3) Bootstrap-Schleife
  boot_coefs <- matrix(NA_real_, nrow = B, ncol = p,
                       dimnames = list(NULL, colnames(X)))
  for (b in seq_len(B)) {
    T_sim <- .simulate_T_given_H0_eta(eta_hat, H0$time, H0$H0)
    C_sim <- censor_sampler(n)
    time_b <- pmin(T_sim, C_sim)
    status_b <- as.integer(T_sim <= C_sim)
    
    yb <- survival::Surv(time_b, status_b)
    
    fit_b <- tryCatch({
      if (use_cv_lambda) {
        cvb <- glmnet::cv.glmnet(
          X, yb, family = "cox",
          alpha = alpha_glm, nfolds = nfolds,
          standardize = FALSE
        )
        lb <- cvb$lambda.min
      } else {
        lb <- lambda_use
      }
      glmnet::glmnet(
        X, yb, family = "cox",
        alpha = alpha_glm,
        lambda = lb,
        standardize = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit_b)) {
      cb <- drop(as.matrix(coef(fit_b)))
      names(cb) <- rownames(coef(fit_b))
      cb_aligned <- cb[colnames(X)]
      boot_coefs[b, ] <- cb_aligned
    }
  }
  boot_coefs
}

.cox_lasso_boot_loughin_inline <- function(X, time, status,
                                           lambda_opt,
                                           B = 200,
                                           seed = NULL,
                                           ties = "efron") {
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- Unpenalized Cox to compute u_i ---
  df <- data.frame(time = time, status = status, X)
  form <- as.formula(
    paste0("Surv(time, status) ~ ", paste(colnames(X), collapse = " + "))
  )
  fit_cox <- survival::coxph(form, data = df, ties = ties)
  
  lp <- as.numeric(predict(fit_cox, type = "lp"))
  g  <- exp(lp)
  
  bh <- survival::basehaz(fit_cox, centered = FALSE)
  H0 <- approx(bh$time, bh$hazard, xout = time,
               method = "constant", rule = 2)$y
  
  u <- exp(-H0 * g)  # Survival probabilities
  
  # --- Bootstrap storage ---
  p <- ncol(X)
  boot_coefs <- matrix(NA_real_, nrow = B, ncol = p,
                       dimnames = list(NULL, colnames(X)))
  
  # --- Loughin bootstrap loop ---
  for (b in seq_len(B)) {
    
    idx <- sample(seq_along(u), replace = TRUE)
    u_star     <- u[idx]
    delta_star <- status[idx]
    
    # probability-scale time transform (stable version)
    y_star <- 1 - exp(log(u_star) / g)
    
    fit_star <- glmnet::glmnet(
      x = X,
      y = survival::Surv(y_star, delta_star),
      family = "cox",
      lambda = lambda_opt,
      standardize = FALSE
    )
    
    boot_coefs[b, ] <- as.numeric(coef(fit_star))
  }
  
  boot_coefs
}


# ===================== main function ====================
#' Analyze a survival dataset with multiple estimators
#'
#' Runs a configurable pipeline for right-censored survival data:
#' optional train/test split, Cox variants, (weighted) bootstrap,
#' (adaptive) penalized GLM via cross-validation, and FLI variants.
#'
#' @param obj List-like dataset (e.g., from a simulator) with at least
#'   `time` (numeric), `status` (0/1), and `X` (matrix/data.frame).
#' @param name Name tag for the dataset/run.
#' @param beta_true Optional true coefficients for oracle/eval.
#' @param full If `TRUE`, compute extra diagnostics/outputs.
#' @param oracle If `TRUE`, fit oracle models using `beta_true`.
#' @param split If `TRUE`, split into train/test (use `splitratio`).
#' @param cox If `TRUE`, fit standard Cox model.
#' @param cox0 If `TRUE`, fit null Cox (baseline only).
#' @param coxwb If `TRUE`, fit Cox with wild bootstrap SEs.
#' @param fli If `TRUE`, fit FLI model.
#' @param fliwb If `TRUE`, fit FLI with wild bootstrap SEs.
#' @param adaptive If `TRUE`, use adaptive penalties where applicable.
#' @param nfolds Number of CV folds (glm-like procedures).
#' @param alpha_glm Elastic-net mixing (1 = lasso, 0 = ridge).
#' @param lambda_choice Rule to pick λ: `"min"`, `"1se"`, `"aic"`,
#'   `"bic"`, or `"negahban"`.
#' @param ties Method for Cox ties (`"efron"`, `"breslow"`, `"exact"`).
#' @param ic_logn Sample size for log-normal IC (`"events"` or `"n"`).
#' @param robust_cox If `TRUE`, use robust (sandwich) variance in Cox.
#' @param splitratio Train fraction if splitting (e.g., 0.5).
#' @param boot If `TRUE`, run bootstrap inference.
#' @param wildboot If `TRUE`, use wild bootstrap scheme.
#' @param B_boot Number of bootstrap replications.
#' @param support_tol Tolerance for support recovery (|β| > tol).
#' @param gamma Additional tuning parameter (method-specific).
#' @param eps Numerical tolerance (convergence/guard).
#' @param seed RNG seed for reproducibility.
#' @param alpha Nominal level for CIs/tests (e.g., 0.10).
#'
#' @details Expects right-censored data. `X` may be numeric or mixed;
#' it will be coerced to a numeric design. Duplicate `split`
#' arguments should be consolidated (keep one).
#'
#' @return A list with fitted objects, CV choices, metrics, (wild) bootstrap
#' results if requested, and the arguments used.
#'
#' @examples
#' # res <- analyze_surv_dataset(obj, cox = TRUE, boot = FALSE)

analyze_surv_dataset <- function(
    obj, name = "dataset",
    beta_true = NULL,
    full = TRUE,
    oracle = TRUE,
    split = TRUE,
    cox = TRUE,
    cox0 = TRUE,
    coxwb = TRUE,
    fli = TRUE,
    fliwb = TRUE,
    debiased = TRUE,                 # <--- NEU
    debiased_wb = TRUE,             # <--- NEU
    debiased_lambda = c("lambda.min","lambda.1se"),
    debiased_nodewise_lambda = NULL,
    debiased_wb_scheme = c("event_mammen","event_rademacher","event_normal",
                           "weird_poisson","weird_exponential",
                           "subject_normal","subject_rademacher"),
    B_boot_db = 500,                 # Replizahlen für debiased-WB
    lasso_boot_residual = TRUE,
    B_lasso_boot_resid = 200,
    a_n_lasso = 0.01,
    lasso_boot_pairs = TRUE,
    B_lasso_boot_pairs = 200,
    refit_lambda_each_boot = FALSE,
    lasso_boot_param = TRUE,
    B_lasso_boot_param = 200,
    adaptive = FALSE,
    nfolds = 5,
    alpha_glm = 1,
    lambda_choice = c("min","1se","aic","bic","negahban","fix"),
    ties = "efron", ic_logn = c("events","n"),
    robust_cox = FALSE,
    splitratio = 0.5,
    boot = FALSE,
    B_boot = 900,
    support_tol = 0,
    gamma = 1,
    eps = 1e-6,
    seed = 123,
    alpha = 0.10,
    lambda_fix = NULL,
    
    randlasso = TRUE,
    randlasso_sampler = c("norejection","adaptMCMC"),
    randlasso_nsample = 5000,
    randlasso_burnin  = 1000,
    
    lasso_boot_loughin = TRUE,
    B_lasso_boot_loughin = 200
    
) {
 
  # lambda from above
  lambda_choice  <- tolower(lambda_choice)
  lambda <- match.arg(lambda_choice,
                             c("min","1se","aic","bic","negahban","fix")
  )
  randlasso_sampler <- match.arg(randlasso_sampler)
  metabric_order <- c(
    "age10",
    "tumor_stage_f0",
    "tumor_stage_f1",
    "tumor_stage_f2",
    "tumor_stage_f3",
    "tumor_stage_f4",
    "er_status_fER_pos",
    "her2_status_fHER2_pos",
    "nottingham_prognostic_index"
  )
  
  X_mm <- to_mm(obj$X)
  tmp  <- l2_normalize(X_mm, return_scale = TRUE)
  X    <- tmp$X
  s_vec <- tmp$scale
  
  # --- FINAL: Names + Order festziehen (genau 1x) ---
  if (is.null(colnames(X))) colnames(X) <- colnames(X_mm)
  if (is.null(colnames(X))) colnames(X) <- names(obj$X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  
  # Metabric order
  if (all(metabric_order %in% colnames(X))) {
    X <- X[, metabric_order, drop = FALSE]
  }
  
  all_coef <- colnames(X)
  names(s_vec) <- all_coef
  
  
  time <- obj$time
  status <- obj$status
  n <- nrow(X); p <- ncol(X)

  
  z_value <- z <- qnorm(1 - alpha/2)   # = 1.959964 für 95%-KI
  
  # glmnet path
  fit <- glmnet(X, Surv(time, status), family = "cox", standardize = FALSE)
  
  lambda_pick   <- NULL
  lambda_1se    <- NULL
  #cvfit         <- NULL
  
  ok_full   <- ok_oracle <- ok_glm   <- ok_adapt     <- ok_split  <- ok_refit  <- NA_integer_
  ok_refit0 <- ok_wb_abs <-ok_wb_std <- ok_wb_normal <- ok_wb_abs <- ok_debiased_wb <- ok_debiased <- NA_integer_
  ok_fli    <- ok_fli_wb  <-kkt_fli   <- kkt_fli_wb <- ok_wb  <- NA_character_
  
  err_refit <- err_refit0<-err_wb_abs<- err_wb_std   <- err_wb_normal <- err_fli <- NA_integer_
  err_oracle<- err_wb_fli <- err_full <- err_debiased_wb <- err_debiased <- err_glm <- err_adaptive <- NA_integer_
  
  ### NEU:
  ok_randlasso  <- NA_integer_
  err_randlasso <- NA_character_
  
  # needed to get the right size for the result dataset below
  ##all_coef <- paste0("X", 1:p)
  
  # returns of split 
  beta_split <- NULL
  bhat_split <- NULL
  m_split    <- NULL
  
  # return of wb
  se_beta_boot  <- setNames(rep(NA_real_, length(all_coef)), all_coef)
  p_boot_std    <- setNames(rep(NA_real_, length(all_coef)), all_coef)
  # tiny helper
  init_df <- function(all_coef) {
    out <- data.frame(
      var     = all_coef,
      bhat    = NA_real_,
      se      = NA_real_,
      lower   = NA_real_,
      upper   = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_
    )
    rownames(out) <- all_coef
    out
  }
  
  all_names <- c(
    "res_full",
    "res_oracle",
    "res_refit",
    "res_refit0",
    "res_split",
    "res_wb_abs",
    "res_wb_std",
    "res_wb_normal",
    "res_fli",
    "res_fli_wb", 
    "res_debiased", 
    "res_debiased_wb", 
    "res_debiased_wb_abs", 
    "res_debiased_wb_std", 
    "res_debiased_wb_normal",
    "res_lasso_boot_resid",
    "res_lasso_boot_resid_abs",
    "res_lasso_boot_resid_std",
    "res_lasso_boot_resid_wald",
    "res_lasso_boot_pairs",
    "res_lasso_boot_pairs_abs",
    "res_lasso_boot_pairs_std",
    "res_lasso_boot_pairs_wald",
    "res_lasso_boot_param",
    "res_lasso_boot_param_abs",
    "res_lasso_boot_param_std",
    "res_lasso_boot_param_wald",
    "res_randlasso",
    "res_lasso_boot_loughin",
    "res_lasso_boot_loughin_abs",
    "res_lasso_boot_loughin_std",
    "res_lasso_boot_loughin_wald" 
  )
  # --- Sanity check: all_coef muss zu X passen ---
  stopifnot(identical(all_coef, colnames(X)))
  
  models <- setNames(
    replicate(length(all_names), init_df(all_coef), simplify = FALSE),
    all_names
  )
  list2env(models, envir = environment())
  rm(models)
  # --- Sanity check: Ergebnis-Container müssen die gleichen Rowname wie X haben ---
  stopifnot(identical(rownames(res_debiased), colnames(X)))
  
  beta <- NULL
  fit <- glmnet(
    X, Surv(time, status), family = "cox", standardize = FALSE, alpha = alpha_glm)
  if (lambda == "min"|| lambda == "1se"){
    cvfit <- cv.glmnet(
      X, Surv(time, status),
      family = "cox", nfolds = nfolds, standardize = FALSE, alpha = alpha_glm
    )
    lambda_pick <- if (lambda == "min")  cvfit$lambda.min else cvfit$lambda.1se
    
  } else if (lambda == "negahban"){
    lambda_pick <- negahban_like_lambda_cox(X, time, status, B = 1000, standardize = TRUE)
    
  } else if (lambda == "aic" || lambda == "bic") {
    dev_path <- deviance(fit, newx = X, y = Surv(time, status))
    df_path  <- fit$df
    lam_path <- fit$lambda
    n_events <- sum(status)
    if(lambda == "aic"){
      AIC_path <- dev_path + 2 * df_path
      lambda_pick <- lam_path[which.min(AIC_path)]
    } else {
      BIC_path <- dev_path + log(n_events) * df_path
      lambda_pick <- lam_path[which.min(BIC_path)]
    }
  } else if (lambda_choice == "fix") {
    if (!is.null(lambda_fix) &&
        length(lambda_fix) == 1 &&
        is.finite(lambda_fix) &&
        lambda_fix > 0) {
      lambda_pick <- lambda_fix
    } else {
      stop("lambda_fix must be a single positive, finite number")
    }
  }
  
  # optionales Refitting, falls lambda_pick außerhalb der Pfadgrenzen
  lambda_grid <- fit$lambda
  lam_min <- min(lambda_grid)
  lam_max <- max(lambda_grid)
  
  if (!(lambda_pick >= lam_min && lambda_pick <= lam_max)) {
    idx_closest <- which.min(abs(lambda_grid - lambda_pick))
    lambda_pick_adj <- lambda_grid[idx_closest]
    warning(sprintf(
      "lambda_pick=%.6g liegt außerhalb des glmnet-Grids [%.6g, %.6g]; verwende nächstliegenden Wert %.6g.",
      lambda_pick, lam_min, lam_max, lambda_pick_adj
    ))
    lambda_pick <- lambda_pick_adj
  }
  
  
  beta <- drop(as.matrix(coef(fit, s = lambda_pick)))
  if (is.null(names(beta)) && !is.null(colnames(X))) names(beta) <- colnames(X)
  if (!is.null(beta)) ok_glm <- 1L else {
    ok_glm <- 0L
    stop("glm didnt work")
  }
  
  lambda_pick_ada <- NULL
  if (adaptive == TRUE){
    w <- abs(beta[beta != 0])
    if (anyNA(w)) w[is.na(w)] <- 0
    penalty.factor <- 1 / (abs(beta) + eps)^gamma
    
    fit_ada <- glmnet(
      X, Surv(obj$time, obj$status), family = "cox",
      alpha = alpha_glm, standardize = FALSE,
      penalty.factor = penalty.factor
    )
    
    lambda_grid_ada <- fit_ada$lambda
    dev_path <- deviance(fit_ada, newx = X, y = Surv(obj$time, obj$status))
    df_path  <- fit_ada$df
    n_events <- sum(obj$status)
    
    if (lambda == "min" || lambda == "1se") {
      # wenn du hier wirklich CV willst, dann cv.glmnet + lambda.min/1se
      cvfit_ada <- cv.glmnet(
        X, Surv(obj$time, obj$status), family = "cox",
        nfolds = nfolds, standardize = FALSE, alpha = alpha_glm,
        penalty.factor = penalty.factor
      )
      lambda_pick_ada <- if (lambda == "min") cvfit_ada$lambda.min else cvfit_ada$lambda.1se
    } else if (lambda == "aic") {
      AIC_path <- dev_path + 2 * df_path
      lambda_pick_ada <- lambda_grid_ada[which.min(AIC_path)]
    } else if (lambda == "bic") {
      BIC_path <- dev_path + log(n_events) * df_path
      lambda_pick_ada <- lambda_grid_ada[which.min(BIC_path)]
    } else {
      lambda_pick_ada <- lambda_pick  # z.B. negahban/fix weiterreichen
    }
    
    # Snap to grid, falls nötig
    lam_min_ada <- min(lambda_grid_ada)
    lam_max_ada <- max(lambda_grid_ada)
    if (!(lambda_pick_ada >= lam_min_ada && lambda_pick_ada <= lam_max_ada)) {
      idx_closest_ada <- which.min(abs(lambda_grid_ada - lambda_pick_ada))
      lambda_pick_ada_adj <- lambda_grid_ada[idx_closest_ada]
      warning(sprintf(
        "adaptive: lambda_pick_ada=%.6g außerhalb des Grids [%.6g, %.6g]; verwende %.6g.",
        lambda_pick_ada, lam_min_ada, lam_max_ada, lambda_pick_ada_adj
      ))
      lambda_pick_ada <- lambda_pick_ada_adj
    }
    
    beta <- drop(as.matrix(coef(fit_ada, s = lambda_pick_ada)))
    if (is.null(names(beta)) && !is.null(colnames(X))) names(beta) <- colnames(X)
    ok_adapt <- if (!is.null(beta)) 1L else 0L
  } else ok_adapt <- 1L
  
  
  # Support & Refit-Design
  m    <- which(beta != 0)
  bhat <- if (length(m)) beta[m] else numeric(0)
  
  Xs <- if (length(m)) X[, m, drop = FALSE] else NULL
  
  if (is.null(colnames(X)))  colnames(X)  <- paste0("X", seq_len(ncol(X)))
  if (!is.null(Xs) && is.null(colnames(Xs))) {
    colnames(Xs) <- paste0("X", seq_len(ncol(Xs)))
  }
  
  df <- if (!is.null(Xs)) data.frame(time = time, status = status, Xs) else NULL
  no_active <- length(m) == 0L
  if (no_active) {
    message("[analyze_surv_dataset] Kein aktiver Support aus glmnet – FLI, Refit, debiased etc. werden übersprungen.")
  }
  
  ##all_coef  <- paste0("X", seq_len(ncol(X)))
  sel_names <- all_coef[m]
  cols      <- c("bhat","se","lower","upper","p_value")
  
  if (isTRUE(full)) {
    tr <- tryCatch({
      df_full <- data.frame(time = time, status = status, X, check.names = FALSE)
      
      full_fit <- survival::coxph(
        Surv(time, status) ~ ., data = df_full,
        ties = ties, x = TRUE, y = TRUE
      )
      
      sm <- summary(full_fit)
      
      tab <- data.frame(
        bhat    = sm$coef[, "coef"],
        se      = sm$coef[, "se(coef)"],
        lower   = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
        upper   = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
        p_value = sm$coef[, "Pr(>|z|)"],
        row.names = rownames(sm$coef),
        check.names = FALSE
      )
      
      # --- Name-Fix: `` `1` `` -> "X1" usw. ---
      rn <- rownames(tab)
      rn <- sub("^`(.*)`$", "\\1", rn)     # Backticks weg
      if (all(grepl("^[0-9]+$", rn))) {
        rn <- paste0("X", as.integer(rn))  # numerische Namen in X1..Xp
      }
      rownames(tab) <- rn
      
      # --- Direkt in res_full schreiben ---
      common <- intersect(rownames(res_full), rownames(tab))
      
      if (length(common)) {
        res_full[common, cols] <- tab[common, cols]
        # z-Werte, falls gewünscht:
        res_full[common, "z_value"] <- res_full[common, "bhat"] / res_full[common, "se"]
      }
      ok_full  <- 1L
      err_full <- NULL
      full_fit
    }, error = function(e) {
      ok_full  <- 0L
      err_full <- conditionMessage(e)
      NULL
    })
  }
  
  if (isTRUE(oracle)) {
    tryCatch({
      tab <- NULL
      sm  <- NULL
      # wahrer Support
      m_oracle <- which(beta_true != 0)
      
      if (!is.null(names(beta_true))) {
        vars_oracle <- names(beta_true)[m_oracle]
        idx <- match(vars_oracle, all_coef)
        idx <- idx[!is.na(idx)]
        m_oracle <- idx
      }
      
      if (length(m_oracle) == 0L) {
        stop("Oracle support is empty (no nonzero entries in beta_true).")
      }
      
      # Design + evtl. init passend auf Subset
      X_or <- X[, m_oracle, drop = FALSE]
      
      # Fit
      df_or <- data.frame(time = time, status = status, X_or, check.names = FALSE)
      fit_oracle <- survival::coxph(
        Surv(time, status) ~ ., data = df_or, ties = ties, x = TRUE, y = TRUE
      )
      
      sm <- summary(fit_oracle)
      #z  <- qnorm(1 - alpha/2)
      
      tab <- data.frame(
        bhat    = sm$coef[, "coef"],
        se      = sm$coef[, "se(coef)"],
        lower   = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
        upper   = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
        p_value = sm$coef[, "Pr(>|z|)"],
        row.names = rownames(sm$coef)
      )
      
      # Row-Namen auf globale ##all_coef mappen (Reihenfolge im df_or = all_coef[m_oracle])
      rn_target <- all_coef[m_oracle]
      if (nrow(tab) == length(rn_target)) rownames(tab) <- rn_target
      
      common_rows <- intersect(rownames(res_oracle), rownames(tab))
      if (length(common_rows)) {
        res_oracle[common_rows, cols] <- tab[common_rows, cols]
      }
      ok_oracle  <- 1L
      err_oracle <- NULL
    }, error = function(e) {
      ok_oracle  <- 0L
      err_oracle <- conditionMessage(e)
    })
  }
  
  if (isTRUE(split)) {
    tab <- NULL
    sm  <- NULL
    set.seed(seed %||% 123)
    n_tr <- max(1L, floor(splitratio * n))
    idx_train <- sort(sample.int(n, n_tr))
    idx_test  <- setdiff(seq_len(n), idx_train)
    
    X_train      <- X[idx_train, , drop = FALSE]
    time_train   <- time[idx_train]
    status_train <- status[idx_train]
    
    X_test      <- X[idx_test, , drop = FALSE]
    time_test   <- time[idx_test]
    status_test <- status[idx_test]
    
    tryCatch({
      # glmnet-Pfad auf dem Trainingsdatensatz
      fit_train <- glmnet(
        X_train, Surv(time_train, status_train),
        family = "cox", standardize = FALSE, alpha = alpha_glm
      )
      lam_grid_tr <- fit_train$lambda
      lam_min_tr  <- min(lam_grid_tr)
      lam_max_tr  <- max(lam_grid_tr)
      
      lambda_pick_split <- NULL
      
      if (lambda == "min" || lambda == "1se") {
        cvfit_tr <- cv.glmnet(
          X_train, Surv(time_train, status_train),
          family = "cox", nfolds = nfolds,
          standardize = FALSE, alpha = alpha_glm
        )
        lambda_pick_split <- if (lambda == "min") cvfit_tr$lambda.min else cvfit_tr$lambda.1se
        
      } else if (lambda == "negahban") {
        lambda_pick_split <- negahban_like_lambda_cox(
          X_train, time_train, status_train,
          B = 1000, standardize = TRUE
        )
        
      } else if (lambda == "aic" || lambda == "bic") {
        dev_path_tr <- deviance(fit_train, newx = X_train, y = Surv(time_train, status_train))
        df_path_tr  <- fit_train$df
        if (lambda == "aic") {
          AIC_path_tr <- dev_path_tr + 2 * df_path_tr
          lambda_pick_split <- lam_grid_tr[which.min(AIC_path_tr)]
        } else {
          n_events_tr <- sum(status_train)
          BIC_path_tr <- dev_path_tr + log(n_events_tr) * df_path_tr
          lambda_pick_split <- lam_grid_tr[which.min(BIC_path_tr)]
        }
      } else if (lambda == "fix") {
        lambda_pick_split <- lambda_fix
      }
      
      # Snap-to-Grid auch für Split
      if (!(lambda_pick_split >= lam_min_tr && lambda_pick_split <= lam_max_tr)) {
        idx_closest_tr <- which.min(abs(lam_grid_tr - lambda_pick_split))
        lambda_pick_split_adj <- lam_grid_tr[idx_closest_tr]
        warning(sprintf(
          "SPLIT: lambda_pick_split=%.6g außerhalb des glmnet-Grids [%.6g, %.6g]; verwende %.6g.",
          lambda_pick_split, lam_min_tr, lam_max_tr, lambda_pick_split_adj
        ))
        lambda_pick_split <- lambda_pick_split_adj
      }
      
      beta_split_full <- drop(as.matrix(coef(fit_train, s = lambda_pick_split)))
      if (is.null(names(beta_split_full))) names(beta_split_full) <- colnames(X_train)
      
      m_split <- which(beta_split_full != 0)
      
      if (length(m_split) == 0L) {
        ok_split  <- 0L
        err_split <- "No active variables in training fit (all betas are zero)."
      } else {
        sel_names_split <- colnames(X_train)[m_split]
        beta_split      <- beta_split_full
        bhat_split      <- beta_split_full[m_split]
        
        df_split <- data.frame(
          time = time_test,
          status = status_test,
          X_test[, m_split, drop = FALSE],
          check.names = FALSE
        )
        
        split_fit <- survival::coxph(
          Surv(time, status) ~ ., data = df_split,
          ties = ties, x = TRUE, y = TRUE
        )
        
        sm <- summary(split_fit)
        tab <- data.frame(
          bhat    = sm$coef[, "coef"],
          se      = sm$coef[, "se(coef)"],
          lower   = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
          upper   = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
          p_value = sm$coef[, "Pr(>|z|)"],
          row.names = rownames(sm$coef),
          check.names = FALSE
        )
        
        # sicherstellen, dass Namen passen
        sel_names_split <- colnames(X_train)[m_split]
        # ggf. Backticks entfernen
        rn <- rownames(tab)
        rn <- sub("^`(.*)`$", "\\1", rn)
        rownames(tab) <- rn
        
        stopifnot(identical(sort(sel_names_split), sort(rownames(tab))))
        
        cols <- c("bhat","se","lower","upper","p_value")
        keep <- intersect(rownames(tab), rownames(res_split))
        res_split[keep, cols] <- tab[keep, cols]
        
        ok_split  <- 1L
        err_split <- NULL
      }
    }, error = function(e) {
      ok_split  <- 0L
      err_split <- conditionMessage(e)
    })
  }
  
  if(isTRUE(cox)&& length(m) > 0L){
    tr <- tryCatch({
      tab <- NULL
      sm  <- NULL
      df_b <- data.frame(time = time, status = status, Xs, check.names = FALSE)
      refit_fit <- survival::coxph(
        Surv(time, status) ~ ., data = df_b, ties = ties,
        init = bhat, x = TRUE, y = TRUE
      )
      sm <- summary(refit_fit)
      
      bh <- sm$coef[, "coef"]
      se <- sm$coef[, "se(coef)"]
      lower   <- sm$coef[, "coef"] - z * sm$coef[, "se(coef)"]
      upper   <- sm$coef[, "coef"] + z * sm$coef[, "se(coef)"]
      p_value <- sm$coef[, "Pr(>|z|)"]
      
      res_refit[sel_names, "bhat"]    <- bhat
      res_refit[sel_names, "se"]      <- se
      res_refit[sel_names, "lower"]   <- lower
      res_refit[sel_names, "upper"]   <- upper
      res_refit[sel_names, "p_value"] <- p_value
      
      ok_refit  <- 1L
      err_refit <- NULL
      refit_fit
    }, error = function(e) {
      ok_refit  <- 0L
      err_refit <- conditionMessage(e)
      NULL
    })
  }
  
  if(isTRUE(cox0)&& length(m) > 0L){
    tr <- tryCatch({
      tab <- NULL
      sm  <- NULL
      df_b <- data.frame(time = time, status = status, Xs, check.names = FALSE)
      refit0_fit <- survival::coxph(
        Surv(time, status) ~ ., data = df_b, ties = ties, iter.max=0,
        init = bhat, x = TRUE, y = TRUE
      )
      sm <- summary(refit0_fit)
   
      bh <- sm$coef[, "coef"]
      se <- sm$coef[, "se(coef)"]
      bh <- sm$coef[, "coef"]
      se <- sm$coef[, "se(coef)"]
      lower   <- sm$coef[, "coef"] - z * sm$coef[, "se(coef)"]
      upper   <- sm$coef[, "coef"] + z * sm$coef[, "se(coef)"]
      p_value <- sm$coef[, "Pr(>|z|)"]
      
      res_refit0[sel_names, "bhat"]    <- bhat
      res_refit0[sel_names, "se"]      <- se
      res_refit0[sel_names, "lower"]   <- lower
      res_refit0[sel_names, "upper"]   <- upper
      res_refit0[sel_names, "p_value"] <- p_value
      
      ok_refit0  <- 1L
      err_refit0 <- NULL
      refit0_fit
    }, error = function(e) {
      ok_refit0  <- 0L
      err_refit0 <- conditionMessage(e)
      NULL
    })
  }
  
  if (isTRUE(fli)&& length(m) > 0L) {
    fli_fit <- try(
      fixedLassoInf(
        x = X, y = time, beta = beta, lambda = lambda_pick,
        family = "cox", alpha = alpha, status = status,
        type = "partial", bits = 200
      ),
      silent = TRUE
    )
    
    if (inherits(fli_fit, "try-error")) {
      ok_fli  <- 0L
      err_fli <- as.character(attr(fli_fit, "condition") %||% "fixedLassoInf failed")
      kkt_fli <- 0L
      # nichts selektiert / kein Ergebnis -> res_fli bleibt einfach NA
    } else {
      ok_fli  <- 1L
      err_fli <- NULL
      kkt_fli <- tryCatch(.extract_kkt_flag(fli_fit), error = function(e) 0L)
      
      idx <- as.integer(fli_fit$vars)
      if (length(idx) > 0L) {
        sel_names <- all_coef[idx]
        
        bhat <- as.numeric(fli_fit$coef0)
        zval <- as.numeric(fli_fit$zscore0)
        pval <- as.numeric(fli_fit$pv)
        lo   <- as.numeric(fli_fit$ci[, 1])
        up   <- as.numeric(fli_fit$ci[, 2])
        se   <- ifelse(is.finite(zval) & zval != 0, abs(bhat) / abs(zval), NA_real_)
        
        fill_df <- data.frame(
          bhat    = bhat,
          se      = se,
          lower   = lo,
          upper   = up,
          z_value = zval,
          p_value = pval,
          row.names = sel_names
        )
        
        cols <- intersect(colnames(res_fli), colnames(fill_df))
        common <- intersect(sel_names, rownames(res_fli))
        if (length(common)) {
          res_fli[common, cols] <- fill_df[common, cols]
        }
      }
    }
  }
  
  sqrtn <- sqrt(n)
  if (isTRUE(coxwb)&& length(m) > 0L) {

    pm <- length(m); n <- length(time)
   
    se_beta_boot  <- setNames(rep(NA_real_, length(all_coef)), all_coef)
    
    # Nichts selektiert? -> leere Tabellen + OK
    if (p == 0L) {
      ok_wb <- 0L
    } else {
      tryCatch({
        # Refit mit fixierten Koeffizienten (iter.max = 0)
        fit <- survival::coxph(Surv(time, status) ~ X[, m],
                               init = bhat, iter.max = 0,
                               robust = FALSE, x = TRUE)
        U_i <- residuals(fit, type = "score"); if (!is.matrix(U_i)) U_i <- matrix(U_i, ncol = pm)
        Iinv <- vcov(fit); if (is.null(Iinv)) stop("vcov is NULL.")
        
        # Multiplier-Draws 
        M <- matrix(rnorm(B_boot * n), B_boot, n)
        
        # Bootstrap-Deltas und SEs
        R_boot <- (M %*% U_i) / sqrtn  # B x p
        Delta  <- R_boot %*% Iinv        # B x p
        se_sel <- apply(Delta, 2, sd)*sqrtn
       
        # Quantile (abs / std)
        q_abs <- apply(abs(Delta)*sqrtn, 2, quantile, 1 - alpha/2, type = 1)
        Delta_std <- sweep(abs(Delta)*sqrtn, 2, se_sel, "/")
        q_std <- apply(Delta_std, 2, quantile, 1 - alpha/2, type = 1)
        
        # CIs 
        ci_abs_lo <- bhat - q_abs;                 ci_abs_up <- bhat + q_abs
        ci_std_lo <- bhat - q_std * se_sel;        ci_std_up <- bhat + q_std * se_sel
        ci_nor_lo <- bhat - z * se_sel;            ci_nor_up <- bhat + z * se_sel
        p_norm    <- 2 * (1 - pnorm(abs(bhat / se_sel)))
        
        # Full-Tables an selektierten Zeilen befüllen
        res_wb_abs[sel_names, "bhat"]    <- bhat
        res_wb_abs[sel_names, "se"]      <- se_sel
        res_wb_abs[sel_names, "lower"]   <- ci_abs_lo
        res_wb_abs[sel_names, "upper"]   <- ci_abs_up
        res_wb_abs[sel_names, "p_value"] <- p_norm
        
        res_wb_std[sel_names, "bhat"]    <- bhat
        res_wb_std[sel_names, "se"]      <- se_sel
        res_wb_std[sel_names, "lower"]   <- ci_std_lo
        res_wb_std[sel_names, "upper"]   <- ci_std_up
        res_wb_std[sel_names, "p_value"] <- p_norm
        
        res_wb_normal[sel_names, "bhat"]    <- bhat
        res_wb_normal[sel_names, "se"]      <- se_sel
        res_wb_normal[sel_names, "lower"]   <- ci_nor_lo
        res_wb_normal[sel_names, "upper"]   <- ci_nor_up
        res_wb_normal[sel_names, "p_value"] <- p_norm
        
        # SE-Vektor in voller Länge
        names(se_sel) <- colnames(X)[m]
        se_beta_boot[sel_names] <- se_sel
        
        
        ok_wb <- 1L
      }, error = function(e) {
        err_wb <- conditionMessage(e)
      })
    }
    
    # -> res_wb_abs, res_wb_std, res_wb_normal, se_beta_boot, ok_wb, err_wb stehen bereit
  }
  
  if (isTRUE(fliwb)&& length(m) > 0L) {
    
    pm <- length(m)
    
    # falls coxwb=FALSE war, brauchst du hier noch einmal die se_sel aus dem
    # Multiplier-Bootstrap
    if (isFALSE(coxwb)) {
      if (p == 0L) {
        ok_wb <- 0L
      } else {
        tryCatch({
          fit <- survival::coxph(Surv(time, status) ~ X[, m],
                                 init = bhat, iter.max = 0,
                                 robust = FALSE, x = TRUE)
          U_i <- residuals(fit, type = "score")
          if (!is.matrix(U_i)) U_i <- matrix(U_i, ncol = pm)
          Iinv <- vcov(fit); if (is.null(Iinv)) stop("vcov is NULL.")
          
          M <- matrix(rnorm(B_boot * n), B_boot, n)
          R_boot <- (M %*% U_i) / sqrt(n)
          Delta  <- R_boot %*% Iinv
          se_sel <- apply(Delta, 2, sd) * sqrt(n)
          names(se_sel) <- colnames(X)[m]
          se_beta_boot[sel_names] <- se_sel
        }, error = function(e) {
          err_fliwb <- conditionMessage(e)
        })
      }
    }
    
    fliwb_fit <- try(
      fixedLassoInf(
        x = X, y = time, beta = beta, lambda = lambda_pick, sigma = se_sel,
        family = "cox", alpha = alpha, status = status, type = "partial", bits = 200
      ),
      silent = TRUE
    )
    
    if (inherits(fliwb_fit, "try-error")) {
      ok_fli_wb  <- 0L
      err_fli_wb <- as.character(attr(fliwb_fit, "condition") %||% "fixedLassoInf failed")
      kkt_fli_wb <- 0L
    } else {
      ok_fli_wb  <- 1L
      err_fli_wb <- NULL
      kkt_fli_wb <- tryCatch(.extract_kkt_flag(fliwb_fit), error = function(e) 0L)
      
      idx <- as.integer(fliwb_fit$vars)
      if (length(idx) > 0L) {
        sel_names <- all_coef[idx]
        
        bhat <- as.numeric(fliwb_fit$coef0)
        zval <- as.numeric(fliwb_fit$zscore0)
        pval <- as.numeric(fliwb_fit$pv)
        lo   <- as.numeric(fliwb_fit$ci[, 1])
        up   <- as.numeric(fliwb_fit$ci[, 2])
        se   <- ifelse(is.finite(zval) & zval != 0, abs(bhat) / abs(zval), NA_real_)
        
        fill_df <- data.frame(
          bhat    = bhat,
          se      = se,
          lower   = lo,
          upper   = up,
          z_value = zval,
          p_value = pval,
          row.names = sel_names
        )
        
        cols <- intersect(colnames(res_fli_wb), colnames(fill_df))
        common <- intersect(sel_names, rownames(res_fli_wb))
        if (length(common)) {
          res_fli_wb[common, cols] <- fill_df[common, cols]
        }
      }
    }
  }
  
  
  if (isTRUE(randlasso)) {
    if (!requireNamespace("selectiveInference", quietly = TRUE)) {
      ok_randlasso  <- 0L
      err_randlasso <- "Package 'selectiveInference' nicht verfügbar; randomizedLasso wird übersprungen."
    } else {
      tryCatch({
        lam_rl <- lambda_pick
        if (!is.finite(lam_rl) || lam_rl <= 0) {
          stop(sprintf("Ungültiges lambda_pick für randomizedLasso: %g", lam_rl))
        }
        
        fit_cox_rl <- selectiveInference::randomizedLasso(
          X      = X,
          y      = survival::Surv(time, status),
          lam    = lam_rl,
          family = "cox"
        )
        
        if (length(fit_cox_rl$active_set) == 0L) {
          ok_randlasso  <- 1L
          err_randlasso <- NULL
        } else {
          targets_rl <- selectiveInference::compute_target(fit_cox_rl, type = "selected")
          
          inf_rl <- selectiveInference::randomizedLassoInf(
            rand_lasso_soln = fit_cox_rl,
            targets         = targets_rl,
            level           = 1 - alpha,
            sampler         = randlasso_sampler,
            nsample         = randlasso_nsample,
            burnin          = randlasso_burnin
          )
          
          beta_hat_rl <- inf_rl$targets$observed_target
          beta_CI_rl  <- inf_rl$ci
          vars_rl     <- names(beta_hat_rl)
          if (is.null(vars_rl)) {
            vars_rl <- all_coef[fit_cox_rl$active_set]
          }
          
          tab_rl <- data.frame(
            bhat    = as.numeric(beta_hat_rl),
            se      = NA_real_,
            lower   = beta_CI_rl[, 1],
            upper   = beta_CI_rl[, 2],
            z_value = NA_real_,
            p_value = inf_rl$pvalues,
            row.names   = vars_rl,
            check.names = FALSE
          )
          
          vars_norm <- sub("^`(.*)`$", "\\1", rownames(tab_rl))
          idx_direct <- match(vars_norm, all_coef)
          if (all(is.na(idx_direct))) {
            num_part <- sub(".*?([0-9]+)$", "\\1", vars_norm)
            is_num   <- grepl("^[0-9]+$", num_part)
            if (all(is_num)) {
              idx <- as.integer(num_part)
              idx <- pmax(pmin(idx, length(all_coef)), 1L)
              vars_norm <- all_coef[idx]
            }
          } else {
            vars_norm[!is.na(idx_direct)] <- all_coef[idx_direct[!is.na(idx_direct)]]
          }
          rownames(tab_rl) <- vars_norm
          
          if (!exists("res_randlasso", inherits = FALSE)) {
            res_randlasso <- tab_rl
          } else {
            common <- intersect(rownames(res_randlasso), rownames(tab_rl))
            if (length(common)) {
              res_randlasso[common, colnames(tab_rl)] <- tab_rl[common, colnames(tab_rl)]
            } else {
              res_randlasso <- tab_rl
            }
          }
          
          ok_randlasso  <- 1L
          err_randlasso <- NULL
        }
      }, error = function(e) {
        ok_randlasso  <- 0L
        err_randlasso <- conditionMessage(e)
      })
    }
  }
  
  
  if (isTRUE(debiased)) {
    tryCatch({
      # sicherstellen, dass lambda_pick im Pfad ist
      if (!(lambda_pick >= min(fit$lambda, na.rm = TRUE) &&
            lambda_pick <= max(fit$lambda, na.rm = TRUE))) {
        fit <- glmnet(
          X, Surv(time, status), family = "cox", standardize = FALSE,
          lambda = c(lambda_pick, fit$lambda), alpha = alpha_glm
        )
      }
      
      beta_hat <- drop(as.matrix(coef(fit, s = lambda_pick)))
      if (is.null(names(beta_hat))) names(beta_hat) <- colnames(X)
      
      # **NEU**: wenn alles 0 -> kein Debiasing, sauber abbrechen
      if (all(abs(beta_hat) < .Machine$double.eps)) {
        ok_debiased  <- 0L
        err_debiased <- "debiased_cox_from_beta: active set empty (all beta_hat == 0); skipping debiasing."
      } else {
        fit_db <- debiased_cox_from_beta(
          X, time, status,
          beta_hat = beta_hat,
          nodewise_lambda = debiased_nodewise_lambda,
          nfolds = nfolds
        )
        
        bh    <- fit_db$beta_debiased
        se    <- fit_db$se
        ci_lo <- bh - z * se
        ci_up <- bh + z * se
        pval  <- 2 * (1 - pnorm(abs(bh / se)))
        
        tab <- data.frame(
          bhat = bh, se = se, lower = ci_lo, upper = ci_up,
          z_value = bh / se, p_value = pval,
          row.names = names(bh), check.names = FALSE
        )
        rn <- rownames(tab)
        rn <- sub("^`(.*)`$", "\\1", rn)
        if (all(grepl("^[0-9]+$", rn))) {
          idx <- as.integer(rn)
          idx <- pmax(pmin(idx, length(all_coef)), 1L)
          rn  <- all_coef[idx]
        }
        rownames(tab) <- rn
        
        common <- intersect(rownames(res_debiased), rownames(tab))
        if (length(common)) {
          res_debiased[common, colnames(tab)] <- tab[common, colnames(tab)]
        }
        
        ok_debiased  <- 1L
        err_debiased <- NULL
        
        assign("._fit_db_internal", fit_db, inherits = FALSE)
      }
    }, error = function(e) {
      ok_debiased  <- 0L
      err_debiased <- conditionMessage(e)
    })
  }
  
  
  if (isTRUE(debiased) && isTRUE(debiased_wb) && isTRUE(ok_debiased == 1L)) {
    tryCatch({
      # Debiasing-Objekt aus dem vorherigen Block holen (falls vorhanden)
      fit_db <- try(get("._fit_db_internal", inherits = FALSE), silent = TRUE)
      if (inherits(fit_db, "try-error") || is.null(fit_db)) {
        fit_db <- debiased_cox_lasso(
          X = X, time = time, status = status,
          lambda_cox     = match.arg(debiased_lambda),
          nodewise_lambda = debiased_nodewise_lambda,
          nfolds         = nfolds
        )
      }
      
      scheme <- match.arg(debiased_wb_scheme)
      beta_star <- wildbootstrap_debiased_cox(
        fit_db = fit_db,
        X      = X,
        time   = time,
        status = status,
        B      = B_boot_db,
        scheme = scheme
      ) # p x B
      
      bh <- fit_db$beta_debiased
      se <- fit_db$se
      z  <- qnorm(1 - alpha / 2)
      
      ## zentrierte Bootstrap-Züge und studentisierte Züge
      t_star         <- sweep(beta_star, 1, bh, "-")      # p x B
      t_star_student <- sweep(t_star, 1, se, "/")         # p x B
      
      ## 1) Absolutes (percentile) Intervall: |beta* - bh|
      q_abs <- apply(
        abs(t_star), 1,
        quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 1
      )
      ci_abs_lo <- bh - q_abs
      ci_abs_up <- bh + q_abs
      
      ## 2) Studentisiertes Intervall: |(beta* - bh)/se|
      q_std <- apply(
        abs(t_star_student), 1,
        quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 1
      )
      ci_std_lo <- bh - q_std * se
      ci_std_up <- bh + q_std * se
      
      ## 3) Normal/Wald
      ci_nor_lo <- bh - z * se
      ci_nor_up <- bh + z * se
      
      ## p-Werte (Wald)
      p_norm <- 2 * (1 - pnorm(abs(bh / se)))
      
      ## Tabellen bauen
      tab_abs <- data.frame(
        bhat    = bh,
        se      = se,
        lower   = ci_abs_lo,
        upper   = ci_abs_up,
        z_value = bh / se,
        p_value = p_norm,
        row.names   = names(bh),
        check.names = FALSE
      )
      tab_std <- data.frame(
        bhat    = bh,
        se      = se,
        lower   = ci_std_lo,
        upper   = ci_std_up,
        z_value = bh / se,
        p_value = p_norm,
        row.names   = names(bh),
        check.names = FALSE
      )
      tab_nor <- data.frame(
        bhat    = bh,
        se      = se,
        lower   = ci_nor_lo,
        upper   = ci_nor_up,
        z_value = bh / se,
        p_value = p_norm,
        row.names   = names(bh),
        check.names = FALSE
      )
      
      normalize_rn <- function(rn, all_coef) {
        rn <- sub("^`(.*)`$", "\\1", rn)
        if (all(grepl("^[0-9]+$", rn))) {
          idx <- as.integer(rn)
          idx <- pmax(pmin(idx, length(all_coef)), 1L)
          rn  <- all_coef[idx]
        }
        rn
      }
      
      rownames(tab_abs) <- normalize_rn(rownames(tab_abs), all_coef)
      rownames(tab_std) <- normalize_rn(rownames(tab_std), all_coef)
      rownames(tab_nor) <- normalize_rn(rownames(tab_nor), all_coef)
      
      ## Falls die Ergebnis-Container existieren, selektiv updaten;
      ## sonst initialisieren
      if (!exists("res_debiased_wb_abs", inherits = FALSE)) {
        res_debiased_wb_abs <- res_debiased[NULL, , drop = FALSE]
      }
      if (!exists("res_debiased_wb_std", inherits = FALSE)) {
        res_debiased_wb_std <- res_debiased[NULL, , drop = FALSE]
      }
      if (!exists("res_debiased_wb_normal", inherits = FALSE)) {
        res_debiased_wb_normal <- res_debiased[NULL, , drop = FALSE]
      }
      
      common_abs <- intersect(rownames(res_debiased_wb_abs), rownames(tab_abs))
      common_std <- intersect(rownames(res_debiased_wb_std), rownames(tab_std))
      common_nor <- intersect(rownames(res_debiased_wb_normal), rownames(tab_nor))
      
      if (length(common_abs)) {
        res_debiased_wb_abs[common_abs, colnames(tab_abs)] <- tab_abs[common_abs, colnames(tab_abs)]
      }
      if (length(common_std)) {
        res_debiased_wb_std[common_std, colnames(tab_std)] <- tab_std[common_std, colnames(tab_std)]
      }
      if (length(common_nor)) {
        res_debiased_wb_normal[common_nor, colnames(tab_nor)] <- tab_nor[common_nor, colnames(tab_nor)]
      }
      
      ## optional: studentisierte Bootstrap-p-Werte
      thr <- abs(bh / se)  # Länge p
      p_boot_std <- sapply(seq_along(thr), function(j) {
        v <- abs(t_star_student[j, ])
        mean(v >= thr[j], na.rm = TRUE)
      })
      names(p_boot_std) <- names(bh)
      # Wenn du willst:
      # tab_std$p_value <- p_boot_std
      
      ok_debiased_wb  <- 1L
      err_debiased_wb <- NULL
    }, error = function(e) {
      ok_debiased_wb  <- 0L
      err_debiased_wb <- conditionMessage(e)
    })
  }
  
  res_lasso_boot_resid <- NULL
  res_lasso_boot_pairs <- NULL
  
  # sauber benannte Referenzschätzer
  beta_hat_named <- beta
  if (is.null(names(beta_hat_named))) names(beta_hat_named) <- colnames(X)
  
  if (isTRUE(lasso_boot_residual)&& length(m) > 0L) {
    bootC_resid <- .cox_lasso_boot_residual_inline(
      X = X, time = time, status = status,
      lambda_opt = lambda_pick,
      B = B_lasso_boot_resid, a_n = a_n_lasso,
      seed = seed
    )
    # drei CI-Varianten + SE/Bias in Tables
    res_resid_all <- .summarize_lasso_boot_all(
      boot_coefs = bootC_resid,
      beta_hat_named = beta_hat_named,
      alpha = alpha
    )
    res_lasso_boot_resid_abs  <- res_resid_all$df_abs
    res_lasso_boot_resid_std  <- res_resid_all$df_std
    res_lasso_boot_resid_wald <- res_resid_all$df_wald
    # (optional: Backwards-Compat)
    res_lasso_boot_resid <- res_lasso_boot_resid_wald
  }
  
  ok_lasso_boot_loughin  <- NA_integer_
  err_lasso_boot_loughin <- NULL
  
  if (isTRUE(lasso_boot_loughin)) {
    
    tryCatch({
      
      bootC_loughin <- .cox_lasso_boot_loughin_inline(
        X = X, time = time, status = status,
        lambda_opt = lambda_pick,
        B = B_lasso_boot_loughin,
        seed = seed
      )
      
      res_loughin_all <- .summarize_lasso_boot_all(
        boot_coefs = bootC_loughin,
        beta_hat_named = beta_hat_named,
        alpha = alpha
      )
      
      res_lasso_boot_loughin_abs  <- res_loughin_all$df_abs
      res_lasso_boot_loughin_std  <- res_loughin_all$df_std
      res_lasso_boot_loughin_wald <- res_loughin_all$df_wald
      
      # Backwards-compatible Wald
      res_lasso_boot_loughin <- res_lasso_boot_loughin_wald
      
      ok_lasso_boot_loughin  <- 1L
      err_lasso_boot_loughin <- NULL
      
    }, error = function(e) {
      
      ok_lasso_boot_loughin  <- 0L
      err_lasso_boot_loughin <- conditionMessage(e)
    })
  }
  
  
  if (isTRUE(lasso_boot_pairs)&& length(m) > 0L) {
    bootC_pairs <- .cox_lasso_boot_pairs_inline(
      X = X, time = time, status = status,
      lambda_base = lambda_pick,
      B = B_lasso_boot_pairs,
      refit_lambda_each_boot = refit_lambda_each_boot,
      alpha_glm = alpha_glm, nfolds = nfolds,
      seed = seed
    )
    res_pairs_all <- .summarize_lasso_boot_all(
      boot_coefs = bootC_pairs,
      beta_hat_named = beta_hat_named,
      alpha = alpha
    )
    res_lasso_boot_pairs_abs  <- res_pairs_all$df_abs
    res_lasso_boot_pairs_std  <- res_pairs_all$df_std
    res_lasso_boot_pairs_wald <- res_pairs_all$df_wald
    # (optional: Backwards-Compat)
    res_lasso_boot_pairs <- res_lasso_boot_pairs_wald
  }
  res_lasso_boot_param <- NULL
  
  if (isTRUE(lasso_boot_param)&& length(m) > 0L) {
    bootC_param <- .cox_lasso_boot_parametric_inline(
      X = X,
      time = time,
      status = status,
      lambda_base = lambda_pick,
      B = B_lasso_boot_param,
      alpha_glm = alpha_glm,
      nfolds = nfolds,
      seed = seed,
      use_cv_lambda = FALSE  # oder TRUE, wenn du im Bootstrap auch CV willst
    )
    
    res_param_all <- .summarize_lasso_boot_all(
      boot_coefs = bootC_param,
      beta_hat_named = beta_hat_named,
      alpha = alpha
    )
    
    res_lasso_boot_param_abs  <- res_param_all$df_abs
    res_lasso_boot_param_std  <- res_param_all$df_std
    res_lasso_boot_param_wald <- res_param_all$df_wald
    res_lasso_boot_param      <- res_lasso_boot_param_wald
  }
  
  
  list(
    name = name,
    s_vec = s_vec,
    lambda_choice = lambda_pick,
    lambda_std = lambda_pick*n,
    beta  = beta,
    active = m,
    bhat = bhat,
    beta_true     = beta_true,
    beta_split    = beta_split,
    active_split  = m_split,
    bhat_split    = bhat_split,
    res_full      = res_full,
    res_oracle    = res_oracle,
    res_refit     = res_refit,
    res_refit0    = res_refit0,
    res_split     = res_split,
    res_fli       = res_fli,
    res_fli_wb    = res_fli_wb,
    res_wb_abs    = res_wb_abs,
    res_wb_std    = res_wb_std,
    res_wb_normal = res_wb_normal,
    res_debiased   = res_debiased,
    res_debiased_wb = res_debiased_wb,
    res_debiased_wb_abs = res_debiased_wb_abs,
    res_debiased_wb_std = res_debiased_wb_std,
    res_debiased_wb_normal = res_debiased_wb_normal,

    res_lasso_boot_resid            = res_lasso_boot_resid,         # (Wald, kompatibel)
    res_lasso_boot_resid_abs        = res_lasso_boot_resid_abs,
    res_lasso_boot_resid_std        = res_lasso_boot_resid_std,
    res_lasso_boot_resid_wald       = res_lasso_boot_resid_wald,
    
    res_lasso_boot_loughin            = res_lasso_boot_loughin_wald,
    res_lasso_boot_loughin_abs        = res_lasso_boot_loughin_abs,
    res_lasso_boot_loughin_std        = res_lasso_boot_loughin_std,
    res_lasso_boot_loughin_wald       = res_lasso_boot_loughin_wald,
    
    ok_lasso_boot_loughin             = ok_lasso_boot_loughin,
    err_lasso_boot_loughin            = err_lasso_boot_loughin,
    
    res_lasso_boot_pairs            = res_lasso_boot_pairs,         # (Wald, kompatibel)
    res_lasso_boot_pairs_abs        = res_lasso_boot_pairs_abs,
    res_lasso_boot_pairs_std        = res_lasso_boot_pairs_std,
    res_lasso_boot_pairs_wald       = res_lasso_boot_pairs_wald,
    ## NEU:
    res_lasso_boot_param            = res_lasso_boot_param, 
    res_lasso_boot_param_abs        = res_lasso_boot_param_abs,
    res_lasso_boot_param_std        = res_lasso_boot_param_std,
    res_lasso_boot_param_wald       = res_lasso_boot_param_wald,
    
    res_randlasso = res_randlasso,
    
    ok_randlasso  = ok_randlasso,
    
    p_boot_std = p_boot_std,
    se_beta_boot  = se_beta_boot,
    ok_glm        = ok_glm,
    ok_adapt      = ok_adapt,
    ok_full       = ok_full,
    ok_oracle     = ok_oracle,
    ok_refit      = ok_refit,
    ok_refit0     = ok_refit0,
    ok_split      = ok_split,
    ok_fli        = ok_fli,
    ok_fli_wb     = ok_fli_wb,
    ok_wb         = ok_wb,
    ok_debiased_wb = ok_debiased_wb,
    ok_debiased   = ok_debiased,
    #err_glm        = err_glm,
    #err_adapt      = err_adapt,
    #err_full       = err_full,
    #err_oracle     = err_oracle,
    #err_refit      = err_refit,
    #err_refit0     = err_refit0,
    #err_split      = err_split,
    #err_fli        = err_fli,
    #err_fli_wb     = err_fli_wb,
    #err_wb         = err_wb,
    #err_debiasedwb = err_debiasedwb,
    #err_debiased   = err_debiased,
    kkt_fli       = kkt_fli,
    kkt_fli_wb    = kkt_fli_wb
  )
}
