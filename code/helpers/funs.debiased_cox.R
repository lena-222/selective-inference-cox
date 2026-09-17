# funs.debiased_cox.R


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
norm_names <- function(nms) {
  nms <- sub("^`(.*)`$", "\\1", nms)
  if (all(grepl("^[0-9]+$", nms))) nms <- paste0("X", as.integer(nms))
  nms
}