# lambda_rules.R

options(encoding = "UTF-8")

# Negahban-like lambda for Cox via multiplier bootstrap of martingale residuals.
# Returns a single scalar lambda compatible with glmnet scaling.
negahban_like_lambda_cox <- function(X, time, status, B = 1000, standardize = TRUE, seed = 1) {
  stopifnot(is.matrix(X))
  if (standardize) X <- scale(X)
  
  n_events <- sum(status)
  if (!is.finite(n_events) || n_events <= 0) stop("No events: cannot compute Negahban-like lambda.")
  
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1)
  mres <- residuals(fit0, type = "martingale")
  
  set.seed(seed)
  G <- matrix(rnorm(length(mres) * B), nrow = length(mres), ncol = B)
  WG <- mres * G
  XT_WG <- t(X) %*% WG
  norms <- apply(XT_WG, 2, function(v) max(abs(v)))
  E_norm <- mean(norms)
  
  as.numeric(2 * E_norm / n_events)
}
