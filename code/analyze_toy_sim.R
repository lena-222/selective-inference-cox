# analyze_toy_sim.R

## Aggregate per-beta selective metrics across simulation replicates
perf_per_beta_split <- function(
    out_list,
    methods = c("full", "oracle", "split", "refit", "refit0", "exact_posi", "debiased"),
    alpha_for_ci = 0.10
) {
  known_methods <- intersect(methods, names(.res_name_map))
  out_list <- Filter(Negate(is.null), out_list)
  stopifnot(length(out_list) >= 1, length(known_methods) >= 1)
  
  idx_bt <- which(vapply(out_list, function(z) !is.null(z$beta_true), logical(1)))
  if (!length(idx_bt)) stop("beta_true missing in all out_list elements.")
  btrue <- out_list[[idx_bt[1]]]$beta_true
  
  terms <- names(btrue)
  if (is.null(terms)) terms <- paste0("X", seq_along(btrue))
  
  res <- vector("list", length(known_methods))
  names(res) <- known_methods
  
  for (meth in known_methods) {
    rec <- lapply(seq_along(out_list), function(s) {
      o <- out_list[[s]]
      if (is.null(o)) return(NULL)
      
      ci <- get_ci_df_present(o, meth, terms, alpha = alpha_for_ci)
      if (!nrow(ci)) return(NULL)
      ci <- normalize_ci_names(ci)
      
      if (!all(c("term","lo","hi") %in% names(ci))) {
        # Skip this replicate for this method if CIs are not available in a usable format
        return(NULL)
      }
      
      
      d <- merge(
        ci,
        data.frame(term = terms, j = seq_along(terms), beta_true = btrue),
        by = "term", sort = FALSE, all.x = TRUE
      )
      
      d$has_ci <- is.finite(d$lo) & is.finite(d$hi)
      if (!any(d$has_ci)) return(NULL)
      
      d$beta_true_aligned <- d$beta_true
      zero_on_scale <- 0
      d$covered_true <- with(d, has_ci & beta_true_aligned >= lo & beta_true_aligned <= hi)
      d$reject0      <- with(d, has_ci & !(zero_on_scale >= lo & zero_on_scale <= hi))
      d$sim <- s
      d
    })
    
    D <- do.call(rbind, rec)
    if (is.null(D) || !nrow(D)) { res[[meth]] <- data.frame(); next }
    
    D <- subset(D, has_ci)
    D_sig  <- subset(D, beta_true != 0)
    D_null <- subset(D, beta_true == 0)
    
    cov_sig  <- agg_bin(D_sig$covered_true,  D_sig$j)
    cov_null <- agg_bin(D_null$covered_true, D_null$j)
    pow_ag   <- agg_bin(D_sig$reject0,       D_sig$j)
    t1_ag    <- agg_bin(D_null$reject0,      D_null$j)
    
    j_all <- sort(unique(c(
      as.integer(names(cov_sig$n)),
      as.integer(names(cov_null$n)),
      as.integer(names(pow_ag$n)),
      as.integer(names(t1_ag$n))
    )))
    
    if (!length(j_all)) { res[[meth]] <- data.frame(); next }
    
    cov_sig_n  <- pick_named(cov_sig$n,  j_all)
    cov_sig_k  <- pick_named(cov_sig$k,  j_all)
    cov_null_n <- pick_named(cov_null$n, j_all)
    cov_null_k <- pick_named(cov_null$k, j_all)
    pow_n      <- pick_named(pow_ag$n,   j_all)
    pow_k      <- pick_named(pow_ag$k,   j_all)
    t1_n       <- pick_named(t1_ag$n,    j_all)
    t1_k       <- pick_named(t1_ag$k,    j_all)
    
    is_sig <- btrue[j_all] != 0
    cov_n  <- ifelse(is_sig, cov_sig_n,  cov_null_n)
    cov_k  <- ifelse(is_sig, cov_sig_k,  cov_null_k)
    
    res[[meth]] <- data.frame(
      method   = meth,
      term     = terms[j_all],
      j        = j_all,
      cov_n    = cov_n,
      cov_k    = cov_k,
      coverage = safe_div(cov_k, cov_n),
      pow_n    = pow_n,
      pow_k    = pow_k,
      power    = safe_div(pow_k, pow_n),
      t1_n     = t1_n,
      t1_k     = t1_k,
      type1    = safe_div(t1_k, t1_n),
      row.names = NULL, check.names = FALSE
    )
  }
  
  do.call(rbind, res)
}
# ===================== main function ====================
#' Analyze a survival dataset with multiple estimators
#'
#' Runs a configurable pipeline for right-censored survival data:
#' optional train/test split, Cox variants, 
#' (adaptive) penalized GLM via cross-validation, and exact_posi variants.
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
#' @param exact_posi If `TRUE`, fit exact_posi model.
#' @param adaptive If `TRUE`, use adaptive penalties where applicable.
#' @param nfolds Number of CV folds (glm-like procedures).
#' @param alpha_glm Elastic-net mixing (1 = lasso, 0 = ridge).
#' @param lambda_choice Rule to pick λ: `"min"`, `"1se"`, `"aic"`,
#'   `"bic"`, or `"negahban"`.
#' @param ties Method for Cox ties (`"efron"`, `"breslow"`, `"exact"`).
#' @param ic_logn Sample size for log-normal IC (`"events"` or `"n"`).
#' @param robust_cox If `TRUE`, use robust (sandwich) variance in Cox.
#' @param splitratio Train fraction if splitting (e.g., 0.5).
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
#' @return A list with fitted objects, CV choices, metrics, 
#' results if requested, and the arguments used.
#'
#' @examples
#' # res <- analyze_surv_dataset(obj, cox = TRUE)



analyze_surv_dataset <- function(
    obj, name = "dataset",
    beta_true = NULL,
    full = FALSE,
    oracle = FALSE,
    split = TRUE,
    cox = TRUE,
    cox0 = TRUE,
    exact_posi = TRUE,
    debiased = TRUE,
    debiased_lambda = c("lambda.min","lambda.1se"),
    debiased_nodewise_lambda = NULL,
    adaptive = FALSE,
    nfolds = 5,
    alpha_glm = 1,
    lambda_choice = c("min","1se","aic","bic","negahban","fix"),
    ties = "efron", ic_logn = c("events","n"),
    robust_cox = FALSE,
    splitratio = 0.5,
    support_tol = 0,
    gamma = 1,
    eps = 1e-6,
    seed = 123,
    alpha = 0.10,
    lambda_fix = NULL
) {
  
  lambda <- match.arg(tolower(lambda_choice), c("min","1se","aic","bic","negahban","fix"))
  
  X_mm <- to_mm(obj$X)
  
  # toy convention: enforce X1..Xp names (prevents name-mismatch everywhere)
  colnames(X_mm) <- paste0("X", seq_len(ncol(X_mm)))
  
  tmp   <- l2_normalize(X_mm, return_scale = TRUE)
  X     <- tmp$X
  s_vec <- tmp$scale
  
  colnames(X) <- colnames(X_mm)
  names(s_vec) <- colnames(X)
  
  time   <- as.numeric(obj$time)
  status <- as.integer(obj$status > 0)
  
  ys <- safe_y_cox(time, status)
  time <- ys$time
  status <- ys$status
  
  n <- nrow(X)
  p <- ncol(X)
  z <- qnorm(1 - alpha/2)
  
  all_coef <- colnames(X)
  
  init_df <- function() {
    out <- data.frame(
      var = all_coef,
      bhat = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      check.names = FALSE
    )
    rownames(out) <- all_coef
    out
  }
  
  all_names <- c("res_full","res_oracle","res_refit","res_refit0","res_split","res_exact_posi","res_debiased")
  models <- setNames(replicate(length(all_names), init_df(), simplify = FALSE), all_names)
  list2env(models, envir = environment())
  rm(models)
  
  ok_full <- ok_oracle <- ok_glm <- ok_adapt <- ok_split <- ok_refit <- ok_refit0 <- ok_debiased <- NA_integer_
  ok_exact_posi <- NA_integer_
  err_full <- err_oracle <- err_glm <- err_adaptive <- err_split <- err_refit <- err_refit0 <- err_exact_posi <- err_debiased <- NULL
  kkt_exact_posi <- 0L
  
  beta_split <- NULL
  bhat_split <- NULL
  m_split <- NULL
  
  fit <- glmnet::glmnet(X, survival::Surv(time, status), family = "cox",
                        standardize = FALSE, alpha = alpha_glm)
  
  lambda_pick <- NA_real_
  
  if (lambda %in% c("min","1se")) {
    cvfit <- glmnet::cv.glmnet(
      X, survival::Surv(time, status),
      family = "cox", nfolds = nfolds,
      standardize = FALSE, alpha = alpha_glm
    )
    lambda_pick <- if (lambda == "min") cvfit$lambda.min else cvfit$lambda.1se
    
  } else if (lambda == "negahban") {
    lambda_pick <- negahban_like_lambda_cox(X, time, status, B = 1000, standardize = TRUE)
    
  } else if (lambda %in% c("aic","bic")) {
    dev_path <- glmnet::deviance(fit, newx = X, y = survival::Surv(time, status))
    df_path  <- fit$df
    lam_path <- fit$lambda
    n_events <- sum(status)
    
    if (lambda == "aic") {
      crit <- dev_path + 2 * df_path
    } else {
      crit <- dev_path + log(n_events) * df_path
    }
    lambda_pick <- lam_path[which.min(crit)]
    
  } else if (lambda == "fix") {
    if (is.numeric(lambda_fix) && length(lambda_fix) == 1L &&
        is.finite(lambda_fix) && lambda_fix > 0) {
      lambda_pick <- lambda_fix
    } else {
      stop("lambda_fix must be a single positive, finite number")
    }
  }
  
  # snap to glmnet grid
  lambda_grid <- fit$lambda
  lam_min <- min(lambda_grid)
  lam_max <- max(lambda_grid)
  if (!(lambda_pick >= lam_min && lambda_pick <= lam_max)) {
    idx_closest <- which.min(abs(lambda_grid - lambda_pick))
    lambda_pick <- lambda_grid[idx_closest]
    warning(sprintf("lambda_pick snapped to closest grid value: %.6g", lambda_pick))
  }
  
  beta <- drop(as.matrix(stats::coef(fit, s = lambda_pick)))
  names(beta) <- all_coef
  ok_glm <- 1L
  
  lambda_pick_ada <- NULL
  if (isTRUE(adaptive)) {
    penalty.factor <- 1 / (abs(beta) + eps)^gamma
    
    fit_ada <- glmnet::glmnet(
      X, survival::Surv(time, status), family = "cox",
      alpha = alpha_glm, standardize = FALSE,
      penalty.factor = penalty.factor
    )
    
    if (lambda %in% c("min","1se")) {
      cvfit_ada <- glmnet::cv.glmnet(
        X, survival::Surv(time, status),
        family = "cox", nfolds = nfolds,
        standardize = FALSE, alpha = alpha_glm,
        penalty.factor = penalty.factor
      )
      lambda_pick_ada <- if (lambda == "min") cvfit_ada$lambda.min else cvfit_ada$lambda.1se
      
    } else if (lambda %in% c("aic","bic")) {
      dev_path <- glmnet::deviance(fit_ada, newx = X, y = survival::Surv(time, status))
      df_path  <- fit_ada$df
      lam_path <- fit_ada$lambda
      n_events <- sum(status)
      
      if (lambda == "aic") crit <- dev_path + 2 * df_path
      else crit <- dev_path + log(n_events) * df_path
      
      lambda_pick_ada <- lam_path[which.min(crit)]
    } else {
      lambda_pick_ada <- lambda_pick
    }
    
    # snap
    g <- fit_ada$lambda
    if (!(lambda_pick_ada >= min(g) && lambda_pick_ada <= max(g))) {
      lambda_pick_ada <- g[which.min(abs(g - lambda_pick_ada))]
      warning(sprintf("adaptive lambda snapped to grid: %.6g", lambda_pick_ada))
    }
    
    beta <- drop(as.matrix(stats::coef(fit_ada, s = lambda_pick_ada)))
    names(beta) <- all_coef
    ok_adapt <- 1L
  } else {
    ok_adapt <- 1L
  }
  
  m <- which(beta != 0)
  bhat <- if (length(m)) beta[m] else numeric(0)
  sel_names <- all_coef[m]
  cols <- c("bhat","se","lower","upper","p_value")
  
  if (length(m) == 0L) {
    return(list(
      name = name,
      s_vec = s_vec,
      lambda_choice = lambda_pick,
      lambda_std = lambda_pick * n,
      beta = beta,
      active = m,
      bhat = bhat,
      beta_true = beta_true,
      beta_split = beta_split,
      active_split = m_split,
      bhat_split = bhat_split,
      res_full = res_full,
      res_oracle = res_oracle,
      res_refit = res_refit,
      res_refit0 = res_refit0,
      res_split = res_split,
      res_exact_posi = res_exact_posi,
      res_debiased = res_debiased,
      ok_glm = ok_glm,
      ok_adapt = ok_adapt,
      ok_full = ok_full,
      ok_oracle = ok_oracle,
      ok_refit = ok_refit,
      ok_refit0 = ok_refit0,
      ok_split = ok_split,
      ok_exact_posi = ok_exact_posi,
      ok_debiased = ok_debiased,
      kkt_exact_posi = kkt_exact_posi
    ))
  }
  
  Xs <- X[, m, drop = FALSE]
  
  if (isTRUE(full)) {
    tryCatch({
      df_full <- data.frame(time = time, status = status, X, check.names = FALSE)
      full_fit <- survival::coxph(survival::Surv(time, status) ~ ., data = df_full,
                                  ties = ties, x = TRUE, y = TRUE)
      sm <- summary(full_fit)
      tab <- data.frame(
        bhat = sm$coef[, "coef"],
        se = sm$coef[, "se(coef)"],
        lower = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
        upper = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
        p_value = sm$coef[, "Pr(>|z|)"],
        row.names = rownames(sm$coef),
        check.names = FALSE
      )
      common <- intersect(rownames(res_full), rownames(tab))
      if (length(common)) {
        res_full[common, cols] <- tab[common, cols]
        res_full[common, "z_value"] <- res_full[common, "bhat"] / res_full[common, "se"]
      }
      ok_full <- 1L
    }, error = function(e) {
      ok_full <<- 0L
      err_full <<- conditionMessage(e)
    })
  }
  
  if (isTRUE(oracle)) {
    tryCatch({
      if (is.null(beta_true)) stop("oracle=TRUE requires beta_true.")
      if (is.null(names(beta_true))) stop("beta_true must be named.")
      bt <- beta_true
      bt <- bt[all_coef]
      bt[is.na(bt)] <- 0
      m_oracle <- which(bt != 0)
      if (length(m_oracle) == 0L) stop("Oracle support empty.")
      
      X_or <- X[, m_oracle, drop = FALSE]
      df_or <- data.frame(time = time, status = status, X_or, check.names = FALSE)
      fit_oracle <- survival::coxph(survival::Surv(time, status) ~ ., data = df_or,
                                    ties = ties, x = TRUE, y = TRUE)
      sm <- summary(fit_oracle)
      tab <- data.frame(
        bhat = sm$coef[, "coef"],
        se = sm$coef[, "se(coef)"],
        lower = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
        upper = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
        p_value = sm$coef[, "Pr(>|z|)"],
        row.names = rownames(sm$coef),
        check.names = FALSE
      )
      
      rn_target <- all_coef[m_oracle]
      if (nrow(tab) == length(rn_target)) rownames(tab) <- rn_target
      common <- intersect(rownames(res_oracle), rownames(tab))
      if (length(common)) res_oracle[common, cols] <- tab[common, cols]
      ok_oracle <- 1L
    }, error = function(e) {
      ok_oracle <<- 0L
      err_oracle <<- conditionMessage(e)
    })
  }
  
  if (isTRUE(split)) {
    tryCatch({
      set.seed(seed %||% 123)
      n_tr <- max(1L, floor(splitratio * n))
      idx_train <- sort(sample.int(n, n_tr))
      idx_test  <- setdiff(seq_len(n), idx_train)
      
      X_train <- X[idx_train, , drop = FALSE]
      time_train <- time[idx_train]
      status_train <- status[idx_train]
      
      X_test <- X[idx_test, , drop = FALSE]
      time_test <- time[idx_test]
      status_test <- status[idx_test]
      
      fit_train <- glmnet::glmnet(X_train, survival::Surv(time_train, status_train),
                                  family = "cox", standardize = FALSE, alpha = alpha_glm)
      
      lambda_pick_split <- NA_real_
      
      if (lambda %in% c("min","1se")) {
        cvfit_tr <- glmnet::cv.glmnet(
          X_train, survival::Surv(time_train, status_train),
          family = "cox", nfolds = nfolds,
          standardize = FALSE, alpha = alpha_glm
        )
        lambda_pick_split <- if (lambda == "min") cvfit_tr$lambda.min else cvfit_tr$lambda.1se
        
      } else if (lambda == "negahban") {
        lambda_pick_split <- negahban_like_lambda_cox(X_train, time_train, status_train,
                                                      B = 1000, standardize = TRUE)
        
      } else if (lambda %in% c("aic","bic")) {
        dev_path <- glmnet::deviance(fit_train, newx = X_train, y = survival::Surv(time_train, status_train))
        df_path <- fit_train$df
        lam_path <- fit_train$lambda
        n_events_tr <- sum(status_train)
        if (lambda == "aic") crit <- dev_path + 2 * df_path
        else crit <- dev_path + log(n_events_tr) * df_path
        lambda_pick_split <- lam_path[which.min(crit)]
        
      } else if (lambda == "fix") {
        lambda_pick_split <- lambda_pick
      }
      
      g <- fit_train$lambda
      if (!(lambda_pick_split >= min(g) && lambda_pick_split <= max(g))) {
        lambda_pick_split <- g[which.min(abs(g - lambda_pick_split))]
        warning(sprintf("split lambda snapped to grid: %.6g", lambda_pick_split))
      }
      
      beta_split_full <- drop(as.matrix(stats::coef(fit_train, s = lambda_pick_split)))
      names(beta_split_full) <- colnames(X_train)
      
      m_split <- which(beta_split_full != 0)
      
      if (length(m_split) == 0L) {
        ok_split <- 0L
        err_split <- "No active variables in training fit."
      } else {
        beta_split <- beta_split_full
        bhat_split <- beta_split_full[m_split]
        
        df_split <- data.frame(
          time = time_test,
          status = status_test,
          X_test[, m_split, drop = FALSE],
          check.names = FALSE
        )
        split_fit <- survival::coxph(survival::Surv(time, status) ~ ., data = df_split,
                                     ties = ties, x = TRUE, y = TRUE)
        sm <- summary(split_fit)
        tab <- data.frame(
          bhat = sm$coef[, "coef"],
          se = sm$coef[, "se(coef)"],
          lower = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
          upper = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
          p_value = sm$coef[, "Pr(>|z|)"],
          row.names = rownames(sm$coef),
          check.names = FALSE
        )
        keep <- intersect(rownames(tab), rownames(res_split))
        if (length(keep)) res_split[keep, cols] <- tab[keep, cols]
        ok_split <- 1L
      }
    }, error = function(e) {
      ok_split <<- 0L
      err_split <<- conditionMessage(e)
    })
  }
  
  if (isTRUE(cox)) {
    tryCatch({
      df_b <- data.frame(time = time, status = status, Xs, check.names = FALSE)
      refit_fit <- survival::coxph(
        survival::Surv(time, status) ~ ., data = df_b, ties = ties,
        init = bhat, x = TRUE, y = TRUE
      )
      sm <- summary(refit_fit)
      tab <- data.frame(
        bhat = sm$coef[, "coef"],
        se = sm$coef[, "se(coef)"],
        lower = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
        upper = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
        p_value = sm$coef[, "Pr(>|z|)"],
        row.names = rownames(sm$coef),
        check.names = FALSE
      )
      keep <- intersect(rownames(tab), rownames(res_refit))
      if (length(keep)) res_refit[keep, cols] <- tab[keep, cols]
      ok_refit <- 1L
    }, error = function(e) {
      ok_refit <<- 0L
      err_refit <<- conditionMessage(e)
    })
  }
  
  if (isTRUE(cox0)) {
    tryCatch({
      df_b <- data.frame(time = time, status = status, Xs, check.names = FALSE)
      refit0_fit <- survival::coxph(
        survival::Surv(time, status) ~ ., data = df_b, ties = ties,
        iter.max = 0, init = bhat, x = TRUE, y = TRUE
      )
      sm <- summary(refit0_fit)
      tab <- data.frame(
        bhat = sm$coef[, "coef"],
        se = sm$coef[, "se(coef)"],
        lower = sm$coef[, "coef"] - z * sm$coef[, "se(coef)"],
        upper = sm$coef[, "coef"] + z * sm$coef[, "se(coef)"],
        p_value = sm$coef[, "Pr(>|z|)"],
        row.names = rownames(sm$coef),
        check.names = FALSE
      )
      keep <- intersect(rownames(tab), rownames(res_refit0))
      if (length(keep)) res_refit0[keep, cols] <- tab[keep, cols]
      ok_refit0 <- 1L
    }, error = function(e) {
      ok_refit0 <<- 0L
      err_refit0 <<- conditionMessage(e)
    })
  }
  
  if (isTRUE(exact_posi)) {
    
    exact_posi_fit <- try(
      fixedLassoInf(
        x = X, y = time, beta = beta, lambda = lambda_pick,
        family = "cox", alpha = alpha, status = status, type = "partial", bits = 200
      ),
      silent = TRUE
    )
    
    ok_exact_posi <- if (!inherits(exact_posi_fit, "try-error")) 1L else 0L
    if (ok_exact_posi == 1L) {
      kkt_exact_posi <- tryCatch(.extract_kkt_flag(exact_posi_fit), error = function(e) 0L)
      idx <- as.integer(exact_posi_fit$vars)
      sel_names_exact_posi <- all_coef[idx]
      
      bh <- as.numeric(exact_posi_fit$coef0)
      zval <- as.numeric(exact_posi_fit$zscore0)
      pval <- as.numeric(exact_posi_fit$pv)
      lo <- as.numeric(exact_posi_fit$ci[,1])
      up <- as.numeric(exact_posi_fit$ci[,2])
      se <- ifelse(is.finite(zval) & zval != 0, abs(bh) / abs(zval), NA_real_)
      
      fill_df <- data.frame(
        bhat = bh, se = se, lower = lo, upper = up,
        z_value = zval, p_value = pval,
        row.names = sel_names_exact_posi,
        check.names = FALSE
      )
      
      common <- intersect(rownames(res_exact_posi), rownames(fill_df))
      if (length(common)) res_exact_posi[common, intersect(colnames(res_exact_posi), colnames(fill_df))] <-
        fill_df[common, intersect(colnames(res_exact_posi), colnames(fill_df))]
    } else {
      err_exact_posi <- as.character(attr(exact_posi_fit, "condition") %||% "fixedLassoInf failed")
      kkt_exact_posi <- 0L
    }
  }
  
  if (isTRUE(debiased)) {
    tryCatch({
      beta_hat <- drop(as.matrix(stats::coef(fit, s = lambda_pick)))
      names(beta_hat) <- all_coef
      
      if (all(abs(beta_hat) < .Machine$double.eps)) {
        ok_debiased <- 0L
        err_debiased <- "debiased: beta_hat all zero; skipping."
      } else {
        fit_db <- debiased_cox_from_beta(
          X, time, status,
          beta_hat = beta_hat,
          nodewise_lambda = debiased_nodewise_lambda,
          nfolds = nfolds
        )
        
        bh <- fit_db$beta_debiased
        se <- fit_db$se
        ci_lo <- bh - z * se
        ci_up <- bh + z * se
        pval <- 2 * (1 - pnorm(abs(bh / se)))
        
        tab <- data.frame(
          bhat = bh, se = se, lower = ci_lo, upper = ci_up,
          z_value = bh / se, p_value = pval,
          row.names = names(bh),
          check.names = FALSE
        )
        
        common <- intersect(rownames(res_debiased), rownames(tab))
        if (length(common)) res_debiased[common, colnames(tab)] <- tab[common, colnames(tab)]
        ok_debiased <- 1L
      }
    }, error = function(e) {
      ok_debiased <<- 0L
      err_debiased <<- conditionMessage(e)
    })
  }
  
  list(
    name = name,
    s_vec = s_vec,
    lambda_choice = lambda_pick,
    lambda_std = lambda_pick * n,
    beta = beta,
    active = m,
    bhat = bhat,
    beta_true = beta_true,
    beta_split = beta_split,
    active_split = m_split,
    bhat_split = bhat_split,
    res_full = res_full,
    res_oracle = res_oracle,
    res_refit = res_refit,
    res_refit0 = res_refit0,
    res_split = res_split,
    res_exact_posi = res_exact_posi,
    res_debiased = res_debiased,
    ok_glm = ok_glm,
    ok_adapt = ok_adapt,
    ok_full = ok_full,
    ok_oracle = ok_oracle,
    ok_refit = ok_refit,
    ok_refit0 = ok_refit0,
    ok_split = ok_split,
    ok_exact_posi = ok_exact_posi,
    ok_debiased = ok_debiased,
    kkt_exact_posi = kkt_exact_posi,
    err_full = err_full,
    err_oracle = err_oracle,
    err_split = err_split,
    err_refit = err_refit,
    err_refit0 = err_refit0,
    err_exact_posi = err_exact_posi,
    err_debiased = err_debiased
  )
}
