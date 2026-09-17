# helpers/ci_metrics.R
# Core selective inference metrics computed from per-simulation "outs" objects.

`%||%` <- function(a, b) if (!is.null(a)) a else b
safe_div <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)

.res_name_map <- list(
  full                   = "res_full",
  oracle                 = "res_oracle",
  refit                  = "res_refit",
  refit0                 = "res_refit0",
  split                  = "res_split",
  exact_posi             = "res_exact_posi",
  debiased               = "res_debiased"
)

normalize_ci_names <- function(ci) {
  nm <- names(ci)
  nm <- sub("(?i)lower|lwr", "lo", nm, perl = TRUE)
  nm <- sub("(?i)upper|upr", "hi", nm, perl = TRUE)
  names(ci) <- nm
  ci
}

pick_named <- function(named_vec, jset) {
  if (!length(named_vec)) return(rep(NA_real_, length(jset)))
  out <- named_vec[as.character(jset)]
  out[is.na(out)] <- NA_real_
  unname(out)
}

agg_bin <- function(success_vec, id_vec) {
  if (!length(success_vec)) return(list(k = numeric(0), n = numeric(0)))
  k <- tapply(success_vec, id_vec, sum, na.rm = TRUE)
  n <- tapply(success_vec, id_vec, function(x) sum(is.finite(x)))
  if (!length(k)) { k <- numeric(0); n <- numeric(0) }
  list(k = k, n = n)
}

# Extract CI limits for those terms that are present for a given method in a given simulation output.
# If the result table has bhat/se but no lower/upper, Wald CIs are constructed at (1-alpha) level.
get_ci_df_present <- function(o, method, terms, alpha = 0.10) {
  tab_name <- .res_name_map[[method]]
  z <- qnorm(1 - alpha / 2)
  
  if (is.null(tab_name)) {
    return(data.frame(term = character(0), lo = numeric(0), hi = numeric(0)))
  }
  
  tab <- o[[as.character(tab_name)]]
  if (is.null(tab)) {
    return(data.frame(term = character(0), lo = numeric(0), hi = numeric(0)))
  }
  
  if (!is.null(tab$var)) rownames(tab) <- as.character(tab$var)
  
  rn <- rownames(tab)
  rn <- sub("^`(.*)`$", "\\1", rn)
  
  if (length(rn) && all(grepl("^V[0-9]+$", rn))) {
    rn <- paste0("X", as.integer(sub("^V", "", rn)))
  }
  nums <- suppressWarnings(as.integer(rn))
  if (length(rn) && all(!is.na(nums))) rn <- paste0("X", nums)
  rownames(tab) <- rn
  
  if (!("bhat" %in% names(tab)) && ("beta_hat" %in% names(tab))) tab$bhat <- tab$beta_hat
  if (!("se"   %in% names(tab)) && ("se_boot" %in% names(tab)))   tab$se   <- tab$se_boot
  
  if (!all(c("lower", "upper") %in% names(tab)) && all(c("bhat", "se") %in% names(tab))) {
    tab$lower <- tab$bhat - z * tab$se
    tab$upper <- tab$bhat + z * tab$se
  }
  
  if (!all(c("lower", "upper") %in% names(tab))) {
    return(data.frame(term = character(0), lo = numeric(0), hi = numeric(0)))
  }
  
  ok <- intersect(terms, rownames(tab))
  if (!length(ok)) {
    return(data.frame(term = character(0), lo = numeric(0), hi = numeric(0)))
  }
  
  data.frame(
    term = ok,
    lo   = suppressWarnings(as.numeric(tab[ok, "lower"])),
    hi   = suppressWarnings(as.numeric(tab[ok, "upper"])),
    row.names = NULL,
    check.names = FALSE
  )
}

# Try to extract a per-simulation C-index for a given method.
# This is intentionally permissive because your result objects differ by method.
extract_cindex <- function(o, method) {
  tab_name <- .res_name_map[[method]]
  if (is.null(tab_name)) return(NA_real_)
  
  tab <- o[[as.character(tab_name)]]
  if (is.null(tab)) return(NA_real_)
  
  if ("cindex" %in% names(tab)) return(as.numeric(tab$cindex[1]))
  if ("Cindex" %in% names(tab)) return(as.numeric(tab$Cindex[1]))
  if ("c_index" %in% names(tab)) return(as.numeric(tab$c_index[1]))
  NA_real_
}

# Variant A: CI-based error metrics (coverage / power / type I error), per term.
compute_ci_error_metrics <- function(out_list, methods, alpha_for_ci = 0.10) {
  known_methods <- intersect(methods, names(.res_name_map))
  out_list <- Filter(Negate(is.null), out_list)
  stopifnot(length(out_list) >= 1, length(known_methods) >= 1)
  
  idx_bt <- which(vapply(out_list, function(z) !is.null(z$beta_true), logical(1)))
  if (!length(idx_bt)) stop("beta_true is missing in all out_list elements.")
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
      
      d <- merge(
        ci,
        data.frame(term = terms, j = seq_along(terms), beta_true = btrue),
        by = "term", sort = FALSE, all.x = TRUE
      )
      
      d$has_ci <- is.finite(d$lo) & is.finite(d$hi)
      if (!any(d$has_ci)) return(NULL)
      
      zero_on_scale <- 0
      d$covered_true <- with(d, has_ci & beta_true >= lo & beta_true <= hi)
      d$reject0      <- with(d, has_ci & !(zero_on_scale >= lo & zero_on_scale <= hi))
      d$width        <- with(d, hi - lo)
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
    
    width_mean <- tapply(D$width, D$j, function(x) mean(x, na.rm = TRUE))
    
    j_all <- sort(unique(c(
      as.integer(names(cov_sig$n)),
      as.integer(names(cov_null$n)),
      as.integer(names(pow_ag$n)),
      as.integer(names(t1_ag$n)),
      as.integer(names(width_mean))
    )))
    if (!length(j_all)) { res[[meth]] <- data.frame(); next }
    
    cov_sig_n <- pick_named(cov_sig$n,  j_all); cov_sig_k <- pick_named(cov_sig$k,  j_all)
    cov_null_n<- pick_named(cov_null$n, j_all); cov_null_k<- pick_named(cov_null$k, j_all)
    pow_n     <- pick_named(pow_ag$n,   j_all); pow_k     <- pick_named(pow_ag$k,   j_all)
    t1_n      <- pick_named(t1_ag$n,    j_all); t1_k      <- pick_named(t1_ag$k,    j_all)
    w_mean    <- pick_named(width_mean, j_all)
    
    is_sig <- btrue[j_all] != 0
    cov_n  <- ifelse(is_sig, cov_sig_n,  cov_null_n)
    cov_k  <- ifelse(is_sig, cov_sig_k,  cov_null_k)
    
    res[[meth]] <- data.frame(
      method             = meth,
      term               = terms[j_all],
      j                  = j_all,
      beta_true          = as.numeric(btrue[j_all]),
      is_signal          = is_sig,
      cov_n              = cov_n,
      cov_k              = cov_k,
      selective_coverage = safe_div(cov_k, cov_n),
      pow_n              = pow_n,
      pow_k              = pow_k,
      selective_power    = safe_div(pow_k, pow_n),
      t1_n               = t1_n,
      t1_k               = t1_k,
      selective_type1    = safe_div(t1_k, t1_n),
      ci_width_mean      = w_mean,
      row.names = NULL,
      check.names = FALSE
    )
  }
  
  do.call(rbind, res)
}

# Variant B: selection-focused summaries (median/IQR CI width, selection probability,
# probability of recovering the exact true model, and mean C-index), per term.
compute_selection_summary_metrics <- function(out_list, methods, alpha_for_ci = 0.10) {
  known_methods <- intersect(methods, names(.res_name_map))
  out_list <- Filter(Negate(is.null), out_list)
  stopifnot(length(out_list) >= 1, length(known_methods) >= 1)
  
  idx_bt <- which(vapply(out_list, function(z) !is.null(z$beta_true), logical(1)))
  if (!length(idx_bt)) stop("beta_true is missing in all out_list elements.")
  btrue <- out_list[[idx_bt[1]]]$beta_true
  terms <- names(btrue)
  if (is.null(terms)) terms <- paste0("X", seq_along(btrue))
  
  true_vars <- terms[as.numeric(btrue) != 0]
  n_rep <- length(out_list)
  
  out_all <- vector("list", length(known_methods))
  names(out_all) <- known_methods
  
  for (meth in known_methods) {
    per_sim <- lapply(seq_along(out_list), function(s) {
      o <- out_list[[s]]
      if (is.null(o)) return(NULL)
      
      ci_df <- get_ci_df_present(o, meth, terms, alpha = alpha_for_ci)
      ci_df <- normalize_ci_names(ci_df)
      
      selected_vars <- if (nrow(ci_df)) ci_df$term else character(0)
      is_true_model <- isTRUE(setequal(selected_vars, true_vars))
      
      c_ind <- extract_cindex(o, meth)
      
      if (!nrow(ci_df)) {
        return(data.frame(
          term = character(0),
          width = numeric(0),
          is_true_model = logical(0),
          c_index = numeric(0),
          sim = integer(0),
          stringsAsFactors = FALSE
        ))
      }
      
      ci_df$width <- ci_df$hi - ci_df$lo
      ci_df$is_true_model <- is_true_model
      ci_df$c_index <- c_ind
      ci_df$sim <- s
      ci_df[, c("term", "width", "is_true_model", "c_index", "sim")]
    })
    
    D <- do.call(rbind, per_sim)
    if (is.null(D) || !nrow(D)) {
      out_all[[meth]] <- data.frame()
      next
    }
    
    # Per-term aggregation:
    # selection_prob = (# times term appears in CI table) / (# repetitions)
    # prob_true_model = P(selected set equals true_vars), conditional on term being selected? (No.)
    # Here: prob_true_model is averaged across simulations where the term was selected (consistent with your earlier code).
    # If you prefer unconditional across all reps, replace mean(...) denominator accordingly.
    res_meth <- suppressMessages(
      dplyr::as_tibble(D) |>
        dplyr::group_by(.data$term) |>
        dplyr::summarise(
          method            = meth,
          beta_true         = as.numeric(btrue[.data$term][1]),
          is_signal         = (beta_true != 0),
          ci_width_median   = median(.data$width, na.rm = TRUE),
          ci_width_iqr      = IQR(.data$width, na.rm = TRUE),
          selection_prob    = dplyr::n() / n_rep,
          prob_true_model   = mean(.data$is_true_model, na.rm = TRUE),
          mean_c_index      = mean(.data$c_index, na.rm = TRUE),
          .groups = "drop"
        )
    )
    
    # Add "j" index so joins are stable and consistent with CI metrics output
    res_meth$j <- match(res_meth$term, terms)
    
    out_all[[meth]] <- as.data.frame(res_meth, check.names = FALSE)
  }
  
  do.call(rbind, out_all)
}

