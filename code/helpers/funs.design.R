# funs.design.R
options(encoding = "UTF-8")

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
    M <- stats::model.matrix(~ . - 1, data = df, na.action = stats::na.pass)
  }
  
  storage.mode(M) <- "double"
  attr(M, "scaled:center") <- NULL
  attr(M, "scaled:scale")  <- NULL
  M
}

l2_normalize <- function(X, return_scale = TRUE) {
  n <- nrow(X)
  s_vec <- sqrt(colSums(X^2)) / sqrt(n)
  s_vec[!is.finite(s_vec) | s_vec == 0] <- 1
  Xn <- sweep(X, 2, s_vec, "/")
  if (return_scale) return(list(X = Xn, scale = s_vec))
  Xn
}


#normalize_ci_names <- function(nms) {
#  nms <- sub("^`(.*)`$", "\\1", nms)
#  nms
#}
normalize_ci_names <- function(ci) {
  if (is.null(ci) || !is.data.frame(ci) || !nrow(ci)) return(ci)
  
  # If 'term' is missing but rownames exist, promote rownames to 'term'
  if (!("term" %in% names(ci))) {
    if ("var" %in% names(ci)) {
      ci$term <- as.character(ci$var)
    } else if (!is.null(rownames(ci)) && any(nzchar(rownames(ci)))) {
      ci$term <- rownames(ci)
    }
  }
  
  # Standardize CI bound column names to 'lo'/'hi'
  nms <- names(ci)
  nms <- sub("(?i)^lower$|^lwr$|^lo$|^ci_lo$|^ci_lower$", "lo", nms, perl = TRUE)
  nms <- sub("(?i)^upper$|^upr$|^hi$|^ci_hi$|^ci_upper$", "hi", nms, perl = TRUE)
  names(ci) <- nms
  
  ci
}


# Normalize coefficient names used across pipelines.
# Removes backticks, trims whitespace, optionally maps pure numeric names -> X<id>
# possibly redundant
norm_names <- function(nms) {
  nms <- sub("^`(.*)`$", "\\1", nms)
  nms <- trimws(nms)
  if (all(grepl("^[0-9]+$", nms))) nms <- paste0("X", as.integer(nms))
  nms
}
