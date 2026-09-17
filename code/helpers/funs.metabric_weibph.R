# funs.metabric_weibph.R
# ---------------------------------------------------------
# Lädt METABRIC, erzeugt Weibull-PH Modell m_os_weibph
# auf einer festen Kovariatenliste (wie im Plot-Script),
# und stellt bereit:
#   - sim_metabric_weibph(n, target_censoring)
#   - get_metabric_beta_true()   # Betas passend zur stabilisierten Dummy-Designmatrix (ohne Intercept)
#
# WICHTIG:
# - Wir fitten flexsurv::flexsurvreg(weibullPH) NICHT direkt auf Faktoren,
#   sondern auf eine explizite Dummy-Designmatrix.
# - Zur Stabilität droppen wir:
#     * konstante Spalten
#     * sehr seltene 0/1-Dummies (min_ones)
#   und skalieren alle Spalten (center/scale).
# - Die für Fit verwendeten Spalten + Skalierung werden global gespeichert:
#     metabric_mm_cols, metabric_mm_center, metabric_mm_scale, metabric_min_ones
#   und in sim_metabric_weibph() identisch angewendet.
# ---------------------------------------------------------

# Pakete
if (!requireNamespace("dplyr", quietly = TRUE))
  stop("Package 'dplyr' is required but not installed.")
if (!requireNamespace("janitor", quietly = TRUE))
  stop("Package 'janitor' is required but not installed.")
if (!requireNamespace("flexsurv", quietly = TRUE))
  stop("Package 'flexsurv' is required but not installed.")
if (!requireNamespace("tidyr", quietly = TRUE))
  stop("Package 'tidyr' is required but not installed.")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is required but not installed.")

# ---------------------------------------------------------
# zentrale Kovariatenliste (wie in deinem Plot-Script)
# ---------------------------------------------------------
.metabric_covars <- c(
  "age_at_diagnosis",
  "tumor_stage_f",
  "er_status_f",
  "her2_status_f",
  "pr_status",
  "chemotherapy",
  "hormone_therapy",
  "radio_therapy",
  "neoplasm_histologic_grade",
  "tumor_size",
  "lymph_nodes_examined_positive",
  "nottingham_prognostic_index"
)

# Designmatrix-Helfer: passt zur analyze_surv_dataset()-Konvention
# (model.matrix(~ ., X_raw)[, -1])
get_metabric_Xmm <- function(X_raw) {
  mm <- stats::model.matrix(~ ., data = X_raw)
  mm <- mm[, -1, drop = FALSE]  # Intercept weg
  mm
}

# ---------------------------------------------------------
# Modell aufbauen (einmalig beim Sourcen)
# ---------------------------------------------------------
if (!exists("m_os_weibph")) {
  
  message("[metabric] Lade Daten und fitte Weibull-PH Modell ...")
  
  project_root <- if (exists("ROOT", inherits = TRUE)) {
    get("ROOT", inherits = TRUE)
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
  fp <- Sys.getenv(
    "METABRIC_CSV",
    unset = file.path(project_root, "data_raw", "cbioportal", "Breast Cancer METABRIC.csv")
  )
  if (!file.exists(fp)) {
    stop(sprintf(
      paste0(
        "[metabric] METABRIC CSV not found at '%s'.\n",
        "Set METABRIC_CSV to the downloaded clinical CSV or place it in ",
        "data_raw/cbioportal/Breast Cancer METABRIC.csv.\n",
        "See the repository README for data instructions."
      ),
      fp
    ))
  }
  
  dt <- read.csv(fp, stringsAsFactors = TRUE, check.names = FALSE)
  dt <- janitor::clean_names(dt)
  
  # Preprocessing (wie im Plot-Script)
  dt_os <- dt %>%
    dplyr::filter(!is.na(overall_survival_months)) %>%
    dplyr::mutate(
      os_time  = overall_survival_months,
      os_event = ifelse(grepl("^Die", patients_vital_status), 1L, 0L),
      
      tumor_stage_f = factor(tumor_stage, levels = sort(unique(tumor_stage))),
      
      er_status_clean = dplyr::na_if(er_status, ""),
      er_status_f = dplyr::case_when(
        er_status_clean == "Positive" ~ "ER_pos",
        er_status_clean == "Negative" ~ "ER_neg",
        TRUE ~ NA_character_
      ) |> factor(),
      
      her2_status_f = dplyr::case_when(
        her2_status %in% c("Positive", "Pos") ~ "HER2_pos",
        her2_status %in% c("Negative", "Neg") ~ "HER2_neg",
        TRUE ~ NA_character_
      ) |> factor()
    ) %>%
    dplyr::filter(os_time > 0)
  
  # Modell-Daten: exakt diese Kovariaten + Outcome
  dt_model <- dt_os %>%
    dplyr::select(os_time, os_event, dplyr::all_of(.metabric_covars)) %>%
    tidyr::drop_na()
  
  if (nrow(dt_model) == 0L) {
    stop("[metabric] dt_model hat 0 Zeilen nach drop_na().")
  }
  
  # global verfügbar machen
  metabric_dt_model <<- dt_model
  
  # --- Designmatrix explizit (ohne Intercept) --------------------------------
  X_raw0 <- dt_model[, .metabric_covars, drop = FALSE]
  X_raw0 <- droplevels(X_raw0)
  
  mm0 <- get_metabric_Xmm(X_raw0)
  
  # 1) konstante Spalten raus
  mm0 <- mm0[, apply(mm0, 2, function(z) stats::sd(z, na.rm = TRUE) > 0), drop = FALSE]
  
  # 2) seltene 0/1-Dummies raus (wichtig für LL/optim Stabilität)
  metabric_min_ones <<- 25L
  is_01 <- apply(mm0, 2, function(z) all(z %in% c(0, 1)))
  ones  <- colSums(mm0 != 0)
  keep  <- (!is_01) | (ones >= metabric_min_ones)
  
  if (!all(keep)) {
    message(sprintf("[metabric] Dropping %d rare dummy cols (<%d ones).",
                    sum(!keep), metabric_min_ones))
  }
  mm0 <- mm0[, keep, drop = FALSE]
  
  if (ncol(mm0) < 2L) {
    stop(sprintf("[metabric] Nach Stabilisierung bleiben zu wenige Spalten: p=%d", ncol(mm0)))
  }
  
  # 3) Zentrieren/Skalieren (damit lp nicht explodiert -> LL nicht definierbar)
  mm0_sc <- scale(mm0)
  mm0_sc <- as.matrix(mm0_sc)
  
  # Center/Scale + Spaltennamen speichern (für Simulation identisch!)
  metabric_mm_cols   <<- colnames(mm0_sc)
  metabric_mm_center <<- attr(mm0_sc, "scaled:center")
  metabric_mm_scale  <<- attr(mm0_sc, "scaled:scale")
  
  dt_fit <- data.frame(
    os_time  = dt_model$os_time,
    os_event = dt_model$os_event,
    mm0_sc,
    check.names = FALSE
  )
  
  # --- Formel ohne "." (damit nix implizit expandiert) ------------------------
  rhs <- paste(colnames(mm0_sc), collapse = " + ")
  fml <- stats::as.formula(paste0("survival::Surv(os_time, os_event) ~ ", rhs))
  
  # --- stabile Startwerte (NATÜRLICHE Skala!) --------------------------------
  init_shape <- 1.0
  init_scale <- stats::median(dt_fit$os_time, na.rm = TRUE)
  if (!is.finite(init_scale) || init_scale <= 0) init_scale <- 1.0
  init_beta  <- rep(0, ncol(mm0_sc))
  
  # Parameterreihenfolge bei flexsurv weibullPH: shape, scale, dann Betas
  inits0 <- c(shape = init_shape, scale = init_scale, init_beta)
  
  # --- Fit mit mehreren Startwerten probieren --------------------------------
  try_inits <- list(
    inits0,
    c(shape = 0.7, scale = init_scale, init_beta),
    c(shape = 1.5, scale = init_scale, init_beta),
    c(shape = 1.0, scale = stats::quantile(dt_fit$os_time, 0.7, na.rm = TRUE), init_beta),
    c(shape = 1.0, scale = stats::quantile(dt_fit$os_time, 0.3, na.rm = TRUE), init_beta)
  )
  
  m_os_weibph <- NULL
  last_err <- NULL
  
  for (k in seq_along(try_inits)) {
    ini <- try_inits[[k]]
    if (!is.finite(ini["shape"]) || ini["shape"] <= 0) next
    if (!is.finite(ini["scale"]) || ini["scale"] <= 0) next
    
    m_os_weibph <- tryCatch(
      flexsurv::flexsurvreg(
        formula = fml,
        data    = dt_fit,
        dist    = "weibullPH",
        inits   = ini,
        method  = "BFGS"
      ),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(m_os_weibph)) break
    
    m_os_weibph <- tryCatch(
      flexsurv::flexsurvreg(
        formula = fml,
        data    = dt_fit,
        dist    = "weibullPH",
        inits   = ini,
        method  = "Nelder-Mead"
      ),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(m_os_weibph)) break
  }
  
  if (is.null(m_os_weibph)) {
    stop("[metabric] flexsurvreg failed for all inits. Last error: ",
         if (!is.null(last_err)) conditionMessage(last_err) else "unknown")
  }
  
  message(sprintf("[metabric] Modell m_os_weibph erfolgreich erzeugt. p=%d (stabilisiert+skaliert)",
                  ncol(mm0_sc)))
  
  # Sanity: stimmen Beta-Namen mit stabilisierter Designmatrix überein?
  cf <- stats::coef(m_os_weibph)
  beta <- cf[!(names(cf) %in% c("shape", "scale"))]
  miss_in_mm <- setdiff(names(beta), metabric_mm_cols)
  if (length(miss_in_mm)) {
    warning(sprintf(
      "[metabric] Beta-Namen nicht in stabilisierter Designmatrix gefunden (n=%d): %s",
      length(miss_in_mm),
      paste(miss_in_mm, collapse = ", ")
    ))
  }
}

# ---------------------------------------------------------
# "Wahre" Betas (log-HR) passend zur stabilisierten Dummy-Designmatrix
# (ohne shape/scale; Namen = metabric_mm_cols)
# ---------------------------------------------------------
get_metabric_beta_true <- function() {
  if (!exists("m_os_weibph"))
    stop("Model m_os_weibph not found! Source funs.metabric_weibph.R first.")
  cf <- stats::coef(m_os_weibph)
  beta <- cf[!(names(cf) %in% c("shape","scale"))]
  beta
}

# ---------------------------------------------------------
# Simulation: ziehe n neue Beobachtungen aus dem Modell
# -> verwendet exakt die beim Fit definierte Dummy-Auswahl + Skalierung
# ---------------------------------------------------------
sim_metabric_weibph <- function(
    n,
    target_censoring = NULL  # optional: 0, 0.1, 0.3
) {
  if (!exists("m_os_weibph"))
    stop("Model m_os_weibph not found!")
  
  if (!exists("metabric_mm_cols", inherits = TRUE) ||
      !exists("metabric_mm_center", inherits = TRUE) ||
      !exists("metabric_mm_scale", inherits = TRUE)) {
    stop("metabric_mm_* not found. Fit must have been created successfully before simulation.")
  }
  
  n_raw <- n
  n <- suppressWarnings(as.integer(n)[1L])
  if (is.na(n) || n <= 0) {
    stop(sprintf("sim_metabric_weibph: n muss positive ganze Zahl sein, bekommen: %s",
                 paste(n_raw, collapse = ",")))
  }
  
  dat <- if (exists("metabric_dt_model", inherits = TRUE)) {
    get("metabric_dt_model", inherits = TRUE)
  } else {
    m_os_weibph$data
  }
  
  Ndat <- tryCatch(nrow(dat), error = function(e) NA_integer_)
  if (is.na(Ndat) || Ndat < 1L) {
    stop("sim_metabric_weibph: Datenobjekt hat keine/undefinierte Zeilen.")
  }
  
  # Parameter (flexsurv speichert shape/scale auf log-Skala)
  cf    <- stats::coef(m_os_weibph)
  shape <- exp(cf["shape"])
  scale <- exp(cf["scale"])
  beta  <- cf[!(names(cf) %in% c("shape", "scale"))]
  
  weibull_par <- c(shape.shape = shape, scale.scale = scale)
  
  # Patienten ziehen
  idx    <- sample.int(Ndat, size = n, replace = TRUE)
  newdat <- dat[idx, , drop = FALSE]
  
  # X_raw wie analyze_surv_dataset erwartet (Rohvariablen)
  X_raw <- newdat[, .metabric_covars, drop = FALSE]
  X_raw <- droplevels(X_raw)
  
  # Designmatrix (ohne Intercept)
  X_mm <- get_metabric_Xmm(X_raw)
  
  # Auf Fit-Spalten bringen: fehlende Spalten als 0 ergänzen
  common <- intersect(metabric_mm_cols, colnames(X_mm))
  X_mm <- X_mm[, common, drop = FALSE]
  
  miss <- setdiff(metabric_mm_cols, colnames(X_mm))
  if (length(miss)) {
    X_mm <- cbind(
      X_mm,
      matrix(0, nrow = nrow(X_mm), ncol = length(miss),
             dimnames = list(NULL, miss))
    )
  }
  X_mm <- X_mm[, metabric_mm_cols, drop = FALSE]
  
  # gleiche Skalierung wie im Fit
  X_mm_sc <- sweep(X_mm, 2, metabric_mm_center[metabric_mm_cols], "-")
  X_mm_sc <- sweep(X_mm_sc, 2, metabric_mm_scale[metabric_mm_cols], "/")
  
  # beta in der richtigen Reihenfolge
  beta_use <- beta[metabric_mm_cols]
  if (anyNA(beta_use)) {
    stop("sim_metabric_weibph: beta_use hat NA (Namen passen nicht zu metabric_mm_cols).")
  }
  
  lp <- as.numeric(X_mm_sc %*% beta_use)
  
  # Ereigniszeiten
  u          <- stats::runif(n)
  time_event <- scale * (-log(u) / exp(lp))^(1 / shape)
  
  # Zensierung (administrativ) wie vorher
  if (is.null(target_censoring)) {
    time   <- time_event
    status <- rep.int(1L, n)
    censoring_target   <- NA_real_
    censoring_observed <- 0
    admin_time         <- Inf
  } else {
    target_censoring_char <- match.arg(
      as.character(target_censoring),
      c("0", "0.1", "0.3")
    )
    censoring_target <- as.numeric(target_censoring_char)
    
    if (censoring_target == 0) {
      time   <- time_event
      status <- rep.int(1L, n)
      admin_time         <- Inf
      censoring_observed <- 0
    } else {
      admin_time <- stats::quantile(
        time_event,
        probs = 1 - censoring_target,
        na.rm = TRUE,
        type = 7
      )
      time   <- pmin(time_event, admin_time)
      status <- as.integer(time_event <= admin_time)
      censoring_observed <- mean(status == 0L)
    }
  }
  
  list(
    time                = time,
    status              = status,
    X                   = X_raw,       # Roh-X (Faktoren/numerisch)
    time_true           = time_event,
    beta_true           = beta_use,     # "Truth" auf stabilisierter+skalierter Designmatrix
    weibull_par         = weibull_par,
    censoring_target    = censoring_target,
    censoring_observed  = censoring_observed,
    admin_time          = admin_time
  )
}

apply_target_censoring <- function(sim_obj, target_censoring = c(0, 0.1, 0.3)) {
  if (is.null(sim_obj$time_true))
    stop("sim_obj muss ein Element 'time_true' enthalten.")
  
  target_censoring_char <- match.arg(
    as.character(target_censoring),
    c("0", "0.1", "0.3")
  )
  censoring_target <- as.numeric(target_censoring_char)
  
  time_event <- sim_obj$time_true
  n <- length(time_event)
  
  if (censoring_target == 0) {
    time   <- time_event
    status <- rep.int(1L, n)
    admin_time         <- Inf
    censoring_observed <- 0
  } else {
    admin_time <- stats::quantile(
      time_event,
      probs = 1 - censoring_target,
      na.rm = TRUE,
      type = 7
    )
    time   <- pmin(time_event, admin_time)
    status <- as.integer(time_event <= admin_time)
    censoring_observed <- mean(status == 0L)
  }
  
  sim_obj$time               <- time
  sim_obj$status             <- status
  sim_obj$censoring_target   <- censoring_target
  sim_obj$censoring_observed <- censoring_observed
  sim_obj$admin_time         <- admin_time
  sim_obj
}
