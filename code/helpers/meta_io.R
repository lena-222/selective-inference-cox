# helpers/meta_io.R
# Utilities for reading meta.txt produced by the simulation runners.

has_ada_token <- function(strings) {
  if (!length(strings)) return(FALSE)
  s <- tolower(paste(strings, collapse = " "))
  grepl("\\bada\\b|adaptive|\\badapt\\b|[_\\.-]ada[_\\.-]|ada-", s, perl = TRUE)
}

# Parse meta.txt written via capture.output(str(meta)) + sessionInfo().
# Only "$key: type value" lines from str(meta) are parsed.
parse_meta_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  
  lines <- readLines(path, warn = FALSE)
  lines_meta <- lines[grepl("^\\s*\\$", lines)]
  if (!length(lines_meta)) return(NULL)
  
  entries <- lapply(lines_meta, function(z) {
    m <- regexec("^\\s*\\$([^:]+):\\s*(.*)$", z)
    r <- regmatches(z, m)[[1]]
    if (length(r) < 3) return(NULL)
    
    name <- trimws(r[2])
    val  <- gsub("\\s+", " ", trimws(r[3]))
    
    parsed <- NA
    if (grepl("^(int|num|dbl)\\s+", val)) {
      parsed <- suppressWarnings(as.numeric(sub("^(int|num|dbl)\\s+", "", val)))
    } else if (grepl("^logi\\s+", val)) {
      lg <- strsplit(sub("^logi\\s+", "", val), " ", fixed = TRUE)[[1]][1]
      parsed <- if (lg %in% c("TRUE", "FALSE")) as.logical(lg) else NA
    } else if (grepl("^chr\\s+\".*\"$", val)) {
      parsed <- sub("^chr\\s+\"(.*)\"$", "\\1", val)
    } else {
      parsed <- sub('^"(.*)"$', "\\1", val)
    }
    setNames(list(parsed), name)
  })
  
  entries <- entries[!vapply(entries, is.null, logical(1))]
  if (!length(entries)) return(NULL)
  
  lst <- list()
  for (e in entries) {
    nm <- names(e)
    if (!nm %in% names(lst)) {
      lst[[nm]] <- e[[1]]
    } else {
      k <- 1L
      newnm <- paste0(nm, ".", k)
      while (newnm %in% names(lst)) {
        k <- k + 1L
        newnm <- paste0(nm, ".", k)
      }
      lst[[newnm]] <- e[[1]]
    }
  }
  
  df <- as.data.frame(lapply(lst, function(x) {
    if (is.logical(x) || is.numeric(x) || is.character(x)) x else as.character(x)
  }), stringsAsFactors = FALSE)
  
  df$meta_path  <- path
  df$run_dir    <- dirname(path)
  df$run_folder <- basename(df$run_dir)
  
  if ("methods" %in% names(df)) {
    methods_raw <- as.character(df$methods[1])
    methods_vec <- strsplit(methods_raw, "\\s*,\\s*")[[1]]
    methods_vec <- methods_vec[nzchar(methods_vec)]
    attr(df, "methods_vec") <- if (length(methods_vec)) methods_vec else NULL
    df$methods <- paste(methods_vec, collapse = ",")
  } else {
    attr(df, "methods_vec") <- NULL
  }
  
  ada_flag <- has_ada_token(c(basename(path), basename(df$run_dir),
                              list.files(df$run_dir, full.names = FALSE)))
  if (!"adaptive" %in% names(df)) {
    df$adaptive <- ada_flag
  } else {
    df$adaptive <- as.logical(df$adaptive)
    if (is.na(df$adaptive)) df$adaptive <- ada_flag
  }
  
  for (nm in c("alpha", "n", "p", "n_sim")) {
    if (nm %in% names(df)) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  
  df
}

