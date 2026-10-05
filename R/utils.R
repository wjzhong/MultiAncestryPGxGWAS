# Utility helpers used across workflow stages.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    return(y)
  }
  if (is.character(x) && length(x) == 1 && !nzchar(x)) {
    return(y)
  }
  x
}

mpgx_ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

mpgx_assert_file <- function(path, arg_name = "path") {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", arg_name, path), call. = FALSE)
  }
  invisible(path)
}

mpgx_assert_scalar_character <- function(x, arg_name) {
  if (!is.character(x) || length(x) != 1 || !nzchar(x)) {
    stop(sprintf("%s must be a non-empty character scalar.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

mpgx_assert_positive_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("%s must be a positive integer.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

mpgx_assert_nonnegative_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x != as.integer(x)) {
    stop(sprintf("%s must be a non-negative integer.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

mpgx_normalize_prefix <- function(prefix) {
  mpgx_assert_scalar_character(prefix, "prefix")
  sub("\\.(bed|bim|fam)$", "", prefix)
}

mpgx_plink_paths <- function(prefix) {
  p <- mpgx_normalize_prefix(prefix)
  list(
    prefix = p,
    bed = paste0(p, ".bed"),
    bim = paste0(p, ".bim"),
    fam = paste0(p, ".fam")
  )
}

mpgx_default_tool_configs <- function() {
  list(
    plink1.9.path = "plink",
    plink2.path = "plink2",
    flashpca.path = "flashpca_x86-64", 
    gcta.path = "gcta", 
    lemma.path = "lemma_1_0_4", 
    gem.path = "GEM_1.4.5_static", 
    bgenix.path = "bgenix", 
    mrmega.path = "MR-MEGA",
    mpirun.path = "mpirun" 
  )
}

mpgx_get_config_value <- function(configs, key) {
  if (!is.null(configs) && !is.list(configs)) {
    stop("configs must be a list.", call. = FALSE)
  }

  defaults <- mpgx_default_tool_configs()
  configs <- configs %||% list()
  keys <- as.character(key)

  config_match <- keys[keys %in% names(configs)]
  if (length(config_match) == 0L) {
    config_match <- NULL
  } else {
    config_match <- config_match[[1]]
  }

  default_match <- keys[keys %in% names(defaults)]
  if (length(default_match) == 0L) {
    default_match <- NULL
  } else {
    default_match <- default_match[[1]]
  }

  match <- config_match %||% default_match

  if (is.null(match) || is.na(match) || !nzchar(match)) {
    return(NULL)
  }

  value <- if (match %in% names(configs)) configs[[match]] else defaults[[match]]
  if (is.null(value) || length(value) == 0) {
    return(NULL)
  }

  mpgx_assert_scalar_character(value, sprintf("configs[['%s']]", match))
  if (!nzchar(value)) {
    return(NULL)
  }

  value
}

mpgx_find_executable <- function(candidate) {
  mpgx_assert_scalar_character(candidate, "executable candidate")

  has_sep <- grepl("[/\\]", candidate)
  if (has_sep) {
    resolved <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    ok <- file.exists(resolved)
  } else {
    resolved <- Sys.which(candidate)
    ok <- nzchar(resolved)
  }

  if (!ok) {
    return("")
  }

  resolved
}

mpgx_resolve_executable <- function(
  tool,
  default = tool,
  explicit_path = NULL,
  configs = NULL,
  config_key = NULL,
  install_hint = NULL
) {
  candidate <- explicit_path %||% mpgx_get_config_value(configs, config_key) %||% default
  mpgx_require_executable(tool = tool, path = candidate, install_hint = install_hint)
}

mpgx_require_executable <- function(tool, path = NULL, install_hint = NULL) {
  mpgx_assert_scalar_character(tool, "tool")

  candidate <- path %||% tool
  resolved <- mpgx_find_executable(candidate)
  ok <- nzchar(resolved)

  if (!ok) {
    msg <- sprintf("Executable for '%s' was not found. Pass an explicit path or install it.", tool)
    if (!is.null(install_hint) && nzchar(install_hint)) {
      msg <- paste(msg, install_hint)
    }
    stop(msg, call. = FALSE)
  }

  resolved
}

mpgx_run_command <- function(command, args = character(), stdout = "", stderr = "") {
  mpgx_assert_scalar_character(command, "command")
  args <- as.character(args)
  cli <- paste(c(shQuote(command), shQuote(args)), collapse = " ")

  status <- system2(command = command, args = args, stdout = stdout, stderr = stderr)
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed with status %s: %s", status, cli), call. = FALSE)
  }

  list(command = command, args = args, cli = cli, status = status)
}

mpgx_safe_scale <- function(x, center = TRUE, scale = TRUE) {
  x <- as.numeric(x)
  if (length(unique(stats::na.omit(x))) <= 1) {
    return(rep(0, length(x)))
  }
  as.numeric(base::scale(x, center = center, scale = scale))
}

mpgx_write_table <- function(x, file, col_names = FALSE, sep = "\t") {
  utils::write.table(
    x,
    file = file,
    quote = FALSE,
    row.names = FALSE,
    col.names = col_names,
    sep = sep
  )
}

mpgx_read_delim_table <- function(file, has_header = TRUE, arg_name = "file") {
  if (!is.logical(has_header) || length(has_header) != 1L || is.na(has_header)) {
    stop("has_header must be a non-missing logical scalar.", call. = FALSE)
  }

  mpgx_assert_scalar_character(file, arg_name)
  mpgx_assert_file(file, arg_name)

  out <- readr::read_table(
    file = file,
    col_names = has_header,
    show_col_types = FALSE,
    progress = FALSE
  )

  as.data.frame(out, stringsAsFactors = FALSE)
}
