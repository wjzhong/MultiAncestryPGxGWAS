#' Run the genotype stage from CLI-style key/value arguments
#'
#' Wrapper entrypoint around `mpgx_prepare_example_genotype()` using
#' `--key=value` arguments for compatibility with script-driven workflows.
#'
#' @param args Character vector of `--key=value` arguments.
#'
#' @return The return value from `mpgx_prepare_example_genotype()`.
#' @export
mpgx_cli_run_genotype_stage <- function(args = commandArgs(trailingOnly = TRUE)) {
  .mpgx_parse_kv_args <- function(args) {
    out <- list()
    if (length(args) == 0) {
      return(out)
    }
    for (x in args) {
      if (!grepl("^--", x)) {
        next
      }
      kv <- strsplit(sub("^--", "", x), "=", fixed = TRUE)[[1]]
      key <- kv[[1]]
      val <- if (length(kv) >= 2) paste(kv[-1], collapse = "=") else ""
      out[[key]] <- val
    }
    out
  }

  .mpgx_arg_or <- function(x, default) {
    if (is.null(x) || !nzchar(x)) {
      return(default)
    }
    x
  }

  .mpgx_parse_bool <- function(x, default = FALSE) {
    if (is.null(x) || !nzchar(x)) {
      return(default)
    }
    x <- tolower(x)
    if (x %in% c("true", "t", "1", "yes", "y")) {
      return(TRUE)
    }
    if (x %in% c("false", "f", "0", "no", "n")) {
      return(FALSE)
    }
    stop(sprintf("Invalid logical value: %s", x), call. = FALSE)
  }

  args <- .mpgx_parse_kv_args(args)

  required <- c("output_dir")
  missing_args <- required[!required %in% names(args)]
  if (length(missing_args) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
  }

  config_keys <- c("plink1.9.path", "plink.path", "plink2.path", "flashpca.path")
  present_config_keys <- intersect(config_keys, names(args))
  configs <- if (length(present_config_keys) == 0L) {
    NULL
  } else {
    vals <- lapply(present_config_keys, function(k) args[[k]])
    names(vals) <- present_config_keys
    vals
  }

  out <- mpgx_prepare_example_genotype(
    output_dir = args$output_dir,
    prefix = .mpgx_arg_or(args$prefix, "toyData"),
    source = .mpgx_arg_or(args$source, "toy_data"),
    configs = configs,
    n_per_ancestry = as.integer(.mpgx_arg_or(args$n_per_ancestry, "50")),
    n_snps = as.integer(.mpgx_arg_or(args$n_snps, "1000")),
    overwrite = .mpgx_parse_bool(args$overwrite, default = FALSE)
  )

  cat("Genotype stage complete for prefix", out$prefix, "\n")
  out
}
