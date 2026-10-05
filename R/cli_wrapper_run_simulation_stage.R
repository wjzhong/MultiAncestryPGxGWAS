#' Run the simulation stage from CLI-style key/value arguments
#'
#' Wrapper entrypoint around `mpgx_simulate_phenotypes()` using
#' `--key=value` arguments for compatibility with script-driven workflows.
#'
#' @param args Character vector of `--key=value` arguments.
#'
#' @return The return value from `mpgx_simulate_phenotypes()`.
#' @export
mpgx_cli_run_simulation_stage <- function(args = commandArgs(trailingOnly = TRUE)) {
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

  args <- .mpgx_parse_kv_args(args)

  required <- c("bfile_prefix", "ancestry_info_file", "output_dir", "result_label", "seed")
  missing_args <- required[!required %in% names(args)]
  if (length(missing_args) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
  }

  out <- mpgx_simulate_phenotypes(
    bfile_prefix = args$bfile_prefix,
    ancestry_info_file = args$ancestry_info_file,
    output_dir = args$output_dir,
    result_label = args$result_label,
    seed = as.integer(args$seed),
    num_pcs = as.integer(.mpgx_arg_or(args$num_pcs, "10")),
    pcs_file = if (!is.null(args$pcs_file) && nzchar(args$pcs_file)) args$pcs_file else NULL,
    variance_proportion = .mpgx_arg_or(args$variance_proportion, "0.25_0.25_0.05_0.05_0.00_0.40"),
    n_snps_maineffects_only = as.integer(.mpgx_arg_or(args$n_snps_maineffects_only, "2")),
    n_snps_interaction_only = as.integer(.mpgx_arg_or(args$n_snps_interaction_only, "2")),
    n_snps_both = as.integer(.mpgx_arg_or(args$n_snps_both, "6")),
    type_of_environment_variable = .mpgx_arg_or(args$type_of_environment_variable, "binaryEnvironmentVariable"),
    type_of_error = .mpgx_arg_or(args$type_of_error, "normal"),
    type_of_phenotype = .mpgx_arg_or(args$type_of_phenotype, "normal")
  )

  cat("Simulation stage complete for result_label", args$result_label, "\n")
  out
}
