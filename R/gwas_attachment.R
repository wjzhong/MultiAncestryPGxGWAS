#' Resolve the first existing GWAS output file path
#'
#' @param candidates Character vector of candidate file paths.
#' @param label Label used in error messages.
#'
#' @return A character scalar path to the first existing GWAS output file.
#' @export
resolve_gwas_file <- function(candidates, label) {
  candidates <- unique(as.character(candidates))
  existing <- candidates[file.exists(candidates)]

  if (!length(existing)) {
    stop(
      "Unable to locate GWAS output for ",
      label,
      ". Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }

  existing[[1]]
}

mpgx_set_attach_metadata <- function(
  result_obj,
  method,
  files,
  label,
  ancestry_labels = NULL,
  output_file = NULL
) {
  info <- list(
    method = method,
    files = files,
    label = label,
    ancestry_labels = ancestry_labels,
    output_file = output_file
  )

  if (is.data.frame(result_obj)) {
    attr(result_obj, "mpgx_attach") <- info
    return(result_obj)
  }

  if (is.list(result_obj)) {
    result_obj$mpgx_attach <- info
    return(result_obj)
  }

  list(raw_result = result_obj, mpgx_attach = info)
}

mpgx_get_attach_metadata <- function(result_obj) {
  if (is.list(result_obj) && !is.data.frame(result_obj)) {
    return(result_obj$mpgx_attach)
  }

  attr(result_obj, "mpgx_attach", exact = TRUE)
}

mpgx_get_attached_gwas_df <- function(result_obj) {
  if (is.list(result_obj) && !is.data.frame(result_obj)) {
    return(result_obj$gwas_df)
  }

  attr(result_obj, "gwas_df", exact = TRUE)
}

#' Extract parsed GWAS data frame from a method result object
#'
#' @param result_obj Method result object.
#' @param label Label used in error messages.
#'
#' @return Parsed GWAS data frame.
#' @export
mpgx_get_method_gwas_df <- function(result_obj, label = "method") {
  if (is.list(result_obj) && !is.data.frame(result_obj) && !is.null(result_obj$gwas_df)) {
    return(result_obj$gwas_df)
  }

  gwas_df_attr <- attr(result_obj, "gwas_df", exact = TRUE)
  if (!is.null(gwas_df_attr)) {
    return(gwas_df_attr)
  }

  if (is.data.frame(result_obj)) {
    return(result_obj)
  }

  stop("No attached gwas_df found for ", label, call. = FALSE)
}

#' Parse and attach GWAS output to a method result object
#'
#' @param result_obj Method result object.
#' @param method Method name understood by [mpgx_parse_method_result_df()].
#' @param files GWAS output file path(s) for `method`.
#' @param label Label used in error messages.
#' @param ancestry_labels Optional ancestry labels for meta methods.
#' @param output_file Optional parsed-output file written by
#'   [mpgx_parse_method_result_df()].
#'
#' @return The input `result_obj` with parsed GWAS results attached.
#'   For list outputs, fields `gwas_file` and `gwas_df` are added.
#'   For data frame outputs, attributes `gwas_file` and `gwas_df` are added.
#' @export
attach_gwas_df <- function(
  result_obj,
  method,
  files,
  label,
  ancestry_labels = NULL,
  output_file = NULL
) {
  mpgx_assert_scalar_character(method, "method")

  if (is.null(label) || !is.character(label) || length(label) != 1 || !nzchar(label)) {
    label <- method
  }

  if (is.character(files)) {
    files <- as.character(files)
    missing_files <- files[!file.exists(files)]
  } else if (is.list(files)) {
    file_values <- as.character(unname(unlist(files, use.names = FALSE)))
    missing_files <- file_values[!file.exists(file_values)]
  } else {
    stop("files must be a character vector or list for ", label, call. = FALSE)
  }

  if (length(missing_files)) {
    stop(
      "Unable to locate GWAS output for ",
      label,
      ": ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  gwas_file <- if (is.null(output_file)) files else output_file
  gwas_df <- mpgx_parse_method_result_df(
    method = method,
    files = files,
    ancestry_labels = ancestry_labels,
    output_file = output_file
  )

  result_obj <- mpgx_set_attach_metadata(
    result_obj = result_obj,
    method = method,
    files = files,
    label = label,
    ancestry_labels = ancestry_labels,
    output_file = output_file
  )

  if (is.data.frame(result_obj)) {
    attr(result_obj, "gwas_file") <- gwas_file
    attr(result_obj, "gwas_df") <- gwas_df
    return(result_obj)
  }

  if (!is.list(result_obj)) {
    result_obj <- list(raw_result = result_obj)
  }

  result_obj$gwas_file <- gwas_file
  result_obj$gwas_df <- gwas_df
  result_obj
}

#' Extract one GWAS file path from each result object
#'
#' @param result_list List of attached result objects.
#'
#' @return Character vector of GWAS file paths.
#' @export
extract_gwas_files <- function(result_list) {
  if (!is.list(result_list)) {
    stop("result_list must be a list.", call. = FALSE)
  }

  vapply(result_list, function(x) {
    gwas_file <- if (is.list(x) && !is.data.frame(x)) x$gwas_file else attr(x, "gwas_file", exact = TRUE)

    if (is.null(gwas_file)) {
      stop("Each element in result_list must contain an attached gwas_file.", call. = FALSE)
    }

    gwas_file <- as.character(unname(unlist(gwas_file, use.names = FALSE)))
    if (length(gwas_file) != 1L) {
      stop("Each element in result_list must resolve to exactly one gwas_file.", call. = FALSE)
    }

    gwas_file
  }, character(1))
}

mpgx_is_multi_file_method <- function(method) {
  method %in% c(
    "meta_StandardLRall",
    "meta_LRallSandwich",
    "meta_StandardLRunrelated",
    "meta_LRunrelatedSandwich",
    "REmeta_LRallSandwich",
    "MRMEGA_LRallSandwich",
    "MRMEGA_GEM"
  )
}

#' Attach GWAS results using stored metadata
#'
#' @param result_obj Method result object.
#' @param method Optional method override.
#' @param files Optional GWAS file path override.
#' @param label Optional label override used in error messages.
#' @param ancestry_labels Optional ancestry labels for meta methods.
#' @param output_file Optional output file path for parsed output.
#' @param force If `TRUE`, force re-attachment even when GWAS data is already
#'   attached.
#'
#' @return The input `result_obj` with attached GWAS outputs.
#' @export
attach_gwas_res <- function(
  result_obj,
  method = NULL,
  files = NULL,
  label = NULL,
  ancestry_labels = NULL,
  output_file = NULL,
  force = FALSE
) {
  if (!isTRUE(force) && !is.null(mpgx_get_attached_gwas_df(result_obj))) {
    return(result_obj)
  }

  attach_info <- mpgx_get_attach_metadata(result_obj)

  method <- method %||% attach_info$method
  files <- files %||% attach_info$files
  label <- label %||% attach_info$label %||% method
  ancestry_labels <- ancestry_labels %||% attach_info$ancestry_labels
  output_file <- output_file %||% attach_info$output_file

  if (is.null(method)) {
    stop(
      "Unable to determine parsing method. Pass method explicitly or run a wrapper that stores attach metadata.",
      call. = FALSE
    )
  }
  if (is.null(files)) {
    stop(
      "Unable to determine GWAS output files. Pass files explicitly or run a wrapper that stores attach metadata.",
      call. = FALSE
    )
  }

  if (is.character(files) && length(files) > 1L && !mpgx_is_multi_file_method(method)) {
    files <- resolve_gwas_file(candidates = files, label = label)
  }

  attach_gwas_df(
    result_obj = result_obj,
    method = method,
    files = files,
    label = label,
    ancestry_labels = ancestry_labels,
    output_file = output_file
  )
}
