#' Locate a toy_data example file
#'
#' @param file File name under `inst/extdata`.
#'
#' @return Absolute file path.
#' @export
mpgx_example_file <- function(file) {
  mpgx_assert_scalar_character(file, "file")
  if (!nzchar(file)) {
    stop("file must be a non-empty character string.", call. = FALSE)
  }

  # If only a basename is provided, search common extdata subfolders first.
  has_dir_component <- grepl("[/\\\\]", file)
  candidates <- if (has_dir_component) {
    file
  } else {
    c(
      file,
      file.path("genotype_data", file),
      file.path("genotype_artifacts", file),
      file.path("simulated_data", file),
      file.path("model_results", file)
    )
  }
  candidates <- unique(candidates)

  # Installed package path.
  for (cand in candidates) {
    path <- system.file("extdata", cand, package = "MultiAncestryPGxGWAS")
    if (nzchar(path) && file.exists(path)) {
      return(path)
    }
  }

  # Development/source-tree fallback.
  rel_roots <- c(".", "..", "../..", "../../..")
  for (root in rel_roots) {
    for (cand in candidates) {
      cand_path <- normalizePath(
        file.path(root, "inst", "extdata", cand),
        winslash = "/",
        mustWork = FALSE
      )
      if (file.exists(cand_path)) {
        return(cand_path)
      }
    }
  }

  stop(sprintf("toy_data example file not found: %s", file), call. = FALSE)
}

#' Stage toy_data PLINK example data
#'
#' Copies toy_data example files into a writable directory.
#'
#' @param output_dir Directory to receive copied files.
#' @param prefix Prefix for output PLINK files.
#' @param overwrite Whether to overwrite existing files.
#' @param pcs_output_dir Optional output directory for `{prefix}.top_pcs.tsv`.
#'
#' @return A list with paths for staged files.
#' @export
mpgx_stage_example_data <- function(
  output_dir,
  prefix = "toyData",
  overwrite = FALSE,
  pcs_output_dir = NULL
) {
  mpgx_assert_scalar_character(output_dir, "output_dir")
  mpgx_assert_scalar_character(prefix, "prefix")

  output_dir <- mpgx_ensure_dir(output_dir)
  if (is.null(pcs_output_dir) || !nzchar(pcs_output_dir)) {
    pcs_output_dir <- file.path(output_dir, "genotype_artifacts")
  }
  genotype_artifacts_dir <- mpgx_ensure_dir(pcs_output_dir)

  resolve_first_existing <- function(candidates) {
    for (cand in candidates) {
      path <- tryCatch(
        mpgx_example_file(cand),
        error = function(e) ""
      )
      if (nzchar(path) && file.exists(path)) {
        return(path)
      }
    }
    stop(
      sprintf(
        "toy_data example file not found. Tried: %s",
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  core_file_candidates <- list(
    "toyData.bed" = c("genotype_data/toyData.bed", "toyData.bed"),
    "toyData.bim" = c("genotype_data/toyData.bim", "toyData.bim"),
    "toyData.fam" = c("genotype_data/toyData.fam", "toyData.fam"),
    "toyData.ancestry.txt" = c("genotype_data/toyData.ancestry.txt", "toyData.ancestry.txt")
  )
  core_src_paths <- vapply(
    core_file_candidates,
    resolve_first_existing,
    FUN.VALUE = character(1)
  )

  core_target_files <- sub("^toyData", prefix, names(core_file_candidates))
  core_dst_paths <- file.path(output_dir, core_target_files)

  pcs_src_path <- resolve_first_existing(c(
    "genotype_artifacts/toyData.top_pcs.tsv",
    "toyData.top_pcs.tsv",
    "toyData_top_pcs.tsv"
  ))
  pcs_dst_path <- file.path(genotype_artifacts_dir, paste0(basename(prefix), ".top_pcs.tsv"))
  dst_paths <- c(core_dst_paths, pcs_dst_path)

  if (!overwrite && any(file.exists(dst_paths))) {
    existing <- dst_paths[file.exists(dst_paths)]
    stop(
      sprintf(
        "Refusing to overwrite existing files in %s: %s",
        output_dir,
        paste(basename(existing), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  copied_core <- file.copy(core_src_paths, core_dst_paths, overwrite = overwrite)
  copied_pcs <- file.copy(pcs_src_path, pcs_dst_path, overwrite = overwrite)
  if (!all(c(copied_core, copied_pcs))) {
    stop("Failed to copy one or more toy_data example files.", call. = FALSE)
  }

  prefix_path <- file.path(output_dir, prefix)
  list(
    prefix = prefix_path,
    ancestry_info_file = paste0(prefix_path, ".ancestry.txt"),
    pcs_file = pcs_dst_path,
    files = dst_paths
  )
}
