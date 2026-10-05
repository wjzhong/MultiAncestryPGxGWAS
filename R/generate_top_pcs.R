#' Generate top principal components for a PLINK dataset
#'
#' Computes principal components with flashpca after LD-pruning with PLINK2.
#' The output table includes sample identifiers (`FID`, `IID`) and the
#' requested principal components (`PC1`, `PC2`, ...).
#'
#' @param plink_prefix PLINK prefix (without extension) for `.bed/.bim/.fam`.
#' @param configs Optional named list of executable overrides such as
#'   `plink2.path` and `flashpca.path`.
#' @param n_pcs_target Number of principal components to generate.
#' @param output_dir Optional output directory. Defaults to
#'   `file.path(dirname(plink_prefix), "genotype_artifacts")`.
#'
#' @return Character scalar path to `{basename(plink_prefix)}.top_pcs.tsv`.
#' @export
mpgx_generate_top_pcs <- function(
  plink_prefix,
  configs = NULL,
  n_pcs_target = 20L,
  output_dir = NULL
) {
  mpgx_assert_scalar_character(plink_prefix, "plink_prefix")
  mpgx_assert_positive_integer(n_pcs_target, "n_pcs_target")

  bed_file <- paste0(plink_prefix, ".bed")
  bim_file <- paste0(plink_prefix, ".bim")
  fam_file <- paste0(plink_prefix, ".fam")
  mpgx_assert_file(bed_file, "plink_prefix (.bed)")
  mpgx_assert_file(bim_file, "plink_prefix (.bim)")
  mpgx_assert_file(fam_file, "plink_prefix (.fam)")

  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- file.path(dirname(plink_prefix), "genotype_artifacts")
  }
  output_dir <- mpgx_ensure_dir(output_dir)

  plink2 <- mpgx_resolve_executable(
    tool = "plink2",
    default = "plink2",
    configs = configs,
    config_key = "plink2.path",
    install_hint = "Install PLINK 2 or pass configs[['plink2.path']] explicitly."
  )

  flashpca <- mpgx_resolve_executable(
    tool = "flashpca",
    default = "flashpca",
    configs = configs,
    config_key = "flashpca.path",
    install_hint = "Install flashpca or pass configs[['flashpca.path']] explicitly."
  )

  fam_all <- utils::read.table(fam_file, stringsAsFactors = FALSE)

  prune_prefix <- paste0(plink_prefix, ".prune")
  pruned_prefix <- paste0(plink_prefix, ".pruned")
  prune_in_file <- paste0(prune_prefix, ".prune.in")
  eigenvalues_file <- paste0(plink_prefix, ".eigenvalues")
  pve_file <- paste0(plink_prefix, ".pve")
  eigenvectors_file <- paste0(plink_prefix, ".eigenvectors")
  pcs_raw_file <- paste0(plink_prefix, ".pcs")

  mpgx_run_command(
    command = plink2,
    args = c(
      "--bfile", plink_prefix,
      "--indep-pairwise", "1000", "50", "0.05",
      "--out", prune_prefix
    )
  )

  mpgx_run_command(
    command = plink2,
    args = c(
      "--bfile", plink_prefix,
      "--extract", prune_in_file,
      "--make-bed",
      "--out", pruned_prefix
    )
  )

  mpgx_run_command(
    command = flashpca,
    args = c(
      "--bfile", pruned_prefix,
      "-d", as.character(n_pcs_target),
      "--numthreads", "8",
      "--outval", eigenvalues_file,
      "--outpve", pve_file,
      "--outvec", eigenvectors_file,
      "--outpc", pcs_raw_file
    )
  )

  pcs_raw <- utils::read.table(
    pcs_raw_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_with_ids <- 2L + n_pcs_target
  if (ncol(pcs_raw) >= required_with_ids) {
    pcs_df <- pcs_raw[, seq_len(required_with_ids), drop = FALSE]
    colnames(pcs_df) <- c("FID", "IID", paste0("PC", seq_len(n_pcs_target)))
  } else if (ncol(pcs_raw) >= n_pcs_target) {
    if (nrow(pcs_raw) == nrow(fam_all) + 1L) {
      pcs_raw <- pcs_raw[-1, , drop = FALSE]
    }
    if (nrow(pcs_raw) != nrow(fam_all)) {
      stop(
        sprintf(
          "flashpca output %s has %d rows; expected %d.",
          pcs_raw_file,
          nrow(pcs_raw),
          nrow(fam_all)
        ),
        call. = FALSE
      )
    }
    pcs_df <- data.frame(
      FID = fam_all$V1,
      IID = fam_all$V2,
      pcs_raw[, seq_len(n_pcs_target), drop = FALSE],
      stringsAsFactors = FALSE
    )
    colnames(pcs_df) <- c("FID", "IID", paste0("PC", seq_len(n_pcs_target)))
  } else {
    stop(
      sprintf(
        "flashpca output %s has %d columns; expected at least %d.",
        pcs_raw_file,
        ncol(pcs_raw),
        n_pcs_target
      ),
      call. = FALSE
    )
  }

  ids_aligned <- nrow(pcs_df) == nrow(fam_all) && all(pcs_df$FID == fam_all$V1 & pcs_df$IID == fam_all$V2)
  if (!ids_aligned) {
    pcs_df <- merge(
      fam_all[, c("V1", "V2"), drop = FALSE],
      pcs_df,
      by.x = c("V1", "V2"),
      by.y = c("FID", "IID"),
      all.x = TRUE,
      sort = FALSE
    )
    colnames(pcs_df)[1:2] <- c("FID", "IID")
  }

  pcs_df <- pcs_df[, c("FID", "IID", paste0("PC", seq_len(n_pcs_target))), drop = FALSE]
  if (anyNA(pcs_df)) {
    stop("Some samples are missing PCs after flashpca processing.", call. = FALSE)
  }

  pcs_file <- file.path(output_dir, paste0(basename(plink_prefix), ".top_pcs.tsv"))
  mpgx_write_table(pcs_df, file = pcs_file, col_names = TRUE, sep = "\t")
  pcs_file
}
