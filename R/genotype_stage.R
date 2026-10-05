#' Prepare example genotype data in PLINK format
#'
#' For local tutorials, this function stages toy_data PLINK example files.
#' Optionally, it can generate dummy PLINK files with `plink --dummy`.
#'
#' @param output_dir Output directory for genotype files.
#' @param prefix Output file prefix.
#' @param source Either `"toy_data"` or `"plink_dummy"`.
#' @param configs Optional named list of executable overrides such as
#'   `plink1.9.path`, `plink2.path`, and `flashpca.path`.
#' @param n_per_ancestry Number of individuals per ancestry for dummy generation.
#' @param n_snps Number of SNPs for dummy generation.
#' @param overwrite Whether to overwrite existing files.
#' @param pcs_output_dir Optional output directory for `{prefix}.top_pcs.tsv`.
#'
#' @return A list with paths needed by downstream workflow steps.
#' @export
mpgx_prepare_example_genotype <- function(
  output_dir,
  prefix = "toyData",
  source = c("toy_data", "plink_dummy"),
  configs = NULL,
  n_per_ancestry = 50,
  n_snps = 3000,
  overwrite = FALSE,
  pcs_output_dir = NULL
) {
  source <- match.arg(source)
  mpgx_assert_scalar_character(output_dir, "output_dir")
  mpgx_assert_scalar_character(prefix, "prefix")
  mpgx_assert_positive_integer(n_per_ancestry, "n_per_ancestry")
  mpgx_assert_positive_integer(n_snps, "n_snps")

  output_dir <- mpgx_ensure_dir(output_dir)
  if (is.null(pcs_output_dir) || !nzchar(pcs_output_dir)) {
    pcs_output_dir <- file.path(output_dir, "genotype_artifacts")
  }
  pcs_output_dir <- mpgx_ensure_dir(pcs_output_dir)

  if (source == "toy_data") {
    return(mpgx_stage_example_data(
      output_dir = output_dir,
      prefix = prefix,
      overwrite = overwrite,
      pcs_output_dir = pcs_output_dir
    ))
  }

  plink <- mpgx_resolve_executable(
    tool = "plink",
    default = "plink",
    configs = configs,
    config_key = c("plink1.9.path", "plink.path"),
    install_hint = "Install PLINK 1.9+ or pass configs[['plink1.9.path']] explicitly."
  )

  ancestry <- c("AFR", "EAS", "EUR", "SAS")
  ancestry_prefixes <- paste0(file.path(output_dir, prefix), ".", ancestry)

  template_prefix <- paste0(file.path(output_dir, prefix), ".template")
  mpgx_run_command(
    command = plink,
    args = c(
      "--dummy", n_per_ancestry, n_snps, "acgt",
      "--seed 12345",
      "--make-bed",
      "--out", template_prefix
    )
  )

  template_bim_file <- paste0(template_prefix, ".bim")
  template_bim <- utils::read.table(template_bim_file, stringsAsFactors = FALSE)
  nsnp <- nrow(template_bim)
  per_chr <- (nsnp + 21L) %/% 22L
  chr <- ((seq_len(nsnp) - 1L) %/% per_chr) + 1L
  chr <- pmin(chr, 22L)
  bp <- seq_len(nsnp) - ((chr - 1L) * per_chr)
  template_bim$V1 <- chr
  template_bim$V4 <- bp
  mpgx_write_table(template_bim, file = template_bim_file, col_names = FALSE, sep = "\t")

  seed <- 23456
  for (i in seq_along(ancestry)) {
    out_prefix <- ancestry_prefixes[[i]]
    mpgx_run_command(
      command = plink,
      args = c(
        "--dummy", n_per_ancestry, n_snps, "acgt",
        "--seed", seed,
        "--make-bed",
        "--out", out_prefix
      )
    )

    file.copy(template_bim_file, paste0(out_prefix, ".bim"), overwrite = TRUE)

    fam <- utils::read.table(paste0(out_prefix, ".fam"), stringsAsFactors = FALSE)
    fam$V1 <- paste0(ancestry[[i]], fam$V1)
    fam$V2 <- paste0(ancestry[[i]], fam$V2)
    fam$V3 <- ifelse(fam$V3 == 0, 0, paste0(ancestry[[i]], fam$V3))
    fam$V4 <- ifelse(fam$V4 == 0, 0, paste0(ancestry[[i]], fam$V4))
    mpgx_write_table(fam, file = paste0(out_prefix, ".fam"), col_names = FALSE, sep = " ")
    seed <- seed + 1
  }

  merge_list <- file.path(output_dir, paste0(prefix, ".merge_list.txt"))
  merge_rows <- vapply(
    ancestry_prefixes[-1],
    function(x) paste0(x, ".bed ", x, ".bim ", x, ".fam"),
    FUN.VALUE = character(1)
  )
  writeLines(merge_rows, con = merge_list)

  merged_prefix <- file.path(output_dir, prefix)
  mpgx_run_command(
    command = plink,
    args = c(
      "--bfile", ancestry_prefixes[[1]],
      "--merge-list", merge_list,
      "--make-bed",
      "--out", merged_prefix
    )
  )

  fam_all <- utils::read.table(paste0(merged_prefix, ".fam"), stringsAsFactors = FALSE)
  id_map <- data.frame(
    FID = fam_all$V1,
    IID = fam_all$V2,
    ancestry = substr(fam_all$V1, 1, 3),
    stringsAsFactors = FALSE
  )
  ancestry_info_file <- paste0(merged_prefix, ".ancestry.txt")
  mpgx_write_table(id_map, file = ancestry_info_file, col_names = FALSE, sep = "\t")

  pcs_file <- mpgx_generate_top_pcs(
    plink_prefix = merged_prefix,
    configs = configs,
    n_pcs_target = 20L,
    output_dir = pcs_output_dir
  )

  list(
    prefix = merged_prefix,
    ancestry_info_file = ancestry_info_file,
    pcs_file = pcs_file,
    files = c(
      paste0(merged_prefix, ".bed"),
      paste0(merged_prefix, ".bim"),
      paste0(merged_prefix, ".fam"),
      ancestry_info_file,
      pcs_file
    )
  )
}

#' Prepare derived genotype artifacts from a PLINK dataset
#'
#' Builds sparse GRM (for fastGWA/GENESIS), BGEN + index (for LEMMA),
#' and GDS (for GENESIS) from an input PLINK prefix.
#'
#' @param genotype_prefix_or_path PLINK prefix (without extension) or `.bed` path.
#' @param out_dir Output directory for generated artifacts.
#' @param configs Optional named list of executable overrides such as
#'   `gcta.path`, `plink2.path`, and `bgenix.path`.
#' @param sparse_grm_cutoff Numeric cutoff passed to GCTA `--make-bK-sparse`.
#'
#' @return A list with `out_dir`, `sparse_grm_prefix`, `bgen_file`, and `gds_file`.
#' @export
mpgx_prepare_genotype_artifacts <- function(
  genotype_prefix_or_path,
  out_dir,
  configs = NULL,
  sparse_grm_cutoff = 0.05
) {
  mpgx_assert_scalar_character(genotype_prefix_or_path, "genotype_prefix_or_path")
  mpgx_assert_scalar_character(out_dir, "out_dir")

  if (!is.numeric(sparse_grm_cutoff) || length(sparse_grm_cutoff) != 1L ||
      is.na(sparse_grm_cutoff) || sparse_grm_cutoff <= 0) {
    stop("sparse_grm_cutoff must be a positive numeric scalar.", call. = FALSE)
  }

  genotype_prefix <- genotype_prefix_or_path
  if (grepl("\\.bed$", genotype_prefix, ignore.case = TRUE)) {
    genotype_prefix <- sub("\\.bed$", "", genotype_prefix, ignore.case = TRUE)
  } else if (grepl("\\.(bim|fam)$", genotype_prefix, ignore.case = TRUE)) {
    stop(
      "genotype_prefix_or_path must be a PLINK prefix or a .bed path.",
      call. = FALSE
    )
  }

  bed_file <- paste0(genotype_prefix, ".bed")
  bim_file <- paste0(genotype_prefix, ".bim")
  fam_file <- paste0(genotype_prefix, ".fam")
  mpgx_assert_file(bed_file, "genotype_prefix_or_path (.bed)")
  mpgx_assert_file(bim_file, "genotype_prefix_or_path (.bim)")
  mpgx_assert_file(fam_file, "genotype_prefix_or_path (.fam)")

  out_dir <- mpgx_ensure_dir(out_dir)
  artifact_stem <- basename(genotype_prefix)

  gcta <- mpgx_resolve_executable(
    tool = "gcta",
    default = "gcta",
    configs = configs,
    config_key = "gcta.path",
    install_hint = "Install GCTA or pass configs[['gcta.path']] explicitly."
  )
  plink2 <- mpgx_resolve_executable(
    tool = "plink2",
    default = "plink2",
    configs = configs,
    config_key = "plink2.path",
    install_hint = "Install PLINK 2 or pass configs[['plink2.path']] explicitly."
  )

  dense_grm_prefix <- file.path(out_dir, paste0(artifact_stem, ".dense_grm"))
  sparse_grm_prefix <- file.path(out_dir, paste0(artifact_stem, ".sp_grm"))
  mpgx_run_external_command(
    executable = gcta,
    args = c("--bfile", genotype_prefix, "--make-grm", "--out", dense_grm_prefix),
    tool_name = "gcta",
    install_hint = "Install GCTA or pass configs[['gcta.path']] explicitly."
  )
  mpgx_run_external_command(
    executable = gcta,
    args = c(
      "--grm", dense_grm_prefix,
      "--make-bK-sparse", as.character(sparse_grm_cutoff),
      "--out", sparse_grm_prefix
    ),
    tool_name = "gcta",
    install_hint = "Install GCTA or pass configs[['gcta.path']] explicitly."
  )

  bgen_prefix <- file.path(out_dir, artifact_stem)
  mpgx_run_external_command(
    executable = plink2,
    args = c(
      "--bfile", genotype_prefix,
      "--max-alleles", "2",
      "--export", "bgen-1.1",
      "--out", bgen_prefix
    ),
    tool_name = "plink2",
    install_hint = "Install PLINK 2 or pass configs[['plink2.path']] explicitly."
  )
  bgen_file <- paste0(bgen_prefix, ".bgen")
  mpgx_run_bgenix_index(bgen_file = bgen_file, configs = configs)

  gds_file <- paste0(bgen_prefix, ".gds")
  SNPRelate::snpgdsBED2GDS(
    bed.fn = bed_file,
    bim.fn = bim_file,
    fam.fn = fam_file,
    out.gdsfn = gds_file
  )

  list(
    out_dir = out_dir,
    sparse_grm_prefix = sparse_grm_prefix,
    bgen_file = bgen_file,
    gds_file = gds_file
  )
}
