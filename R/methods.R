#' Validate configured executable paths
#'
#' Validates tool paths provided in `configs` and optionally falls back to PATH for
#' requested tools. Recognized config keys are `plink1.9.path` (or legacy `plink.path`),
#' `plink2.path`, `flashpca.path`, `gcta.path`, `lemma.path`, `gem.path`, `bgenix.path`, and
#' `mrmega.path`. `mpirun.path` is optional and only used when you explicitly run
#' LEMMA through `mpirun` (for example in SLURM submissions).
#'
#' @param configs Optional named list of executable overrides.
#' @param tools Optional character vector of tool names to validate.
#'
#' @return Named character vector of resolved executable paths.
#' @export
mpgx_validate_tool_configs <- function(configs = NULL, tools = NULL) {
  specs <- list(
    plink = list(command = "plink", config_key = c("plink1.9.path", "plink.path"), label = "plink"),
    plink2 = list(command = "plink2", config_key = "plink2.path", label = "plink2"),
    flashpca = list(command = "flashpca", config_key = "flashpca.path", label = "flashpca"),
    gcta = list(command = "gcta", config_key = "gcta.path", label = "gcta"),
    lemma = list(command = "lemma_1_0_4", config_key = "lemma.path", label = "lemma_1_0_4"),
    gem = list(command = "GEM", config_key = "gem.path", label = "GEM"),
    bgenix = list(command = "bgenix", config_key = "bgenix.path", label = "bgenix"),
    mrmega = list(command = "MR-MEGA", config_key = "mrmega.path", label = "MR-MEGA"),
    mpirun = list(command = "mpirun", config_key = "mpirun.path", label = "mpirun")
  )

  if (!is.null(configs) && !is.list(configs)) {
    stop("configs must be a list.", call. = FALSE)
  }

  if (is.null(tools)) {
    if (is.null(configs)) {
      selected <- character(0)
    } else {
      selected <- names(specs)[vapply(specs, function(spec) {
        any(as.character(spec$config_key) %in% names(configs))
      }, FUN.VALUE = logical(1))]
    }
  } else {
    tools <- as.character(tools)
    unknown <- setdiff(tools, names(specs))
    if (length(unknown) > 0) {
      stop(sprintf("Unknown tool names: %s", paste(unknown, collapse = ", ")), call. = FALSE)
    }
    selected <- tools
  }

  resolved <- vapply(
    selected,
    function(name) {
      spec <- specs[[name]]
      config_key_label <- if (length(spec$config_key) > 1L) spec$config_key[[1]] else spec$config_key
      legacy_hint <- if (identical(spec$label, "plink")) {
        " (legacy configs[['plink.path']] also supported)"
      } else {
        ""
      }
      mpgx_resolve_executable(
        tool = spec$label,
        default = spec$command,
        configs = configs,
        config_key = spec$config_key,
        install_hint = sprintf(
          "Provide configs[['%s']]%s or add %s to PATH.",
          config_key_label,
          legacy_hint,
          spec$command
        )
      )
    },
    FUN.VALUE = character(1)
  )

  names(resolved) <- selected
  invisible(resolved)
  cat("OK", "\n")
}

#' Execute an external command with validation
#'
#' @param executable Executable path or command name.
#' @param args Command arguments.
#' @param tool_name Optional human-readable tool label for error messages.
#' @param install_hint Optional installation hint appended to errors.
#'
#' @return A list containing command execution metadata.
#' @export
mpgx_run_external_command <- function(
  executable,
  args = character(),
  tool_name = NULL,
  install_hint = NULL
) {
  label <- tool_name %||% executable
  exe <- mpgx_require_executable(label, path = executable, install_hint = install_hint)
  mpgx_run_command(command = exe, args = args)
}

#' Run GCTA linear regression (G + GxE)
#'
#' Wraps the `--fastGWA-lr` command path from the original scripts.
#'
#' @param bfile_prefix PLINK prefix.
#' @param phenotype_file Simulated phenotype file.
#' @param qcovar_file Quantitative covariate file.
#' @param covar_file Binary covariate file.
#' @param envir_file Environment file.
#' @param output_prefix Output prefix for GCTA results.
#' @param gcta_path Optional path to `gcta`.
#' @param keep_file Optional keep file.
#' @param sandwich If `FALSE`, adds `--noSandwich`.
#' @param threads Number of threads.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#' @param configs Optional named list of executable overrides such as `gcta.path`.
#'
#' @return Command metadata list, optionally enriched with `gwas_file` and `gwas_df`.
#' @export
mpgx_run_lr_sandwich <- function(
  bfile_prefix,
  phenotype_file,
  qcovar_file,
  covar_file,
  envir_file,
  output_prefix,
  gcta_path = NULL,
  keep_file = NULL,
  sandwich = TRUE,
  threads = 2,
  attach_result = TRUE,
  configs = NULL
) {
  mpgx_assert_scalar_character(bfile_prefix, "bfile_prefix")
  mpgx_assert_scalar_character(output_prefix, "output_prefix")
  mpgx_assert_positive_integer(threads, "threads")

  mpgx_assert_file(paste0(bfile_prefix, ".bed"), "bfile_prefix .bed")
  mpgx_assert_file(phenotype_file, "phenotype_file")
  mpgx_assert_file(qcovar_file, "qcovar_file")
  mpgx_assert_file(covar_file, "covar_file")
  mpgx_assert_file(envir_file, "envir_file")
  if (!is.null(keep_file)) {
    mpgx_assert_file(keep_file, "keep_file")
  }

  args <- c(
    "--bfile", bfile_prefix,
    "--fastGWA-lr",
    "--pheno", phenotype_file,
    "--qcovar", qcovar_file,
    "--covar", covar_file,
    "--envir", envir_file,
    "--threads", as.integer(threads),
    "--out", output_prefix
  )

  if (!isTRUE(sandwich)) {
    args <- c(args, "--noSandwich")
  }
  if (!is.null(keep_file)) {
    args <- c(args, "--keep", keep_file)
  }

  gcta <- mpgx_resolve_executable(
    tool = "gcta",
    default = "gcta",
    explicit_path = gcta_path,
    configs = configs,
    config_key = "gcta.path",
    install_hint = "Install GCTA >= 1.95 or provide gcta_path or configs[['gcta.path']]."
  )

  run_result <- mpgx_run_external_command(
    executable = gcta,
    args = args,
    tool_name = "gcta",
    install_hint = "Install GCTA >= 1.95 or provide gcta_path or configs[['gcta.path']]."
  )

  run_result <- mpgx_set_attach_metadata(
    result_obj = run_result,
    method = "LRallSandwich",
    files = paste0(output_prefix, ".fastGWA"),
    label = output_prefix
  )

  if (!isTRUE(attach_result)) {
    return(run_result)
  }

  attach_gwas_res(run_result)
}

#' Run fastGWA-GE mixed model analysis
#'
#' Wraps `gcta --fastGWA-mlm` used in the original workflow.
#'
#' @param bfile_prefix PLINK prefix.
#' @param grm_sparse_prefix Prefix of sparse GRM generated by GCTA.
#' @param phenotype_file Phenotype file.
#' @param qcovar_file Quantitative covariate file.
#' @param covar_file Binary covariate file.
#' @param envir_file Environment file.
#' @param output_prefix Output prefix.
#' @param gcta_path Optional path to `gcta`.
#' @param grid_size Grid size used by fastGWA.
#' @param sandwich If `FALSE`, adds `--noSandwich`.
#' @param threads Number of threads.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#' @param configs Optional named list of executable overrides such as `gcta.path`.
#'
#' @return Command metadata list, optionally enriched with `gwas_file` and `gwas_df`.
#' @export
mpgx_run_fastgwa_ge <- function(
  bfile_prefix,
  grm_sparse_prefix,
  phenotype_file,
  qcovar_file,
  covar_file,
  envir_file,
  output_prefix,
  gcta_path = NULL,
  grid_size = 11,
  sandwich = TRUE,
  threads = 8,
  attach_result = TRUE,
  configs = NULL
) {
  mpgx_assert_scalar_character(grm_sparse_prefix, "grm_sparse_prefix")
  mpgx_assert_positive_integer(grid_size, "grid_size")
  mpgx_assert_positive_integer(threads, "threads")

  mpgx_assert_file(paste0(grm_sparse_prefix, ".grm.sp"), "grm sparse file (.grm.sp)")

  args <- c(
    "--bfile", bfile_prefix,
    "--grm-sparse", grm_sparse_prefix,
    "--fastGWA-mlm",
    "--grid-size", as.integer(grid_size),
    "--pheno", phenotype_file,
    "--qcovar", qcovar_file,
    "--covar", covar_file,
    "--envir", envir_file,
    "--threads", as.integer(threads),
    "--optimal-rho",
    "--force-gwa",
    "--h2-limit", 30,
    "--out", output_prefix
  )

  if (!isTRUE(sandwich)) {
    args <- c(args, "--noSandwich")
  }

  gcta <- mpgx_resolve_executable(
    tool = "gcta",
    default = "gcta",
    explicit_path = gcta_path,
    configs = configs,
    config_key = "gcta.path",
    install_hint = "Install GCTA >= 1.95 or provide gcta_path or configs[['gcta.path']]."
  )

  run_result <- mpgx_run_external_command(
    executable = gcta,
    args = args,
    tool_name = "gcta",
    install_hint = "Install GCTA >= 1.95 or provide gcta_path or configs[['gcta.path']]."
  )

  run_result <- mpgx_set_attach_metadata(
    result_obj = run_result,
    method = "fastGWA_GE",
    files = paste0(output_prefix, ".fastGWA"),
    label = output_prefix
  )

  if (!isTRUE(attach_result)) {
    return(run_result)
  }

  attach_gwas_res(run_result)
}

#' Run LEMMA analysis
#'
#' @param bgen_file BGEN file path.
#' @param phenotype_file LEMMA phenotype file.
#' @param envir_file LEMMA environment file.
#' @param covar_file LEMMA covariates file.
#' @param output_prefix Output prefix.
#' @param lemma_path Optional path to `lemma_1_0_4`.
#' @param mpirun_path Optional path to `mpirun`. Only needed when launching
#'   LEMMA with MPI (for example in SLURM submissions).
#' @param n_tasks Number of MPI tasks when `mpirun_path` or `configs[['mpirun.path']]`
#'   is provided.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#' @param configs Optional named list of executable overrides such as `lemma.path`.
#'   Include `mpirun.path` only when you want MPI launch behavior.
#'
#' @return Command metadata list, optionally enriched with `gwas_file` and `gwas_df`.
#' @export
mpgx_run_lemma <- function(
  bgen_file,
  phenotype_file,
  envir_file,
  covar_file,
  output_prefix,
  lemma_path = NULL,
  mpirun_path = NULL,
  n_tasks = 1,
  attach_result = TRUE,
  configs = NULL
) {
  mpgx_assert_file(bgen_file, "bgen_file")
  mpgx_assert_file(phenotype_file, "phenotype_file")
  mpgx_assert_file(envir_file, "envir_file")
  mpgx_assert_file(covar_file, "covar_file")
  mpgx_assert_positive_integer(n_tasks, "n_tasks")

  lemma <- mpgx_resolve_executable(
    tool = "lemma_1_0_4",
    default = "lemma_1_0_4",
    explicit_path = lemma_path,
    configs = configs,
    config_key = "lemma.path",
    install_hint = "Install LEMMA and pass lemma_path or configs[['lemma.path']] when not on PATH."
  )

  use_mpirun <- !is.null(mpirun_path) || !is.null(mpgx_get_config_value(configs, "mpirun.path"))

  lemma_args <- c(
    "--pheno", phenotype_file,
    "--environment", envir_file,
    "--covar", covar_file,
    "--VB",
    "--bgen", bgen_file,
    "--singleSnpStats",
    "--RHEreg",
    "--random-seed", 1,
    "--out", output_prefix
  )

  run_lemma_once <- function() {
    if (isTRUE(use_mpirun)) {
      mpirun <- mpgx_resolve_executable(
        tool = "mpirun",
        default = "mpirun",
        explicit_path = mpirun_path,
        configs = configs,
        config_key = "mpirun.path",
        install_hint = "Install OpenMPI or pass mpirun_path or configs[['mpirun.path']]."
      )

      return(mpgx_run_external_command(
        executable = mpirun,
        args = c("-n", as.integer(n_tasks), lemma, lemma_args),
        tool_name = "mpirun/LEMMA"
      ))
    }

    mpgx_run_external_command(
      executable = lemma,
      args = lemma_args,
      tool_name = "LEMMA"
    )
  }

  run_result <- tryCatch(
      run_lemma_once(),
      error = function(e) {
        mpi_tmpdir <- file.path(dirname(output_prefix), "mpi_tmp")
        dir.create(mpi_tmpdir, recursive = TRUE, showWarnings = FALSE)
        Sys.setenv(OMPI_MCA_orte_tmpdir_base = mpi_tmpdir)
        run_lemma_once()
      }
    )

  run_result <- mpgx_set_attach_metadata(
    result_obj = run_result,
    method = "LEMMA",
    files = sub("^(.+?)\\.LEMMA(.*)$", "\\1_loco_pvals.LEMMA\\2", output_prefix),
    label = output_prefix
  )

  if (!isTRUE(attach_result)) {
    return(run_result)
  }

  attach_gwas_res(run_result)
}

#' Run GEM analysis
#'
#' @param bfile_prefix PLINK prefix.
#' @param phenotype_file GEM phenotype file.
#' @param output_prefix Output prefix.
#' @param gem_path Optional path to GEM executable.
#' @param covar_names Covariate column names in phenotype file.
#' @param exposure_name Exposure variable name.
#' @param sampleid_name Sample ID column.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#' @param configs Optional named list of executable overrides such as `gem.path`.
#'
#' @return Command metadata list, optionally enriched with `gwas_file` and `gwas_df`.
#' @export
mpgx_run_gem <- function(
  bfile_prefix,
  phenotype_file,
  output_prefix,
  gem_path = NULL,
  covar_names = c(
    "binaryCovariate", "continuousCovariate",
    paste0("PC", 1:10)
  ),
  exposure_name = "envir",
  sampleid_name = "IID",
  attach_result = TRUE,
  configs = NULL
) {
  mpgx_assert_file(paste0(bfile_prefix, ".bed"), "bfile_prefix .bed")
  mpgx_assert_file(phenotype_file, "phenotype_file")
  mpgx_assert_scalar_character(output_prefix, "output_prefix")

  gem <- mpgx_resolve_executable(
    tool = "GEM",
    default = "GEM",
    explicit_path = gem_path,
    configs = configs,
    config_key = "gem.path",
    install_hint = "Install GEM 1.4.5+ and pass gem_path or configs[['gem.path']] when not on PATH."
  )

  args <- c(
    "--bfile", bfile_prefix,
    "--pheno-file", phenotype_file,
    "--sampleid-name", sampleid_name,
    "--pheno-name", "phenotype",
    "--covar-names", covar_names,
    "--exposure-names", exposure_name,
    "--robust", 1,
    "--center", 0,
    "--out", output_prefix
  )

  run_result <- mpgx_run_external_command(
    executable = gem,
    args = args,
    tool_name = "GEM"
  )

  run_result <- mpgx_set_attach_metadata(
    result_obj = run_result,
    method = "GEM",
    files = c(
      paste0(output_prefix, ".txt"),
      paste0(output_prefix, ".assoc.txt"),
      output_prefix
    ),
    label = output_prefix
  )

  if (!isTRUE(attach_result)) {
    return(run_result)
  }

  attach_gwas_res(run_result)
}

#' Prepare per-ancestry GEM outputs for MR-MEGA
#'
#' Creates per-ancestry PLINK subsets and GEM phenotype CSV files, runs GEM
#' for each ancestry, and writes MR-MEGA-ready OR-scale files.
#'
#' @param bfile_prefix PLINK prefix for the full cohort genotype data.
#' @param gem_phenotype_file GEM phenotype file containing at least `FID` and
#'   `IID` columns.
#' @param keep_files Named list or character vector of keep files by ancestry.
#' @param output_dir Output directory for per-ancestry temporary and result
#'   files.
#' @param plink_path Optional path to `plink` executable.
#' @param configs Optional named list of executable overrides such as
#'   `plink1.9.path`, `plink.path`, and `gem.path`.
#'
#' @return Named character vector of per-ancestry MR-MEGA input file paths.
#' @export
mpgx_prepare_mrmega_gem_inputs <- function(
  bfile_prefix,
  gem_phenotype_file,
  keep_files,
  output_dir,
  plink_path = NULL,
  configs = NULL
) {
  mpgx_assert_file(paste0(bfile_prefix, ".bed"), "bfile_prefix .bed")
  mpgx_assert_file(gem_phenotype_file, "gem_phenotype_file")
  mpgx_assert_scalar_character(output_dir, "output_dir")

  keep_list <- as.list(keep_files)
  keep_names <- names(keep_list)

  if (!length(keep_list)) {
    stop("keep_files must contain at least one ancestry keep file.", call. = FALSE)
  }
  if (is.null(keep_names) || any(!nzchar(keep_names))) {
    stop("keep_files must be named with ancestry labels.", call. = FALSE)
  }

  keep_paths <- as.character(unname(unlist(keep_list, use.names = FALSE)))
  missing_keep <- keep_paths[!file.exists(keep_paths)]
  if (length(missing_keep)) {
    stop(
      "Missing keep_files: ",
      paste(missing_keep, collapse = ", "),
      call. = FALSE
    )
  }

  mpgx_ensure_dir(output_dir)

  plink <- mpgx_resolve_executable(
    tool = "plink",
    default = "plink",
    explicit_path = plink_path,
    configs = configs,
    config_key = c("plink1.9.path", "plink.path"),
    install_hint = "Install PLINK 1.9 or provide plink_path or configs[['plink1.9.path']]."
  )

  sim_gem <- utils::read.csv(gem_phenotype_file, stringsAsFactors = FALSE)
  required_pheno_cols <- c("FID", "IID")
  missing_pheno_cols <- setdiff(required_pheno_cols, colnames(sim_gem))
  if (length(missing_pheno_cols)) {
    stop(
      "gem_phenotype_file is missing required columns: ",
      paste(missing_pheno_cols, collapse = ", "),
      call. = FALSE
    )
  }

  base_stem <- basename(bfile_prefix)

  out_files <- vapply(
    keep_names,
    function(ances) {
      keep_file <- as.character(keep_list[[ances]])[[1]]

      keep_df <- utils::read.table(keep_file, stringsAsFactors = FALSE)
      if (ncol(keep_df) < 2) {
        stop("keep_file must contain at least FID and IID columns: ", keep_file, call. = FALSE)
      }
      colnames(keep_df)[1:2] <- c("FID", "IID")
      keep_df <- keep_df[, c("FID", "IID"), drop = FALSE]

      eth_bfile <- file.path(output_dir, paste0(base_stem, ".", ances))
      mpgx_run_external_command(
        executable = plink,
        args = c(
          "--bfile", bfile_prefix,
          "--keep", keep_file,
          "--make-bed",
          "--out", eth_bfile
        ),
        tool_name = "plink",
        install_hint = "Install PLINK 1.9 or provide plink_path or configs[['plink1.9.path']]."
      )

      gem_eth <- merge(keep_df, sim_gem, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
      gem_eth_file <- file.path(output_dir, paste0("binary.", ances, ".forGEM.csv"))
      utils::write.csv(gem_eth, gem_eth_file, row.names = FALSE, quote = FALSE)

      gem_out_prefix <- file.path(output_dir, paste0("binary.GEM.", ances))
      gem_run <- mpgx_run_gem(
        bfile_prefix = eth_bfile,
        phenotype_file = gem_eth_file,
        output_prefix = gem_out_prefix,
        configs = configs
      )

      gem_result <- gem_run$gwas_file
      if (is.null(gem_result) || length(gem_result) != 1L || !file.exists(gem_result)) {
        stop("Unable to resolve GEM output file for ancestry '", ances, "'.", call. = FALSE)
      }

      gem_tab <- mpgx_read_delim_table(gem_result, has_header = TRUE, arg_name = "gem_result")
      required_gem_cols <- c(
        "Beta_Marginal", "robust_SE_Beta_Marginal",
        "Beta_G-envir", "robust_SE_Beta_G-envir"
      )
      missing_gem_cols <- setdiff(required_gem_cols, colnames(gem_tab))
      if (length(missing_gem_cols)) {
        stop(
          "Missing expected GEM output columns for ancestry '",
          ances,
          "': ",
          paste(missing_gem_cols, collapse = ", "),
          call. = FALSE
        )
      }

      gem_tab$OR_G <- exp(gem_tab$Beta_Marginal)
      gem_tab$OR_95L_G <- exp(gem_tab$Beta_Marginal - 1.96 * gem_tab$robust_SE_Beta_Marginal)
      gem_tab$OR_95U_G <- exp(gem_tab$Beta_Marginal + 1.96 * gem_tab$robust_SE_Beta_Marginal)
      gem_tab$OR_G_by_E <- exp(gem_tab[["Beta_G-envir"]])
      gem_tab$OR_95L_G_by_E <- exp(gem_tab[["Beta_G-envir"]] - 1.96 * gem_tab[["robust_SE_Beta_G-envir"]])
      gem_tab$OR_95U_G_by_E <- exp(gem_tab[["Beta_G-envir"]] + 1.96 * gem_tab[["robust_SE_Beta_G-envir"]])

      out_for_mrmega <- file.path(output_dir, paste0("binary.GEM.", ances, ".forMRMEGA.txt"))
      utils::write.table(gem_tab, out_for_mrmega, quote = FALSE, row.names = FALSE, sep = "\t")
      out_for_mrmega
    },
    FUN.VALUE = character(1)
  )

  names(out_files) <- keep_names
  out_files
}

mpgx_read_two_id_file <- function(path, value_prefix) {
  dat <- utils::read.table(path, stringsAsFactors = FALSE)
  if (ncol(dat) < 2) {
    stop(sprintf("File must contain at least FID and IID columns: %s", path), call. = FALSE)
  }
  extra <- ncol(dat) - 2
  colnames(dat) <- c("FID", "IID", if (extra > 0) paste0(value_prefix, seq_len(extra)) else NULL)
  dat
}

mpgx_extract_beta_se <- function(coef_vec, vcov_mat, term) {
  if (!(term %in% names(coef_vec)) || !(term %in% rownames(vcov_mat))) {
    return(c(beta = NA_real_, se = NA_real_))
  }
  beta <- as.numeric(coef_vec[[term]])
  se <- sqrt(as.numeric(vcov_mat[term, term]))
  c(beta = beta, se = se)
}

#' Run GENESIS association analysis
#'
#' Refactors the original `runGENESIS.R` logic into a callable package function.
#'
#' @param gds_file Path to genotype GDS file.
#' @param sparse_grm_prefix Prefix for sparse GRM files (`.grm.id`, `.grm.sp`).
#' @param phenotype_file Phenotype file with FID, IID, phenotype.
#' @param qcovar_file Quantitative covariate file.
#' @param covar_file Binary covariate file.
#' @param envir_file Environment file.
#' @param output_file Optional output file path.
#' @param snp_block SNP block size for iterator.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#'
#' @return Data frame of GENESIS association results with optional attached
#'   `gwas_file` and `gwas_df` attributes.
#' @export
mpgx_run_genesis <- function(
  gds_file,
  sparse_grm_prefix,
  phenotype_file,
  qcovar_file,
  covar_file,
  envir_file,
  output_file = NULL,
  snp_block = 5000,
  attach_result = TRUE
) {
  mpgx_assert_file(gds_file, "gds_file")
  mpgx_assert_file(paste0(sparse_grm_prefix, ".grm.id"), "sparse_grm_prefix .grm.id")
  mpgx_assert_file(paste0(sparse_grm_prefix, ".grm.sp"), "sparse_grm_prefix .grm.sp")
  mpgx_assert_file(phenotype_file, "phenotype_file")
  mpgx_assert_file(qcovar_file, "qcovar_file")
  mpgx_assert_file(covar_file, "covar_file")
  mpgx_assert_file(envir_file, "envir_file")
  mpgx_assert_positive_integer(snp_block, "snp_block")

  geno <- GWASTools::GdsGenotypeReader(filename = gds_file)
  on.exit(GWASTools::close(geno), add = TRUE)

  pheno <- mpgx_read_two_id_file(phenotype_file, "pheno")
  qcovar <- mpgx_read_two_id_file(qcovar_file, "qcovar")
  covar <- mpgx_read_two_id_file(covar_file, "covar")
  envir <- mpgx_read_two_id_file(envir_file, "envir")

  if (!("pheno1" %in% colnames(pheno))) {
    stop("phenotype_file must include one phenotype column after FID/IID.", call. = FALSE)
  }
  if (!("envir1" %in% colnames(envir))) {
    stop("envir_file must include one environment column after FID/IID.", call. = FALSE)
  }

  mydat <- Reduce(
    function(x, y) merge(x, y, by = c("FID", "IID"), all = TRUE, sort = FALSE),
    list(pheno, qcovar, covar, envir)
  )
  mydat$scanID <- mydat$IID

  scan_ids <- data.frame(scanID = GWASTools::getScanID(geno), stringsAsFactors = FALSE)
  mydat <- merge(scan_ids, mydat, by = "scanID", all.x = TRUE, sort = FALSE)

  scan_annot <- GWASTools::ScanAnnotationDataFrame(mydat)
  geno_data <- GWASTools::GenotypeData(geno, scanAnnot = scan_annot)

  grm_id <- utils::read.table(paste0(sparse_grm_prefix, ".grm.id"), stringsAsFactors = FALSE)
  grm_sp <- utils::read.table(paste0(sparse_grm_prefix, ".grm.sp"), stringsAsFactors = FALSE)

  sparse_grm <- Matrix::forceSymmetric(
    stats::xtabs(V3 ~ V1 + V2, data = grm_sp, sparse = TRUE),
    uplo = "L"
  )
  rownames(sparse_grm) <- colnames(sparse_grm) <- grm_id$V2

  covars <- c(grep("^(qcovar|covar)", colnames(mydat), value = TRUE), "envir1")
  binary_outcome <- nlevels(factor(mydat$pheno1)) == 2

  null_model <- if (binary_outcome) {
    GENESIS::fitNullModel(
      scan_annot,
      outcome = "pheno1",
      covars = covars,
      cov.mat = sparse_grm,
      family = stats::binomial
    )
  } else {
    GENESIS::fitNullModel(
      scan_annot,
      outcome = "pheno1",
      covars = covars,
      cov.mat = sparse_grm,
      family = stats::gaussian
    )
  }

  geno_iter <- GWASTools::GenotypeBlockIterator(geno_data, snpBlock = as.integer(snp_block))
  assoc <- if (binary_outcome) {
    GENESIS::assocTestSingle(geno_iter, GxE = "envir1", null.model = null_model, test = "Score.SPA")
  } else {
    GENESIS::assocTestSingle(geno_iter, GxE = "envir1", null.model = null_model)
  }

  assoc <- as.data.frame(assoc)

  if (!is.null(output_file)) {
    mpgx_ensure_dir(dirname(output_file))
    mpgx_write_table(assoc, file = output_file, col_names = TRUE, sep = "\t")
  }

  assoc <- mpgx_set_attach_metadata(
    result_obj = assoc,
    method = "GENESIS",
    files = output_file,
    label = output_file %||% "GENESIS"
  )

  if (!isTRUE(attach_result)) {
    return(assoc)
  }

  if (is.null(output_file)) {
    warning(
      "attach_result = TRUE but output_file is NULL in mpgx_run_genesis(); returning GENESIS results without attached gwas_df.",
      call. = FALSE
    )
    return(assoc)
  }

  attach_gwas_res(assoc)
}

#' Index a BGEN file with bgenix
#'
#' @param bgen_file Path to BGEN file.
#' @param bgenix_path Optional path to `bgenix`.
#' @param configs Optional named list of executable overrides such as `bgenix.path`.
#'
#' @return Command metadata list.
#' @export
mpgx_run_bgenix_index <- function(
  bgen_file,
  bgenix_path = NULL,
  configs = NULL
) {
  mpgx_assert_file(bgen_file, "bgen_file")

  bgenix <- mpgx_resolve_executable(
    tool = "bgenix",
    default = "bgenix",
    explicit_path = bgenix_path,
    configs = configs,
    config_key = "bgenix.path",
    install_hint = "Install bgenix or provide bgenix_path or configs[['bgenix.path']]."
  )

  mpgx_run_external_command(
    executable = bgenix,
    args = c("-g", bgen_file, "-index -clobber"),
    tool_name = "bgenix",
    install_hint = "Install bgenix or provide bgenix_path or configs[['bgenix.path']]."
  )
}

#' Run MR-MEGA meta-analysis command
#'
#' @param input_file_list File listing cohort result files.
#' @param output_prefix Output prefix.
#' @param mrmega_path Optional path to `MR-MEGA` executable.
#' @param extra_args Additional command-line arguments passed to MR-MEGA.
#' @param configs Optional named list of executable overrides such as `mrmega.path`.
#'
#' @return Command metadata list.
#' @export
mpgx_run_mrmega <- function(
  input_file_list,
  output_prefix,
  mrmega_path = NULL,
  extra_args = character(),
  configs = NULL
) {
  mpgx_assert_file(input_file_list, "input_file_list")
  mpgx_assert_scalar_character(output_prefix, "output_prefix")

  args <- c("-i", input_file_list, "-o", output_prefix, as.character(extra_args))

  mrmega <- mpgx_resolve_executable(
    tool = "MR-MEGA",
    default = "MR-MEGA",
    explicit_path = mrmega_path,
    configs = configs,
    config_key = "mrmega.path",
    install_hint = "Install MR-MEGA or provide mrmega_path or configs[['mrmega.path']]."
  )

  mpgx_run_external_command(
    executable = mrmega,
    args = args,
    tool_name = "MR-MEGA",
    install_hint = "Install MR-MEGA or provide mrmega_path or configs[['mrmega.path']]."
  )
}

#' Write an MR-MEGA input file list from cohort result files
#'
#' @param files Character vector or list of cohort result file paths.
#' @param output_file Output `.in` file path.
#' @param label Label used in error messages.
#'
#' @return Character scalar path to the written input list.
#' @export
mpgx_write_mrmega_input_list <- function(files, output_file, label = "MR-MEGA input list") {
  mpgx_assert_scalar_character(output_file, "output_file")

  file_list <- as.character(unname(unlist(files, use.names = FALSE)))
  if (!length(file_list)) {
    stop("No cohort result files were provided for ", label, ".", call. = FALSE)
  }

  missing_files <- file_list[!file.exists(file_list)]
  if (length(missing_files)) {
    stop(
      "Missing cohort result file(s) for ",
      label,
      ": ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  mpgx_ensure_dir(dirname(output_file))
  writeLines(file_list, output_file)
  if (!file.exists(output_file)) {
    stop("Failed to create ", label, ": ", output_file, call. = FALSE)
  }

  output_file
}

#' Run MR-MEGA for LR-Sandwich per-ancestry results
#'
#' @param ancestry_result_files Character vector or list of per-ancestry LR
#'   result files.
#' @param output_prefix_base Output prefix base for MR-MEGA outputs.
#' @param mrmega_path Optional path to `MR-MEGA` executable.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#' @param configs Optional named list of executable overrides such as
#'   `mrmega.path`.
#'
#' @return Result list with MR-MEGA file paths and command metadata, optionally
#'   enriched with `gwas_file` and `gwas_df`.
#' @export
mpgx_run_mrmega_for_lr <- function(
  ancestry_result_files,
  output_prefix_base,
  mrmega_path = NULL,
  attach_result = TRUE,
  configs = NULL
) {
  mpgx_assert_scalar_character(output_prefix_base, "output_prefix_base")

  input_file_list <- paste0(output_prefix_base, ".in")
  mpgx_write_mrmega_input_list(
    files = ancestry_result_files,
    output_file = input_file_list,
    label = "MR-MEGA input list"
  )

  n_cohorts <- length(as.character(unname(unlist(ancestry_result_files, use.names = FALSE))))
  pc_num <- max(n_cohorts - 3, 0)

  out_g <- paste0(output_prefix_base, ".test_G")
  out_gxe <- paste0(output_prefix_base, ".test_G_by_E")

  run_test_G <- mpgx_run_mrmega(
    input_file_list = input_file_list,
    output_prefix = out_g,
    mrmega_path = mrmega_path,
    extra_args = c(
      "--name_chr", "CHR", "--name_marker", "SNP", "--name_pos", "POS",
      "--name_n", "N", "--name_ea", "A1", "--name_nea", "A2",
      "--name_eaf", "AF1", "--qt", "--name_beta", "BETA_G",
      "--name_se", "SE_G", "--pc", as.character(pc_num)
    ),
    configs = configs
  )

  run_test_G_by_E <- mpgx_run_mrmega(
    input_file_list = input_file_list,
    output_prefix = out_gxe,
    mrmega_path = mrmega_path,
    extra_args = c(
      "--name_chr", "CHR", "--name_marker", "SNP", "--name_pos", "POS",
      "--name_n", "N", "--name_ea", "A1", "--name_nea", "A2",
      "--name_eaf", "AF1", "--qt", "--name_beta", "BETA_G_by_E",
      "--name_se", "SE_G_by_E", "--pc", as.character(pc_num)
    ),
    configs = configs
  )

  result_obj <- list(
    input_file_list = input_file_list,
    test_G = paste0(out_g, ".result"),
    test_G_by_E = paste0(out_gxe, ".result"),
    run_test_G = run_test_G,
    run_test_G_by_E = run_test_G_by_E
  )

  if (!isTRUE(attach_result)) {
    return(result_obj)
  }

  result_obj$test_G <- resolve_gwas_file(
    candidates = result_obj$test_G,
    label = paste0(output_prefix_base, " test_G")
  )
  result_obj$test_G_by_E <- resolve_gwas_file(
    candidates = result_obj$test_G_by_E,
    label = paste0(output_prefix_base, " test_G_by_E")
  )

  attach_gwas_df(
    result_obj = result_obj,
    method = "MRMEGA_LRallSandwich",
    files = list(
      test_G = result_obj$test_G,
      test_G_by_E = result_obj$test_G_by_E
    ),
    label = output_prefix_base
  )
}

#' Run MR-MEGA for GEM per-ancestry results
#'
#' @param ancestry_result_files Character vector or list of per-ancestry GEM
#'   result files prepared for MR-MEGA input.
#' @param output_prefix_base Output prefix base for MR-MEGA outputs.
#' @param mrmega_path Optional path to `MR-MEGA` executable.
#' @param attach_result If `TRUE`, parse and attach GWAS results.
#' @param configs Optional named list of executable overrides such as
#'   `mrmega.path`.
#'
#' @return Result list with MR-MEGA file paths and command metadata, optionally
#'   enriched with `gwas_file` and `gwas_df`.
#' @export
mpgx_run_mrmega_for_gem <- function(
  ancestry_result_files,
  output_prefix_base,
  mrmega_path = NULL,
  attach_result = TRUE,
  configs = NULL
) {
  mpgx_assert_scalar_character(output_prefix_base, "output_prefix_base")

  input_file_list <- paste0(output_prefix_base, ".in")
  mpgx_write_mrmega_input_list(
    files = ancestry_result_files,
    output_file = input_file_list,
    label = "MR-MEGA input list"
  )

  n_cohorts <- length(as.character(unname(unlist(ancestry_result_files, use.names = FALSE))))
  pc_num <- max(n_cohorts - 3, 0)

  out_g <- paste0(output_prefix_base, ".test_G")
  out_gxe <- paste0(output_prefix_base, ".test_G_by_E")

  run_test_G <- mpgx_run_mrmega(
    input_file_list = input_file_list,
    output_prefix = out_g,
    mrmega_path = mrmega_path,
    extra_args = c(
      "--name_chr", "CHR", "--name_marker", "SNPID", "--name_pos", "POS",
      "--name_n", "N_Samples", "--name_ea", "Effect_Allele", "--name_nea", "Non_Effect_Allele",
      "--name_eaf", "AF", "--name_or", "OR_G", "--name_or_95l", "OR_95L_G",
      "--name_or_95u", "OR_95U_G", "--pc", as.character(pc_num)
    ),
    configs = configs
  )

  run_test_G_by_E <- mpgx_run_mrmega(
    input_file_list = input_file_list,
    output_prefix = out_gxe,
    mrmega_path = mrmega_path,
    extra_args = c(
      "--name_chr", "CHR", "--name_marker", "SNPID", "--name_pos", "POS",
      "--name_n", "N_Samples", "--name_ea", "Effect_Allele", "--name_nea", "Non_Effect_Allele",
      "--name_eaf", "AF", "--name_or", "OR_G_by_E", "--name_or_95l", "OR_95L_G_by_E",
      "--name_or_95u", "OR_95U_G_by_E", "--pc", as.character(pc_num)
    ),
    configs = configs
  )

  result_obj <- list(
    input_file_list = input_file_list,
    test_G = paste0(out_g, ".result"),
    test_G_by_E = paste0(out_gxe, ".result"),
    run_test_G = run_test_G,
    run_test_G_by_E = run_test_G_by_E
  )

  if (!isTRUE(attach_result)) {
    return(result_obj)
  }

  result_obj$test_G <- resolve_gwas_file(
    candidates = result_obj$test_G,
    label = paste0(output_prefix_base, " test_G")
  )
  result_obj$test_G_by_E <- resolve_gwas_file(
    candidates = result_obj$test_G_by_E,
    label = paste0(output_prefix_base, " test_G_by_E")
  )

  attach_gwas_df(
    result_obj = result_obj,
    method = "MRMEGA_GEM",
    files = list(
      test_G = result_obj$test_G,
      test_G_by_E = result_obj$test_G_by_E
    ),
    label = output_prefix_base
  )
}

#' Run fixed-effects meta-analysis for LR-Sandwich and attach parsed output
#'
#' @param ancestry_result_files Character vector (or named list) of per-
#'   ancestry LR output files.
#' @param ancestry_labels Optional character vector of ancestry labels.
#' @param output_file Output file path written by the fixed-effects meta-analysis.
#' @param attach_result If `TRUE`, parse and attach GWAS results from
#'   `output_file`.
#'
#' @return Data frame from [mpgx_fixed_effect_meta_lr()], optionally enriched
#'   with attached `gwas_file` and `gwas_df` attributes.
#' @export
mpgx_run_fe_meta_for_lr <- function(
  ancestry_result_files,
  ancestry_labels = NULL,
  output_file,
  attach_result = TRUE
) {
  mpgx_assert_scalar_character(output_file, "output_file")

  run_result <- mpgx_fixed_effect_meta_lr(
    ancestry_result_files = ancestry_result_files,
    ancestry_labels = ancestry_labels,
    output_file = output_file
  )

  if (!isTRUE(attach_result)) {
    return(run_result)
  }

  attach_gwas_df(
    result_obj = run_result,
    method = "LRallSandwich",
    files = output_file,
    label = output_file
  )
}
#' Fixed-effects meta-analysis for LR-Sandwich cohort files
#'
#' Re-implements the fixed-effect branch from `scripts/summarize_results.R`
#' for `meta_*` methods, while retaining legacy reference snippets in an
#' `if (F)` block for traceability.
#'
#' @param ancestry_result_files Character vector (or named list) of per-ancestry
#'   LR output files (e.g. `*.fastGWA`).
#' @param ancestry_labels Optional character vector of ancestry labels. If
#'   `NULL`, uses names from `ancestry_result_files` when available, otherwise
#'   uses `cohort1`, `cohort2`, ...
#' @param output_file Optional output file path to write a tab-delimited table.
#'
#' @return Data frame with one row per SNP and meta-analysis statistics.
#' @export
mpgx_fixed_effect_meta_lr <- function(
  ancestry_result_files,
  ancestry_labels = NULL,
  output_file = NULL
) {
  file_list <- as.character(unname(unlist(ancestry_result_files, use.names = FALSE)))
  if (!length(file_list)) {
    stop("ancestry_result_files must contain at least one LR result file.", call. = FALSE)
  }
  if (!all(file.exists(file_list))) {
    missing_files <- file_list[!file.exists(file_list)]
    stop(
      "Missing ancestry_result_files: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(ancestry_labels)) {
    name_guess <- names(ancestry_result_files)
    if (!is.null(name_guess) && length(name_guess) == length(file_list) && all(nzchar(name_guess))) {
      ancestry_labels <- name_guess
    } else {
      ancestry_labels <- paste0("cohort", seq_along(file_list))
    }
  }
  ancestry_labels <- as.character(ancestry_labels)
  if (length(ancestry_labels) != length(file_list)) {
    stop("ancestry_labels length must match ancestry_result_files length.", call. = FALSE)
  }

  required_cols <- c(
    "CHR", "SNP", "BETA_G", "BETA_G_by_E",
    "SE_G", "SE_G_by_E", "Cov_BETA_G_and_G_by_E"
  )

  result_list <- lapply(seq_along(file_list), function(i) {
    dat <- mpgx_read_delim_table(file_list[[i]], has_header = TRUE, arg_name = "ancestry_result_files")
    missing_cols <- setdiff(required_cols, colnames(dat))
    if (length(missing_cols)) {
      stop(
        "Missing expected columns in LR result file: ",
        file_list[[i]],
        " (missing: ",
        paste(missing_cols, collapse = ", "),
        ")",
        call. = FALSE
      )
    }

    dat <- dat[, required_cols, drop = FALSE]
    colnames(dat) <- c(
      "CHR", "SNP",
      paste0(ancestry_labels[[i]], "_", c(
        "BETA_G", "BETA_G_by_E", "SE_G", "SE_G_by_E", "Cov_BETA_G_and_G_by_E"
      ))
    )
    dat
  })

  result <- Reduce(function(x, y) dplyr::full_join(x, y, by = c("CHR", "SNP")), result_list)
  eth_list <- ancestry_labels

  # Legacy reference implementation kept intentionally for traceability.
  if (F) {
    meta_analysis <- t(apply(dplyr::select(result[1:20, , drop = FALSE], -CHR, -SNP), 1, function(x) {
      #print(x[1])
      #b=NULL
      Sigma_inv_list <- list()
      #W_list=list()
      tW_Sigma_inv_b_list <- list()
      for (ancestry in eth_list) {
        if (!is.na(x[paste(ancestry, "BETA_G", sep = "_")])) {
          temp_beta <- t(x[paste(ancestry, c("BETA_G", "BETA_G_by_E"), sep = "_")])
          #b=c(b, temp_beta )
          temp_mat <- matrix(NA, 2, 2)
          temp_mat[1, 1] <- x[paste(ancestry, "SE_G", sep = "_")]^2
          temp_mat[1, 2] <- temp_mat[2, 1] <- as.numeric(x[paste(ancestry, "Cov_BETA_G_and_G_by_E", sep = "_")])
          temp_mat[2, 2] <- x[paste(ancestry, "SE_G_by_E", sep = "_")]^2
          temp_mat_inv <- solve(temp_mat)
          Sigma_inv_list[[ancestry]] <- temp_mat_inv
          tW_Sigma_inv_b_list[[ancestry]] <- temp_mat_inv %*% as.vector(temp_beta)
          #W_list[[ancestry]]=diag(2)
        }
      }
      #Sigma_inv=Reduce(magic::adiag,Sigma_inv_list)
      #W=do.call(rbind,W_list)
      #Cov_BETA_meta=MASS::ginv(t(W)%*%Sigma_inv%*%W)
      Cov_BETA_meta <- solve(Reduce(`+`, Sigma_inv_list))

      tW_Sigma_inv_b <- Reduce(`+`, tW_Sigma_inv_b_list)
      #BETA_meta=Cov_BETA_meta%*%(t(W)%*%(Sigma_inv%*%b))
      BETA_meta <- Cov_BETA_meta %*% tW_Sigma_inv_b

      chisq_2df <- t(BETA_meta) %*% solve(Cov_BETA_meta) %*% BETA_meta
      P_2df <- 1 - stats::pchisq(chisq_2df, df = 2)

      c("chisq_2df" = chisq_2df, "P_2df" = P_2df)
    }))
  }

  manual_inverse_2by2 <- function(mat_2by2) {
    out <- matrix(NA_real_, 2, 2)
    mat_det <- mat_2by2[1, 1] * mat_2by2[2, 2] - mat_2by2[1, 2] * mat_2by2[2, 1]
    out[1, 1] <- mat_2by2[2, 2]
    out[2, 2] <- mat_2by2[1, 1]
    out[1, 2] <- -mat_2by2[1, 2]
    out[2, 1] <- -mat_2by2[2, 1]
    out <- out / mat_det
    out
  }

  manual_inverse_2by2_vector <- function(a, b, d) {
    out <- matrix(NA_real_, 2, 2)
    mat_det <- a * d - b^2
    out[1, 1] <- d
    out[2, 2] <- a
    out[1, 2] <- out[2, 1] <- -b
    out <- out / mat_det
    out
  }

  meta_analysis <- t(apply(dplyr::select(result, -CHR, -SNP), 1, function(x) {
    Cov_BETA_meta <- matrix(0, 2, 2)
    tW_Sigma_inv_b <- matrix(0, 2, 1)

    sum_invvar_G <- 0
    sum_invvar_G_by_E <- 0
    sum_weighted_beta_G <- 0
    sum_weighted_beta_G_by_E <- 0

    for (ancestry in eth_list) {
      if (!is.na(x[paste(ancestry, "BETA_G", sep = "_")])) {
        temp_beta <- t(x[paste(ancestry, c("BETA_G", "BETA_G_by_E"), sep = "_")])

        temp_beta_G <- as.numeric(x[paste(ancestry, c("BETA_G"), sep = "_")])
        temp_beta_G_by_E <- as.numeric(x[paste(ancestry, c("BETA_G_by_E"), sep = "_")])
        temp_invvar_G <- 1 / as.numeric(x[paste(ancestry, "SE_G", sep = "_")]^2)
        temp_invvar_G_by_E <- 1 / as.numeric(x[paste(ancestry, "SE_G_by_E", sep = "_")]^2)

        sum_weighted_beta_G <- sum_weighted_beta_G + temp_invvar_G * temp_beta_G
        sum_weighted_beta_G_by_E <- sum_weighted_beta_G_by_E + temp_invvar_G_by_E * temp_beta_G_by_E
        sum_invvar_G <- sum_invvar_G + temp_invvar_G
        sum_invvar_G_by_E <- sum_invvar_G_by_E + temp_invvar_G_by_E

        #temp_mat=matrix(NA,2,2)
        #temp_mat[1,1]=x[paste(ancestry,"SE_G",sep="_")]^2
        #temp_mat[1,2]=temp_mat[2,1]=as.numeric(x[paste(ancestry,"Cov_BETA_G_and_G_by_E",sep="_")])
        #temp_mat[2,2]=x[paste(ancestry,"SE_G_by_E",sep="_")]^2
        # temp_mat_inv=manual_inverse_2by2(temp_mat)
        temp_mat_inv <- manual_inverse_2by2_vector(
          as.numeric(x[paste(ancestry, "SE_G", sep = "_")]^2),
          as.numeric(x[paste(ancestry, "Cov_BETA_G_and_G_by_E", sep = "_")]),
          as.numeric(x[paste(ancestry, "SE_G_by_E", sep = "_")]^2)
        )
        Cov_BETA_meta <- Cov_BETA_meta + temp_mat_inv
        tW_Sigma_inv_b <- tW_Sigma_inv_b + temp_mat_inv %*% t(temp_beta)
      }
    }

    Cov_BETA_meta <- manual_inverse_2by2(Cov_BETA_meta)
    BETA_meta <- Cov_BETA_meta %*% tW_Sigma_inv_b
    chisq_2df <- t(BETA_meta) %*% manual_inverse_2by2(Cov_BETA_meta) %*% BETA_meta
    P_2df <- 1 - stats::pchisq(chisq_2df, df = 2)

    BETA_G <- sum_weighted_beta_G / sum_invvar_G
    SE_G <- sqrt(1 / sum_invvar_G)
    chisq_G <- BETA_G^2 * sum_invvar_G
    P_G <- 1 - stats::pchisq(chisq_G, df = 1)

    BETA_G_by_E <- sum_weighted_beta_G_by_E / sum_invvar_G_by_E
    SE_G_by_E <- sqrt(1 / sum_invvar_G_by_E)
    chisq_G_by_E <- BETA_G_by_E^2 * sum_invvar_G_by_E
    P_G_by_E <- 1 - stats::pchisq(chisq_G_by_E, df = 1)

    c(
      "BETA_G" = BETA_G,
      "SE_G" = SE_G,
      "BETA_G_by_E" = BETA_G_by_E,
      "SE_G_by_E" = SE_G_by_E,
      "chisq_G" = chisq_G,
      "chisq_G_by_E" = chisq_G_by_E,
      "chisq_2df" = chisq_2df,
      "P_G" = P_G,
      "P_G_by_E" = P_G_by_E,
      "P_2df" = P_2df
    )
  }))

  meta_analysis <- as.data.frame(meta_analysis, stringsAsFactors = FALSE)
  for (nm in colnames(meta_analysis)) {
    meta_analysis[[nm]] <- as.numeric(meta_analysis[[nm]])
  }

  out <- cbind(result[, c("CHR", "SNP"), drop = FALSE], meta_analysis)

  if (!is.null(output_file)) {
    mpgx_ensure_dir(dirname(output_file))
    utils::write.table(out, output_file, quote = FALSE, row.names = FALSE, sep = "\t")
  }

  out
}

