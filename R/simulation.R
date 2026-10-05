#' Parse variance proportion settings
#'
#' Accepts either a numeric vector of length 6 or a script-style underscore
#' delimited string in the order:
#' `main`, `interaction`, `covariates`, `stratification`, `shared`, `error`.
#'
#' @param x Numeric vector or single character string.
#'
#' @return Named numeric vector of length 6.
#' @export
mpgx_parse_variance_proportion <- function(x) {
  if (is.character(x) && length(x) == 1) {
    x <- as.numeric(strsplit(x, split = "_", fixed = TRUE)[[1]])
  }

  if (!is.numeric(x) || length(x) != 6 || any(is.na(x))) {
    stop(
      "variance_proportion must contain 6 numeric values in the order main_interaction_covariates_stratification_shared_error.",
      call. = FALSE
    )
  }

  names(x) <- c("main", "interaction", "covariates", "stratification", "shared", "error")
  x
}

#' Select causal SNPs from a BIM file
#'
#' Mirrors the logic of the original `random_select_causal_SNPs.R` script.
#'
#' @param bim_file Path to PLINK BIM file.
#' @param n_snps_maineffects_only Number of SNPs with G-only effects.
#' @param n_snps_interaction_only Number of SNPs with GxE-only effects.
#' @param n_snps_both Number of SNPs with both G and GxE effects.
#' @param seed Seed for sampling SNP identities.
#' @param ethnicity Ethnicity label appended to output file names.
#' @param seed_effect_size Seed for sampling effect sizes.
#' @param effect_size_mean Mean of effect-size distribution.
#' @param effect_size_sd SD of effect-size distribution.
#' @param out_dir Optional output directory for writing script-compatible files.
#'
#' @return List with `main`, `interaction`, and optional file paths.
#' @export
mpgx_select_causal_snps <- function(
  bim_file,
  n_snps_maineffects_only,
  n_snps_interaction_only,
  n_snps_both,
  seed = 5555,
  ethnicity = "EUR",
  seed_effect_size = 101,
  effect_size_mean = 0,
  effect_size_sd = 1,
  out_dir = NULL
) {
  mpgx_assert_file(bim_file, "bim_file")
  mpgx_assert_nonnegative_integer(n_snps_maineffects_only, "n_snps_maineffects_only")
  mpgx_assert_nonnegative_integer(n_snps_interaction_only, "n_snps_interaction_only")
  mpgx_assert_nonnegative_integer(n_snps_both, "n_snps_both")

  n_snps_maineffects_only <- as.integer(n_snps_maineffects_only)
  n_snps_interaction_only <- as.integer(n_snps_interaction_only)
  n_snps_both <- as.integer(n_snps_both)

  bim <- readr::read_tsv(
    bim_file,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )

  odd <- bim[bim$X1 %% 2 == 1 & bim$X1 <= 22, , drop = FALSE]
  num_causal <- n_snps_maineffects_only + n_snps_interaction_only + n_snps_both

  if (num_causal < 0) {
    stop("Number of causal SNPs cannot be negative.", call. = FALSE)
  }
  if (num_causal > nrow(odd)) {
    stop(
      sprintf("Requested %d causal SNPs but only %d odd-chromosome SNPs are available.", num_causal, nrow(odd)),
      call. = FALSE
    )
  }

  idx <- integer(0)
  if (num_causal > 0) {
    set.seed(seed)
    idx <- sample(seq_len(nrow(odd)), num_causal)
  }

  n_main <- n_snps_maineffects_only + n_snps_both
  i_main <- if (n_main > 0) idx[seq_len(n_main)] else integer(0)

  i_inter <- integer(0)
  inter_start <- n_snps_maineffects_only + 1
  if (inter_start <= length(idx)) {
    i_inter <- idx[inter_start:length(idx)]
  }

  set.seed(seed_effect_size)

  main_df <- data.frame(
    SNP = odd$X2[i_main],
    POS = odd$X6[i_main],
    BETA = stats::rnorm(length(i_main), mean = effect_size_mean, sd = effect_size_sd),
    stringsAsFactors = FALSE
  )

  interaction_df <- data.frame(
    SNP = odd$X2[i_inter],
    POS = odd$X6[i_inter],
    BETA = stats::rnorm(length(i_inter), mean = effect_size_mean, sd = effect_size_sd),
    stringsAsFactors = FALSE
  )

  out_files <- NULL
  if (!is.null(out_dir)) {
    out_dir <- mpgx_ensure_dir(out_dir)

    main_file <- file.path(
      out_dir,
      sprintf("SNPs_maineffects.%d.%s.txt", nrow(main_df), ethnicity)
    )
    interaction_file <- file.path(
      out_dir,
      sprintf("SNPs_interaction.%d.%s.txt", nrow(interaction_df), ethnicity)
    )

    mpgx_write_table(main_df, file = main_file, col_names = FALSE, sep = "\t")
    mpgx_write_table(interaction_df, file = interaction_file, col_names = FALSE, sep = "\t")
    out_files <- list(main = main_file, interaction = interaction_file)
  }

  list(main = main_df, interaction = interaction_df, files = out_files)
}

mpgx_read_plink_bundle <- function(bfile_prefix) {
  if (!requireNamespace("BEDMatrix", quietly = TRUE)) {
    stop("The BEDMatrix package is required to read PLINK BED files. Install it via install.packages(\"BEDMatrix\").", call. = FALSE)
  }

  paths <- mpgx_plink_paths(bfile_prefix)
  mpgx_assert_file(paths$bed, "bed file")
  mpgx_assert_file(paths$bim, "bim file")
  mpgx_assert_file(paths$fam, "fam file")

  geno <- BEDMatrix::BEDMatrix(paths$bed)
  bim <- utils::read.table(paths$bim, stringsAsFactors = FALSE)
  fam <- utils::read.table(paths$fam, stringsAsFactors = FALSE)

  colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
  colnames(fam)[1:2] <- c("FID", "IID")

  list(geno = geno, bim = bim, fam = fam, paths = paths)
}

mpgx_read_sscore <- function(sscore_file) {
  if (!file.exists(sscore_file) || file.info(sscore_file)$size == 0) {
    return(data.frame(FID = character(), IID = character(), score = numeric(), stringsAsFactors = FALSE))
  }

  sscore <- tryCatch(
    # Keep legacy script behavior: PLINK2 '#FID' header is skipped as a comment,
    # and column 5 is used as the score component.
    utils::read.table(sscore_file, stringsAsFactors = FALSE),
    error = function(e) {
      utils::read.table(sscore_file, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    }
  )

  if (nrow(sscore) == 0 || ncol(sscore) < 2) {
    return(data.frame(FID = character(), IID = character(), score = numeric(), stringsAsFactors = FALSE))
  }

  score_col <- if (ncol(sscore) >= 5) 5 else ncol(sscore)

  data.frame(
    FID = as.character(sscore[[1]]),
    IID = as.character(sscore[[2]]),
    score = as.numeric(sscore[[score_col]]),
    stringsAsFactors = FALSE
  )
}

mpgx_run_plink2_score <- function(plink2, bfile_prefix, score_file, out_prefix, keep_file = NULL) {
  sscore_file <- paste0(out_prefix, ".sscore")

  if (!file.exists(score_file) || file.info(score_file)$size == 0 || !file.exists(paste0(bfile_prefix, ".fam"))) {
    writeLines(character(0), con = sscore_file)
    return(sscore_file)
  }

  args <- c("--bfile", bfile_prefix, "--score", score_file)
  if (!is.null(keep_file) && file.exists(keep_file)) {
    args <- c(args, "--keep", keep_file)
  }
  args <- c(args, "--out", out_prefix, "--threads", "1")

  mpgx_run_command(command = plink2, args = args)
  sscore_file
}

mpgx_zise <- function(x, only = NULL, by = NULL) {
  if (length(x) < 1) {
    return(x)
  }

  if (is.null(only)) {
    only <- which(!is.na(x))
  } else {
    stopifnot(length(only) == length(x))
    only <- which(only & !is.na(x))
  }

  zx <- rep(NA_real_, length(x))
  if (length(only) >= 1) {
    if (is.null(by)) {
      zx[only] <- stats::qnorm((rank(x[only]) - 0.5) / length(only))
    } else {
      stopifnot(length(by) == length(x))
      by <- factor(by)
      for (by1 in levels(by)) {
        onlyby <- intersect(only, which(by == by1))
        if (length(onlyby) >= 1) {
          zx[onlyby] <- stats::qnorm((rank(x[onlyby]) - 0.5) / length(onlyby))
        }
      }
    }
  }
  zx
}

mpgx_write_simulation_parameters <- function(
  result_dir,
  params
) {
  out <- data.frame(
    parameter = names(params),
    meaning = unname(vapply(params, function(x) x[["meaning"]], FUN.VALUE = character(1))),
    value = unname(vapply(params, function(x) as.character(x[["value"]]), FUN.VALUE = character(1))),
    stringsAsFactors = FALSE
  )
  mpgx_write_table(out, file = file.path(result_dir, "simulationParameters.csv"), col_names = TRUE, sep = ",")
}

#' Prepare ancestry-specific keep files from ancestry metadata
#'
#' @param ancestry_info_file File with three columns: FID, IID, ancestry/ethnicity.
#' @param out_dir Output directory for generated keep files.
#'
#' @return A list with `keep_files`, keyed by ancestry label.
#' @export
mpgx_prepare_ancestry_keep_files <- function(
  ancestry_info_file,
  out_dir
) {
  if (missing(out_dir) || is.null(out_dir) || !nzchar(out_dir)) {
    stop("out_dir must be provided.", call. = FALSE)
  }
  if (missing(ancestry_info_file) || is.null(ancestry_info_file) || !nzchar(ancestry_info_file)) {
    stop("ancestry_info_file must be provided.", call. = FALSE)
  }
  if (!file.exists(ancestry_info_file)) {
    stop(sprintf("ancestry_info_file does not exist: %s", ancestry_info_file), call. = FALSE)
  }

  out_dir <- mpgx_ensure_dir(out_dir)
  ancestry_map <- utils::read.table(ancestry_info_file, stringsAsFactors = FALSE)
  colnames(ancestry_map)[1:3] <- c("FID", "IID", "ancestry")

  keep_files <- lapply(split(ancestry_map, ancestry_map$ancestry), function(df) {
    ances <- unique(as.character(df$ancestry))[1]
    keep_path <- file.path(out_dir, paste0("keep_", ances, ".txt"))
    if (!file.exists(keep_path)) {
      utils::write.table(
        df[, c("FID", "IID")],
        keep_path,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
    }
    keep_path
  })

  list(keep_files = keep_files)
}

#' Simulate phenotypes from PLINK genotype data
#'
#' Refactors and modularizes the original simulation script workflow while
#' preserving key parameter semantics and output file contracts.
#'
#' @param bfile_prefix PLINK prefix (without extension).
#' @param ancestry_info_file File with three columns: FID, IID, ancestry.
#' @param output_dir Parent output directory where `result_label` folder is created.
#' @param result_label Simulation result label.
#' @param seed Seed controlling the single simulation replicate.
#' @param num_pcs Number of PCs to include in qcovar.
#' @param pcs_file Optional PC file (FID, IID, PC1..).
#' @param variance_proportion String or vector for component variances.
#' @param n_snps_maineffects_only Number of G-only causal SNPs.
#' @param n_snps_interaction_only Number of GxE-only causal SNPs.
#' @param n_snps_both Number of SNPs with both effects.
#' @param type_of_environment_variable Environment type.
#' @param contributed_related_individuals Kept for script compatibility.
#' @param type_of_error Error distribution.
#' @param type_of_phenotype Phenotype scale transformation.
#' @param type_of_effect_size Whether effect sizes differ by ethnicity.
#' @param effect_size_mean Mean effect size.
#' @param effect_size_sd SD of effect sizes.
#' @param causal_snps_only_in_eur If `TRUE`, only EUR carries causal signal.
#' @param seed_causal Seed for selecting causal variants.
#' @param configs Optional named list of executable overrides (reserved for compatibility).
#'
#' @return Structured list with result path, single-run seed metadata, output file paths, and causal SNP metadata.
#' @export
mpgx_simulate_phenotypes <- function(
  bfile_prefix,
  ancestry_info_file,
  output_dir,
  result_label = "forTutorialLocal",
  seed = 100001,
  num_pcs = 10,
  pcs_file = NULL,
  variance_proportion = "0.25_0.25_0.05_0.05_0.00_0.40",
  n_snps_maineffects_only = 2,
  n_snps_interaction_only = 2,
  n_snps_both = 6,
  type_of_environment_variable = c("binaryEnvironmentVariable", "continuousEnvironmentVariable"),
  contributed_related_individuals = "none",
  type_of_error = c("normal", "chisquare"),
  type_of_phenotype = c(
    "normal",
    "exponentiallyTransformed",
    "binary",
    "inverseNormalTransformed",
    "INTexponentiallyTransformed"
  ),
  type_of_effect_size = c("differentAcrossEthnicity", "sameAcrossEthnicity"),
  effect_size_mean = 0,
  effect_size_sd = 1,
  causal_snps_only_in_eur = FALSE,
  seed_causal = 5555,
  configs = NULL
) {
  type_of_environment_variable <- match.arg(type_of_environment_variable)
  type_of_error <- match.arg(type_of_error)
  type_of_phenotype <- match.arg(type_of_phenotype)
  type_of_effect_size <- match.arg(type_of_effect_size)

  mpgx_assert_scalar_character(result_label, "result_label")
  mpgx_assert_positive_integer(seed, "seed")
  seed <- as.integer(seed)

  variance_proportion <- mpgx_parse_variance_proportion(variance_proportion)
  weights <- sqrt(variance_proportion)

  paths <- mpgx_plink_paths(bfile_prefix)
  mpgx_assert_file(paths$bed, "bed file")
  mpgx_assert_file(paths$bim, "bim file")
  mpgx_assert_file(paths$fam, "fam file")

  fam <- utils::read.table(paths$fam, stringsAsFactors = FALSE)
  colnames(fam)[1:2] <- c("FID", "IID")

  plink2 <- mpgx_resolve_executable(
    tool = "plink2",
    default = "plink2",
    configs = configs,
    config_key = "plink2.path",
    install_hint = "Install PLINK 2 or pass configs[['plink2.path']] explicitly."
  )

  mpgx_assert_file(ancestry_info_file, "ancestry_info_file")
  id_map <- utils::read.table(ancestry_info_file, stringsAsFactors = FALSE)
  colnames(id_map)[1:3] <- c("FID", "IID", "ethnicity")

  fam_ids <- fam[, c("FID", "IID"), drop = FALSE]
  all_fam <- merge(fam_ids, id_map, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
  if (any(is.na(all_fam$ethnicity))) {
    all_fam$ethnicity[is.na(all_fam$ethnicity)] <- substr(all_fam$FID[is.na(all_fam$ethnicity)], 1, 3)
  }

  if (is.null(pcs_file)) {
    stop("pcs_file cannot be NULL. Please provide a valid PC file path.", call. = FALSE)
  }

  pcs_raw <- utils::read.table(
    pcs_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_cols <- 2 + max(0, num_pcs)
  if (ncol(pcs_raw) < required_cols) {
    stop(
      sprintf(
        "PC file %s has %d columns, but at least %d are required for num_pcs=%d.",
        pcs_file,
        ncol(pcs_raw),
        required_cols,
        num_pcs
      ),
      call. = FALSE
    )
  }

  pcs_df <- pcs_raw[, seq_len(required_cols), drop = FALSE]
  colnames(pcs_df)[1:2] <- c("FID", "IID")
  if (num_pcs > 0) {
    colnames(pcs_df)[3:ncol(pcs_df)] <- paste0("PC", seq_len(num_pcs))
  }
  pcs_df <- merge(all_fam[, c("FID", "IID"), drop = FALSE], pcs_df, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
  pcs_df[is.na(pcs_df)] <- 0

  result_dir <- mpgx_ensure_dir(file.path(output_dir, result_label))

  n_maineffects <- n_snps_maineffects_only + n_snps_both
  n_interaction <- n_snps_interaction_only + n_snps_both

  parameter_meta <- list(
    result_label = list(meaning = "output folder name", value = result_label),
    bfile_prefix = list(meaning = "plink prefix", value = bfile_prefix),
    plink2 = list(meaning = "plink2 executable path", value = plink2),
    num_pcs = list(meaning = "number of principal components", value = num_pcs),
    seed = list(meaning = "seed for the simulation replicate", value = seed),
    variance_proportion = list(meaning = "main,interaction,covariates,stratification,shared,error", value = paste(variance_proportion, collapse = "_")),
    n_snps_maineffects_only = list(meaning = "SNPs with G-only effects", value = n_snps_maineffects_only),
    n_snps_interaction_only = list(meaning = "SNPs with GxE-only effects", value = n_snps_interaction_only),
    n_snps_both = list(meaning = "SNPs with both effects", value = n_snps_both),
    type_of_environment_variable = list(meaning = "environment type", value = type_of_environment_variable),
    contributed_related_individuals = list(meaning = "kept for compatibility", value = contributed_related_individuals),
    type_of_error = list(meaning = "error distribution", value = type_of_error),
    type_of_phenotype = list(meaning = "phenotype type", value = type_of_phenotype),
    type_of_effect_size = list(meaning = "effect-size sharing across ethnicities", value = type_of_effect_size),
    effect_size_mean = list(meaning = "effect-size mean", value = effect_size_mean),
    effect_size_sd = list(meaning = "effect-size SD", value = effect_size_sd),
    causal_snps_only_in_eur = list(meaning = "signal only in EUR", value = causal_snps_only_in_eur),
    pcs_file = list(meaning = "principal components file", value = pcs_file)
  )
  mpgx_write_simulation_parameters(result_dir, parameter_meta)

  ethnicity_levels <- unique(as.character(all_fam$ethnicity))
  ethnicity_levels <- ethnicity_levels[!is.na(ethnicity_levels)]

  causal_by_eth <- vector("list", length(ethnicity_levels))
  names(causal_by_eth) <- ethnicity_levels

  seed_for_effect <- 100
  for (eth in ethnicity_levels) {
    if (type_of_effect_size == "differentAcrossEthnicity") {
      seed_for_effect <- seed_for_effect + 1
    }

    causal_by_eth[[eth]] <- mpgx_select_causal_snps(
      bim_file = paths$bim,
      n_snps_maineffects_only = n_snps_maineffects_only,
      n_snps_interaction_only = n_snps_interaction_only,
      n_snps_both = n_snps_both,
      seed = seed_causal,
      ethnicity = eth,
      seed_effect_size = seed_for_effect,
      effect_size_mean = effect_size_mean,
      effect_size_sd = effect_size_sd,
      out_dir = result_dir
    )
  }

  n_samples <- nrow(all_fam)
  main_scores <- rep(0, n_samples)
  interaction_scores <- rep(0, n_samples)
  all_keys <- paste(all_fam$FID, all_fam$IID, sep = "	")
  bfile_name <- basename(paths$prefix)

  for (eth in ethnicity_levels) {
    if (isTRUE(causal_snps_only_in_eur) && eth != "EUR") {
      next
    }

    idx <- which(all_fam$ethnicity == eth)
    if (length(idx) == 0) {
      next
    }

    bfile_eth <- paste0(paths$prefix, ".", eth)
    keep_file <- NULL

    if (!file.exists(paste0(bfile_eth, ".fam"))) {
      bfile_eth <- paths$prefix
      keep_file <- file.path(result_dir, sprintf("tmp.keep.%s.txt", eth))
      mpgx_write_table(
        all_fam[idx, c("FID", "IID"), drop = FALSE],
        file = keep_file,
        col_names = FALSE,
        sep = "	"
      )
    }

    main_out_prefix <- file.path(result_dir, sprintf("%s.main.%d.%s", bfile_name, n_maineffects, eth))
    interaction_out_prefix <- file.path(result_dir, sprintf("%s.interaction.%d.%s", bfile_name, n_interaction, eth))

    main_sscore <- mpgx_run_plink2_score(
      plink2 = plink2,
      bfile_prefix = bfile_eth,
      score_file = causal_by_eth[[eth]]$files$main,
      out_prefix = main_out_prefix,
      keep_file = keep_file
    )

    interaction_sscore <- mpgx_run_plink2_score(
      plink2 = plink2,
      bfile_prefix = bfile_eth,
      score_file = causal_by_eth[[eth]]$files$interaction,
      out_prefix = interaction_out_prefix,
      keep_file = keep_file
    )

    if (!is.null(keep_file) && file.exists(keep_file)) {
      unlink(keep_file)
    }

    main_df <- mpgx_read_sscore(main_sscore)
    interaction_df <- mpgx_read_sscore(interaction_sscore)

    if (nrow(main_df) > 0) {
      main_match <- match(all_keys[idx], paste(main_df$FID, main_df$IID, sep = "	"))
      main_values <- main_df$score[main_match]
      main_values[is.na(main_values)] <- 0
      main_scores[idx] <- main_values
    }

    if (nrow(interaction_df) > 0) {
      interaction_match <- match(all_keys[idx], paste(interaction_df$FID, interaction_df$IID, sep = "	"))
      interaction_values <- interaction_df$score[interaction_match]
      interaction_values[is.na(interaction_values)] <- 0
      interaction_scores[idx] <- interaction_values
    }
  }


  set.seed(666)
  continuous_covariate <- sample(20:60, n_samples, replace = TRUE)
  binary_covariate <- sample(c(0, 1), n_samples, replace = TRUE)

  if (type_of_environment_variable == "continuousEnvironmentVariable") {
    envir <- stats::rnorm(n_samples, mean = 1, sd = 2)
  } else {
    envir <- sample(c(0, 1), n_samples, replace = TRUE)
  }

  all_components <- data.frame(
    FID = all_fam$FID,
    IID = all_fam$IID,
    ethnicity = all_fam$ethnicity,
    main = main_scores,
    interaction = interaction_scores,
    stringsAsFactors = FALSE
  )

  if (length(unique(all_components$ethnicity)) > 1) {
    map <- c(AFR = 0, EAS = 0.4, EUR = -1, SAS = 2)
    strat <- unname(map[as.character(all_components$ethnicity)])
    strat[is.na(strat)] <- 0
    all_components$scaled_stratification <- mpgx_safe_scale(strat)
  } else {
    all_components$scaled_stratification <- 0
  }

  all_components$continuousCovariate <- continuous_covariate
  all_components$binaryCovariate <- binary_covariate
  all_components$envir <- envir

  covariates_effects <- cbind(
    mpgx_safe_scale(continuous_covariate),
    mpgx_safe_scale(binary_covariate),
    mpgx_safe_scale(envir)
  ) %*% c(-0.5, -0.6, 0.4)

  all_components$scaled_covariates_effects <- mpgx_safe_scale(as.numeric(covariates_effects))
  all_components$scaled_shared_effects <- 0
  all_components$scaled_main <- mpgx_safe_scale(all_components$main)

  interaction_centered <- mpgx_safe_scale(all_components$interaction, scale = FALSE)
  envir_centered <- mpgx_safe_scale(all_components$envir, scale = FALSE)
  interaction_term <- interaction_centered * envir_centered
  all_components$scaled_interaction <- mpgx_safe_scale(interaction_term)

  set.seed(seed)
  all_components$error <- if (type_of_error == "normal") {
    mpgx_safe_scale(stats::rnorm(n_samples))
  } else {
    mpgx_safe_scale(stats::rchisq(n = n_samples, df = 1))
  }

  design <- cbind(
    main = all_components$scaled_main,
    interaction = all_components$scaled_interaction,
    covariates = all_components$scaled_covariates_effects,
    stratification = all_components$scaled_stratification,
    shared = all_components$scaled_shared_effects,
    error = all_components$error
  )

  if (n_maineffects == 0) {
    design <- design[, setdiff(colnames(design), "main"), drop = FALSE]
  }
  if (n_interaction == 0) {
    design <- design[, setdiff(colnames(design), "interaction"), drop = FALSE]
  }

  weight_subset <- weights[colnames(design)]
  temp_phenotype <- as.numeric(design %*% weight_subset)

  phenotype <- switch(
    type_of_phenotype,
    normal = temp_phenotype,
    exponentiallyTransformed = exp(temp_phenotype),
    binary = as.numeric(temp_phenotype > stats::median(temp_phenotype)),
    inverseNormalTransformed = mpgx_zise(temp_phenotype),
    INTexponentiallyTransformed = mpgx_zise(exp(temp_phenotype))
  )

  sim_covar <- data.frame(FID = all_components$FID, IID = all_components$IID, binaryCovariate = binary_covariate)
  sim_qcovar <- data.frame(
    FID = all_components$FID,
    IID = all_components$IID,
    continuousCovariate = continuous_covariate,
    pcs_df[, setdiff(colnames(pcs_df), c("FID", "IID")), drop = FALSE],
    stringsAsFactors = FALSE
  )
  sim_envir <- data.frame(FID = all_components$FID, IID = all_components$IID, envir = envir)
  sim_phenotype <- data.frame(FID = all_components$FID, IID = all_components$IID, phenotype = phenotype)

  if (num_pcs > 0) {
    pc_cols <- paste0("PC", seq_len(num_pcs))
    fit <- stats::lm(phenotype ~ ., data = data.frame(phenotype = phenotype, sim_qcovar[, pc_cols, drop = FALSE]))
    adj_r2 <- summary(fit)$adj.r.squared
    writeLines(as.character(adj_r2), con = file.path(result_dir, sprintf("sim_stratification_adj_r2_seed_%d.txt", seed)))
  }

  covar_file <- file.path(result_dir, "sim_covar.txt")
  qcovar_file <- file.path(result_dir, "sim_qcovar.txt")
  envir_file <- file.path(result_dir, "sim_envir.txt")

  mpgx_write_table(sim_covar, file = covar_file, col_names = FALSE, sep = "\t")
  mpgx_write_table(sim_qcovar, file = qcovar_file, col_names = FALSE, sep = "\t")
  mpgx_write_table(sim_envir, file = envir_file, col_names = FALSE, sep = "\t")

  phenotype_file <- file.path(result_dir, sprintf("sim_phenotype_seed_%d.txt", seed))
  mpgx_write_table(sim_phenotype, file = phenotype_file, col_names = FALSE, sep = "\t")

  sim_for_gem <- Reduce(
    function(x, y) merge(x, y, by = c("FID", "IID"), all = TRUE, sort = FALSE),
    list(sim_covar, sim_qcovar, sim_envir, sim_phenotype)
  )
  mpgx_write_table(
    sim_for_gem,
    file = file.path(result_dir, sprintf("sim_phenotype_forGEM_seed_%d.txt", seed)),
    col_names = TRUE,
    sep = ","
  )

  sample_ids <- fam[, c("FID", "IID"), drop = FALSE]

  covar_for_lemma <- merge(sample_ids, sim_covar, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
  qcovar_for_lemma <- merge(sample_ids, sim_qcovar, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
  envir_for_lemma <- merge(sample_ids, sim_envir, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
  phenotype_for_lemma <- merge(sample_ids, sim_phenotype, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)

  covar_qcovar_lemma <- cbind(
    covar_for_lemma[, setdiff(colnames(covar_for_lemma), c("FID", "IID")), drop = FALSE],
    qcovar_for_lemma[, setdiff(colnames(qcovar_for_lemma), c("FID", "IID")), drop = FALSE]
  )

  mpgx_write_table(covar_qcovar_lemma, file = file.path(result_dir, "sim_covar_qcovar_forLEMMA.txt"), col_names = TRUE, sep = "\t")
  mpgx_write_table(
    envir_for_lemma[, setdiff(colnames(envir_for_lemma), c("FID", "IID")), drop = FALSE],
    file = file.path(result_dir, "sim_envir_forLEMMA.txt"),
    col_names = TRUE,
    sep = "\t"
  )
  mpgx_write_table(
    phenotype_for_lemma[, setdiff(colnames(phenotype_for_lemma), c("FID", "IID")), drop = FALSE],
    file = file.path(result_dir, sprintf("sim_phenotype_forLEMMA_seed_%d.txt", seed)),
    col_names = TRUE,
    sep = "\t"
  )

  causal_main <- unique(unlist(lapply(causal_by_eth, function(x) x$main$SNP), use.names = FALSE))
  causal_interaction <- unique(unlist(lapply(causal_by_eth, function(x) x$interaction$SNP), use.names = FALSE))

  list(
    result_label = result_label,
    result_dir = result_dir,
    seed = seed,
    phenotype_file = phenotype_file,
    covar_file = covar_file,
    qcovar_file = qcovar_file,
    envir_file = envir_file,
    causal_by_ethnicity = causal_by_eth,
    causal_main_snps = causal_main,
    causal_interaction_snps = causal_interaction,
    n_maineffects = n_maineffects,
    n_interaction = n_interaction
  )
}
