mpgx_test_prepare_runtime_env <- function() {
  tmp_dir <- file.path(tempdir(), "mpgx_mpi_tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  Sys.setenv(TMPDIR = tmp_dir)
  Sys.setenv(OMPI_MCA_orte_tmpdir_base = tmp_dir)
  Sys.setenv(OMPI_MCA_tmpdir_base = tmp_dir)

  invisible(tmp_dir)
}

mpgx_test_toy_inputs <- function() {
  ancestries <- c("AFR", "EAS", "EUR", "SAS")

  strip_suffix <- function(path, suffix) {
    if (endsWith(path, suffix)) {
      return(substr(path, 1, nchar(path) - nchar(suffix)))
    }
    path
  }

  toy_prefix <- strip_suffix(
    mpgx_example_file("genotype_data/toyData.bed"),
    ".bed"
  )
  sparse_grm_prefix <- strip_suffix(
    mpgx_example_file("genotype_artifacts/toyData.sp_grm.grm.sp"),
    ".grm.sp"
  )

  lr_by_eth_files <- setNames(
    vapply(
      ancestries,
      function(ances) {
        mpgx_example_file(
          file.path(
            "model_results",
            "continuous",
            "lr_sandwich",
            paste0("continuous.LRallSandwich.", ances, ".fastGWA")
          )
        )
      },
      FUN.VALUE = character(1)
    ),
    ancestries
  )

  keep_files_binary <- setNames(
    vapply(
      ancestries,
      function(ances) {
        mpgx_example_file(file.path("model_results", "binary", paste0("keep_", ances, ".txt")))
      },
      FUN.VALUE = character(1)
    ),
    ancestries
  )

  list(
    toy_prefix = toy_prefix,
    sparse_grm_prefix = sparse_grm_prefix,
    bgen_file = mpgx_example_file("genotype_artifacts/toyData.bgen"),
    gds_file = mpgx_example_file("genotype_artifacts/toyData.gds"),
    continuous = list(
      phenotype_file = mpgx_example_file("simulated_data/continuous/sim_phenotype_seed_100001.txt"),
      qcovar_file = mpgx_example_file("simulated_data/continuous/sim_qcovar.txt"),
      covar_file = mpgx_example_file("simulated_data/continuous/sim_covar.txt"),
      envir_file = mpgx_example_file("simulated_data/continuous/sim_envir.txt"),
      lemma_phenotype_file = mpgx_example_file("simulated_data/continuous/sim_phenotype_forLEMMA_seed_100001.txt"),
      lemma_covar_file = mpgx_example_file("simulated_data/continuous/sim_covar_qcovar_forLEMMA.txt"),
      lemma_envir_file = mpgx_example_file("simulated_data/continuous/sim_envir_forLEMMA.txt")
    ),
    binary = list(
      phenotype_file = mpgx_example_file("simulated_data/binary/sim_phenotype_seed_100001.txt"),
      qcovar_file = mpgx_example_file("simulated_data/binary/sim_qcovar.txt"),
      covar_file = mpgx_example_file("simulated_data/binary/sim_covar.txt"),
      envir_file = mpgx_example_file("simulated_data/binary/sim_envir.txt"),
      gem_phenotype_file = mpgx_example_file("simulated_data/binary/sim_phenotype_forGEM_seed_100001.txt"),
      lemma_phenotype_file = mpgx_example_file("simulated_data/binary/sim_phenotype_forLEMMA_seed_100001.txt"),
      lemma_covar_file = mpgx_example_file("simulated_data/binary/sim_covar_qcovar_forLEMMA.txt"),
      lemma_envir_file = mpgx_example_file("simulated_data/binary/sim_envir_forLEMMA.txt")
    ),
    lr_by_eth_files = lr_by_eth_files,
    keep_files_binary = keep_files_binary
  )
}

mpgx_test_collect_causal_snps <- function(outcome = c("continuous", "binary")) {
  outcome <- match.arg(outcome)
  ancestries <- c("AFR", "EAS", "EUR", "SAS")

  main_files <- file.path(
    "simulated_data",
    outcome,
    paste0("SNPs_maineffects.8.", ancestries, ".txt")
  )
  interaction_files <- file.path(
    "simulated_data",
    outcome,
    paste0("SNPs_interaction.8.", ancestries, ".txt")
  )

  read_snp_col <- function(rel_path) {
    abs_path <- mpgx_example_file(rel_path)
    tab <- utils::read.table(abs_path, stringsAsFactors = FALSE)
    as.character(tab[[1]])
  }

  unique(unlist(lapply(c(main_files, interaction_files), read_snp_col), use.names = FALSE))
}

mpgx_test_wrap_gwas_df <- function(gwas_df) {
  list(gwas_df = as.data.frame(gwas_df, stringsAsFactors = FALSE))
}

mpgx_test_reference_methods_continuous <- function() {
  list(
    LRallSandwich = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "LRallSandwich",
        files = mpgx_example_file("model_results/continuous/lr_sandwich/continuous.LRallSandwich.fastGWA")
      )
    ),
    meta_LRallSandwich = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "LRallSandwich",
        files = mpgx_example_file("model_results/continuous/fe_meta_lr/continuous.meta_LRallSandwich.txt")
      )
    ),
    MRMEGA_LRallSandwich = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "MRMEGA_LRallSandwich",
        files = list(
          test_G = mpgx_example_file("model_results/continuous/mrmega_lr/continuous.MRMEGA_LRallSandwich.test_G.result"),
          test_G_by_E = mpgx_example_file("model_results/continuous/mrmega_lr/continuous.MRMEGA_LRallSandwich.test_G_by_E.result")
        )
      )
    ),
    fastGWA_GE = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "fastGWA_GE",
        files = mpgx_example_file("model_results/continuous/fastgwa_ge/continuous.fastGWA_GE.fastGWA")
      )
    ),
    GENESIS = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "GENESIS",
        files = mpgx_example_file("model_results/continuous/genesis/continuous.GENESIS.txt")
      )
    ),
    LEMMA = mpgx_test_wrap_gwas_df(
      suppressWarnings(
        mpgx_parse_method_result_df(
          method = "LEMMA",
          files = mpgx_example_file("model_results/continuous/lemma/continuous_loco_pvals.LEMMA")
        )
      )
    )
  )
}

mpgx_test_reference_methods_binary <- function() {
  list(
    fastGWA_GE_naive = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "fastGWA_GE",
        files = mpgx_example_file("model_results/binary/fastgwa_ge_naive/binary.fastGWA_GE_naive.fastGWA")
      )
    ),
    LEMMA_naive = mpgx_test_wrap_gwas_df(
      suppressWarnings(
        mpgx_parse_method_result_df(
          method = "LEMMA",
          files = mpgx_example_file("model_results/binary/lemma_naive/binary_loco_pvals.LEMMA_naive")
        )
      )
    ),
    GENESIS = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "GENESIS",
        files = mpgx_example_file("model_results/binary/genesis/binary.GENESIS.txt")
      )
    ),
    GEM = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "GEM",
        files = mpgx_example_file("model_results/binary/gem/binary.GEM")
      )
    ),
    MRMEGA_GEM = mpgx_test_wrap_gwas_df(
      mpgx_parse_method_result_df(
        method = "MRMEGA_GEM",
        files = list(
          test_G = mpgx_example_file("model_results/binary/mrmega_gem/binary.MRMEGA_GEM.test_G.result"),
          test_G_by_E = mpgx_example_file("model_results/binary/mrmega_gem/binary.MRMEGA_GEM.test_G_by_E.result")
        )
      )
    )
  )
}

mpgx_test_summarize_method_list <- function(
  method_results,
  causal_snps,
  scenario,
  seed,
  mrmega_methods = character(),
  mrmega_pc_num = NULL
) {
  summaries <- lapply(names(method_results), function(method_name) {
    method_df <- mpgx_get_method_gwas_df(method_results[[method_name]], label = method_name)
    mrmega_pc <- if (method_name %in% mrmega_methods) mrmega_pc_num else NULL
    has_non_missing_2df <- "P_2df" %in% colnames(method_df) && any(!is.na(method_df$P_2df))
    summary_df <- mpgx_summarize_association(
      result = method_df,
      causal_snps = causal_snps,
      rm_na_2df = has_non_missing_2df,
      mrmega_pc_num = mrmega_pc
    )

    cbind(
      data.frame(
        method = method_name,
        scenario = scenario,
        seed = seed,
        stringsAsFactors = FALSE
      ),
      summary_df
    )
  })

  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  out
}

mpgx_test_order_summary <- function(df) {
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  out <- out[order(out$method), , drop = FALSE]
  rownames(out) <- NULL
  out
}

mpgx_test_expect_summary_equal <- function(actual, expected, tolerance = 1e-4) {
  actual <- mpgx_test_order_summary(actual)
  expected <- mpgx_test_order_summary(expected)

  testthat::expect_equal(colnames(actual), colnames(expected))
  testthat::expect_equal(nrow(actual), nrow(expected))

  char_cols <- names(actual)[vapply(actual, function(x) is.character(x) || is.factor(x), logical(1))]
  for (nm in char_cols) {
    testthat::expect_equal(
      as.character(actual[[nm]]),
      as.character(expected[[nm]]),
      info = paste("Mismatch in", nm)
    )
  }

  num_cols <- names(actual)[vapply(actual, is.numeric, logical(1))]
  for (nm in num_cols) {
    testthat::expect_equal(
      unname(actual[[nm]]),
      unname(expected[[nm]]),
      tolerance = tolerance,
      info = paste("Mismatch in", nm)
    )
  }
}
