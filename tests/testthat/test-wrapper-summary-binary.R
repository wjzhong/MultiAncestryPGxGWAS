testthat::test_that("binary wrapper summaries match bundled toy-data references", {
  mpgx_test_prepare_runtime_env()

  configs <- mpgx_default_tool_configs()
  inputs <- mpgx_test_toy_inputs()

  out_dir <- file.path(tempdir(), paste0("mpgx_test_binary_", Sys.getpid()))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fastgwa_bin <- mpgx_run_fastgwa_ge(
    bfile_prefix = inputs$toy_prefix,
    grm_sparse_prefix = inputs$sparse_grm_prefix,
    phenotype_file = inputs$binary$phenotype_file,
    qcovar_file = inputs$binary$qcovar_file,
    covar_file = inputs$binary$covar_file,
    envir_file = inputs$binary$envir_file,
    output_prefix = file.path(out_dir, "binary.fastGWA_GE_naive"),
    threads = 1,
    configs = configs
  )

  lemma_bin <- suppressWarnings(
    mpgx_run_lemma(
      bgen_file = inputs$bgen_file,
      phenotype_file = inputs$binary$lemma_phenotype_file,
      envir_file = inputs$binary$lemma_envir_file,
      covar_file = inputs$binary$lemma_covar_file,
      output_prefix = file.path(out_dir, "binary.LEMMA_naive"),
      n_tasks = 1,
      configs = configs
    )
  )

  genesis_file <- file.path(out_dir, "binary.GENESIS.txt")
  genesis_bin <- mpgx_run_genesis(
    gds_file = inputs$gds_file,
    sparse_grm_prefix = inputs$sparse_grm_prefix,
    phenotype_file = inputs$binary$phenotype_file,
    qcovar_file = inputs$binary$qcovar_file,
    covar_file = inputs$binary$covar_file,
    envir_file = inputs$binary$envir_file,
    output_file = genesis_file,
    attach_result = TRUE
  )

  gem_bin <- mpgx_run_gem(
    bfile_prefix = inputs$toy_prefix,
    phenotype_file = inputs$binary$gem_phenotype_file,
    output_prefix = file.path(out_dir, "binary.GEM"),
    configs = configs
  )

  mrmega_inputs <- mpgx_prepare_mrmega_gem_inputs(
    bfile_prefix = inputs$toy_prefix,
    gem_phenotype_file = inputs$binary$gem_phenotype_file,
    keep_files = inputs$keep_files_binary,
    output_dir = file.path(out_dir, "binary_gem_by_ancestry"),
    configs = configs
  )

  mrmega_gem <- mpgx_run_mrmega_for_gem(
    ancestry_result_files = mrmega_inputs,
    output_prefix_base = file.path(out_dir, "binary.MRMEGA_GEM"),
    configs = configs
  )

  testthat::expect_true(file.exists(fastgwa_bin$gwas_file))
  testthat::expect_true(file.exists(lemma_bin$gwas_file))
  testthat::expect_true(file.exists(genesis_file))
  testthat::expect_true(file.exists(gem_bin$gwas_file))
  testthat::expect_true(file.exists(mrmega_gem$test_G))
  testthat::expect_true(file.exists(mrmega_gem$test_G_by_E))

  method_results <- list(
    fastGWA_GE_naive = fastgwa_bin,
    LEMMA_naive = lemma_bin,
    GENESIS = genesis_bin,
    GEM = gem_bin,
    MRMEGA_GEM = mrmega_gem
  )

  causal_snps <- mpgx_test_collect_causal_snps("binary")
  pc_num <- max(length(mrmega_inputs) - 3, 0)

  actual_summary <- mpgx_test_summarize_method_list(
    method_results = method_results,
    causal_snps = causal_snps,
    scenario = "binary",
    seed = 100001,
    mrmega_methods = "MRMEGA_GEM",
    mrmega_pc_num = pc_num
  )

  expected_summary <- mpgx_test_summarize_method_list(
    method_results = mpgx_test_reference_methods_binary(),
    causal_snps = causal_snps,
    scenario = "binary",
    seed = 100001,
    mrmega_methods = "MRMEGA_GEM",
    mrmega_pc_num = pc_num
  )

  mpgx_test_expect_summary_equal(actual_summary, expected_summary, tolerance = 1e-4)
})
