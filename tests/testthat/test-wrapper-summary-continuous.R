testthat::test_that("continuous wrapper summaries match bundled toy-data references", {
  mpgx_test_prepare_runtime_env()

  configs <- mpgx_default_tool_configs()
  inputs <- mpgx_test_toy_inputs()

  out_dir <- file.path(tempdir(), paste0("mpgx_test_continuous_", Sys.getpid()))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  lr <- mpgx_run_lr_sandwich(
    bfile_prefix = inputs$toy_prefix,
    phenotype_file = inputs$continuous$phenotype_file,
    qcovar_file = inputs$continuous$qcovar_file,
    covar_file = inputs$continuous$covar_file,
    envir_file = inputs$continuous$envir_file,
    output_prefix = file.path(out_dir, "continuous.LRallSandwich"),
    threads = 1,
    configs = configs
  )

  fastgwa <- mpgx_run_fastgwa_ge(
    bfile_prefix = inputs$toy_prefix,
    grm_sparse_prefix = inputs$sparse_grm_prefix,
    phenotype_file = inputs$continuous$phenotype_file,
    qcovar_file = inputs$continuous$qcovar_file,
    covar_file = inputs$continuous$covar_file,
    envir_file = inputs$continuous$envir_file,
    output_prefix = file.path(out_dir, "continuous.fastGWA_GE"),
    threads = 1,
    configs = configs
  )

  lemma <- suppressWarnings(
    mpgx_run_lemma(
      bgen_file = inputs$bgen_file,
      phenotype_file = inputs$continuous$lemma_phenotype_file,
      envir_file = inputs$continuous$lemma_envir_file,
      covar_file = inputs$continuous$lemma_covar_file,
      output_prefix = file.path(out_dir, "continuous.LEMMA"),
      n_tasks = 1,
      configs = configs
    )
  )

  genesis_file <- file.path(out_dir, "continuous.GENESIS.txt")
  genesis <- mpgx_run_genesis(
    gds_file = inputs$gds_file,
    sparse_grm_prefix = inputs$sparse_grm_prefix,
    phenotype_file = inputs$continuous$phenotype_file,
    qcovar_file = inputs$continuous$qcovar_file,
    covar_file = inputs$continuous$covar_file,
    envir_file = inputs$continuous$envir_file,
    output_file = genesis_file,
    attach_result = TRUE
  )

  mrmega_lr <- mpgx_run_mrmega_for_lr(
    ancestry_result_files = inputs$lr_by_eth_files,
    output_prefix_base = file.path(out_dir, "continuous.MRMEGA_LRallSandwich"),
    configs = configs
  )

  meta_file <- file.path(out_dir, "continuous.meta_LRallSandwich.txt")
  meta_lr <- mpgx_run_fe_meta_for_lr(
    ancestry_result_files = inputs$lr_by_eth_files,
    ancestry_labels = names(inputs$lr_by_eth_files),
    output_file = meta_file,
    attach_result = TRUE
  )

  testthat::expect_true(file.exists(lr$gwas_file))
  testthat::expect_true(file.exists(fastgwa$gwas_file))
  testthat::expect_true(file.exists(lemma$gwas_file))
  testthat::expect_true(file.exists(genesis_file))
  testthat::expect_true(file.exists(mrmega_lr$test_G))
  testthat::expect_true(file.exists(mrmega_lr$test_G_by_E))
  testthat::expect_true(file.exists(meta_file))

  method_results <- list(
    LRallSandwich = lr,
    meta_LRallSandwich = meta_lr,
    MRMEGA_LRallSandwich = mrmega_lr,
    fastGWA_GE = fastgwa,
    GENESIS = genesis,
    LEMMA = lemma
  )

  causal_snps <- mpgx_test_collect_causal_snps("continuous")
  pc_num <- max(length(inputs$lr_by_eth_files) - 3, 0)

  actual_summary <- mpgx_test_summarize_method_list(
    method_results = method_results,
    causal_snps = causal_snps,
    scenario = "continuous",
    seed = 100001,
    mrmega_methods = "MRMEGA_LRallSandwich",
    mrmega_pc_num = pc_num
  )

  expected_summary <- mpgx_test_summarize_method_list(
    method_results = mpgx_test_reference_methods_continuous(),
    causal_snps = causal_snps,
    scenario = "continuous",
    seed = 100001,
    mrmega_methods = "MRMEGA_LRallSandwich",
    mrmega_pc_num = pc_num
  )

  mpgx_test_expect_summary_equal(actual_summary, expected_summary, tolerance = 1e-4)
})
