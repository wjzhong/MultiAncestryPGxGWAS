test_that("example-data lookup belongs to the package namespace", {
  expect_identical(
    unname(getNamespaceName(environment(mpgx_example_file))),
    "MultiAncestryPGxGWAS"
  )

  check_from_unrelated_directory <- function() {
    work_dir <- tempfile("mpgx-installed-data-")
    dir.create(work_dir)
    previous_dir <- getwd()
    on.exit({
      setwd(previous_dir)
      unlink(work_dir, recursive = TRUE)
    }, add = TRUE)
    setwd(work_dir)

    expected <- system.file(
      "extdata", "genotype_data", "toyData.bed",
      package = "MultiAncestryPGxGWAS", mustWork = TRUE
    )
    expect_identical(
      normalizePath(mpgx_example_file("genotype_data/toyData.bed")),
      normalizePath(expected)
    )
    expect_identical(
      normalizePath(mpgx_example_file("toyData.bed")),
      normalizePath(expected)
    )
    expect_error(
      mpgx_example_file("genotype_data/nonexistent-file-test.bed"),
      "toy_data example file not found", fixed = TRUE
    )
  }

  check_from_unrelated_directory()
})

test_that("bundled example path lists use the package vignette directory", {
  package_name <- unname(getNamespaceName(environment(mpgx_example_file)))
  directory_component <- paste0("/", package_name, "_vignette/")
  path_lists <- list(
    list(
      file = "genotype_data/toyData.merge_list.txt",
      n_rows = 3L, n_columns = 3L
    ),
    list(
      file = "model_results/binary/mrmega_gem/binary.MRMEGA_GEM.in",
      n_rows = 4L, n_columns = 1L
    ),
    list(
      file = "model_results/continuous/mrmega_lr/continuous.MRMEGA_LRallSandwich.in",
      n_rows = 4L, n_columns = 1L
    )
  )

  for (path_list in path_lists) {
    lines <- readLines(mpgx_example_file(path_list$file), warn = FALSE)
    paths <- strsplit(trimws(lines), "[[:space:]]+")
    expect_length(paths, path_list$n_rows)
    expect_true(
      all(lengths(paths) == path_list$n_columns), info = path_list$file
    )
    expect_true(
      all(grepl(
        directory_component, unlist(paths, use.names = FALSE), fixed = TRUE
      )),
      info = path_list$file
    )
  }
})
