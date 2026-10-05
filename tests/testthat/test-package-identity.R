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
