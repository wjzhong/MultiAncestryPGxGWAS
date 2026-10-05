# MultiAncestryPGxGWAS

R package for multi-ancestry genome-wide genotype-by-environment (GxE) and genotype-by-treatment (GxT) interaction analysis.

License: GPL (>= 3)

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [External Software and Tool Configuration](#external-software-and-tool-configuration)
- [Tutorials](#tutorials)
- [Available Analysis Wrappers](#available-analysis-wrappers)
- [End-to-End Example with Toy Data](#end-to-end-example-with-toy-data)
- [Script Wrappers](#script-wrappers)
- [Common Issues](#common-issues)
- [References](#references)

## Overview

MultiAncestryPGxGWAS provides wrapper functions to run and compare multiple GxE/GxT GWAS methods in a consistent workflow. It includes packaged toy data under `inst/extdata` so users can run end-to-end examples quickly.

## Installation

```r
# Install devtools if it is not already available
install.packages("devtools")

# Install the R package MultiAncestryPGxGWAS
devtools::install_github("wjzhong/MultiAncestryPGxGWAS")
```

## Prerequisites

- An R environment with support for building packages.
- External command-line tools listed in [External Software and Tool Configuration](#external-software-and-tool-configuration).
- For LEMMA parallel runs, OpenMPI `mpirun` is optional but recommended.

## Quick Start

Use this minimal example to validate tool paths and run one method on packaged toy data.

```r
library(MultiAncestryPGxGWAS)

configs <- list(
  plink1.9.path = "plink",
  plink2.path = "plink2",
  flashpca.path = "flashpca_x86-64",
  gcta.path = "gcta",
  lemma.path = "lemma_1_0_4",
  gem.path = "GEM_1.4.5_static",
  bgenix.path = "bgenix",
  mrmega.path = "MR-MEGA",
  mpirun.path = "mpirun"
)
# Optionally override executable locations when tools are not on PATH.
# configs[["gcta.path"]] <- "/abs/path/to/gcta"

mpgx_validate_tool_configs(configs)

out_dir <- file.path(tempdir(), "mpgx_quickstart")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

toy_prefix <- sub("\\.bed$", "", mpgx_example_file("genotype_data/toyData.bed"))
pheno <- mpgx_example_file("simulated_data/continuous/sim_phenotype_seed_100001.txt")
qcovar <- mpgx_example_file("simulated_data/continuous/sim_qcovar.txt")
covar <- mpgx_example_file("simulated_data/continuous/sim_covar.txt")
envir <- mpgx_example_file("simulated_data/continuous/sim_envir.txt")

lr <- mpgx_run_lr_sandwich(
  bfile_prefix = toy_prefix,
  phenotype_file = pheno,
  qcovar_file = qcovar,
  covar_file = covar,
  envir_file = envir,
  output_prefix = file.path(out_dir, "quickstart.LRallSandwich"),
  configs = configs
)

str(lr, max.level = 1)
```

## External Software and Tool Configuration

The following versions are known to work with this package.

| Tool | Known Version | Notes |
| --- | --- | --- |
| PLINK1.9 | v1.90b5.1 | Used by LR-Sandwich and related PLINK-based steps |
| PLINK2 | v2.00a2LM | Used by simulation/scoring workflows |
| GCTA | v1.95.3 | Used by LR-Sandwich and fastGWA-GE wrappers |
| LEMMA | 1.0.2+ | Used by `mpgx_run_lemma()` |
| GEM | 1.4.5+ | Used by `mpgx_run_gem()` |
| bgenix | 1.1.7 | Used for BGEN indexing and related steps |
| MR-MEGA | 0.2 | Used by MR-MEGA wrappers |
| mpirun (OpenMPI, optional) | 4.1.7 | Optional helper for launching LEMMA with MPI |

Download links:
- PLINK1.9: https://www.cog-genomics.org/plink/1.9/
- PLINK2: https://www.cog-genomics.org/plink/2.0/
- GCTA: https://yanglab.westlake.edu.cn/software/gcta/#Download
- LEMMA: https://github.com/mkerin/LEMMA.git
- GEM: https://github.com/large-scale-gxe-methods/GEM.git
- bgenix: https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk
- MR-MEGA: https://genomics.ut.ee/en/tools
- OpenMPI: https://docs.open-mpi.org/en/main/

Users can either keep these tools on `PATH` or provide explicit locations in a `configs` list.

Path resolution order used by wrappers:
1. Function-specific `*_path` argument
2. `configs` list value
3. `PATH`

Example explicit configuration:

```r
configs <- list(
  plink1.9.path = "/path/to/plink",
  plink2.path = "/path/to/plink2",
  flashpca.path = "/path/to/flashpca_x86-64",
  gcta.path = "/path/to/gcta",
  lemma.path = "/path/to/lemma_1_0_4",
  gem.path = "/path/to/GEM",
  bgenix.path = "/path/to/bgenix",
  mrmega.path = "/path/to/MR-MEGA",
  mpirun.path = "/path/to/mpirun"
)

mpgx_validate_tool_configs(configs)
```

Default executable names used by the wrappers:

```r
configs <- list(
  plink1.9.path = "plink",
  plink2.path = "plink2",
  flashpca.path = "flashpca_x86-64",
  gcta.path = "gcta",
  lemma.path = "lemma_1_0_4",
  gem.path = "GEM_1.4.5_static",
  bgenix.path = "bgenix",
  mrmega.path = "MR-MEGA",
  mpirun.path = "mpirun"
)
```

Pass `configs` to any method wrapper, for example:

```r
mpgx_run_lr_sandwich(
  bfile_prefix = "/path/to/data/geno",
  phenotype_file = "/path/to/pheno.tsv",
  qcovar_file = "/path/to/qcovar.tsv",
  covar_file = "/path/to/covar.tsv",
  envir_file = "/path/to/envir.tsv",
  output_prefix = "/path/to/out/gcta_lr",
  configs = configs
)
```

## Tutorials

Tutorials are available in both rendered and source formats:

- `vignettes/MultiAncestryPGxGWAS-simulate-apply.html`
- `vignettes/MultiAncestryPGxGWAS-simulate-apply.pdf`
- `vignettes/MultiAncestryPGxGWAS-simulate-apply.Rmd`
- `vignettes/MultiAncestryPGxGWAS-quick-start.html`
- `vignettes/MultiAncestryPGxGWAS-quick-start.pdf`
- `vignettes/MultiAncestryPGxGWAS-quick-start.Rmd`

## Available Analysis Wrappers

| Method | Wrapper | Typical Use |
| --- | --- | --- |
| Linear regression with robust variance | `mpgx_run_lr_sandwich()` | Baseline GxE/GxT test using GCTA |
| fastGWA-GE | `mpgx_run_fastgwa_ge()` | Sparse-GRM accelerated mixed model |
| LEMMA | `mpgx_run_lemma()` | Bayesian whole-genome interaction modeling |
| GENESIS | `mpgx_run_genesis()` | GDS-based mixed model framework |
| GEM | `mpgx_run_gem()` | Efficient interaction testing for large cohorts |
| MR-MEGA (generic) | `mpgx_run_mrmega()` | Direct MR-MEGA execution wrapper |
| MR-MEGA for LR outputs | `mpgx_run_mrmega_for_lr()` | Prepare/run MR-MEGA from per-ancestry LR outputs |
| MR-MEGA for GEM outputs | `mpgx_run_mrmega_for_gem()` | Prepare/run MR-MEGA from per-ancestry GEM outputs |
| Fixed-effect meta-analysis | `mpgx_run_fe_meta_for_lr()` | Meta-analyze per-ancestry LR outputs |

## End-to-End Example with Toy Data

The full toy-data example below runs multiple wrappers and meta-analysis helpers.

<details>
<summary>Show full example script</summary>

```r
library(MultiAncestryPGxGWAS)

# ----- Configuration -----
configs <- list(
  plink1.9.path = "plink",
  plink2.path = "plink2",
  flashpca.path = "flashpca_x86-64",
  gcta.path = "gcta",
  lemma.path = "lemma_1_0_4",
  gem.path = "GEM_1.4.5_static",
  bgenix.path = "bgenix",
  mrmega.path = "MR-MEGA",
  mpirun.path = "mpirun"
)
# Optionally override executable locations when tools are not on PATH.
# configs[["gcta.path"]] <- "/abs/path/to/gcta"

mpgx_validate_tool_configs(configs)

out_dir <- file.path(tempdir(), "mpgx_wrapper_examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Packaged toy data -----
toy_prefix <- sub("\\.bed$", "", mpgx_example_file("genotype_data/toyData.bed"))
sparse_grm_prefix <- sub(
  "\\.grm\\.sp$",
  "",
  mpgx_example_file("genotype_artifacts/toyData.sp_grm.grm.sp")
)
bgen_file <- mpgx_example_file("genotype_artifacts/toyData.bgen")
gds_file <- mpgx_example_file("genotype_artifacts/toyData.gds")

# continuous outcome
pheno_cont <- mpgx_example_file("simulated_data/continuous/sim_phenotype_seed_100001.txt")
qcovar_cont <- mpgx_example_file("simulated_data/continuous/sim_qcovar.txt")
covar_cont <- mpgx_example_file("simulated_data/continuous/sim_covar.txt")
envir_cont <- mpgx_example_file("simulated_data/continuous/sim_envir.txt")

pheno_lemma_cont <- mpgx_example_file("simulated_data/continuous/sim_phenotype_forLEMMA_seed_100001.txt")
covar_lemma_cont <- mpgx_example_file("simulated_data/continuous/sim_covar_qcovar_forLEMMA.txt")
envir_lemma_cont <- mpgx_example_file("simulated_data/continuous/sim_envir_forLEMMA.txt")

# binary outcome
pheno_bin <- mpgx_example_file("simulated_data/binary/sim_phenotype_seed_100001.txt")
qcovar_bin <- mpgx_example_file("simulated_data/binary/sim_qcovar.txt")
covar_bin <- mpgx_example_file("simulated_data/binary/sim_covar.txt")
envir_bin <- mpgx_example_file("simulated_data/binary/sim_envir.txt")
pheno_gem_bin <- mpgx_example_file("simulated_data/binary/sim_phenotype_forGEM_seed_100001.txt")

# ----- Wrapper examples for this toy-data run -----
# Methods are shown with one set of example inputs; applicability can differ by method and model setup.
# 1) GCTA linear regression with robust (sandwich) variance for a continuous outcome.
lr <- mpgx_run_lr_sandwich(
  bfile_prefix = toy_prefix,
  phenotype_file = pheno_cont,
  qcovar_file = qcovar_cont,
  covar_file = covar_cont,
  envir_file = envir_cont,
  output_prefix = file.path(out_dir, "continuous.LRallSandwich"),
  configs = configs
)

# 2) fastGWA-GE for a continuous outcome.
fastgwa <- mpgx_run_fastgwa_ge(
  bfile_prefix = toy_prefix,
  grm_sparse_prefix = sparse_grm_prefix,
  phenotype_file = pheno_cont,
  qcovar_file = qcovar_cont,
  covar_file = covar_cont,
  envir_file = envir_cont,
  output_prefix = file.path(out_dir, "continuous.fastGWA_GE"),
  configs = configs
)

# 3) LEMMA for a continuous outcome.
lemma <- mpgx_run_lemma(
  bgen_file = bgen_file,
  phenotype_file = pheno_lemma_cont,
  envir_file = envir_lemma_cont,
  covar_file = covar_lemma_cont,
  output_prefix = file.path(out_dir, "continuous.LEMMA"),
  n_tasks = 1,
  configs = configs
)

# 4) GENESIS for a continuous outcome.
genesis <- mpgx_run_genesis(
  gds_file = gds_file,
  sparse_grm_prefix = sparse_grm_prefix,
  phenotype_file = pheno_cont,
  qcovar_file = qcovar_cont,
  covar_file = covar_cont,
  envir_file = envir_cont,
  output_file = file.path(out_dir, "continuous.GENESIS.txt")
)

# 5) GEM for a binary outcome
gem <- mpgx_run_gem(
  bfile_prefix = toy_prefix,
  phenotype_file = pheno_gem_bin,
  output_prefix = file.path(out_dir, "binary.GEM"),
  configs = configs
)

# 6) MR-MEGA wrapper for continuous per-ancestry LR-Sandwich outputs.
lr_by_eth_files <- c(
  AFR = mpgx_example_file("model_results/continuous/lr_sandwich/continuous.LRallSandwich.AFR.fastGWA"),
  EAS = mpgx_example_file("model_results/continuous/lr_sandwich/continuous.LRallSandwich.EAS.fastGWA"),
  EUR = mpgx_example_file("model_results/continuous/lr_sandwich/continuous.LRallSandwich.EUR.fastGWA"),
  SAS = mpgx_example_file("model_results/continuous/lr_sandwich/continuous.LRallSandwich.SAS.fastGWA")
)

mrmega_lr <- mpgx_run_mrmega_for_lr(
  ancestry_result_files = lr_by_eth_files,
  output_prefix_base = file.path(out_dir, "continuous.MRMEGA_LRallSandwich"),
  configs = configs
)

# 7) Fixed-effect meta-analysis for binary per-ancestry LR-Sandwich outputs.
ancestry_labels <- c("AFR", "EAS", "EUR", "SAS")
lr_bin_by_eth <- lapply(ancestry_labels, function(ances) {
  mpgx_run_lr_sandwich(
    bfile_prefix = toy_prefix,
    phenotype_file = pheno_bin,
    qcovar_file = qcovar_bin,
    covar_file = covar_bin,
    envir_file = envir_bin,
    keep_file = mpgx_example_file(paste0("model_results/binary/keep_", ances, ".txt")),
    output_prefix = file.path(out_dir, paste0("binary.LRallSandwich.", ances)),
    configs = configs
  )
})
names(lr_bin_by_eth) <- ancestry_labels

fe_meta_binary <- mpgx_run_fe_meta_for_lr(
  ancestry_result_files = extract_gwas_files(lr_bin_by_eth),
  ancestry_labels = ancestry_labels,
  output_file = file.path(out_dir, "binary.meta_LRallSandwich.txt"),
  attach_result = TRUE
)
```

</details>

## Script Wrappers

CLI-style wrappers are provided in the package R source:

- `R/cli_wrapper_run_simulation_stage.R` (`mpgx_cli_run_simulation_stage()`)
- `R/cli_wrapper_run_genotype_stage.R` (`mpgx_cli_run_genotype_stage()`)


## References

Wujuan Zhong, Judong Shen*. Review and simulation evaluation of statistical methods for discovering predictive biomarkers in multi-ancestry pharmacogenomics GWAS.
