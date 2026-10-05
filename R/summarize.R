#' Summarize association results
#'
#' Computes genomic inflation, FPR, and power metrics consistent with the
#' original `summarize_results.R` logic.
#'
#' @param result Data frame with association results.
#' @param causal_snps Character vector of causal SNP IDs.
#' @param significance_level Significance level for FPR.
#' @param genomewide_level Genome-wide threshold for power from p-values.
#' @param rm_na_2df If `TRUE`, rows with missing `P_2df` are dropped before
#'   metric calculation.
#' @param mrmega_pc_num Optional MR-MEGA PC count for inflation denominator.
#'
#' @return One-row data frame with summary metrics.
#' @export
mpgx_summarize_association <- function(
  result,
  causal_snps,
  significance_level = 0.05,
  genomewide_level = 5e-8,
  rm_na_2df = TRUE,
  mrmega_pc_num = NULL
) {
  required_cols <- c(
    "CHR", "SNP",
    "chisq_G", "chisq_G_by_E", "chisq_2df",
    "P_G", "P_G_by_E", "P_2df"
  )

  res <- as.data.frame(result, stringsAsFactors = FALSE)
  for (nm in required_cols) {
    if (!(nm %in% colnames(res))) {
      res[[nm]] <- NA_real_
    }
  }

  res$CHR <- as.numeric(res$CHR)
  res$SNP <- as.character(res$SNP)

  if (isTRUE(rm_na_2df)) {
    res <- res[!is.na(res$P_2df), , drop = FALSE]
  }

  if (nrow(res) == 0) {
    empty <- data.frame(
      genomic_inflation_factor_G = NA_real_,
      genomic_inflation_factor_G_by_E = NA_real_,
      genomic_inflation_factor_2df = NA_real_,
      FPR_G = NA_real_,
      FPR_G_by_E = NA_real_,
      FPR_2df = NA_real_,
      power_G = NA_real_,
      power_G_by_E = NA_real_,
      power_2df = NA_real_,
      power_G_pval = NA_real_,
      power_G_by_E_pval = NA_real_,
      power_2df_pval = NA_real_,
      stringsAsFactors = FALSE
    )
    return(empty)
  }

  even <- res[res$CHR %% 2 == 0 & res$CHR <= 22, , drop = FALSE]
  if (nrow(even) == 0) {
    even <- res
  }

  median1 <- stats::qchisq(0.5, 1)
  median2 <- stats::qchisq(0.5, 2)

  if (is.null(mrmega_pc_num)) {
    gif_g <- stats::median(even$chisq_G, na.rm = TRUE) / median1
    gif_ge <- stats::median(even$chisq_G_by_E, na.rm = TRUE) / median1
    gif_2df <- stats::median(even$chisq_2df, na.rm = TRUE) / median2
  } else {
    denom <- stats::qchisq(0.5, mrmega_pc_num + 1)
    gif_g <- stats::median(even$chisq_G, na.rm = TRUE) / denom
    gif_ge <- stats::median(even$chisq_G_by_E, na.rm = TRUE) / denom
    gif_2df <- NA_real_
  }

  prop_sig <- function(x, cutoff) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    mean(x < cutoff)
  }

  fpr_g <- prop_sig(even$P_G, significance_level)
  fpr_ge <- prop_sig(even$P_G_by_E, significance_level)
  fpr_2df <- prop_sig(even$P_2df, significance_level)

  causal <- res[res$SNP %in% causal_snps, , drop = FALSE]

  mean_scaled <- function(x, inflation) {
    if (nrow(causal) == 0 || !is.finite(inflation) || inflation == 0) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE) / inflation
  }

  power_g <- mean_scaled(causal$chisq_G, gif_g)
  power_ge <- mean_scaled(causal$chisq_G_by_E, gif_ge)
  power_2df <- mean_scaled(causal$chisq_2df, gif_2df)

  power_g_p <- prop_sig(causal$P_G, genomewide_level)
  power_ge_p <- prop_sig(causal$P_G_by_E, genomewide_level)
  power_2df_p <- prop_sig(causal$P_2df, genomewide_level)

  data.frame(
    genomic_inflation_factor_G = gif_g,
    genomic_inflation_factor_G_by_E = gif_ge,
    genomic_inflation_factor_2df = gif_2df,
    FPR_G = fpr_g,
    FPR_G_by_E = fpr_ge,
    FPR_2df = fpr_2df,
    power_G = power_g,
    power_G_by_E = power_ge,
    power_2df = power_2df,
    power_G_pval = power_g_p,
    power_G_by_E_pval = power_ge_p,
    power_2df_pval = power_2df_p,
    stringsAsFactors = FALSE
  )
}

