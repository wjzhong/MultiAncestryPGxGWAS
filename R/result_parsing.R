#' Parse and normalize method-specific GWAS result files
#'
#' Standardizes result files from methods used in `scripts/summarize_results.R`
#' into a common schema with columns `CHR`, `SNP`, `chisq_G`, `chisq_G_by_E`,
#' `chisq_2df`, `P_G`, `P_G_by_E`, and `P_2df` (plus method-specific optional
#' beta/se columns when available).
#'
#' @param method Method name from the original workflow.
#' @param files Input files required by `method`.
#'   For single-file methods (`GENESIS`, `LEMMA`, `GEM`, and default fastGWA-like
#'   methods), pass a length-1 character path.
#'   For fixed-effect and random-effect meta methods, pass one file per ancestry
#'   (vector or named list).
#'   For `MRMEGA_LRallSandwich` / `MRMEGA_GEM`, pass either a named list with
#'   entries `test_G` and `test_G_by_E`, or a character vector of length 2 in
#'   that order.
#' @param ancestry_labels Optional ancestry labels used by meta methods.
#' @param output_file Optional output file path to write a tab-delimited table.
#'
#' @return Parsed result data frame in a common schema.
#' @export
mpgx_parse_method_result_df <- function(
  method,
  files,
  ancestry_labels = NULL,
  output_file = NULL
) {
  mpgx_assert_scalar_character(method, "method")

  normalize_out <- function(df) {
    out <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!("CHR" %in% colnames(out))) {
      stop("Parsed result does not contain CHR.", call. = FALSE)
    }
    if (!("SNP" %in% colnames(out))) {
      stop("Parsed result does not contain SNP.", call. = FALSE)
    }

    required_stat_cols <- c("chisq_G", "chisq_G_by_E", "chisq_2df", "P_G", "P_G_by_E", "P_2df")
    for (nm in required_stat_cols) {
      if (!(nm %in% colnames(out))) {
        out[[nm]] <- NA_real_
      }
    }

    out$CHR <- suppressWarnings(as.numeric(out$CHR))
    out$SNP <- as.character(out$SNP)
    for (nm in setdiff(colnames(out), c("SNP"))) {
      if (is.character(out[[nm]]) || is.factor(out[[nm]])) {
        next
      }
      out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
    }

    out
  }

  read_single <- function(x, label) {
    mpgx_read_delim_table(x, has_header = TRUE, arg_name = label)
  }

  meta_methods <- c(
    "meta_StandardLRall", "meta_LRallSandwich",
    "meta_StandardLRunrelated", "meta_LRunrelatedSandwich"
  )

  if (method %in% meta_methods) {
    out <- mpgx_fixed_effect_meta_lr(
      ancestry_result_files = files,
      ancestry_labels = ancestry_labels,
      output_file = NULL
    )
    out <- normalize_out(out)
  } else if (identical(method, "REmeta_LRallSandwich")) {
    if (!requireNamespace("metafor", quietly = TRUE)) {
      stop(
        "Method 'REmeta_LRallSandwich' requires package 'metafor'. Install it to parse this method.",
        call. = FALSE
      )
    }

    file_list <- as.character(unname(unlist(files, use.names = FALSE)))
    if (!length(file_list)) {
      stop("files must contain per-ancestry LR result files for REmeta_LRallSandwich.", call. = FALSE)
    }
    if (is.null(ancestry_labels)) {
      name_guess <- names(files)
      if (!is.null(name_guess) && length(name_guess) == length(file_list) && all(nzchar(name_guess))) {
        ancestry_labels <- name_guess
      } else {
        ancestry_labels <- paste0("cohort", seq_along(file_list))
      }
    }

    required_cols <- c("CHR", "SNP", "BETA_G", "BETA_G_by_E", "SE_G", "SE_G_by_E")
    result_list <- lapply(seq_along(file_list), function(i) {
      dat <- read_single(file_list[[i]], "files")
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
        paste0(ancestry_labels[[i]], "_", c("BETA_G", "BETA_G_by_E", "SE_G", "SE_G_by_E"))
      )
      dat
    })
    result <- Reduce(function(x, y) dplyr::full_join(x, y, by = c("CHR", "SNP")), result_list)
    eth_list <- as.character(ancestry_labels)

    reMetaAnalysis <- function(vec_BETA, vec_SE) {
      ans <- try(metafor::rma(yi = vec_BETA, sei = vec_SE, ni = NA, method = "REML"), silent = TRUE)
      if (inherits(ans, "try-error")) {
        beta <- NA_real_
        se <- NA_real_
        chisq <- NA_real_
        pval <- NA_real_
      } else {
        beta <- as.numeric(ans$b)
        se <- as.numeric(ans$se)
        chisq <- as.numeric(ans$zval)^2
        pval <- as.numeric(ans$pval)
      }
      data.frame(beta = beta, se = se, chisq = chisq, pval = pval, stringsAsFactors = FALSE)
    }

    REmeta_analysis <- t(apply(dplyr::select(result, -CHR, -SNP), 1, function(x) {
      temp_BETA_G <- temp_SE_G <- temp_BETA_G_by_E <- temp_SE_G_by_E <- NULL

      for (ancestry in eth_list) {
        if (!is.na(x[paste(ancestry, "BETA_G", sep = "_")])) {
          temp_BETA_G[ancestry] <- as.numeric(x[paste(ancestry, c("BETA_G"), sep = "_")])
          temp_SE_G[ancestry] <- as.numeric(x[paste(ancestry, c("SE_G"), sep = "_")])
          temp_BETA_G_by_E[ancestry] <- as.numeric(x[paste(ancestry, c("BETA_G_by_E"), sep = "_")])
          temp_SE_G_by_E[ancestry] <- as.numeric(x[paste(ancestry, c("SE_G_by_E"), sep = "_")])
        }
      }

      output_G <- reMetaAnalysis(temp_BETA_G, temp_SE_G)
      output_G_by_E <- reMetaAnalysis(temp_BETA_G_by_E, temp_SE_G_by_E)

      c(
        "BETA_G" = as.numeric(output_G[["beta"]]),
        "SE_G" = as.numeric(output_G[["se"]]),
        "BETA_G_by_E" = as.numeric(output_G_by_E[["beta"]]),
        "SE_G_by_E" = as.numeric(output_G_by_E[["se"]]),
        "chisq_G" = as.numeric(output_G[["chisq"]]),
        "chisq_G_by_E" = as.numeric(output_G_by_E[["chisq"]]),
        "chisq_2df" = NA_real_,
        "P_G" = as.numeric(output_G[["pval"]]),
        "P_G_by_E" = as.numeric(output_G_by_E[["pval"]]),
        "P_2df" = NA_real_
      )
    }))

    out <- cbind(result[, c("CHR", "SNP"), drop = FALSE], as.data.frame(REmeta_analysis, stringsAsFactors = FALSE))
    out <- normalize_out(out)
  } else if (method %in% c("MRMEGA_LRallSandwich", "MRMEGA_GEM")) {
    if (is.list(files) && !is.null(names(files))) {
      g_file <- files[["test_G"]]
      gxe_file <- files[["test_G_by_E"]]
    } else {
      file_list <- as.character(unname(unlist(files, use.names = FALSE)))
      if (length(file_list) != 2L) {
        stop("For MRMEGA methods, files must provide test_G and test_G_by_E result files.", call. = FALSE)
      }
      g_file <- file_list[[1]]
      gxe_file <- file_list[[2]]
    }

    result_G <- read_single(g_file, "files[['test_G']]")
    result_G_by_E <- read_single(gxe_file, "files[['test_G_by_E']]")

    out <- dplyr::full_join(
      dplyr::rename(
        dplyr::select(result_G, MarkerName, Chromosome, Position, `P-value_association`, chisq_association),
        P_G = `P-value_association`,
        chisq_G = chisq_association
      ),
      dplyr::rename(
        dplyr::select(result_G_by_E, MarkerName, Chromosome, Position, `P-value_association`, chisq_association),
        P_G_by_E = `P-value_association`,
        chisq_G_by_E = chisq_association
      ),
      by = c("MarkerName", "Chromosome", "Position")
    )

    out <- dplyr::rename(out, SNP = MarkerName, CHR = Chromosome)
    out$chisq_2df <- NA_real_
    out$P_2df <- NA_real_
    out <- normalize_out(out)
  } else if (identical(method, "GENESIS")) {
    result <- read_single(as.character(unname(unlist(files, use.names = FALSE)))[[1]], "files")
    out <- result
    out$chisq_G <- (out$Est.G / out$SE.G)^2
    out$P_G <- 1 - stats::pchisq(out$chisq_G, df = 1)
    out$chisq_G_by_E <- out$GxE.Stat^2
    out$P_G_by_E <- out$GxE.pval
    out$chisq_2df <- out$Joint.Stat^2
    out$P_2df <- out$Joint.pval
    out$CHR <- as.numeric(out$chr)
    out$SNP <- as.character(out$variant.id)
    out <- normalize_out(out)
  } else if (identical(method, "LEMMA")) {
    result <- read_single(as.character(unname(unlist(files, use.names = FALSE)))[[1]], "files")
    colnames(result) <- c(colnames(result)[-1], "NA")

    out <- result
    out$chisq_G <- out$tStat_main^2
    out$P_G <- 10^(-out$neglogP_main)
    out$chisq_G_by_E <- out$chiSqStat_gxe_robust
    out$P_G_by_E <- 10^(-out$neglogP_gxe_robust)
    out$chisq_2df <- NA_real_
    out$P_2df <- NA_real_
    out$CHR <- as.numeric(out$chr)
    out$SNP <- as.character(out$rsid)
    out <- normalize_out(out)
  } else if (identical(method, "GEM")) {
    result <- read_single(as.character(unname(unlist(files, use.names = FALSE)))[[1]], "files")

    out <- result
    out$SNP <- as.character(out$SNPID)
    out$chisq_G <- (out$Beta_Marginal / out$robust_SE_Beta_Marginal)^2
    out$P_G <- out$robust_P_Value_Marginal
    out$chisq_G_by_E <- (out$`Beta_G-envir` / out$`robust_SE_Beta_G-envir`)^2
    out$P_G_by_E <- out$robust_P_Value_Interaction
    out$chisq_2df <- stats::qchisq(p = out$robust_P_Value_Joint, df = 2, lower.tail = FALSE)
    out$P_2df <- out$robust_P_Value_Joint
    out <- normalize_out(out)
  } else {
    result <- read_single(as.character(unname(unlist(files, use.names = FALSE)))[[1]], "files")
    out <- normalize_out(result)
  }

  if (!is.null(output_file)) {
    mpgx_ensure_dir(dirname(output_file))
    utils::write.table(out, output_file, quote = FALSE, row.names = FALSE, sep = "\t")
  }

  out
}
