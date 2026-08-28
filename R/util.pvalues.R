correction.multiplicity <- function(
  pvals,
  correction = c(
    "ADDIS",
    "online_fallback",
    "holm",
    "hochberg",
    "hommel",
    "bonferroni",
    "BH",
    "BY"
  ),
  alpha = 0.05,
  verbose = TRUE
) {
  correction <- match.arg(correction)

  if (!is.numeric(pvals)) {
    stop("'pvals' must be numeric.")
  }

  if (
    length(alpha) != 1L ||
      !is.numeric(alpha) ||
      is.na(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop("'alpha' must be in (0, 1).")
  }

  if (verbose) {
    message(
      paste0(
        "Correct multiple hypothesis testing bias using ",
        correction,
        " with a threshold of ",
        alpha
      )
    )
  }

  non_na_idx <- which(!is.na(pvals))

  if (length(non_na_idx) == 0L) {
    return(pvals)
  }

  pvals_non_na <- pvals[non_na_idx]

  if (correction == "ADDIS") {
    pvals_corrected <- onlineFDR::ADDIS_spending(
      pvals_non_na,
      alpha = alpha
    ) |>
      dplyr::mutate(
        p_adjusted = pmin(
          1,
          alpha / alphai * pval
        )
      ) |>
      dplyr::pull(
        p_adjusted
      )
  } else if (correction == "online_fallback") {
    pvals_corrected <- onlineFDR::online_fallback(
      pvals_non_na,
      alpha = alpha
    ) |>
      dplyr::mutate(
        p_adjusted = pmin(
          1,
          alpha / alphai * pval
        )
      ) |>
      dplyr::pull(
        p_adjusted
      )
  } else {
    pvals_corrected <- stats::p.adjust(
      pvals_non_na,
      method = correction
    )
  }

  result <- pvals

  result[non_na_idx] <- pvals_corrected

  result
}
