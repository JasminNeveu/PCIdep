#' Compute selective hierarchical-clustering p-values
#'
#' @param X Data matrix.
#' @param Y Optional test data.
#' @param hcl Hierarchical clustering object.
#' @param kmax Maximum number of clusters to evaluate.
#' @param dismat Optional distance matrix.
#' @param parallel_config Parallel configuration.
#' @param U Optional nuisance matrix.
#' @param Sigma Optional covariance matrix.
#' @param UY Optional nuisance matrix for Y.
#' @param precUY Optional precision matrix.
#' @param linkage Linkage method.
#' @param ndraws Number of Monte Carlo draws.
#' @param sample_split Whether to sample split.
#' @param nY Number of observations in Y.
#' @param return_Sigma Return Sigma.
#' @param return_X_clus Return clustered X.
#' @param early_stop Enable branch-wise early stopping.
#' @param alpha Significance level.
#' @param show_progress Show progress.
#' @param correction Multiple-testing correction.
#' @param min_pts Minimum cluster size.
#' @param verbose Verbose output.
#'
#' @return Named numeric vector of p-values, indexed by k.
#' @keywords internal
compute.pvals <- function(
  X,
  Y,
  hcl,
  dismat,
  linkage,
  kmax = nrow(X) - 5L,
  parallel_config = NULL,
  U = NULL,
  Sigma = NULL,
  UY = NULL,
  precUY = NULL,
  ndraws = 2000L,
  early_stop = FALSE,
  alpha = 0.05,
  show_progress = FALSE,
  correction = NULL,
  min_pts = 1L,
  verbose = FALSE
) {
  if (is.null(hcl)) {
    stop("'hcl' object is required.")
  }

  if (!inherits(hcl, "hclust")) {
    stop("'hcl' must be an object inheriting from 'hclust'.")
  }

  n <- length(hcl$order)

  if (n < 2L) {
    stop("'hcl' must contain at least two observations.")
  }

  if (
    length(kmax) != 1L ||
      !is.numeric(kmax) ||
      is.na(kmax) ||
      kmax != as.integer(kmax)
  ) {
    stop("'kmax' must be a single integer.")
  }

  kmax <- as.integer(kmax)

  if (kmax < 2L || kmax > n - 1L) {
    stop("'kmax' must be between 2 and n - 1.")
  }

  online_methods <- c("ADDIS", "online_fallback")

  is_online_correction <- !is.null(correction) &&
    correction %in% online_methods

  if (
    early_stop &&
      !is.null(parallel_config) &&
      isTRUE(parallel_config$enabled)
  ) {
    stop(
      "Can't use parallelization and early stopping at the same time."
    )
  }

  if (
    early_stop &&
      !is.null(correction) &&
      !is_online_correction
  ) {
    stop(
      "With early_stop = TRUE, only an online correction (",
      paste(online_methods, collapse = ", "),
      ") can be applied."
    )
  }

  size_text <- if (min_pts > 1L) {
    paste0(", min cluster size = ", min_pts)
  } else {
    ""
  }

  corr_text <- if (is.null(correction)) {
    ""
  } else if (is_online_correction) {
    paste0(", ", correction, " correction (alpha = ", alpha, ")")
  } else {
    paste0(", ", correction, " correction")
  }

  if (!early_stop) {
    cat(
      sprintf(
        "Computing p-values for k = 2..%d%s%s.\n",
        kmax,
        corr_text,
        size_text
      )
    )
  } else {
    cat(
      sprintf(
        paste0(
          "Computing p-values with branch-wise early stopping ",
          "(alpha = %s, kmax = %d)%s%s.\n"
        ),
        alpha,
        kmax,
        corr_text,
        size_text
      )
    )
  }

  # Full vector for k = 2,...,n.
  pvals <- stats::setNames(
    rep(NA_real_, n - 1L),
    as.character(seq.int(2L, n))
  )

  k_seq <- seq.int(2L, kmax)

  node_sizes <- get.node.sizes(hcl)

  node.for.k <- function(k) {
    n - k + 1L
  }

  is_node_blocked <- rep(FALSE, n - 1L)

  if (!early_stop) {
    run_one_k <- function(k) {
      node <- node.for.k(k)

      c1 <- hcl$merge[node, 1]
      c2 <- hcl$merge[node, 2]

      size1 <- if (c1 < 0) 1L else node_sizes[c1]
      size2 <- if (c2 < 0) 1L else node_sizes[c2]

      if (size1 >= min_pts && size2 >= min_pts) {
        clusters <- get.merged.clusters(
          hcl,
          k
        )

        hc_test <- PCIdep::test.clusters.hc(
          X = X,
          U = U,
          Sigma = Sigma,
          Y = Y,
          UY = UY,
          precUY = precUY,
          linkage = linkage,
          NC = k,
          clusters = clusters,
          ndraws = ndraws,
          sample_split = FALSE,
          nY = nY,
          return_Sigma = return_Sigma,
          return_X_clus = return_X_clus,
          hcl = hcl,
          dismat = dismat
        )

        return(
          max(
            hc_test$pval,
            2.2e-16,
            na.rm = TRUE
          )
        )
      }

      NA_real_
    }

    if (
      !is.null(parallel_config) &&
        isTRUE(parallel_config$enabled)
    ) {
      workers <- if (!is.null(parallel_config$workers)) {
        parallel_config$workers
      } else {
        max(1L, parallel::detectCores() - 1L)
      }

      future::plan(
        future::multisession,
        workers = workers
      )

      on.exit(
        future::plan(future::sequential),
        add = TRUE
      )

      pvals_raw <- progressr::with_progress(
        {
          p <- progressr::progressor(along = k_seq)

          unlist(
            furrr::future_map(
              k_seq,
              function(k) {
                result <- run_one_k(k)
                p()
                result
              }
            )
          )
        },
        enable = show_progress
      )
    } else {
      pvals_raw <- sapply(
        k_seq,
        run_one_k
      )
    }

    pvals_raw <- as.numeric(pvals_raw)

    names(pvals_raw) <- as.character(k_seq)

    if (!is.null(correction)) {
      pvals_corrected <- correction.multiplicity(
        pvals_raw,
        correction = correction,
        alpha = alpha,
        verbose = FALSE
      )
    } else {
      pvals_corrected <- pvals_raw
    }

    pvals[names(pvals_corrected)] <- pvals_corrected

    # Propagate blocked nodes top-down.
    for (i in seq_along(k_seq)) {
      k <- k_seq[i]
      k_str <- as.character(k)
      node <- node.for.k(k)

      c1 <- hcl$merge[node, 1]
      c2 <- hcl$merge[node, 2]

      size1 <- if (c1 < 0) 1L else node_sizes[c1]
      size2 <- if (c2 < 0) 1L else node_sizes[c2]

      if (is_node_blocked[node]) {
        if (verbose) {
          cat(
            sprintf(
              "k=%d: Ignored (branch already blocked by ancestor)\n",
              k
            )
          )
        }

        pvals[k_str] <- NA_real_

        if (c1 > 0) {
          is_node_blocked[c1] <- TRUE
        }

        if (c2 > 0) {
          is_node_blocked[c2] <- TRUE
        }

        next
      }

      if (size1 >= min_pts && size2 >= min_pts) {
        pval_adj <- pvals_corrected[k_str]

        pvals[k_str] <- pval_adj

        if (
          !is.na(pval_adj) &&
            pval_adj >= alpha
        ) {
          if (c1 > 0) {
            is_node_blocked[c1] <- TRUE
          }

          if (c2 > 0) {
            is_node_blocked[c2] <- TRUE
          }
        }
      } else if (
        size1 < min_pts &&
          size2 >= min_pts
      ) {
        pvals[k_str] <- NA_real_

        if (c1 > 0) {
          is_node_blocked[c1] <- TRUE
        }
      } else if (
        size1 >= min_pts &&
          size2 < min_pts
      ) {
        pvals[k_str] <- NA_real_

        if (c2 > 0) {
          is_node_blocked[c2] <- TRUE
        }
      } else {
        pvals[k_str] <- NA_real_

        if (c1 > 0) {
          is_node_blocked[c1] <- TRUE
        }

        if (c2 > 0) {
          is_node_blocked[c2] <- TRUE
        }
      }
    }

    return(pvals)
  }

  # ---------------------------------------------------------
  # Early stopping
  # ---------------------------------------------------------

  raw_history <- numeric(0)
  active_count <- 1L

  for (i in seq_along(k_seq)) {
    k <- k_seq[i]
    k_str <- as.character(k)
    node <- node.for.k(k)

    c1 <- hcl$merge[node, 1]
    c2 <- hcl$merge[node, 2]

    if (is_node_blocked[node]) {
      pvals[k_str] <- NA_real_

      if (c1 > 0) {
        is_node_blocked[c1] <- TRUE
      }

      if (c2 > 0) {
        is_node_blocked[c2] <- TRUE
      }

      next
    }

    size1 <- if (c1 < 0) 1L else node_sizes[c1]
    size2 <- if (c2 < 0) 1L else node_sizes[c2]

    if (
      size1 >= min_pts &&
        size2 >= min_pts
    ) {
      clusters <- get.merged.clusters(
        hcl,
        k
      )

      hc_test <- PCIdep::test.clusters.hc(
        X = X,
        U = U,
        Sigma = Sigma,
        Y = Y,
        UY = UY,
        precUY = precUY,
        linkage = linkage,
        NC = k,
        clusters = clusters,
        ndraws = ndraws,
        sample_split = FALSE,
        nY = nY,
        return_Sigma = return_Sigma,
        return_X_clus = return_X_clus,
        hcl = hcl,
        dismat = dismat
      )

      raw_pval <- max(
        hc_test$pval,
        2.2e-16,
        na.rm = TRUE
      )

      raw_history <- c(
        raw_history,
        raw_pval
      )

      if (!is.null(correction)) {
        corrected_history <- correction.multiplicity(
          raw_history,
          correction = correction,
          alpha = alpha,
          verbose = FALSE
        )

        decision_pval <- corrected_history[
          length(corrected_history)
        ]
      } else {
        decision_pval <- raw_pval
      }

      pvals[k_str] <- decision_pval

      if (decision_pval < alpha) {
        active_count <- active_count + 1L

        if (verbose) {
          cat(
            sprintf(
              paste0(
                "k=%d: Split VALIDATED ",
                "(p=%.4g < %g) -> continue\n"
              ),
              k,
              decision_pval,
              alpha
            )
          )
        }
      } else {
        if (c1 > 0) {
          is_node_blocked[c1] <- TRUE
        }

        if (c2 > 0) {
          is_node_blocked[c2] <- TRUE
        }

        active_count <- active_count - 1L

        if (verbose) {
          cat(
            sprintf(
              "k=%d: Early stop (p=%.4f)\n",
              k,
              decision_pval
            )
          )
        }
      }
    } else if (
      size1 < min_pts &&
        size2 >= min_pts
    ) {
      pvals[k_str] <- NA_real_

      if (c1 > 0) {
        is_node_blocked[c1] <- TRUE
      }
    } else if (
      size1 >= min_pts &&
        size2 < min_pts
    ) {
      pvals[k_str] <- NA_real_

      if (c2 > 0) {
        is_node_blocked[c2] <- TRUE
      }
    } else {
      pvals[k_str] <- NA_real_

      if (c1 > 0) {
        is_node_blocked[c1] <- TRUE
      }

      if (c2 > 0) {
        is_node_blocked[c2] <- TRUE
      }

      active_count <- active_count - 1L
    }

    if (active_count <= 0L) {
      if (verbose) {
        cat(
          sprintf(
            paste0(
              "All branches blocked at k=%d. ",
              "Stopping before kmax=%d.\n"
            ),
            k,
            kmax
          )
        )
      }

      break
    }
  }

  pvals
}
