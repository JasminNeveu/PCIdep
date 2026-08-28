#' Compute a selective p-dendrogram
#'
#' @param X Data matrix used for selective inference.
#' @param Y Optional test data for covariance estimation. If NULL, sample splitting is performed on X.
#' @param hcl Optional hclust object (requires Y to be provided).
#' @param dismat Optional distance matrix (requires Y to be provided).
#' @param split_prop Proportion of samples from X used for clustering when Y is NULL (default: 0.8).
#' @param min_pts Minimum size for a valid cluster.
#' @param alpha Significance level.
#' @param correction Multiple-testing correction.
#' @param linkage Method for hierarchical clustering when hcl is NULL (default: "ward.D").
#' @param ... Additional arguments passed to [compute.pvals()].
#'
#' @return A p.hclust object.
#' @export
p.dendrogram <- function(
  X,
  Y = NULL,
  hcl = NULL,
  dismat = NULL,
  split_prop = 0.8,
  min_pts = 1L,
  alpha = 0.05,
  correction = NULL,
  linkage = "ward.D",
  ...
) {
  # --- 1. Validation des arguments d'entrée de base ---
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("'X' must be a matrix or data.frame.")
  }

  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop("'alpha' must be a single value in (0, 1).")
  }

  if (
    !is.numeric(min_pts) ||
      length(min_pts) != 1L ||
      is.na(min_pts) ||
      min_pts < 1 ||
      min_pts != as.integer(min_pts)
  ) {
    stop("'min_pts' must be a positive integer.")
  }

  min_pts <- as.integer(min_pts)

  # --- 2. Détection des erreurs de combinaisons incompatibles ---
  if (is.null(Y) && (!is.null(hcl) || !is.null(dismat))) {
    stop(
      "Incompatible arguments: When 'Y' is NULL, sample splitting is automatically performed on 'X'.\n",
      "You cannot pass pre-computed 'hcl' or 'dismat' objects without providing 'Y'."
    )
  }

  # --- 3. Gestion du Sample Splitting vs Y fourni ---
  if (is.null(Y)) {
    if (
      !is.numeric(split_prop) ||
        length(split_prop) != 1L ||
        is.na(split_prop) ||
        split_prop <= 0 ||
        split_prop >= 1
    ) {
      stop(
        "'split_prop' must be a single numeric value strictly between 0 and 1."
      )
    }

    n <- nrow(X)
    n_cluster <- round(n * split_prop)
    idx_cluster <- sort(sample.int(n, size = n_cluster))
    idx_cov <- setdiff(seq_len(n), idx_cluster)

    message(sprintf(
      "Info: 'Y' was not provided. Sample splitting applied: %d observations for clustering (%g%%) and %d for covariance estimation (%g%%).",
      length(idx_cluster),
      round(split_prop * 100, 1),
      length(idx_cov),
      round((1 - split_prop) * 100, 1)
    ))

    X_work <- X[idx_cluster, , drop = FALSE]
    Y_work <- X[idx_cov, , drop = FALSE]

    split_indices <- list(
      idx_clustering = idx_cluster,
      idx_covariance = idx_cov
    )
  } else {
    if (!is.matrix(Y) && !is.data.frame(Y)) {
      stop("'Y' must be a matrix or data.frame.")
    }
    X_work <- X
    Y_work <- Y
    split_indices <- NULL
  }

  # --- 4. Construction ou vérification de dismat ---
  if (is.null(dismat)) {
    dismat <- stats::dist(X_work, method = "euclidean")^2
  } else {
    if (!inherits(dismat, "dist") && !is.matrix(dismat)) {
      stop("'dismat' must be a 'dist' object or a matrix.")
    }
  }

  # --- 5. Construction ou vérification de hcl ---
  if (is.null(hcl)) {
    hcl <- fastcluster::hclust(dismat, method = linkage)
  } else {
    if (!inherits(hcl, "hclust")) {
      stop("'hcl' must be an object inheriting from 'hclust'.")
    }
    if (length(hcl$order) != nrow(X_work)) {
      stop(sprintf(
        "Dimension mismatch: 'hcl' contains %d elements, but 'X_work' has %d observations.",
        length(hcl$order),
        nrow(X_work)
      ))
    }
  }

  # --- 6. Calcul des p-values ---
  pvalues <- compute.pvals(
    X = X_work,
    Y = Y_work,
    hcl = hcl,
    dismat = dismat,
    min_pts = min_pts,
    alpha = alpha,
    correction = correction,
    ...
  )

  # --- 7. Assemblage de l'objet p.hclust final ---
  new_p_hclust(
    hcl = hcl,
    pvalues = pvalues,
    alpha = alpha,
    min_pts = min_pts,
    correction = correction,
    split_indices = split_indices
  )
}
