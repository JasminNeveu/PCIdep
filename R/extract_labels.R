#' Extract the selective clustering partition
#'
#' @param x Object from which labels should be extracted.
#' @param ... Additional arguments passed to the class method.
#'
#' @return A vector of cluster labels.
#' @export
extract_labels <- function(x, ...) {
  UseMethod("extract_labels")
}


#' Extract labels from a p.hclust object
#'
#' @param x A p.hclust object.
#' @param min_pts Optional override for the minimum cluster size.
#' @param alpha Optional override for the significance level.
#' @param ... Unused arguments.
#'
#' @return Integer vector of cluster labels.
#' @export
extract_labels.p.hclust <- function(
  x,
  min_pts = NULL,
  alpha = NULL,
  ...
) {
  validate_p_hclust(x)

  config <- p_hclust_config(x)

  if (is.null(min_pts)) {
    min_pts <- config$min_pts
  }

  if (is.null(alpha)) {
    alpha <- config$alpha
  }

  selective.cutree(
    hcl = x,
    pvals = x$pvalues,
    min_pts = min_pts,
    alpha = alpha
  )
}


# Compatibility alias with singular naming.
#' @rdname extract_labels
#' @export
extract.label <- function(x, ...) {
  extract_labels(x, ...)
}


selective.cutree <- function(
  hcl,
  pvals,
  min_pts = 10L,
  alpha = 0.05
) {
  if (!inherits(hcl, "hclust")) {
    stop("'hcl' must inherit from 'hclust'.")
  }

  n_obs <- length(hcl$order)

  expected_names <- as.character(
    seq.int(2L, n_obs)
  )

  if (
    !is.numeric(pvals) ||
      length(pvals) != n_obs - 1L
  ) {
    stop(
      "'pvals' must contain one value for each k = 2, ..., n."
    )
  }

  # Accept unnamed vectors for internal compatibility.
  if (is.null(names(pvals))) {
    names(pvals) <- expected_names
  }

  node_sizes <- get.node.sizes(hcl)

  labels <- rep(
    0L,
    n_obs
  )

  cluster_id <- 1L

  visit <- function(row, peeling_root) {
    if (row < 0) {
      return(invisible(NULL))
    }

    c1 <- hcl$merge[row, 1]
    c2 <- hcl$merge[row, 2]

    s1 <- if (c1 < 0) 1L else node_sizes[c1]
    s2 <- if (c2 < 0) 1L else node_sizes[c2]

    if (
      s1 < min_pts &&
        s2 < min_pts
    ) {
      leafs <- get.node.leafs(
        hcl,
        peeling_root,
        node_sizes
      )

      labels[leafs] <<- cluster_id

      cluster_id <<- cluster_id + 1L

      return(invisible(NULL))
    }

    if (
      s1 < min_pts &&
        s2 >= min_pts
    ) {
      if (c2 > 0) {
        visit(
          c2,
          peeling_root
        )
      }

      return(invisible(NULL))
    }

    if (
      s1 >= min_pts &&
        s2 < min_pts
    ) {
      if (c1 > 0) {
        visit(
          c1,
          peeling_root
        )
      }

      return(invisible(NULL))
    }

    if (
      s1 >= min_pts &&
        s2 >= min_pts
    ) {
      k <- n_obs - row + 1L

      p_val <- pvals[
        as.character(k)
      ]

      if (
        !is.na(p_val) &&
          p_val < alpha
      ) {
        if (c1 > 0) {
          visit(
            c1,
            peeling_root = c1
          )
        }

        if (c2 > 0) {
          visit(
            c2,
            peeling_root = c2
          )
        }
      } else {
        leafs <- get.node.leafs(
          hcl,
          row,
          node_sizes
        )

        labels[leafs] <<- cluster_id

        cluster_id <<- cluster_id + 1L
      }

      return(invisible(NULL))
    }

    invisible(NULL)
  }

  root_node <- nrow(
    hcl$merge
  )

  visit(
    root_node,
    peeling_root = root_node
  )

  labels
}
