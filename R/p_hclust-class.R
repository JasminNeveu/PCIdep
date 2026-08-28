#' Constructor for p.hclust objects
#'
#' @param hcl A valid object inheriting from "hclust".
#' @param pvalues Numeric vector of p-values indexed by number of clusters.
#' @param alpha Significance level used to construct the object.
#' @param min_pts Minimum cluster size used during selective testing.
#' @param correction Multiple-testing correction used.
#' @param split_indices Optional list containing sample splitting indices.
#'
#' @return An object of class "p.hclust", "hclust".
#' @keywords internal
new_p_hclust <- function(
  hcl,
  pvalues,
  alpha = 0.05,
  min_pts = 1L,
  correction = NULL,
  split_indices = NULL
) {
  if (!inherits(hcl, "hclust")) {
    stop("'hcl' must be an object inheriting from 'hclust'.")
  }

  if (!is.numeric(pvalues)) {
    stop("'pvalues' must be a numeric vector.")
  }

  n <- length(hcl$order)
  expected_length <- n - 1L

  if (length(pvalues) != expected_length) {
    stop(
      "'pvalues' must have length n - 1, where n is the number ",
      "of observations in 'hcl'."
    )
  }

  if (any(!is.na(pvalues) & (pvalues < 0 | pvalues > 1))) {
    stop("'pvalues' must lie in [0, 1] or be NA.")
  }

  if (
    length(alpha) != 1L ||
      !is.numeric(alpha) ||
      is.na(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop("'alpha' must be a single numeric value strictly between 0 and 1.")
  }

  if (
    length(min_pts) != 1L ||
      !is.numeric(min_pts) ||
      is.na(min_pts) ||
      min_pts < 1 ||
      min_pts != as.integer(min_pts)
  ) {
    stop("'min_pts' must be a positive integer.")
  }

  if (
    !is.null(correction) &&
      (!is.character(correction) || length(correction) != 1L)
  ) {
    stop("'correction' must be NULL or a single character string.")
  }

  # Copy the hclust object.
  x <- hcl

  # Canonical p-value indexing: k = 2, ..., n
  names(pvalues) <- as.character(seq.int(2L, n))

  x$pvalues <- pvalues

  # Champ optionnel conservant les indices issus du sample splitting
  if (!is.null(split_indices)) {
    x$split_indices <- split_indices
  }

  # Keep analysis parameters outside the public list structure.
  attr(
    x,
    "p.hclust_config"
  ) <- list(
    alpha = alpha,
    min_pts = as.integer(min_pts),
    correction = correction
  )

  class(x) <- c("p.hclust", "hclust")

  validate_p_hclust(x)

  x
}


#' Validate a p.hclust object
#'
#' @param x Object to validate.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
validate_p_hclust <- function(x) {
  if (!inherits(x, "p.hclust")) {
    stop("Object must inherit from class 'p.hclust'.")
  }

  if (!inherits(x, "hclust")) {
    stop("A 'p.hclust' object must also inherit from 'hclust'.")
  }

  required_hclust_fields <- c(
    "merge",
    "height",
    "order",
    "labels",
    "method",
    "call",
    "dist.method"
  )

  missing_fields <- setdiff(
    required_hclust_fields,
    names(x)
  )

  if (length(missing_fields) > 0L) {
    stop(
      "Invalid 'p.hclust' object. Missing hclust fields: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (!"pvalues" %in% names(x)) {
    stop("Invalid 'p.hclust' object: missing 'pvalues'.")
  }

  n <- length(x$order)

  if (
    !is.matrix(x$merge) ||
      nrow(x$merge) != n - 1L ||
      ncol(x$merge) != 2L
  ) {
    stop("Invalid 'merge' component.")
  }

  if (length(x$height) != n - 1L) {
    stop("Invalid 'height' component.")
  }

  if (length(x$pvalues) != n - 1L) {
    stop(
      "'pvalues' must contain exactly one value for each possible ",
      "number of clusters k = 2, ..., n."
    )
  }

  if (
    any(
      !is.na(x$pvalues) &
        (x$pvalues < 0 | x$pvalues > 1)
    )
  ) {
    stop("'pvalues' must lie in [0, 1] or be NA.")
  }

  expected_names <- as.character(seq.int(2L, n))

  if (!identical(names(x$pvalues), expected_names)) {
    stop(
      "'pvalues' must be named with the number of clusters ",
      "from 2 to n."
    )
  }

  # Validation optionnelle de split_indices si le champ existe
  if ("split_indices" %in% names(x) && !is.null(x$split_indices)) {
    if (
      !is.list(x$split_indices) ||
        !all(c("idx_clustering", "idx_covariance") %in% names(x$split_indices))
    ) {
      stop(
        "Invalid 'split_indices' component: must be a list containing 'idx_clustering' and 'idx_covariance'."
      )
    }
  }

  config <- attr(x, "p.hclust_config")

  if (is.null(config) || !is.list(config)) {
    stop(
      "Invalid 'p.hclust' object: missing 'p.hclust_config' metadata."
    )
  }

  invisible(TRUE)
}


#' Test whether an object is a p.hclust
#'
#' @param x Object to test.
#'
#' @return TRUE or FALSE.
#' @export
is.p.hclust <- function(x) {
  inherits(x, "p.hclust")
}


#' Extract p.hclust configuration
#'
#' @param x A p.hclust object.
#'
#' @return A named list.
#' @keywords internal
p_hclust_config <- function(x) {
  if (!is.p.hclust(x)) {
    stop("'x' must be a 'p.hclust' object.")
  }

  config <- attr(x, "p.hclust_config")

  if (is.null(config)) {
    stop("Missing 'p.hclust_config' metadata.")
  }

  config
}
