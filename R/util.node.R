get.node.sizes <- function(hcl) {
  n <- length(hcl$order)

  sizes <- integer(n - 1L)

  for (i in seq_len(n - 1L)) {
    c1 <- hcl$merge[i, 1]
    c2 <- hcl$merge[i, 2]

    s1 <- if (c1 < 0) 1L else sizes[c1]
    s2 <- if (c2 < 0) 1L else sizes[c2]

    sizes[i] <- s1 + s2
  }

  sizes
}


get.merged.clusters <- function(hcl, k) {
  n <- length(hcl$order)

  if (
    length(k) != 1L ||
      is.na(k) ||
      k != as.integer(k) ||
      k < 2L ||
      k > n
  ) {
    stop("'k' must be an integer between 2 and n.")
  }

  k <- as.integer(k)

  merge_step <- n - k + 1L

  c1_internal <- hcl$merge[
    merge_step,
    1
  ]

  c2_internal <- hcl$merge[
    merge_step,
    2
  ]

  clusters_k <- stats::cutree(
    hcl,
    k = k
  )

  get_leaf_index <- function(id) {
    if (id < 0) {
      return(-id)
    }

    get_leaf_index(
      hcl$merge[id, 1]
    )
  }

  leaf1 <- get_leaf_index(c1_internal)
  leaf2 <- get_leaf_index(c2_internal)

  c(
    unname(clusters_k[leaf1]),
    unname(clusters_k[leaf2])
  )
}


get.individuals.merged.clusters <- function(hcl, k) {
  n <- length(hcl$order)

  if (
    length(k) != 1L ||
      is.na(k) ||
      k != as.integer(k) ||
      k < 2L ||
      k > n
  ) {
    stop("'k' must be an integer between 2 and n.")
  }

  k <- as.integer(k)

  merge_step <- n - k + 1L

  c1_internal <- hcl$merge[
    merge_step,
    1
  ]

  c2_internal <- hcl$merge[
    merge_step,
    2
  ]

  get_individuals <- function(id) {
    if (id < 0) {
      return(-id)
    }

    left <- get_individuals(
      hcl$merge[id, 1]
    )

    right <- get_individuals(
      hcl$merge[id, 2]
    )

    c(left, right)
  }

  individuals_1 <- get_individuals(
    c1_internal
  )

  individuals_2 <- get_individuals(
    c2_internal
  )

  list(
    unname(individuals_1),
    unname(individuals_2)
  )
}


get.node.leafs <- function(
  hcl,
  row,
  node_sizes = get.node.sizes(hcl)
) {
  if (
    length(row) != 1L ||
      is.na(row) ||
      row < 1L ||
      row > nrow(hcl$merge)
  ) {
    stop("'row' must be a valid internal node index.")
  }

  n_obs <- length(hcl$order)

  leafs <- integer(
    node_sizes[row]
  )

  idx <- 1L

  stack <- integer(n_obs)
  stack_ptr <- 1L

  stack[stack_ptr] <- row

  while (stack_ptr > 0L) {
    curr <- stack[stack_ptr]

    stack_ptr <- stack_ptr - 1L

    c1 <- hcl$merge[curr, 1]
    c2 <- hcl$merge[curr, 2]

    if (c1 < 0) {
      leafs[idx] <- -c1

      idx <- idx + 1L
    } else {
      stack_ptr <- stack_ptr + 1L

      stack[stack_ptr] <- c1
    }

    if (c2 < 0) {
      leafs[idx] <- -c2

      idx <- idx + 1L
    } else {
      stack_ptr <- stack_ptr + 1L

      stack[stack_ptr] <- c2
    }
  }

  leafs
}
