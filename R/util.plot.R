data.pdendrogram <- function(
  pvals,
  hcl
) {
  if (!inherits(hcl, "hclust")) {
    stop("'hcl' must inherit from 'hclust'.")
  }

  n <- length(hcl$order)

  expected_length <- n - 1L

  if (length(pvals) != expected_length) {
    stop(
      "'pvals' must have length n - 1."
    )
  }

  hcldata <- ggdendro::dendro_data(
    hcl,
    type = "rectangle"
  )

  hcldata$segments <- hcldata$segments[
    order(
      -hcldata$segments$y
    ),
  ]

  pvals_rep <- rep(
    pvals,
    each = 4L
  )

  if (length(pvals_rep) != nrow(hcldata$segments)) {
    stop(
      "Unexpected number of dendrogram segments returned by ",
      "'ggdendro::dendro_data()'."
    )
  }

  hcldata$segments$pval <- pvals_rep

  merge_points <- hcldata$segments |>
    dplyr::filter(
      y == yend
    ) |>
    dplyr::group_by(
      y
    ) |>
    dplyr::summarise(
      x_min = min(x),
      x_max = max(x),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(y)
    )

  n_merges <- nrow(merge_points)

  merge_points$x_mid <- (merge_points$x_min +
    merge_points$x_max) /
    2

  merge_points$k <- seq.int(
    from = 2L,
    length.out = n_merges
  )

  merge_points$label <- NA_character_

  n_labels <- min(
    length(pvals),
    n_merges
  )

  if (n_labels > 0L) {
    merge_points$label[
      seq_len(n_labels)
    ] <- signif(
      pvals[
        seq_len(n_labels)
      ],
      3
    ) |>
      as.character()
  }

  hcldata$merge_points <- merge_points

  hcldata
}


validate_groups_labels <- function(
  groups
) {
  if (
    !is.list(groups) ||
      length(groups) == 0L
  ) {
    stop("'groups' must be a non-empty list.")
  }

  n_groups <- length(groups)
  group_names <- names(groups)

  # Si la liste n'a pas de noms du tout
  if (is.null(group_names)) {
    group_names <- rep("", n_groups)
  }

  # Remplacer les noms manquants ou vides par "group_X"
  missing_names <- is.na(group_names) | !nzchar(group_names)
  if (any(missing_names)) {
    group_names[missing_names] <- paste0(
      "group_",
      which(missing_names)
    )
  }

  if (anyDuplicated(group_names)) {
    stop(
      "The names in 'groups' list must be unique."
    )
  }

  group_names
}


build_leaves <- function(
  groups,
  ord,
  labels = NULL
) {
  if (is.null(labels)) {
    labels <- validate_groups_labels(
      groups
    )
  }

  rows <- lapply(
    seq_along(groups),
    function(i) {
      group <- groups[[i]]

      if (length(group) != length(ord)) {
        stop(
          "Each element of 'groups' must have length ",
          "equal to the number of leaves in the dendrogram."
        )
      }

      data.frame(
        x = seq_along(ord),
        group_value = as.character(
          group[ord]
        ),
        group_id = i,
        group_name = labels[i],
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    rows
  )
}


assign_strip_positions <- function(
  leaves_long,
  strip_height_factor,
  base_y = 0,
  gap = 0.00001
) {
  n_indiv <- max(
    leaves_long$x,
    na.rm = TRUE
  )

  strip_height <- strip_height_factor *
    n_indiv

  gap_size <- gap *
    n_indiv

  leaves_long$strip_y <- -((leaves_long$group_id - 0.5) *
    strip_height +
    gap_size) +
    base_y

  list(
    data = leaves_long,
    strip_height = strip_height
  )
}


add_cluster_strips <- function(
  p,
  leaves_long,
  strip_height,
  colorblind_target
) {
  group_ids <- sort(
    unique(
      leaves_long$group_id
    )
  )

  group_info <- lapply(
    group_ids,
    function(gid) {
      df <- leaves_long[
        leaves_long$group_id == gid,
        ,
        drop = FALSE
      ]

      levels_val <- sort(
        unique(
          as.character(
            df$group_value
          )
        )
      )

      df$group_value <- factor(
        as.character(
          df$group_value
        ),
        levels = levels_val
      )

      list(
        gid = gid,
        data = df,
        levels = levels_val,
        name = unique(
          df$group_name
        )[1]
      )
    }
  )

  n_total_colors <- sum(
    vapply(
      group_info,
      function(x) length(x$levels),
      integer(1)
    )
  )

  all_colors <- Polychrome::createPalette(
    N = n_total_colors,
    seedcolors = "#f3756c",
    target = colorblind_target,
    range = c(30, 80),
    M = 50000
  )

  all_colors <- unname(
    all_colors
  )

  all_keys <- unlist(
    lapply(
      group_info,
      function(x) {
        paste0(
          "cluster_",
          x$gid,
          "::",
          x$levels
        )
      }
    )
  )

  names(all_colors) <- all_keys

  if (length(unique(all_colors)) != length(all_colors)) {
    stop(
      "The generated color palette contains duplicates."
    )
  }

  for (i in seq_along(group_info)) {
    cluster_i <- group_info[[i]]

    df <- cluster_i$data

    levels_val <- cluster_i$levels

    cluster_keys <- paste0(
      "cluster_",
      cluster_i$gid,
      "::",
      levels_val
    )

    colors_i <- all_colors[
      cluster_keys
    ]

    names(colors_i) <- levels_val

    p <- p +
      ggplot2::geom_tile(
        data = df,
        ggplot2::aes(
          x = x,
          y = strip_y,
          fill = group_value
        ),
        width = 1,
        height = strip_height
      ) +
      ggplot2::scale_fill_manual(
        values = colors_i,
        breaks = levels_val,
        drop = FALSE,
        name = cluster_i$name
      )

    if (i < length(group_info)) {
      p <- p +
        ggnewscale::new_scale_fill()
    }
  }

  p
}

#' @export
labels.groups <- function(
  data_pvalues,
  hcl,
  groups,
  group_names
) {
  clusters_list <- lapply(
    data_pvalues$merge_points$k,
    function(k) {
      cl <- get.individuals.merged.clusters(
        hcl,
        k
      )

      unique(
        c(
          cl[[1]],
          cl[[2]]
        )
      )
    }
  )

  cluster_sizes <- vapply(
    clusters_list,
    length,
    numeric(1)
  )

  for (i in seq_along(groups)) {
    group <- groups[[i]]
    group_name <- group_names[i]

    tmp <- lapply(
      clusters_list,
      function(indiv) {
        values <- group[
          indiv
        ]

        values <- values[
          !is.na(values)
        ]

        if (length(values) == 0L) {
          return(
            list(
              class = NA_character_,
              prop = NA_real_
            )
          )
        }

        tab <- table(
          values
        )

        idx <- which.max(
          tab
        )

        list(
          class = names(tab)[idx],
          prop = as.numeric(
            tab[idx] /
              sum(tab)
          )
        )
      }
    )

    data_pvalues$merge_points[[paste0(
      group_name,
      "_class"
    )]] <- vapply(
      tmp,
      `[[`,
      character(1),
      "class"
    )

    data_pvalues$merge_points[[paste0(
      group_name,
      "_prop"
    )]] <- vapply(
      tmp,
      `[[`,
      numeric(1),
      "prop"
    )
  }

  data_pvalues$merge_points$cluster_size <-
    cluster_sizes

  data_pvalues
}
