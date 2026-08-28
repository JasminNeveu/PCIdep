#' Plot a p.hclust object
#'
#' @param x A p.hclust object.
#' @param plot_config Named list controlling graphical parameters.
#' @param groups Optional named list of group annotations.
#' @param ... Additional graphical arguments.
#'
#' @return A ggplot object.
#' @export
plot.p.hclust <- function(
  x,
  plot_config = list(),
  groups = NULL,
  ...
) {
  validate_p_hclust(x)
  n_clust <- length(x$order)

  # --- Filtrage automatique des groupes en cas de sample splitting ---
  if (!is.null(groups) && is.list(groups) && !is.null(x$split_indices)) {
    idx_clust <- x$split_indices$idx_clustering

    groups <- lapply(groups, function(g) {
      if (is.null(g)) {
        return(NULL)
      }

      # Cas 1 : Le vecteur a la taille initiale N (supérieure à n_clust) -> On filtre sur idx_clustering
      if (length(g) > n_clust && length(g) >= max(idx_clust)) {
        return(g[idx_clust])
      } else if (length(g) == n_clust) {
        # Cas 2 : Le vecteur a déjà été pré-filtré à la taille n_clust -> On le conserve tel quel
        return(g)
      } else {
        # Cas 3 : Dimension incompatible
        warning(sprintf(
          "La taille de l'annotation (%d) ne correspond ni au nombre total d'individus, ni au sous-ensemble de clustering (%d).",
          length(g),
          n_clust
        ))
        return(g)
      }
    })
  }

  if (!is.list(plot_config)) {
    stop("'plot_config' must be a list.")
  }

  config <- list(
    labels_pvalues = FALSE,
    strip_height_factor = 1,
    log_height = FALSE,
    colorblind_target = "normal",
    break_limits = NULL,
    low = NULL,
    mid = NULL,
    high = NULL
  )

  unknown_config <- setdiff(
    names(plot_config),
    names(config)
  )

  if (length(unknown_config) > 0L) {
    stop(
      "Unknown plot_config entries: ",
      paste(unknown_config, collapse = ", ")
    )
  }

  config <- utils::modifyList(
    config,
    plot_config
  )

  labels_pvalues <- config$labels_pvalues
  strip_height_factor <- config$strip_height_factor
  log_height <- config$log_height
  colorblind_target <- config$colorblind_target
  break_limits <- config$break_limits

  low <- config$low
  mid <- config$mid
  high <- config$high

  if (
    !is.logical(labels_pvalues) ||
      length(labels_pvalues) != 1L ||
      is.na(labels_pvalues)
  ) {
    stop("'labels_pvalues' must be TRUE or FALSE.")
  }

  if (
    !is.numeric(strip_height_factor) ||
      length(strip_height_factor) != 1L ||
      is.na(strip_height_factor) ||
      strip_height_factor <= 0
  ) {
    stop("'strip_height_factor' must be positive.")
  }

  if (
    !is.logical(log_height) ||
      length(log_height) != 1L ||
      is.na(log_height)
  ) {
    stop("'log_height' must be TRUE or FALSE.")
  }

  valid_targets <- c(
    "normal",
    "protanope",
    "deuteranope",
    "tritanope"
  )

  if (!colorblind_target %in% valid_targets) {
    stop(
      "'colorblind_target' must be one of: ",
      paste(valid_targets, collapse = ", ")
    )
  }

  filled <- c(
    !is.null(low),
    !is.null(mid),
    !is.null(high)
  )

  partial_defined <- sum(filled) > 0L &&
    sum(filled) < 3L

  if (partial_defined) {
    stop(
      "'low', 'mid' and 'high' must either all be provided ",
      "or all be NULL."
    )
  }

  if (all(filled)) {
    colours <- c(
      low,
      mid,
      high
    )
  } else {
    colours <- c(
      "#B2182B",
      "#D9C27A",
      "#2166AC"
    )
  }

  if (!is.null(break_limits)) {
    if (length(break_limits) == 1L) {
      if (
        !is.numeric(break_limits) ||
          is.na(break_limits) ||
          break_limits <= 0
      ) {
        stop(
          "'break_limits' must be NULL, a positive scalar, ",
          "or a numeric vector of length 2."
        )
      }

      break_limits <- c(
        0,
        break_limits
      )
    } else if (length(break_limits) == 2L) {
      if (
        !is.numeric(break_limits) ||
          anyNA(break_limits) ||
          break_limits[1] >= break_limits[2]
      ) {
        stop(
          "When length 2, 'break_limits' must contain ",
          "increasing numeric limits."
        )
      }
    } else {
      stop(
        "'break_limits' must be NULL, a scalar, ",
        "or a numeric vector of length 2."
      )
    }
  }

  p_hcl <- x

  if (log_height) {
    p_hcl$height <- log10(
      p_hcl$height + 1
    )
  }

  analysis_config <- p_hclust_config(x)

  alpha <- analysis_config$alpha

  p_floor <- 2.2e-16

  breaks <- c(
    p_floor,
    alpha,
    1
  )

  vals <- scales::rescale(
    log10(breaks)
  )

  data_pvalues <- data.pdendrogram(
    pvals = x$pvalues,
    hcl = p_hcl
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = data_pvalues$segments,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = pval
      )
    ) +
    ggplot2::labs(
      title = "p-dendrogram",
      x = "Observations",
      y = "Merging distances"
    ) +
    ggplot2::scale_color_gradientn(
      colours = colours,
      values = vals,
      trans = "log10",
      limits = c(
        p_floor,
        1
      ),
      breaks = breaks,
      labels = c(
        "< 2.2e-16",
        as.character(alpha),
        "1"
      ),
      na.value = "grey",
      name = "p-value"
    ) +
    ggplot2::theme_bw()

  if (log_height) {
    y_breaks <- pretty(
      p_hcl$height,
      n = 6
    )

    p <- p +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        labels = function(z) {
          parse(
            text = paste0(
              "10^",
              z
            )
          )
        }
      )
  }

  if (!is.null(groups)) {
    resolved_labels <- validate_groups_labels(
      groups
    )

    leaves_long <- build_leaves(
      groups = groups,
      ord = x$order,
      labels = resolved_labels
    )

    tmp <- assign_strip_positions(
      leaves_long,
      strip_height_factor
    )

    leaves_long <- tmp$data
    strip_height <- tmp$strip_height

    p <- add_cluster_strips(
      p,
      leaves_long,
      strip_height,
      colorblind_target
    )

    # Calcul des proportions pour le hover_text (uniquement si groups est présent)
    data_pvalues <- labels.groups(
      data_pvalues = data_pvalues,
      hcl = x,
      groups = groups,
      group_names = resolved_labels
    )

    df <- data_pvalues$merge_points

    data_pvalues$merge_points$hover_text <- vapply(
      seq_len(nrow(df)),
      function(i) {
        row <- df[
          i,
          ,
          drop = TRUE
        ]

        txt <- paste0(
          "\nsize = ",
          row[["cluster_size"]]
        )

        for (group_name in resolved_labels) {
          class_col <- paste0(
            group_name,
            "_class"
          )

          prop_col <- paste0(
            group_name,
            "_prop"
          )

          txt <- paste0(
            txt,
            "\n",
            group_name,
            ": ",
            row[[class_col]],
            " (",
            round(
              100 *
                as.numeric(
                  row[[prop_col]]
                ),
              1
            ),
            "%)"
          )
        }

        txt
      },
      character(1)
    )
  }

  # 2. Affichage des p-values sur l'arbre (ne nécessite plus groups)
  if (labels_pvalues) {
    p <- p +
      ggplot2::geom_text(
        data = data_pvalues$merge_points |>
          dplyr::filter(
            !is.na(label)
          ),
        ggplot2::aes(
          x = x_mid,
          y = y,
          label = label
        ),
        vjust = "top",
        size = 3
      )
  }

  # 3. Cassure de l'axe Y si demandée
  if (!is.null(break_limits)) {
    p <- p +
      ggbreak::scale_y_break(
        break_limits
      ) +
      ggplot2::theme(
        axis.text.y.right = ggplot2::element_blank(),
        axis.ticks.y.right = ggplot2::element_blank(),
        axis.line.y.right = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_blank()
      )
  }

  p
}
