#' Plot results of ancestral state reconstruction model
#'
#' @param data Tibble of D-PLACE data
#' @param fit Results of fitted model
#' @param tree Tree object of class multiPhylo
#' @param tree_id Indexes for trees
#' @param family (optional) Character. Resulting plot will summarise
#'   only ancestral nodes for taxa in a particular language family.
#' @param start_time Start of the time window (thousands of years before
#'   present). Defaults to -20.
#' @param end_time End of the time window (thousands of years before
#'   present). Defaults to -0.25.
#' @param time_slice How frequently to take time slices (in thousands of years
#'   before present). Defaults to 0.25.
#'
#' @returns A ggplot object
#'
plot_model <- function(data, fit, tree, tree_id, family = NULL,
                       start_time = -20, end_time = -0.25, time_slice = 0.25) {

  # get sequence of time slices
  times_seq <- seq(start_time, end_time, by = time_slice)

  # loop over trees
  results <- lapply(tree_id, function(id) {

    # get tree
    tree_obj <- tree[[id]]

    # extract edges and timings
    edges <- tree_obj$edge
    times <- ape::node.depth.edgelength(tree_obj)
    times <- times - max(times)
    parent_node <- edges[, 1]
    child_node <- edges[, 2]
    time_start <- times[parent_node]
    time_end <- times[child_node]

    # get posterior samples for this tree
    prob_matrix <-
      fit |>
      filter(tree_id == id) |>
      dplyr::select(ends_with("P(3)") & !starts_with("Root")) |>
      as.matrix()

    colnames(prob_matrix) <-
      parse_number(str_sub(colnames(prob_matrix), 2, 5))

    # family filtering (once per tree)
    if (!is.null(family)) {

      # get taxa in language family
      taxa <-
        data |>
        filter(!is.na(language_family) & language_family == family) |>
        pull(xd_id)

      # get most recent common ancestor
      mrca <- ape::getMRCA(tree_obj, taxa)

      # get ancestors
      ancestors <- unique(unlist(
        phangorn::Ancestors(tree_obj, node = taxa, type = "all")
      ))

      # retain ancestors younger than mrca
      ancestors <- ancestors[ancestors >= mrca]

      # filter to family
      keep <- parent_node %in% ancestors
      parent_node <- parent_node[keep]
      child_node <- child_node[keep]
      time_start <- time_start[keep]
      time_end <- time_end[keep]
    }

    # time slicing (vectorised)
    tree_res <- lapply(times_seq, function(t) {

      # is lineage active at this time?
      active <- (time_start <= t) & (time_end > t)

      # number of lineages active
      n_lineages <- sum(active)

      # if no lineages, return empty tibble
      if (n_lineages == 0) {
        return(
          tibble(
            id = id,
            time = t,
            maximums = NA,
            n_lineages = 0
          )
        )
      }

      # filter posterior samples to active lineages
      mat <- prob_matrix[, as.character(parent_node)[active], drop = FALSE]

      # get maximum probability of inequality across all active lineages
      maximums <- apply(mat, 1, max)

      # return tibble
      tibble(
        id = id,
        time = t,
        maximums = list(maximums),
        n_lineages = n_lineages
      )
    })

    # bind results
    do.call(rbind, tree_res)

  })

  out <- do.call(rbind, results)

  # aggregate across trees
  out_summary <-
    out |>
    unnest(maximums) |>
    group_by(time) |>
    dplyr::summarise(
      median  = median(maximums, na.rm = TRUE),
      lower95 = quantile(maximums, 0.025, na.rm = TRUE),
      upper95 = quantile(maximums, 0.975, na.rm = TRUE),
      lower50 = quantile(maximums, 0.375, na.rm = TRUE),
      upper50 = quantile(maximums, 0.625, na.rm = TRUE),
      lower25 = quantile(maximums, 0.025, na.rm = TRUE),
      upper25 = quantile(maximums, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  # plot time slices
  p <-
    ggplot(data = out_summary) +
    geom_hline(
      yintercept = 1/3,
      linewidth = 0.1,
      linetype = "dashed"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower95,
        ymax = upper95
      ),
      fill = "grey90",
      alpha = 0.5
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower50,
        ymax = upper50
      ),
      fill = "grey80",
      alpha = 0.5
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower25,
        ymax = upper25
      ),
      fill = "grey70",
      alpha = 0.5
    ) +
    geom_line(
      aes(
        x = time,
        y = median
      )
    ) +
    scale_x_continuous(name = "Time before present (thousands of years)") +
    scale_y_continuous(
      name = "Maximum probability of inequality",
      limits = c(0, 1)
    ) +
    ggtitle(ifelse(!is.null(family), family, "Global")) +
    theme_classic() +
    theme(plot.title = element_text(size = 9))
  # cleanup
  rm(data, fit, out, out_summary, results, tree, start_time, end_time,
     time_slice, times_seq, tree_id, family)
  # return
  p
}
