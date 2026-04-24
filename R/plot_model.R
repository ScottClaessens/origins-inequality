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
  # get ancestral states
  d <-
    fit |>
    dplyr::select(c(tree_id, ends_with("P(3)") & !starts_with("Root"))) |>
    pivot_longer(
      cols = !tree_id,
      names_to = "node",
      values_to = "prob_inequality"
    ) |>
    mutate(node = parse_number(str_sub(node, 2, 5))) |>
    group_by(tree_id, node) |>
    summarise(
      prob_inequality = list(prob_inequality),
      .groups = "drop"
    ) |>
    rename(id = tree_id) |>
    ungroup()
  # for each tree, get internal nodes, times, and trait values
  edges <-
    tibble(id = tree_id) |>
    mutate(
      # get tree
      tree = map(id, function(id) tree[[id]]),
      # get all edges in the tree
      edges = map(tree, function(tree) tree$edge),
      # get timings of the splits in the tree
      times = map(tree, function(tree) {
        node.depth.edgelength(tree) - max(node.depth.edgelength(tree))
      }),
      # get parent and child nodes
      parent_node = map(edges, function(edges) edges[, 1]),
      child_node = map(edges, function(edges) edges[, 2]),
      # get start and end times
      time_start = map2(times, edges, function(times, edges) times[edges[, 1]]),
      time_end = map2(times, edges, function(times, edges) times[edges[, 2]])
    ) |>
    dplyr::select(-c(tree, edges, times)) |>
    unnest(c(parent_node, child_node, time_start, time_end)) |>
    left_join(
      transmute(
        d,
        id = id,
        parent_node = node,
        prob_start = prob_inequality
      ),
      by = c("id", "parent_node")
    ) |>
    left_join(
      transmute(
        d,
        id = id,
        child_node = node,
        prob_end = prob_inequality
      ),
      by = c("id", "child_node")
    )
  # if family specified, filter to ancestors of taxa in the language family
  if (!is.null(family)) {
    edges <-
      edges |>
      mutate(
        # is this node an ancestral node for taxa in this language family?
        is_ancestor =
          map2(id, parent_node, function(id, parent_node) {
            # get taxa in language family according to glottolog
            taxa <-
              data |>
              filter(!is.na(language_family) & language_family == family) |>
              pull(xd_id)
            # get most recent common ancestor
            mrca <-
              getMRCA(
                phy = tree[[id]],
                tip = taxa
              )
            # get all ancestral nodes
            ancestors <-
              Ancestors(
                x = tree[[id]],
                node = taxa,
                type = "all"
              )
            # keep only mrca and younger
            ancestors <- unlist(ancestors)
            ancestors <- unique(ancestors[ancestors >= mrca])
            # is the current node in the list of ancestors?
            parent_node %in% ancestors
          })
      ) |>
      # filter to ancestors only
      unnest(is_ancestor) |>
      filter(is_ancestor)
  }
  # summarise lineages at time slices
  out <-
    map(seq(start_time, end_time, by = time_slice), function(t) {
      # get number of posterior samples
      n_samples <- length(edges$prob_start[[1]])
      # wrangle data
      edges |>
        dplyr::select(!prob_end) |>
        # filter to lineages alive at time t
        filter(time_start <= t & time_end > t) |>
        # get number of lineages alive at time t
        mutate(n_lineages = n()) |>
        # unnest posterior samples
        unnest(c(prob_start)) |>
        # add index for each posterior sample
        mutate(iter = rep_len(1:n_samples, length.out = n())) |>
        # get average across lineages
        group_by(id, iter) |>
        summarise(
          prob_start = median(prob_start),
          n_lineages = unique(n_lineages),
          .groups = "drop"
        ) |>
        # summary of posterior samples
        summarise(
          time = t,
          median = median(prob_start),
          lower95 = quantile(prob_start, 0.975),
          upper95 = quantile(prob_start, 0.025),
          lower50 = quantile(prob_start, 0.25),
          upper50 = quantile(prob_start, 0.75),
          lower25 = quantile(prob_start, 0.375),
          upper25 = quantile(prob_start, 0.625),
          n_lineages = median(n_lineages)
        ) |>
        mutate(n_lineages = ifelse(is.na(n_lineages), 0, n_lineages))
    }) |>
    list_rbind() |>
    # plot time slices
    ggplot() +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower95,
        ymax = upper95
      ),
      fill = "grey90"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower50,
        ymax = upper50
      ),
      fill = "grey80"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower25,
        ymax = upper25
      ),
      fill = "grey70"
    ) +
    geom_line(
      aes(
        x = time,
        y = median
      )
    ) +
    scale_x_continuous(name = "Time before present (thousands of years)") +
    scale_y_continuous(
      name = "Probability of inequality",
      limits = c(0, 1)
    ) +
    ggtitle(ifelse(!is.null(family), family, "Global")) +
    theme_classic() +
    theme(plot.title = element_text(size = 9))
  # cleanup
  rm(d, data, edges, fit, tree, family, tree_id)
  # return
  out
}
