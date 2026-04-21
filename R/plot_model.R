#' Plot results of ancestral state reconstruction model
#'
#' @param ancestral_states Tibble of ancestral states from the model
#' @param tree Tree object of class multiPhylo
#' @param tree_ids Indexes for trees
#' @param family (optional) Character. Resulting plot will summarise
#'   only ancestral nodes for taxa in a particular language family.
#'
#' @returns A ggplot object
#'
plot_model <- function(data, ancestral_states, tree, tree_ids,
                       family = NULL) {
  # use subset of trees
  tree <- tree[tree_ids]
  # for each tree, get internal nodes, times, and trait values
  edges <-
    tibble(tree_id = 1:length(tree_ids)) |>
    mutate(
      # get tree
      tree = map(tree_id, function(tree_id) tree[[tree_id]]),
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
      ancestral_states |>
        rowwise() |>
        transmute(
          tree_id = tree,
          parent_node = node,
          prob_start = list(prob3 + prob4 + prob5)
        ),
      by = c("tree_id", "parent_node")
    ) |>
    left_join(
      ancestral_states |>
        rowwise() |>
        transmute(
          tree_id = tree,
          child_node = node,
          prob_end = list(prob3 + prob4 + prob5)
        ),
      by = c("tree_id", "child_node")
    )
  # if family specified, filter to ancestors of taxa in the language family
  if (!is.null(family)) {
    edges <-
      edges |>
      mutate(
        # is this node an ancestral node for taxa in this language family?
        is_ancestor =
          map2(tree_id, parent_node, function(tree_id, parent_node) {
            # get taxa in language family according to glottolog
            taxa <-
              data |>
              filter(!is.na(language_family) & language_family == family) |>
              pull(xd_id)
            # get all ancestor nodes
            ancestors <- Ancestors(
              x = tree[[tree_id]],
              node = taxa,
              type = "all"
            )
            # is the current node in the list of ancestors?
            parent_node %in% unlist(ancestors)
          })
      ) |>
      # filter to ancestors only
      unnest(is_ancestor) |>
      filter(is_ancestor)
  }
  # summarise lineages at time slices
  out <-
    map(seq(-20, -0.25, by = 0.25), function(t) {
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
        group_by(iter, tree_id) |>
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
    theme(plot.title = element_text(size = 7))
  # cleanup
  rm(data, ancestral_states, edges, tree, tree_ids)
  # return
  out
}
