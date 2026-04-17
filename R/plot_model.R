#' Plot results of ancestral state reconstruction model
#'
#' @param ancestral_states Tibble of ancestral states from the model
#' @param tree Tree object of class multiPhylo
#' @param tree_ids Indexes for trees
#'
#' @returns A ggplot object
#'
plot_model <- function(ancestral_states, tree, tree_ids) {
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
  # summarise lineages at time slices
  out <-
    map(seq(-20, -0.25, by = 0.25), function(t) {
      edges |>
        dplyr::select(!prob_end) |>
        filter(time_start <= t & time_end > t) |>
        mutate(n = n()) |>
        unnest(c(prob_start)) |>
        summarise(
          time = t,
          median = median(prob_start),
          lower50 = quantile(prob_start, 0.25),
          upper50 = quantile(prob_start, 0.75),
          lower25 = quantile(prob_start, 0.375),
          upper25 = quantile(prob_start, 0.625),
          n = unique(n)
        )
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
      fill = "grey50"
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
    theme_classic()
  # save plot and return
  ggsave(
    filename = "plots/results.pdf",
    plot = out,
    height = 4,
    width = 5
  )
  # cleanup
  rm(ancestral_states, edges, tree, tree_ids)
  # return
  out
}
