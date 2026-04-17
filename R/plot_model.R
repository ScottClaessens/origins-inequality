#' Plot results of ancestral state reconstruction model
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param fit Fitted coevolve model
#' @param tree_ids Indexes for trees
#'
#' @returns A ggplot object
#'
plot_model <- function(data, tree, fit, tree_ids) {
  # use subset of trees
  tree <- tree[tree_ids]
  # extract samples
  post <- extract_samples(fit)
  rm(fit)
  # get ancestral states
  ancestral_states <- tibble()
  for (t in 1:length(tree_ids)) {
    for (n in 1:Nnode(tree[t], internal.only = FALSE)) {
      # get latent trait value
      eta <- post$eta[, t, n, 1]
      # get cumulative probabilities of each ordinal category
      p <- plogis(post$c1 - eta)
      # add row
      ancestral_states <-
        bind_rows(
          ancestral_states,
          tibble(
            tree = t,
            node = n,
            eta = list(eta),
            prob1 = list(p[, 1]),
            prob2 = list(p[, 2] - p[, 1]),
            prob3 = list(p[, 3] - p[, 2]),
            prob4 = list(p[, 4] - p[, 3]),
            prob5 = list(1 - p[, 4])
          )
        )
    }
  }
  # cleanup
  rm(post, t, n, p, eta)
  # get every language family in D-PLACE with at least ten associated taxon
  families <-
    data |>
    filter(!is.na(language_family)) |>
    group_by(language_family) |>
    summarise(n = n()) |>
    filter(n >= 10) |>
    pull(language_family)
  # reconstruct most recent common ancestors for language families
  family_reconstruct <- tibble()
  # for all language families:
  for (family in families) {
    # get taxa in language family according to D-PLACE
    taxa <- data$xd_id[data$language_family == family]
    taxa <- taxa[!is.na(taxa)]
    # then, for every tree:
    for (t in 1:length(tree_ids)) {
      # 1. get most recent common ancestor of taxa in language family
      node_id <- getMRCA(tree[[t]], taxa)
      # 2. record the node time
      node_times <-
        node.depth.edgelength(tree[[t]]) -
        max(node.depth.edgelength(tree[[t]]))
      node_time <- node_times[node_id]
      # 3. record the posterior probability of inequality
      post_prob <-
        ancestral_states |>
        filter(tree == t & node == node_id) |>
        rowwise() |>
        mutate(post = list(prob3 + prob4 + prob5)) |>
        pull(post) |>
        unlist()
      # add row
      family_reconstruct <-
        bind_rows(
          family_reconstruct,
          tibble(
            language_family = family,
            tree = t,
            node = node_id,
            time = node_time,
            prob = list(post_prob)
          )
        )
    }
  }
  # cleanup
  rm(families, family, node_id, node_time, node_times, taxa, post_prob)
  # plot language family reconstructions
  family_reconstruct |>
    unnest(c(prob)) |>
    ggplot(
      mapping = aes(
        x = time,
        y = log(prob) - log(1 - prob),
        group = language_family
      )
    ) +
    stat_ellipse(
      geom = "polygon",
      type = "norm",
      fill = "skyblue1",
      level = 0.99
    ) +
    stat_ellipse(
      geom = "polygon",
      type = "norm",
      fill = "skyblue2",
      level = 0.90
    ) +
    stat_ellipse(
      geom = "polygon",
      type = "norm",
      fill = "skyblue3",
      level = 0.80
    ) +
    stat_ellipse(
      geom = "polygon",
      type = "norm",
      fill = "skyblue4",
      level = 0.50
    ) +
    scale_x_continuous(limits = c(-25, 1)) +
    facet_wrap(. ~ language_family) +
    theme_classic()
}
