#' Get ancestral states estimated by the model
#'
#' @param tree Tree object of class multiPhylo
#' @param fit Fitted coevolve model
#' @param tree_ids Indexes for trees
#'
#' @returns A tibble
#'
get_ancestral_states <- function(tree, fit, tree_ids) {
  # use subset of trees
  tree <- tree[tree_ids]
  # extract samples
  post <- extract_samples(fit)
  # get numbers of nodes and tips
  n_nodes <- Nnode(tree[[1]], internal.only = FALSE)
  n_tips <- length(tree[[1]]$tip.label)
  # get ancestral states across trees and nodes
  out <-
    expand_grid(
      tree = 1:length(tree_ids),
      node = 1:n_nodes
    ) |>
    mutate(
      # is node a tip?
      is_tip = node %in% 1:n_tips,
      # get latent trait value
      eta = map2(tree, node, function(tree, node) post$eta[, tree, node, 1]),
      # get cumulative probabilities
      p = map(eta, function(eta) plogis(post$c1 - eta)),
      # get probabilities
      prob1 = map(p, function(p) p[, 1]),
      prob2 = map(p, function(p) p[, 2] - p[, 1]),
      prob3 = map(p, function(p) p[, 3] - p[, 2]),
      prob4 = map(p, function(p) p[, 4] - p[, 3]),
      prob5 = map(p, function(p) 1 - p[, 4])
    ) |>
    dplyr::select(-p)
  # cleanup
  rm(tree, fit, post, tree_ids)
  # return
  out
}
