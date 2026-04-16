#' Fit ancestral state reconstruction model
#'
#' Fit the ancestral state reconstruction model using the coevolve package
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_ids Indexes for trees to use
#' @param taxa_ids Tip labels
#'
#' @returns A coevfit object
#'
fit_model <- function(data, tree, tree_ids, taxa_ids) {
  # use subset of trees
  tree <- tree[tree_ids]
  # use subset of taxa
  tree <- keep.tip.multiPhylo(tree, taxa_ids)
  data <- filter(data, xd_id %in% taxa_ids)
  # fit model
  coev_fit(
    data = data,
    variables = list(class_differentiation = "ordered_logistic"),
    id = "xd_id",
    tree = tree,
    prior = list(eta_anc = "normal(0, 2)"),
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 1234
  )
}
