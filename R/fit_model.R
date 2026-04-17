#' Fit ancestral state reconstruction model
#'
#' Fit the ancestral state reconstruction model using the coevolve package
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_ids Indexes for trees to use
#'
#' @returns A coevfit object
#'
fit_model <- function(data, tree, tree_ids) {
  # use subset of trees
  tree <- tree[tree_ids]
  # fit model
  coev_fit(
    data = data,
    variables = list(class_differentiation = "ordered_logistic"),
    id = "xd_id",
    tree = tree,
    parallel_chains = 4,
    iter_warmup = 400,
    iter_sampling = 400,
    seed = 1234
  )
}
