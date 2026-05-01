#' Simulate data and tree for model comparison
#'
#' Simulate a case where we know that all nodes before 11.7 kya were
#' egalitarian to test whether the model comparison works as expected.
#'
#' @returns Numeric. Log Bayes factor for model comparison.
#'
simulate_model_comparison <- function(data, mcc_tree) {
  # slice the tree into subtrees at 11.7 kya
  times <- ape::node.depth.edgelength(mcc_tree)
  tree_slice <- phytools::treeSlice(
    tree = mcc_tree,
    slice = max(times) - 11.7,
    trivial = TRUE
  )
  # simulate data on subtrees younger than 11.7 kya using rates from multistate
  sim_list <- lapply(tree_slice, function(phy) {
    ape::rTraitDisc(
      phy = phy,
      model = matrix(
        c(0.00, 0.15, 0.04,
          0.29, 0.00, 0.40,
          0.16, 0.20, 0.00),
        nrow = 3,
        byrow = TRUE
      ),
      # ancestor was egalitarian
      root.value = 1
    )
  })
  # add simulated values to dataset
  sim_states <- unlist(sim_list)
  sim_data <-
    data |>
    left_join(
      tibble(
        xd_id = names(sim_states),
        sim = as.character(sim_states)
      ),
      by = "xd_id"
    ) |>
    mutate(
      sim = ifelse(sim == "A", "Absence of distinctions", sim),
      sim = ifelse(sim == "B", "Wealth distinctions", sim),
      sim = ifelse(sim == "C", "Elite stratification", sim),
      class_differentiation =
        ordered(sim, levels = levels(data$class_differentiation))
    ) |>
    dplyr::select(!sim)
  # return result of model comparison
  run_model_comparison(sim_data, mcc_tree)
}
