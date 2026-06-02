#' Simulate data and tree for model comparison
#'
#' Simulate a case where we know that all nodes before 11.7 kya were
#' a particular value to test whether the model comparison works as expected.
#'
#' @param data Tibble of D-PLACE data
#' @param mcc_tree Maximum clade credibility tree object of class phylo
#' @param fossil Value to fossilise
#'
#' @returns Numeric. Log Bayes factor for model comparison.
#'
simulate_model_comparison <- function(data, mcc_tree, fossil = "1") {

  # slice the tree into subtrees at 11.7 kya
  times <- ape::node.depth.edgelength(mcc_tree)
  tree_slice <- phytools::treeSlice(
    tree = mcc_tree,
    slice = max(times) - 11.7,
    trivial = TRUE
  )

  # simulate data on subtrees younger than 11.7 kya
  # under equal-rates relaxed unilinear model
  sim_list <- lapply(tree_slice, function(phy) {
    ape::rTraitDisc(
      phy = phy,
      model = matrix(
        c(0.00, 0.10, 0.00,
          0.10, 0.00, 0.10,
          0.10, 0.10, 0.00),
        nrow = 3,
        byrow = TRUE
      ),
      # set ancestor
      root.value = as.numeric(fossil)
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

  # run model comparison
  log_lik_1 <- get_log_lik_fossilised(sim_data, mcc_tree, "1")
  log_lik_2 <- get_log_lik_fossilised(sim_data, mcc_tree, "2")
  log_lik_3 <- get_log_lik_fossilised(sim_data, mcc_tree, "3")

  # return table of log bayes factors
  get_table_fossilised(log_lik_1, log_lik_2, log_lik_3)

}
