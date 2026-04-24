#' Calculate diagnostics for the ancestral state reconstruction model
#'
#' @param fit Tibble of posterior samples from the fitted model
#'
#' @returns A tibble of model diagnostics
#'
calculate_model_diagnostics <- function(fit) {
  fit |>
    dplyr::select(-c(Iteration, `Tree No`)) |>
    # collect posterior samples
    group_by(tree_id, chain) |>
    summarise(across(everything(), list), .groups = "drop") |>
    # get parameters in rows and mcmc chains in columns
    pivot_longer(
      cols = !c(tree_id, chain),
      names_to = "parameter"
    ) |>
    pivot_wider(
      names_from = chain,
      values_from = value
    ) |>
    # calculate model diagnostics with the rstan package
    rowwise() |>
    mutate(
      post = list(as.matrix(cbind(`1`, `2`, `3`, `4`))),
      rhat = rstan::Rhat(post),
      ess_bulk = rstan::ess_bulk(post),
      ess_tail = rstan::ess_tail(post)
    ) |>
    ungroup() |>
    dplyr::select(!c(`1`, `2`, `3`, `4`, `post`))
}
