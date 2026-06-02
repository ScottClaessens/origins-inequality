#' Plot MCMC trace for rate parameters in BayesTraits model
#'
#' @param fit Tibble of posterior samples from the fitted model
#' @param model Character of length 1. Which model was fitted. One of: "full",
#'   "rectilinear", "unilinear", "relaxed_unilinear", "alternative", or
#'   "alternative_reversible"
#'
#' @returns A ggplot object
#'
plot_mcmc_trace <- function(fit, model = "full") {

  # plot trace
  out <-
    fit |>
    dplyr::select(c(chain, Iteration, q12, q13, q21, q23, q31, q32)) |>
    pivot_longer(
      cols = !c(chain, Iteration),
      names_to = "parameter"
    ) |>
    ggplot(
      aes(
        x = Iteration / 1000,
        y = value,
        colour = as.factor(chain)
      )
    ) +
    geom_line() +
    facet_wrap(
      . ~ parameter,
      scales = "free_y"
    ) +
    labs(
      x = "Iterations (in thousands)",
      y = NULL
    ) +
    scale_colour_discrete(
      name = "Chain",
      palette = "BuPu"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7))

  # save
  ggsave(
    filename = paste0("plots/traceplots/traceplot_", model, ".pdf"),
    plot = out,
    height = 4,
    width = 6
  )

  # cleanup
  rm(fit, tree_id, id)

  # return
  out

}
