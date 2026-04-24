#' Plot MCMC trace for rate parameters in BayesTraits model
#'
#' @param fit Tibble of posterior samples from the fitted model
#' @param tree_id Index for tree to plot
#'
#' @returns A ggplot object
#'
plot_mcmc_trace <- function(fit, tree_id) {
  # get tree id
  id <- tree_id
  # plot trace
  out <-
    fit |>
    filter(tree_id == id) |>
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
      title = paste0("Tree ID: ", id),
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
    filename = "plots/traceplot.pdf",
    plot = out,
    height = 4,
    width = 6
  )
  # return
  out
}
