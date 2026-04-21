#' Combine and save plots of results
#'
#' @param plot_list List of ggplots objects to combine
#'
#' @returns A patchwork composition of ggplot objects
#'
combine_plots <- function(plot_list) {
  # combine plots
  out <-
    wrap_plots(plot_list) +
    plot_layout(
      nrow = 4,
      ncol = 4,
      axes = "collect_y",
      axis_titles = "collect"
    )
  # save
  ggsave(
    filename = "plots/results.pdf",
    plot = out,
    width = 7,
    height = 7
  )
  # return
  out
}
