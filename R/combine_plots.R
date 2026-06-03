#' Combine and save plots of ancestral state recontruction time slices
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
      nrow = 2,
      ncol = 3,
      axes = "collect_y",
      axis_titles = "collect"
    )

  # save
  ggsave(
    filename = "plots/time_slices.pdf",
    plot = out,
    width = 6,
    height = 4
  )

  # return
  out

}
