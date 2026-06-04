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
      nrow = 5,
      ncol = 5,
      axes = "collect_y",
      axis_titles = "collect"
    )

  # save
  ggsave(
    filename = "plots/time_slices.pdf",
    plot = out,
    width = 7.5,
    height = 7.5
  )

  # return
  out

}
