#' Combine and save plots of ancestral state recontruction time slices
#'
#' @param plot_list List of ggplots objects to combine
#'
#' @returns A patchwork composition of ggplot objects
#'
combine_plots <- function(plot_list) {

  # get global plot
  left_plot <-
    plot_list[1][[1]] +
    theme(legend.position = "none")

  # combine family plots
  right_plot <-
    wrap_plots(plot_list[-1]) +
    plot_layout(
      nrow = 5,
      ncol = 5,
      axes = "collect",
      axis_titles = "collect",
      guides = "collect"
    ) &
    scale_x_continuous(
      name = "Time before present (thousands of years)",
      breaks = c(-20, -10, 0)
    ) &
    scale_y_continuous(
      name = NULL,
      limits = c(0, 1),
      breaks = c(0, 1)
    ) &
    theme(
      plot.title = element_text(size = 6),
      legend.key.height = unit(15, "pt"),
      legend.spacing.y = unit(0, "pt")
    )

  # combine all plots
  out <-
    left_plot + right_plot +
    plot_layout(axis_titles = "collect")

  # save
  ggsave(
    filename = "plots/time_slices.pdf",
    plot = out,
    width = 8.5,
    height = 4
  )

  # return
  out

}
