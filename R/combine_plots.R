#' Combine and save plots of results
#'
#' @param plot_Global Global plot of results
#' @param plot_Africa Plot of results in Africa only
#' @param plot_Asia Plot of results in Asia only
#' @param plot_Europe Plot of results in Europe only
#' @param plot_North.America Plot of results in North America only
#' @param plot_Oceania Plot of results in Oceania only
#' @param plot_South.America Plot of results in South America only
#'
#' @returns A patchwork composition of ggplot objects
#'
combine_plots <- function(plot_Global, plot_Africa, plot_Asia, plot_Europe,
                          plot_North.America, plot_Oceania,
                          plot_South.America) {
  # combine plots
  out <-
    (
      plot_Global + plot_Africa + plot_Asia + plot_Europe +
        plot_North.America + plot_Oceania + plot_South.America +
        plot_spacer()
    ) +
    plot_layout(
      nrow = 2,
      ncol = 4,
      axes = "collect_y",
      axis_titles = "collect"
    )
  # save
  ggsave(
    filename = "plots/results.pdf",
    plot = out,
    width = 7,
    height = 5
  )
  # return
  out
}
