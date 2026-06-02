#' Produce table of log Bayes Factors for model comparison
#'
#' @param fits Model outputs combined with \code{bind_rows()}
#'
#' @returns Tibble
#'
get_table_model_comparison <- function(fits) {

  # get table of log marginal likelihoods
  log_liks <-
    fits |>
    filter(chain == 1) |>
    group_by(model) |>
    summarise(log_lik = unique(log_lik)) |>
    arrange(desc(log_lik))

  # calculate log bayes factors compared to full model
  log_liks$log_bf <-
    2 * (log_liks$log_lik - log_liks$log_lik[log_liks$model == "full"])

  # return table
  log_liks |>
    transmute(
      Model = str_to_sentence(str_replace(model, "_", " ")),
      Rank = 1:n(),
      `Log marginal likelihood` = format(round(log_lik, 2), nsmall = 2),
      `Log BF vs. full model` = ifelse(
        log_bf == 0,
        " ",
        format(round(log_bf, 2), trim = TRUE, nsmall = 2)
      )
    )

}
