#' Produce table of log Bayes Factors for model comparison
#'
#' @param log_lik_1 Log marginal likelihood for model with internal nodes prior
#'   to 11.7 kya fossilised to "1"
#' @param log_lik_2 Log marginal likelihood for model with internal nodes prior
#'   to 11.7 kya fossilised to "2"
#' @param log_lik_3 Log marginal likelihood for model with internal nodes prior
#'   to 11.7 kya fossilised to "3"
#'
#' @returns Tibble
#'
get_table_log_bayes_factors <- function(log_lik_1, log_lik_2, log_lik_3) {

  tibble(
    Model = as.character(1:3),
    `Fossilised state` = c("Egalitarian", "Wealth distinctions",
                           "Stratification"),
    `Log likelihood` = c(log_lik_1, log_lik_2, log_lik_3),
    `Log BF vs. Model 1` = c(
      2 * (log_lik_1 - log_lik_1),
      2 * (log_lik_2 - log_lik_1),
      2 * (log_lik_3 - log_lik_1)
    ),
    `Log BF vs. Model 2` = c(
      2 * (log_lik_1 - log_lik_2),
      2 * (log_lik_2 - log_lik_2),
      2 * (log_lik_3 - log_lik_2)
    ),
    `Log BF vs. Model 3` = c(
      2 * (log_lik_1 - log_lik_3),
      2 * (log_lik_2 - log_lik_3),
      2 * (log_lik_3 - log_lik_3)
    )
  ) |>
  # round and convert to character
  mutate(
    across(
      where(is.numeric), function(x) {
        as.character(format(round(x, 2), nsmall = 2))
      }
    )
  ) |>
  # replace zeros with dashes
  mutate(across(everything(), function(x) ifelse(str_ends(x, "0.00"), "-", x)))

}
