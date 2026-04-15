#' Load and wrangle data from D-PLACE
#'
#' Load and wrangle data on inequality from the Ethnographic Atlas via the
#' online database D-PLACE (release v3.3.0). The dataset is filtered to 1258
#' societies that can be linked to the phylogenetic tree.
#'
#' @details The dataset produced by this function is a tibble with 1258
#'   observations and 10 variables:
#' \describe{
#'  \item{soc_id}{Character, society ID}
#'  \item{xd_id}{Character, cross-dataset ID (see
#'    https://d-place.org/glossary#q9)}
#'  \item{society}{Name of the society}
#'  \item{glottocode}{Glottocode for the language or dialect of the society}
#'  \item{language_family}{Language family for the language or dialect of the
#'    society}
#'  \item{region}{Region of the society}
#'  \item{focal_year}{Principal year to which data refer}
#'  \item{latitude}{Latitude of the society}
#'  \item{longitude}{Longitude of the society}
#'  \item{class_differentiation}{Ordered factor (five levels), extent of class
#'    differentiation; coded from EA066}
#' }
#'
#' @param dplace_data_url URL to access cldf/data.csv from D-PLACE v3.3.0
#' @param dplace_societies_url URL to access cldf/societies.csv from D-PLACE
#'   v3.3.0
#' @param glottolog_languages_url URL to access cldf/languages.csv from
#'   Glottolog v5.3
#' @param mcc_tree Maximum clade credibility tree of D-PLACE societies used to
#'   filter the dataset
#'
#' @returns A tibble
#'
load_dplace_data <- function(dplace_data_url, dplace_societies_url,
                             glottolog_languages_url, mcc_tree) {
  # load csv files
  data <- read.csv(file = dplace_data_url)
  societies <- read.csv(file = dplace_societies_url)
  languages <- read.csv(file = glottolog_languages_url)
  # ordered levels
  levels_EA066 <- c("Absence of distinctions", "Wealth distinctions",
                    "Elite stratification", "Dual stratification",
                    "Complex stratification")
  # wrangle ethnographic atlas data
  data |>
    # filter to ethnographic atlas data only
    left_join(societies, by = c("Soc_ID" = "ID")) |>
    filter(Contribution_ID == "dplace-dataset-ea") |>
    # pivot wider
    pivot_wider(
      id_cols = c(Soc_ID, xd_id, Name, Latitude, Longitude, region,
                  Glottocode, main_focal_year),
      names_from = Var_ID,
      values_from = Value
    ) |>
    # match duplicates to tree
    mutate(
      xd_id = ifelse(Soc_ID == "Nd53", "xd1189a", xd_id),
      xd_id = ifelse(Soc_ID == "Nd55", "xd1189b", xd_id)
    ) |>
    # retain variables
    transmute(
      soc_id                = Soc_ID,
      xd_id                 = xd_id,
      society               = Name,
      glottocode            = Glottocode,
      region                = region,
      focal_year            = main_focal_year,
      latitude              = Latitude,
      longitude             = Longitude,
      class_differentiation = ordered(EA066, levels = levels_EA066)
    ) |>
    # filter to societies in phylogenetic tree (n = 1258)
    filter(xd_id %in% mcc_tree$tip.label) |>
    # add data on language family affiliation
    left_join(
      dplyr::select(languages, Glottocode, Family_ID),
      by = c("glottocode" = "Glottocode")
    ) |>
    mutate(Family_ID = ifelse(Family_ID == "", NA, Family_ID)) |>
    left_join(
      dplyr::select(languages, Glottocode, Name),
      by = c("Family_ID" = "Glottocode")
    ) |>
    dplyr::select(soc_id:glottocode, Name, region:class_differentiation) |>
    rename(language_family = Name)
}
