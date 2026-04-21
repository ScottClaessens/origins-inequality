options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)

tar_option_set(
  packages = c("ape", "coevolve", "ggtree", "patchwork",
               "phangorn", "tidyverse")
)
tar_source()

# pipeline
list(
  # get data urls
  tar_target(
    dplace_data_url,
    paste0(
      "https://raw.githubusercontent.com/D-PLACE/dplace-cldf/",
      "6c2008c187a297d1955b41d8ae80d8e31d404f6c/cldf/data.csv"
    ),
    format = "url"
  ),
  tar_target(
    dplace_societies_url,
    paste0(
      "https://raw.githubusercontent.com/D-PLACE/dplace-cldf/",
      "6c2008c187a297d1955b41d8ae80d8e31d404f6c/cldf/societies.csv"
    ),
    format = "url"
  ),
  tar_target(
    glottolog_languages_url,
    paste0(
      "https://raw.githubusercontent.com/glottolog/glottolog-cldf/",
      "072ca0d0410039fb8b779be8fc165bac575d2cda/cldf/languages.csv"
    ),
    format = "url"
  ),
  # get data file paths
  tar_target(tree_file, "data/tree/dplace.nxs", format = "file"),
  # load tree
  tar_target(tree, read.nexus(tree_file)),
  # compute maximum clade credibility tree
  tar_target(mcc_tree, phangorn::mcc(tree)),
  # load dplace data
  tar_target(
    data,
    load_dplace_data(
      dplace_data_url, dplace_societies_url,
      glottolog_languages_url, mcc_tree
    )
  ),
  # get random tree ids
  tar_target(tree_ids, sample(1:length(tree), size = 20)),
  # fit ancestral state reconstruction model
  tar_target(fit, fit_model(data, tree, tree_ids)),
  # get ancestral states
  tar_target(ancestral_states, get_ancestral_states(tree, fit, tree_ids)),
  # plot results globally and by language family
  tar_target(plot_Global, plot_model(data, ancestral_states, tree, tree_ids)),
  tar_map(
    values = tibble(
      family = c("Atlantic-Congo", "Austronesian", "Afro-Asiatic",
                 "Uto-Aztecan", "Indo-European", "Nilotic", "Algic",
                 "Athabaskan-Eyak-Tlingit", "Sino-Tibetan", "Mande", "Salishan",
                 "Uralic", "Eskimo-Aleut", "Austroasiatic", "Dravidian")
    ),
    tar_target(plot, plot_model(data, ancestral_states, tree, tree_ids,
                                family = family))
  ),
  # combine plots
  tar_target(
    combined_plots,
    combine_plots(
      list(
        plot_Global, plot_Atlantic.Congo, plot_Austronesian, plot_Afro.Asiatic,
        plot_Uto.Aztecan, plot_Indo.European, plot_Nilotic, plot_Algic,
        plot_Athabaskan.Eyak.Tlingit, plot_Sino.Tibetan, plot_Mande,
        plot_Salishan, plot_Uralic, plot_Eskimo.Aleut, plot_Austroasiatic,
        plot_Dravidian
      )
    )
  ),
  # print session info
  tar_target(
    sessionInfo,
    writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  )
)
