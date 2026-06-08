options(tidyverse.quiet = TRUE)
library(crew)
library(targets)
library(tarchetypes)
library(tidyverse)

tar_option_set(
  packages = c("ape", "deeptime", "ggtree", "patchwork", "phangorn",
               "phytools", "rstan", "tidyverse", "withr"),
  controller = crew_controller_local(workers = 8),
  deployment = "main"
)
tar_source()

# pipeline
list(

  # ─────────────────────────────────────────
  # Load D-PLACE data and phylogeny
  # ─────────────────────────────────────────

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

  # get tree file path
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

  # ─────────────────────────────────────────
  # Compare models of evolution
  # ─────────────────────────────────────────

  # get independent mcmc chains
  tar_target(chain, 1:4),

  # loop over models
  tar_map(

    values = tibble(
      model = c(
        "full", "rectilinear", "unilinear", "relaxed_unilinear",
        "alternative", "alternative_reversible"
      )
    ),

    # fit model
    tar_target(
      fit,
      fit_model(data, tree, chain, model),
      pattern = map(chain),
      deployment = "worker",
      storage = "worker",
      retrieval = "worker"
    ),

    # get diagnostics
    tar_target(diagnostics, calculate_model_diagnostics(fit)),

    # plot MCMC trace
    tar_target(plot_trace, plot_mcmc_trace(fit, model))

  ),

  # model comparison table
  tar_target(
    table_model_comparison,
    get_table_model_comparison(
      bind_rows(
        fit_full, fit_rectilinear, fit_unilinear, fit_relaxed_unilinear,
        fit_alternative, fit_alternative_reversible
      )
    )
  ),

  # ─────────────────────────────────────────
  # Estimate ancestral states
  # ─────────────────────────────────────────

  # get tree ids
  tar_target(tree_id, sample(1:length(tree), size = 100, replace = FALSE)),

  # fit ancestral state reconstruction model
  tar_target(
    fit_asr,
    fit_model(data, tree, chain, model = "relaxed_unilinear",
              stones = FALSE, asr = TRUE, tree_id = tree_id),
    pattern = cross(chain, tree_id),
    # run in parallel
    deployment = "worker",
    storage = "worker",
    retrieval = "worker"
  ),

  # calculate model diagnostics
  tar_target(diagnostics_asr, calculate_model_diagnostics(fit_asr)),

  # plot tree
  tar_target(
    plot_tree_states,
    plot_tree(data, tree, tree_id, fit_asr)
  ),

  # plot results globally and by language family
  tar_target(plot_Global, plot_model(data, fit_asr, tree, tree_id,
                                     end_time = -0.2, time_slice = 0.2)),
  tar_map(
    values = tibble(
      family = c(
        "Atlantic-Congo", "Austronesian", "Afro-Asiatic", "Uto-Aztecan",
        "Indo-European", "Nilotic", "Algic", "Athabaskan-Eyak-Tlingit",
        "Sino-Tibetan", "Mande", "Salishan", "Uralic", "Eskimo-Aleut",
        "Austroasiatic", "Arawakan", "Cariban", "Central Sudanic", "Turkic",
        "Cochimi-Yuman", "Dravidian", "Nuclear Trans New Guinea", "Siouan",
        "Tupian", "Mayan"
      )
    ),
    tar_target(plot, plot_model(data, fit_asr, tree, tree_id,
                                family = family,
                                end_time = -0.2, time_slice = 0.2))
  ),

  # combine plots
  tar_target(
    combined_plots,
    combine_plots(
      list(
        plot_Global, plot_Atlantic.Congo, plot_Afro.Asiatic, plot_Nilotic,
        plot_Mande, plot_Central.Sudanic, plot_Indo.European, plot_Sino.Tibetan,
        plot_Uralic, plot_Austroasiatic, plot_Turkic, plot_Dravidian,
        plot_Uto.Aztecan, plot_Algic, plot_Athabaskan.Eyak.Tlingit,
        plot_Salishan, plot_Eskimo.Aleut, plot_Cochimi.Yuman, plot_Siouan,
        plot_Mayan, plot_Austronesian, plot_Nuclear.Trans.New.Guinea,
        plot_Arawakan, plot_Cariban, plot_Tupian
      )
    )
  ),

  # ─────────────────────────────────────────
  # Fossilising nodes
  # ─────────────────────────────────────────

  # run comparison of different fossilisations
  tar_target(log_lik_1, get_log_lik_fossilised(data, mcc_tree, "1")),
  tar_target(log_lik_2, get_log_lik_fossilised(data, mcc_tree, "2")),
  tar_target(log_lik_3, get_log_lik_fossilised(data, mcc_tree, "3")),

  # summarise log bayes factors
  tar_target(
    table_fossilised_log_bfs,
    get_table_fossilised(log_lik_1, log_lik_2, log_lik_3)
  ),

  # simulate model comparison
  tar_target(sim_1, simulate_model_comparison(data, mcc_tree, fossil = "1")),
  tar_target(sim_2, simulate_model_comparison(data, mcc_tree, fossil = "2")),
  tar_target(sim_3, simulate_model_comparison(data, mcc_tree, fossil = "3")),

  # ─────────────────────────────────────────
  # Produce manuscript
  # ─────────────────────────────────────────

  # knit quarto manuscript
  tar_quarto(manuscript, "quarto/manuscript.qmd", quiet = FALSE),

  # write session info
  tar_target(
    sessionInfo,
    writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  )

)
