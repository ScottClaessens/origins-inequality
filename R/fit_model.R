#' Fit ancestral state reconstruction model
#'
#' Fit the ancestral state reconstruction model using BayesTraits
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_id Index for tree
#' @param chain Index of the independent chain for the MCMC sampler. Used as a
#'   random seed for reproducibility.
#'
#' @returns A tibble of posterior samples
#'
fit_model <- function(data, tree, tree_id, chain) {
  # subset to tree
  tree <- tree[[tree_id]]
  # get file names for data, tree, and commands
  data_file <- paste0("data_", tree_id, "_", chain, ".txt")
  tree_file <- paste0("tree_", tree_id, "_", chain, ".txt")
  commands_file <- paste0("commands_", tree_id, "_", chain, ".txt")
  # create tab-separated data file in /bayestraits directory
  data |>
    transmute(
      xd_id = xd_id,
      trait = as.character(class_differentiation),
      trait = ifelse(trait == "Absence of distinctions", "1", trait),
      trait = ifelse(trait == "Wealth distinctions",     "2", trait),
      trait = ifelse(str_ends(trait, "stratification"),  "3", trait),
      trait = ifelse(is.na(trait), "-", trait)
    ) |>
    write_tsv(
      file = paste0("bayestraits/", data_file),
      col_names = FALSE
    )
  # create tree file in /bayestraits directory
  write.nexus(
    phy = tree,
    file = paste0("bayestraits/", tree_file)
  )
  # create command file in /bayestraits directory
  commands <-
    c(
      "1",                       # multistate
      "2",                       # mcmc
      "HyperPriorAll exp 0 10",  # priors
      "Iterations 200000",       # number of iterations
      paste("Seed", chain)       # seed
    )
  # declare ancestral states for all internal nodes in the tree
  internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  for (i in internal_nodes) {
    # get all descendant taxa from current node
    taxa <- extract.clade(tree, node = i)$tip.label
    # get tag command
    tag_command <-
      paste(
        "AddTag",
        paste0("T", i),
        paste(taxa, collapse = " ")
      )
    # get node command
    node_command <-
      paste(
        "AddNode",
        paste0("N", i),
        paste0("T", i)
      )
    # add to commands
    commands <- c(commands, tag_command, node_command)
  }
  commands <- c(commands, "Run")
  # write command file
  writeLines(
    text = commands,
    con = paste0("bayestraits/", commands_file)
  )
  # run bayes traits in command line
  dir <- paste0(here::here(), "/bayestraits/")
  system(
    paste(
      paste0(dir, "BayesTraitsV5"),
      paste0(dir, tree_file),
      paste0(dir, data_file),
      paste0("< ", dir, commands_file)
    )
  )
  # read log file
  model_log <-
    read_tsv(
      file = paste0("bayestraits/", data_file, ".Log.txt"),
      skip = 2 * length(internal_nodes) + 50,
      show_col_types = FALSE
    ) |>
    # remove final empty column
    dplyr::select(-last_col()) |>
    # can save storage size by removing one of the redundant probabilities
    # it is possible to calculate with 1 - sum(other probabilities)
    dplyr::select(-ends_with("P(1)"))
  # cleanup
  invisible(
    file.remove(
      c(
        paste0("bayestraits/", commands_file),
        paste0("bayestraits/", data_file),
        paste0("bayestraits/", data_file, ".Log.txt"),
        paste0("bayestraits/", data_file, ".Schedule.txt"),
        paste0("bayestraits/", tree_file)
      )
    )
  )
  # return
  tibble(
    tree_id = rep(tree_id, times = nrow(model_log)),
    chain = rep(chain, times = nrow(model_log))
  ) |>
    bind_cols(model_log)
}
