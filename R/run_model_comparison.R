#' Run comparison of unfossilised and fossilised models
#'
#' Fit and compare two Multistate models, one without any fossilised nodes and
#' another where all nodes prior to 11.7kya are set to egalitarian.
#'
#' @param data Tibble of D-PLACE data
#' @param mcc_tree Maximum clade credibility tree object of class phylo
#'
#' @returns Numeric. Log Bayes factor for model comparison.
#'
run_model_comparison <- function(data, mcc_tree) {
  # get log marginal likelihoods for both models
  log_lik_fossilised   <- get_log_lik(data, mcc_tree, fossilised = TRUE)
  log_lik_unfossilised <- get_log_lik(data, mcc_tree, fossilised = FALSE)
  # return log bayes factor in favour of fossilised model
  2 * (log_lik_fossilised - log_lik_unfossilised)
}

# function to fit model and return log likelihood
get_log_lik <- function(data, mcc_tree, fossilised = TRUE) {
  # get file names for data, tree, and commands
  suffix <- ifelse(fossilised, "fossilised", "unfossilised")
  data_file <- paste0("data_", suffix, ".txt")
  tree_file <- paste0("tree_", suffix, ".txt")
  commands_file <- paste0("commands_", suffix, ".txt")
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
  mcc_tree$node.label <- NULL
  write.nexus(
    phy = mcc_tree,
    file = paste0("bayestraits/", tree_file)
  )
  # create command file in /bayestraits directory
  commands <-
    c(
      "1",                       # multistate
      "2",                       # mcmc
      "HyperPriorAll exp 0 10",  # priors
      "Iterations 500000",       # number of iterations
      "Stones 100 10000",        # stepping stone sampler for model comparison
      "Seed 123"                 # seed
    )
  # declare internal nodes in the tree
  internal_nodes <- (Ntip(mcc_tree) + 1):(Ntip(mcc_tree) + Nnode(mcc_tree))
  # get timings of nodes
  times <- ape::node.depth.edgelength(mcc_tree)
  times <- times - max(times)
  # loop over internal nodes
  for (i in internal_nodes) {
    # if internal node is older than 11.7 kya, fossilise
    if (fossilised & times[i] <= -11.7) {
      # get all descendant taxa from current node
      taxa <- extract.clade(mcc_tree, node = i)$tip.label
      # get tag command
      tag_command <-
        paste(
          "AddTag",
          paste0("T", i),
          paste(taxa, collapse = " ")
        )
      # get fossil command
      fossil_command <-
        paste(
          "Fossil",
          paste0("N", i),
          paste0("T", i),
          "12"
        )
      # add to commands
      commands <- c(commands, tag_command, fossil_command)
    }
  }
  commands <- c(commands, "Run")
  # write command file
  writeLines(
    text = commands,
    con = paste0("bayestraits/", commands_file)
  )
  # run bayes traits in command line
  dir <- paste0(here::here(), "/bayestraits/")
  invisible(
    system(
      paste(
        paste0(dir, "BayesTraitsV3"),
        paste0(dir, tree_file),
        paste0(dir, data_file),
        paste0("< ", dir, commands_file)
      ),
      intern = TRUE
    )
  )
  # read log marginal likelihood
  stones <- readLines(paste0("bayestraits/", data_file, ".Stones.txt"))
  log_marginal_likelihood <- parse_number(stones[length(stones)])
  # cleanup
  invisible(
    file.remove(
      c(
        paste0("bayestraits/", commands_file),
        paste0("bayestraits/", data_file),
        paste0("bayestraits/", data_file, ".Log.txt"),
        paste0("bayestraits/", data_file, ".Schedule.txt"),
        paste0("bayestraits/", data_file, ".Stones.txt"),
        paste0("bayestraits/", tree_file)
      )
    )
  )
  # return log marginal likelihood
  log_marginal_likelihood
}
