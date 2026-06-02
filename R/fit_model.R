#' Fit the MultiState model using BayesTraits
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param chain Index of the independent chain for the MCMC sampler. Used as a
#'   random seed for reproducibility.
#' @param model Character of length 1. Which model to fit. One of: "full",
#'   "rectilinear", "unilinear", "relaxed_unilinear", "alternative", or
#'   "alternative_reversible"
#' @param iter Numeric. Number of MCMC sampling iterations.
#' @param burnin Numeric. Number of MCMC burn-in iterations.
#' @param stones Logical. If \code{TRUE} (default), include stepping stone
#'   sampler.
#' @param asr Logical. If \code{FALSE} (default), the model does not estimate
#'   ancestral states. If \code{TRUE}, the model estimates ancestral states for
#'   all internal nodes of the tree.
#' @param tree_id Numeric. Index for tree. If \code{NULL} (default), the full
#'   posterior treeset is used. Must be set if \code{asr = TRUE}.
#'
#' @returns A tibble of posterior samples
#'
fit_model <- function(data, tree, chain, model = "full", iter = 550000,
                      burnin = 50000, stones = TRUE, asr = FALSE,
                      tree_id = NULL) {

  # tree_id must be set if asr is true
  if (asr & is.null(tree_id)) {
    stop(
      "Tree ID must be set when estimating ancestral states.",
      call. = FALSE
    )
  }

  # subset to tree
  if (is.null(tree_id)) {
    tree_id <- 0
  } else {
    tree <- tree[[tree_id]]
  }

  # get file names for data, tree, and commands
  data_file     <- paste0("data_",     model, "_", tree_id, "_", chain, ".txt")
  tree_file     <- paste0("tree_",     model, "_", tree_id, "_", chain, ".txt")
  commands_file <- paste0("commands_", model, "_", tree_id, "_", chain, ".txt")

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
  withr::with_options(list(scipen = 999), {
    commands <-
      c(
        "1",                       # multistate
        "2",                       # mcmc
        "HyperPriorAll exp 0 10",  # priors
        paste("Burnin", burnin),   # number of burnin iterations
        paste("Iterations", iter), # number of sampling iterations
        paste("Seed", chain)       # seed
      )
  })

  # include stepping stone sampler?
  if (stones) {
    commands <- c(commands, "Stones 100 10000")
  }

  # get restricted parameters for different models
  restricted_pars <- list(
    "rectilinear"            = c("q13", "q21", "q31", "q32"),
    "unilinear"              = c("q13", "q31"),
    "relaxed_unilinear"      = c("q13"),
    "alternative"            = c("q21", "q23", "q31", "q32"),
    "alternative_reversible" = c("q23", "q32")
  )

  # add command for model restrictions? no restrictions for full model
  if (model != "full") {
    commands <- c(
      commands,
      paste("Res", paste(restricted_pars[[model]], collapse = " "), "0")
    )
  }

  # estimate ancestral states?
  if (asr) {

    # declare all internal nodes in the tree
    internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))

    # loop over internal nodes
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

  }

  # add run command
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

  # on which row does the log output start?
  row_start <- which(
    str_starts(
      readLines(paste0("bayestraits/", data_file, ".Log.txt")),
      "Iteration\t"
    )
  )

  # read log file
  model_log <-
    read_tsv(
      file = paste0("bayestraits/", data_file, ".Log.txt"),
      skip = row_start - 1,
      show_col_types = FALSE
    ) |>
    # remove final empty column
    dplyr::select(-last_col()) |>
    # can save storage size by removing one of the redundant probabilities
    # it is possible to calculate with 1 - sum(other probabilities)
    dplyr::select(-ends_with("P(1)"))

  # add log marginal likelihood to model log?
  if (stones) {
    stones_output <- readLines(paste0("bayestraits/", data_file, ".Stones.txt"))
    model_log$log_lik <- parse_number(stones_output[length(stones_output)])
  }

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

  if (stones) {
    invisible(file.remove(paste0("bayestraits/", data_file, ".Stones.txt")))
  }

  # return
  tibble(
    tree_id = rep(tree_id, times = nrow(model_log)),
    model   = rep(model,   times = nrow(model_log)),
    chain   = rep(chain,   times = nrow(model_log))
  ) |>
    bind_cols(model_log)

}
