#' Plot single tree with ancestral state reconstructions and trait data
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_id Indexes for trees
#' @param fit Results of fitted model
#'
#' @returns A ggplot object
#'
plot_tree <- function(data, tree, tree_id, fit) {
  # use first tree only
  id <- tree_id[1]
  tree <- tree[[id]]
  # wrangle data
  d <-
    data |>
    transmute(
      var = as.numeric(class_differentiation),
      var = ifelse(var %in% 3:5, 3, var),
      var = as.factor(var)
    ) |>
    as.data.frame()
  rownames(d) <- data$xd_id
  # summarise ancestral states
  dd <-
    fit |>
    filter(tree_id == id) |>
    # get probability of inequality at internal nodes
    dplyr::select(ends_with("P(3)") & !starts_with("Root")) |>
    pivot_longer(
      cols = everything(),
      names_to = "node",
      values_to = "prob_inequality"
    ) |>
    mutate(node = parse_number(str_sub(node, 2, 5))) |>
    # get posterior median for each node
    group_by(node) |>
    summarise(
      prob_inequality = median(prob_inequality),
      .groups = "drop"
    ) |>
    # add data for tips
    bind_rows(
      tibble(
        node = match(data$xd_id, tree$tip.label),
        prob_inequality =
          ifelse(as.numeric(data$class_differentiation) %in% 3:5, 1, 0)
      )
    )
  # plot tree
  out <-
    ggtree(
      tr = tree,
      layout = "circular",
      colour = "grey",
      linewidth = 0.15
    ) %<+% dd +
    geom_nodepoint(
      mapping = aes(colour = prob_inequality),
      size = 0.05
    ) +
    scale_colour_gradient(
      name = "Probability of\nstratification",
      limits = c(0, 1),
      low = "grey95",
      high = "black"
    )
  # add trait data
  out <-
    gheatmap(
      p = out,
      data = d,
      offset = -0.5,
      width = 0.05,
      colnames = FALSE,
      color = NA
    ) +
    scale_fill_brewer(
      name = "EA066",
      labels = function(x) {
        labs <- c("Absence of distinctions",
                  "Wealth distinctions",
                  "Stratification")
        labs[as.numeric(x)]
      },
      type = "seq",
      palette = 7,
      na.value = "grey95"
    )
  # get taxa bookends for language families
  taxa_bookends <- list(
    "Atlantic-Congo"          = c("xd10",   "xd236"),
    "Mande"                   = c("xd188",  "xd292"),
    #"Athabaskan-Eyak-Tlingit" = c("xd1026", "xd1082"),
    "Algic"                   = c("xd1047", "xd1097"),
    "Uto-Aztecan"             = c("xd1118", "xd1292"),
    "Afro-Asiatic"            = c("xd305",  "xd590"),
    "Indo-European"           = c("xd528",  "xd604"),
    "Dravidian"               = c("xd668",  "xd680"),
    "Uralic"                  = c("xd544",  "xd632"),
    "Nilotic"                 = c("xd2",    "xd391"),
    "Austronesian"            = c("xd687",  "xd732"),
    "Sino-Tibetan"            = c("xd639",  "xd658"),
    "Austroasiatic"           = c("xd661",  "xd724"),
    "Salishan"                = c("xd1071", "xd1106"),
    "Eskimo-Aleut"            = c("xd1022", "xd1067")
  )
  # add clade labels for language families
  for (family in names(taxa_bookends)) {
    # node number for most recent common ancestor
    node <- getMRCA(tree, taxa_bookends[[family]])
    # add clade label
    out <-
      out +
      geom_cladelab(
        node = node,
        label = family,
        offset = 15,
        offset.text = 3,
        barsize = 0.2,
        fontsize = 2.5,
        hjust = ifelse(
          family %in% c("Afro-Asiatic", "Indo-European", "Uralic",
                        "Dravidian", "Austronesian", "Austroasiatic",
                        "Sino-Tibetan", "Nilotic"),
          1, 0
        )
      )
  }
  # save
  ggsave(
    filename = "plots/tree.pdf",
    plot = out,
    height = 10,
    width = 10
  )
  # cleanup
  rm(data, tree, tree_id, fit, id, d, dd, taxa_bookends, node)
  # return
  out
}
