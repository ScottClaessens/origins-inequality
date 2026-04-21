#' Plot single tree with ancestral state reconstructions and trait data
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_ids Indexes for trees
#' @param ancestral_states Tibble of ancestral states from the model
#'
#' @returns A ggplot object
#'
plot_tree <- function(data, tree, tree_ids, ancestral_states) {
  # use subset of trees
  tree <- tree[tree_ids]
  # wrangle data
  d <-
    data |>
    transmute(var = as.factor(as.numeric(class_differentiation))) |>
    as.data.frame()
  rownames(d) <- data$xd_id
  # summarise ancestral states
  dd <-
    ancestral_states |>
    # first tree only
    filter(tree == 1) |>
    # probability of categories 3, 4, or 5
    rowwise() |>
    mutate(prob_inequality = list(prob3 + prob4 + prob5)) |>
    dplyr::select(node, prob_inequality) |>
    unnest(c(prob_inequality)) |>
    group_by(node) |>
    # get posterior median
    summarise(
      prob_inequality = median(prob_inequality),
      .groups = "drop"
    )
  # plot tree
  out <-
    ggtree(
      tr = tree[[1]],
      layout = "circular",
      mapping = aes(colour = prob_inequality),
      linewidth = 0.2
    ) %<+% dd +
    scale_colour_gradient(
      name = "Probability of\nstratification",
      limits = c(0, 1),
      low = "grey85",
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
      labels = function(x) levels(data$class_differentiation)[as.numeric(x)],
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
    "Afro-Asiatic"            = c("xd305",  "xd526"),
    "Indo-European"           = c("xd528",  "xd604"),
    "Dravidian"               = c("xd668",  "xd680"),
    "Uralic"                  = c("xd544",  "xd632"),
    "Nilotic"                 = c("xd2",    "xd406"),
    "Austronesian"            = c("xd687",  "xd749"),
    "Sino-Tibetan"            = c("xd639",  "xd704"),
    "Austroasiatic"           = c("xd661",  "xd727"),
    "Salishan"                = c("xd1071", "xd1146"),
    "Eskimo-Aleut"            = c("xd1022", "xd1067")
  )
  # add clade labels for language families
  for (family in names(taxa_bookends)) {
    # node number for most recent common ancestor
    node <- getMRCA(tree[[1]], taxa_bookends[[family]])
    # add clade label
    out <-
      out +
      geom_cladelab(
        node = node,
        label = family,
        offset = 10,
        offset.text = 2,
        barsize = 0.2,
        fontsize = 2.5,
        hjust = ifelse(
          family %in% c("Afro-Asiatic", "Indo-European", "Uralic", "Dravidian",
                        "Austronesian", "Austroasiatic", "Sino-Tibetan"),
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
  # return
  out
}
