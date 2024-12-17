
#' @title Detect doublets in a single cytometry sample
#' @description Augments data with simulated doublets, computes
#' nearest neighbors for augmented dataset, identifies doublets as those
#' events with a high share of simulated doublets among nearest neighbors.
#' @param df A data frame containing protein expression data.
#' @param cols Columns to use in analysis.
#' @param cofactor Parameter of arcsinh transformation, applied before
#' computing nearest neighbors. Recommended values are 5 for mass cytometry
#' and 500-1000 for flow cytometry.
#' @param thresh Among the 15 nearest neighbors, how many should be simulated
#' doublets in order for the event to be classified as doublet?
#' @param is_debris Optional, binary array with length matching
#' the number of rows in df. TRUE for debris events, FALSE for everything else.
#' This package includes helper functions to compute this for flow or
#' mass cytometry data.
#' @returns A list with multiple elements, among them the singlet/doublet
#' status of each event.
#' @examples
#' path <- system.file("extdata", "df_mdipa.csv", package="Cleanet")
#' df_mdipa <- read.csv(path, check.names=FALSE)
#' cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
#'           "CD56", "CD294", "CD14", "CD3", "CD20",
#'           "CD66b", "CD38", "HLA-DR", "CD45RA",
#'           "DNA1", "DNA2")
#' cleanet_res <- cleanet(df_mdipa, cols, cofactor=5)
#' @export
cleanet <- function(df, cols, cofactor,
                    thresh=5, is_debris=NULL) {
  message("Starting Cleanet...")

  if (!is.null(is_debris)) {
    if (length(is_debris) != nrow(df))
      stop("Length mismatch between debris array and data frame!")

    idx_keep <- which(!is_debris)
    n_keep <- length(idx_keep)
    perc_keep <- round(100*n_keep/nrow(df),1)
    message(paste0("Doublet simulation using ", n_keep, " events, ", perc_keep, "% of total..."))
  } else {
    idx_keep <- seq(nrow(df))
    message(paste0("Doublet simulation using all ", nrow(df),  " events..."))
  }

  n_doub <- floor(length(idx_keep)/2)
  sel1 <- sample(idx_keep, n_doub)
  sel2 <- sample(idx_keep, n_doub)
  x_aug <- augment_data(df, cols, sel1, sel2, cofactor)

  all_knn <- hnsw_knn(x_aug, k=15, distance= 'l2', n_threads=1, M=48)
  n_nb_doub <- apply(all_knn$idx, 1, function(row) length(which(row<n_doub)))
  status <- rep("Filtered", nrow(df))
  status <- if_else(n_nb_doub[-seq(n_doub)] > thresh, "Doublet", "Singlet")

  n_found <- length(which(status=="Doublet"))
  perc_doublets <- round(100*n_found/nrow(df),1)
  message(paste0("Done! Found ", n_found, " doublets, ", perc_doublets, "% of total."))

  return(list(n_nb_doub = n_nb_doub[-seq(n_doub)],
              status = status,
              nbrs = all_knn$idx[-seq(n_doub),],
              sel1 = sel1,
              sel2 = sel2,
              n_doub = n_doub,
              thresh = thresh,
              n_nb_sim = n_nb_doub[seq(n_doub)],
              sensitivity = length(which(n_nb_doub[seq(n_doub)] > 5))/n_doub,
              is_debris=is_debris))
}


augment_data <- function(df, cols, sel1, sel2, cofactor) {
  x0 <- as.matrix(df[,cols])
  x_doub <- x0[sel1,] + x0[sel2,]
  x_aug <- rbind(x_doub, x0)
  x_aug <- asinh(x_aug/cofactor)
  return(x_aug)
}


#' @title Flag debris in mass cytometry data.
#' @description Detect events with low distance from 0 in protein space.
#' This function aims for high specificity, but not high sensitivity: for
#' Cleanet's purposes, it suffices to deplete debris, even if not all of
#' it is eliminated.
#' @param df A data frame containing protein expression data.
#' @param cols Columns to use in analysis. It is recommended to use the same
#' ones in the call to cleanet.
#' @param cols_plot Two columns that are used for visual feedback.
#' @param cofactor Parameter for arcsinh transformation used before computing
#' distances. 5 is a good default for mass cytometry data.
#' @param threshold Number between 0 and 1; distances are scaled between
#' 0 and 1 and events whose distance to the origin is smaller than the
#' threshold are flagged.
#' @returns A binary array with the same length as the number of rows in df.
#' TRUE for debris, FALSE for everything else.
#' @examples
#' path <- system.file("extdata", "df_mdipa.csv", package="Cleanet")
#' df_mdipa <- read.csv(path, check.names=FALSE)
#' cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
#'           "CD56", "CD294", "CD14", "CD3", "CD20",
#'           "CD66b", "CD38", "HLA-DR", "CD45RA",
#'           "DNA1", "DNA2")
#' is_debris <- filter_debris_cytof(df_mdipa, cols)
#' @export
filter_debris_cytof <- function(df, cols, cols_plot=c("DNA1", "CD45"),
                                cofactor=5, threshold=0.3) {
  x0 <- asinh(as.matrix(df[cols])/cofactor)
  ncol <- length(cols)

  norm <- sqrt(apply(apply(x0, 2, function(col) col/max(col))^2, 1, sum))/sqrt(ncol)

  sub <- sample(nrow(df), min(nrow(df),1e4))
  df_sub <- df[sub,] %>%
    mutate(Debris = norm[sub]<threshold) %>%
    mutate(across(cols_plot, ~asinh(.x/cofactor)))
  p <- ggplot(df_sub, aes(x=.data[[cols_plot[1]]],
                          y=.data[[cols_plot[2]]],
                          color=.data[["Debris"]])) +
    geom_point(size=0.5) +
    theme_bw()
  show(p)

  return(norm < threshold)
}


#' @title Flag debris in flow cytometry data.
#' @description Detect events in the lower left corner of FSC-A/SSC-A plots.
#' This function aims for high specificity, but not high sensitivity: for
#' Cleanet's purposes, it suffices to deplete debris, even if not all of
#' it is eliminated.
#' @param df A data frame containing scattering channels and
#' protein expression data.
#' @param sum_max Numeric; events whose sum of FSC-A and SSC-A is smaller than
#' this value are flagged.
#' @param cols Names of columns to use. This function is intended for use
#' with the area channel of forward and side scatter.
#' @returns A binary array with the same length as the number of rows in df.
#' TRUE for debris, FALSE for everything else.
#' @export
filter_debris_flow <- function(df, sum_max=5e4, cols=c("FSC-A", "SSC-A")) {
  is_debris <- df[[cols[1]]] + df[[cols[2]]] < sum_max &
    df[[cols[1]]] > 0 & df[[cols[2]]] > 0

  pal <- setNames(c("#66A61E", "#666666"), c("FALSE", "TRUE"))

  sub <- sample(nrow(df), min(nrow(df),1e4))
  df_sub <- df[sub,] %>%
    mutate(Debris = is_debris[sub])
  p <- ggplot(df_sub, aes(x=.data[[cols[1]]], y=.data[[cols[2]]], color=.data[["Debris"]])) +
    geom_point(size=0.5) +
    scale_color_manual(values=pal) +
    guides(color=guide_legend(override.aes = list(size=4))) +
    theme_bw(base_size=16) +
    theme(legend.position = "top", legend.box.spacing = margin(0))
  show(p)

  return(is_debris)
}


#' @title Classify doublets (or multiplets) based on component singlets.
#' @description Extends a classification of singlets into a classification
#' of doublets.
#' @param cleanet_res The output of a call to the cleanet function.
#' @param singlet_clas An array giving a classification of the singlets, whose
#' length must match the number of singlet events returned in cleanet_res.
#' @param max_multi The highest cardinality of a multiplet to be considered.
#' @returns An array with the same length as the number of doublets found
#' in cleanet_res, specifying the composition of each doublet.
#' @examples
#' path <- system.file("extdata", "df_mdipa.csv", package="Cleanet")
#' df_mdipa <- read.csv(path, check.names=FALSE)
#' cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
#'           "CD56", "CD294", "CD14", "CD3", "CD20",
#'           "CD66b", "CD38", "HLA-DR", "CD45RA",
#'           "DNA1", "DNA2")
#' cleanet_res <- cleanet(df_mdipa, cols, cofactor=5)
#' singlet_clas <- df_mdipa$label[which(cleanet_res$status!="Doublet")]
#' doublet_clas <- classify_doublets(cleanet_res, singlet_clas)
#' @export
classify_doublets <- function(cleanet_res, singlet_clas, max_multi=4) {

  sel1 <- cleanet_res[["sel1"]]
  sel2 <- cleanet_res[["sel2"]]
  nbrs <- cleanet_res[["nbrs"]]
  n_doub <- cleanet_res[["n_doub"]]
  thresh <- cleanet_res[["thresh"]]
  n_nb_doub <- cleanet_res[["n_nb_doub"]]

  lab_with_doub <- rep("Doublet", nrow(nbrs))
  lab_with_doub[which(cleanet_res$status!="Doublet")] <- singlet_clas

  for (iter in seq(max_multi)) {
    sim_lab <- if_else(lab_with_doub[sel1] < lab_with_doub[sel2],
                       paste0(lab_with_doub[sel1], "_", lab_with_doub[sel2]),
                       paste0(lab_with_doub[sel2], "_", lab_with_doub[sel1]))
    sim_lab_fac <- as.factor(sim_lab)

    doub_c <- doub_lab_c(nbrs, sim_lab_fac, n_doub, length(levels(sim_lab_fac)))
    doub_lab <- c("",levels(sim_lab_fac))[doub_c+1]

    to_replace <- which(lab_with_doub=="Doublet" & !grepl("Doublet", doub_lab))
    lab_with_doub[to_replace] <- doub_lab[to_replace]
  }

  return(lab_with_doub[which(n_nb_doub>thresh)])
}


#' @title Tabulate expected and observed proportions of doublet types.
#' @description Given compatible classifications of singlets and doublets,
#' this function computes expected proportions of doublets as the product
#' of the proportions of their components.
#' @param cleanet_res The output of a call to the cleanet function.
#' @param singlet_clas An array giving a classification of the singlets, whose
#' length must match the number of singlet events returned in cleanet_res.
#' @param doublet_clas An array giving a classification of the doublets, whose
#' length must match the number of doublet events returned in cleanet_res.
#' @returns A data frame tabulating expected and observed proportions for
#' each unique doublet type.
#' @examples
#' path <- system.file("extdata", "df_mdipa.csv", package="Cleanet")
#' df_mdipa <- read.csv(path, check.names=FALSE)
#' cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
#'           "CD56", "CD294", "CD14", "CD3", "CD20",
#'           "CD66b", "CD38", "HLA-DR", "CD45RA",
#'           "DNA1", "DNA2")
#' cleanet_res <- cleanet(df_mdipa, cols, cofactor=5)
#' singlet_clas <- df_mdipa$label[which(cleanet_res$status!="Doublet")]
#' doublet_clas <- classify_doublets(cleanet_res, singlet_clas)
#' df_exp_obs <- compare_doublets_exp_obs(doublet_clas, singlet_clas, cleanet_res)
#' @export
compare_doublets_exp_obs <- function(doublet_clas, singlet_clas, cleanet_res) {
  if (is.null(cleanet_res$is_debris)) {
    tab_sing <- table(singlet_clas)
  } else {
    tab_sing <- table(singlet_clas[which(!cleanet_res$is_debris[which(cleanet_res$status!="Doublet")])])
  }

  out <- outer(tab_sing, tab_sing) / (sum(tab_sing))^2
  df_exp <- melt(out)
  names(df_exp) <- c("Var1", "Var2", "value")

  df_exp <- df_exp %>%
    as_tibble() %>%
    filter(as.character(.data[["Var1"]]) <= as.character(.data[["Var2"]])) %>%
    mutate(value = if_else(.data[["Var1"]] == .data[["Var2"]],
                           .data[["value"]],
                           2*.data[["value"]])) %>%
    mutate(Type = if_else(as.character(.data[["Var1"]]) < as.character(.data[["Var2"]]),
                          paste0(.data[["Var1"]], "_", .data[["Var2"]]),
                          paste0(.data[["Var2"]], "_", .data[["Var1"]])),
           Expected = .data[["value"]],
           .keep="unused")

  tab_doub <- table(doublet_clas)
  df_obs <- tibble(Type = names(tab_doub),
                   Observed = as.numeric(tab_doub) / length(doublet_clas))

  df_exp_obs <- full_join(df_exp, df_obs) %>%
    replace_na(list("Observed"=0, "Expected"=0))
  return(df_exp_obs)
}


