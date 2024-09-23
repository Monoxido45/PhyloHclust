#' @title Phylogenetical Feature Importance Score for Hierarchical clustering methods
#' @description
#' Phylogenetical Feature Importance Score of each feature with respect to each hierarchical clustering method
#' and distance combination.
#' @param data Data frame used to build hierarchical clustering
#' @param cl.list Named list of hierarchical clustering methods: each entry is composed by the name of the hierarchical
#' clustering algorithm (e.g. hclust) followed by a character vector of all it's underlying agglomeration methods selected
#' for comparison (e.g. "complete", "ward.D2").
#' @param dists Chosen distance (or list of distances) to compute dissimilarity matrix for each hierarchical clustering method.
#' Default is "euclidean".
#' @param mixed_dist Chosen mixed distance (or list of mixed distances) to compute dissimilarity matrix for each hierarchical
#' clustering method. Mixed distances are preferable to use when there are categorical features present in the dataset.
#' Default is "gower".
#' @param seed Fixed random seed for SIMMAP algorithm. Default is 99.
#' @param ncores Number of cores to parallelize the computing process. Default is 2.
#' @param plot Set whether the feature importance of each feature for all selected hierarchical clustering methods and distances
#' should be plotted. Default is TRUE.
#' @return If plot is set to FALSE, a feature importance data frame is returned, with hierarchical clustering method and distance
#' combinations placed in rows and features and types of distances in columns. If plot is set to TRUE, a list with the feature
#' importance data frame and the plot object is returned instead.
#' @export
PFIS <- function(data,
                 cl.list,
                 dists = "euclidean",
                 mixed_dist = "gower",
                 seed = 99,
                 ncores = 2,
                 plot = TRUE) {
  list.names <- names(cl.list)
  p <- ncol(data)
  l <- length(list.names)
  all_results <- list()
  types <- sapply(data, class)
  bool <- (types == "integer" | types == "numeric")
  row.names(data) <- 1:nrow(data)
  no_dists <- TRUE

  if (length(bool[bool != TRUE]) == 0) {
    if (length(dists) > 1) no_dists <- FALSE
    used.dists <- dists
  } else {
    if (length(mixed_dist) > 1) no_dists <- FALSE
    used.dists <- mixed_dist
  }

  for (k in (1:l)) {
    hc.func <- list.names[k]
    methods <- cl.list[[hc.func]]
    m <- length(methods)
    if (no_dists == T) {
      # testing conditions
      if (length(bool[bool != TRUE]) == 0) {
        d <- stats::dist(scale(data), method = dists)
      } else {
        # specifying numerical and categorical indexes
        id_num <- which(types == "integer" | types == "numeric")
        id_cat <- which(types == "character" | types == "factor")

        d <- kmed::distmix(
          data,
          method = mixed_dist,
          idnum = id_num,
          idcat = id_cat
        ) %>% stats::as.dist()
      }

      if (anyNA(methods) == FALSE) {
        for (c in (1:m)) {
          clust <- get(hc.func)(d, method = methods[c])
          tree <- convert_to_phylo(clust)
          results <- L_score(tree,
            data,
            score = TRUE,
            ncores = 2
          )
          all_results[[paste0(hc.func, ".", methods[c])]] <- 1 - results
        }
      } else {
        clust <- get(hc.func)(d)
        tree <- convert_to_phylo(clust)
        results <- L_score(tree,
          data,
          score = TRUE,
          ncores = 2
        )
        all_results[[hc.func]] <- 1 - results
      }
    } else {
      dists.length <- length(used.dists)
      for (h in 1:dists.length) {
        # testing conditions
        if (length(bool[bool != TRUE]) == 0) {
          d <- stats::dist(scale(data), method = dists[h])
        } else {
          # specifying numerical and categorical indexes
          id_num <- which(types == "integer" | types == "numeric")
          id_cat <- which(types == "character" | types == "factor")

          d <- kmed::distmix(
            data,
            method = mixed_dist,
            idnum = id_num,
            idcat = id_cat
          ) %>% stats::as.dist()
        }

        if (anyNA(methods) == FALSE) {
          for (c in (1:m)) {
            clust <- get(hc.func)(d, method = methods[c])
            tree <- convert_to_phylo(clust)
            results <- L_score(tree,
              data,
              score = TRUE,
              ncores = 2
            )
            all_results[[paste0(
              hc.func, ".",
              methods[c], ".",
              used.dists[h]
            )]] <- c(
              1 - results,
              used.dists[h]
            )
          }
        } else {
          clust <- get(hc.func)(d)
          tree <- convert_to_phylo(clust)
          results <- L_score(tree,
            data,
            score = TRUE,
            ncores = 2
          )
          all_results[[paste0(
            hc.func, ".",
            used.dists[h]
          )]] <- c(
            1 - results,
            used.dists[h]
          )
        }
      }
    }
  }
  if (no_dists == FALSE) {
    all_cvs <- do.call(rbind, all_results)
    cvs_data.frame <- as.data.frame(all_cvs)
    cvs_data.frame[, 1:p] <- sapply(cvs_data.frame, as.numeric)
    colnames(cvs_data.frame)[1:p] <- colnames(data)
    if (plot) {
      p1 <- PFIS_plot(cvs_data.frame, no_dists = FALSE)
      methods::show(p1)
      return(list(
        importance_data = cvs_data.frame,
        plot = p1
      ))
    }
    return(cvs_data.frame)
  }
  all_cvs <- do.call(rbind, all_results)
  cvs_data.frame <- as.data.frame(all_cvs)
  colnames(cvs_data.frame) <- colnames(data)
  if (plot) {
    p1 <- PFIS_plot(cvs_data.frame, no_dists = TRUE)
    methods::show(p1)
    return(list(
      importance_data = cvs_data.frame,
      plot = p1
    ))
  }
  return(cvs_data.frame)
}

#' @title Plot PFIS for each feature of the data
#' @description
#' Phylogenetical Feature Importance Score plot of each feature with respect to each hierarchical clustering method
#' and distance combination.
#' @param importance_data Feature importance data frame returned by PFIS function. Rows must be each of the hierarchical clustering
#' methods and distance combination and columns the features from the dataset of interest.
#' @param no_dists Set whether only one distance is used to compute dissimilarity matrix (TRUE) or not (FALSE). Default is TRUE.
#' @return plot object with PFIS plot
#' @export
PFIS_plot <- function(importance_data, no_dists) {
  if (no_dists == TRUE) {
    if (nrow(importance_data) > 1) {
      p1 <- importance_data %>%
        tibble::rownames_to_column(var = "methods") %>%
        tidyr::pivot_longer(!methods, names_to = "variable", values_to = "value") %>%
        dplyr::mutate(
          variable =
            as.factor(variable)
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = variable, y = value, group = methods)) +
        ggplot2::geom_line(ggplot2::aes(colour = methods)) +
        ggplot2::geom_point(ggplot2::aes(colour = methods)) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = "Features",
          y = "PFIS"
        )
    } else {
      p1 <- importance_data %>%
        tidyr::pivot_longer(
          cols = tidyr::everything(),
          names_to = "variable",
          values_to = "value"
        ) %>%
        dplyr::mutate(
          variable =
            as.factor(variable)
        ) %>%
        dplyr::mutate(
          variable =
            forcats::fct_reorder(variable,
              value,
              .fun = "median",
              .desc = T
            )
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = variable, y = value, group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = "Features",
          y = "PFIS"
        )
    }
    return(p1)
  }
  p <- ncol(importance_data)
  n <- nrow(importance_data)

  p1 <- importance_data %>%
    dplyr::rename(distances = names(.)[p]) %>%
    tibble::rownames_to_column(var = "methods") %>%
    dplyr::mutate(methods = sub("^(.*)[.].*", "\\1", methods)) %>%
    tidyr::pivot_longer(cols = 2:p, names_to = "variable", values_to = "value") %>%
    dplyr::mutate(
      distances = as.factor(paste0("distance = ", distances)),
      variable = as.factor(variable)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = variable,
      y = value,
      group = methods
    )) +
    ggplot2::geom_line(ggplot2::aes(color = methods)) +
    ggplot2::geom_point(ggplot2::aes(color = methods)) +
    ggplot2::theme_minimal() +
    ggplot2::facet_grid(rows = ggplot2::vars(distances)) +
    ggplot2::labs(
      x = "Variables names",
      y = "PFIS",
      colour = "Methods"
    )
  return(p1)
}
