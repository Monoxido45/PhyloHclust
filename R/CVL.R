#' @title Cross validation Loss for hierarchical clustering methods
#' @description
#' Leave one feature out Cross validation loss for comparing a list of hierarchical clustering methods.
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
#' @param median Set whether to compute the median of the leave one feature out loss associated to every feature (TRUE) or the
#'  mean (FALSE). Default is FALSE.
#' @return A data frame comparing the CVL of each hierarchical clustering method and distance combinations.
#' @export
CVL <- function(data,
                cl.list,
                dists = "euclidean",
                mixed_dist = "gower",
                seed = 99,
                ncores = 2,
                median = FALSE) {
  list.names <- names(cl.list)
  p <- ncol(data)
  l <- length(list.names)
  all_results <- list()
  types <- sapply(data, class)
  # specifying numerical and categorical indexes
  id_num <- which(types == "integer" | types == "numeric")
  id_cat <- which(types == "character" | types == "factor")
  bool <- (types == "integer" | types == "numeric")
  no_dists <- TRUE
  row.names(data) <- 1:nrow(data)

  if (length(bool[bool != TRUE]) == 0) {
    if (length(dists) > 1) no_dists <- F
    used.dists <- dists
  } else {
    if (length(mixed_dist) > 1) no_dists <- F
    used.dists <- mixed_dist
  }

  for (k in (1:l)) {
    hc.func <- list.names[k]
    methods <- cl.list[[hc.func]]
    m <- length(methods)
    if (no_dists == TRUE) {
      if (anyNA(methods) == F) {
        for (c in (1:m)) {
          results <- numeric(p)
          for (j in 1:p) {
            test_data <- data %>%
              dplyr::select(j) %>%
              as.data.frame()

            training_data <- data %>%
              dplyr::select(-j) %>%
              as.data.frame()
            if (length(bool[bool != TRUE]) == 0) {
              d <- stats::dist(scale(training_data), method = dists)
            } else {
              d <- kmed::distmix(
                training_data,
                method = mixed_dist,
                idnum = id_num,
                idcat = id_cat
              ) %>% stats::as.dist()
            }
            clust <- get(hc.func)(d, method = methods[c])
            tree <- convert_to_phylo(clust)
            results[j] <- L_score(tree, test_data, ncores = ncores)
          }
          if (median == F) {
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[1])]] <-
              c(
                mean(results, na.rm = TRUE),
                used.dists[1]
              )
          } else {
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[1])]] <-
              c(
                median(results, na.rm = TRUE),
                used.dists[1]
              )
          }
        }
      } else {
        results <- numeric(p)
        for (j in 1:p) {
          test_data <- data %>%
            dplyr::select(j) %>%
            as.data.frame()

          training_data <- data %>%
            dplyr::select(-j) %>%
            as.data.frame()

          if (length(bool[bool != TRUE]) == 0) {
            d <- stats::dist(scale(training_data), method = dists)
          } else {
            d <- kmed::distmix(
              training_data,
              method = mixed_dist,
              idnum = id_num,
              idcat = id_cat
            ) %>% stats::as.dist()
          }

          clust <- get(hc.func)(d)
          tree <- convert_to_phylo(clust)
          results[j] <- L_score(tree, test_data, ncores = ncores)
        }
        if (median == F) {
          all_results[[paste0(hc.func, ".", used.dists[1])]] <-
            c(mean(results, na.rm = TRUE), used.dists[1])
        } else {
          all_results[[paste0(hc.func, ".", used.dists[1])]] <-
            c(median(results, na.rm = TRUE), used.dists[1])
        }
      }
    } else {
      dists.length <- length(used.dists)
      for (h in 1:dists.length) {
        if (anyNA(methods) == F) {
          for (c in (1:m)) {
            results <- numeric(p)
            for (j in 1:p) {
              test_data <- data %>%
                dplyr::select(j) %>%
                as.data.frame()

              training_data <- data %>%
                dplyr::select(-j) %>%
                as.data.frame()
              if (length(bool[bool != TRUE]) == 0) {
                d <- stats::dist(scale(training_data), method = dists[h])
              } else {
                d <- kmed::distmix(
                  training_data,
                  method = mixed_dist[h],
                  idnum = id_num,
                  idcat = id_cat
                ) %>% stats::as.dist()
              }
              clust <- get(hc.func)(d, method = methods[c])
              tree <- convert_to_phylo(clust)
              results[j] <- L_score(tree, test_data, ncores = ncores)
            }
            if (median == F) {
              all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] <-
                c(
                  mean(results, na.rm = TRUE),
                  used.dists[h]
                )
            } else {
              all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] <-
                c(
                  median(results, na.rm = TRUE),
                  used.dists[h]
                )
            }
          }
        } else {
          results <- numeric(p)
          for (j in 1:p) {
            test_data <- data %>%
              dplyr::select(j) %>%
              as.data.frame()

            training_data <- data %>%
              dplyr::select(-j) %>%
              as.data.frame()

            if (length(bool[bool != TRUE]) == 0) {
              d <- stats::dist(scale(training_data), method = dists[h])
            } else {
              d <- kmed::distmix(
                training_data,
                method = mixed_dist[h],
                idnum = id_num,
                idcat = id_cat
              ) %>% stats::as.dist()
            }

            clust <- get(hc.func)(d)
            tree <- convert_to_phylo(clust)
            results[j] <- L_score(tree, test_data, ncores = ncores)
          }
          if (median == F) {
            all_results[[paste0(hc.func, ".", used.dists[h])]] <-
              c(mean(results, na.rm = TRUE), used.dists[h])
          } else {
            all_results[[paste0(hc.func, ".", used.dists[h])]] <-
              c(median(results, na.rm = TRUE), used.dists[h])
          }
        }
      }
    }
  }
  if (no_dists == F) {
    all_cvs <- do.call(rbind, all_results)
    cvs_data.frame <- as.data.frame(all_cvs)
    cvs_data.frame$V1 <- as.numeric(cvs_data.frame$V1)
    colnames(cvs_data.frame) <- c("loss", "distance")
    return(cvs_data.frame)
  } else {
    all_cvs <- do.call(rbind, all_results)
    cvs_data.frame <- as.data.frame(all_cvs)
    return(cvs_data.frame)
  }
}
