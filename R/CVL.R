#' @title Cross validation Loss for hierarchical clustering methods
#' @export
CVL <- function(data,
                cl.list,
                dists = "euclidean",
                mixed_dist = "gower",
                tol = 1e-20,
                seed = 99,
                median = F){

  list.names <- names(cl.list)
  p <- ncol(data)
  l <- length(list.names)
  all_results <- list()
  types <- sapply(data, class)
  bool <- (types == "integer" | types == "numeric")
  no_dists <- T
  row.names(data) = 1:nrow(data)

  if(length(bool[bool != T]) == 0){
    if(length(dists) > 1) no_dists <- F
    used.dists <- dists
  }else{
    if(length(mixed_dist) > 1) no_dists <- F
    used.dists <- mixed_dist
  }

  for(k in (1:l)){
    hc.func <- list.names[k]
    methods <- cl.list[[hc.func]]
    m <- length(methods)
    if(no_dists == T){
      if(anyNA(methods) == F){
        for(c in (1:m)){
          results <- numeric(p)
          for(j in 1:p){
            test_data <- data %>%
              dplyr::select(j) %>%
              as.data.frame()

            training_data <- data %>%
              dplyr::select(-j) %>%
              as.data.frame()
            if(length(bool[bool != T]) == 0){
              d <- stats::dist(scale(training_data), method = dists)
            }else{
              d <- kmed::distmix(training_data, method = mixed_dist)
            }
            clust <- get(hc.func)(d, method = methods[c])
            tree <- convert_to_phylo(clust)
            results[j] <- L_score(tree, test_data, ncores = 1)
          }
          if(median == F){
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[1])]] <-
              c(mean(results, na.rm = T),
                used.dists[1])
          }else{
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[1])]] <-
              c(median(results, na.rm = T),
                used.dists[1])
          }
        }}else{
          results <- numeric(p)
          for(j in 1:p){
            test_data <- data %>%
              dplyr::select(j) %>%
              as.data.frame()

            training_data <- data %>%
              dplyr::select(-j) %>%
              as.data.frame()

            if(length(bool[bool != T]) == 0){
              d <- stats::dist(scale(training_data), method = dists)
            }else{
              d <- kmed::distmix(training_data, method = mixed_dist)
            }

            clust <- get(hc.func)(d)
            tree <- convert_to_phylo(clust)
            results[j] <- L_score(tree, test_data, ncores = 1)
          }
          if(median == F){
            all_results[[paste0(hc.func, ".", used.dists[1])]] <-
              c(mean(results, na.rm = T), used.dists[1])
          }else{
            all_results[[paste0(hc.func, ".", used.dists[1])]] <-
              c(median(results, na.rm = T), used.dists[1])
          }
        }
    }
    else{
      dists.length <- length(used.dists)
      for(h in 1:dists.length){
        if(anyNA(methods) == F){
          for(c in (1:m)){
            results <- numeric(p)
            for(j in 1:p){
              test_data <- data %>%
                dplyr::select(j) %>%
                as.data.frame()

              training_data <- data %>%
                dplyr::select(-j) %>%
                as.data.frame()
              if(length(bool[bool != T]) == 0){
                d <- stats::dist(scale(training_data), method = dists[h])
              }else{
                d <- kmed::distmix(training_data, method = mixed_dist[h])
              }
              clust <- get(hc.func)(d, method = methods[c])
              tree <- convert_to_phylo(clust)
              results[j] <- L_score(tree, test_data, ncores = 1)
            }
            if(median == F){
              all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] <-
                c(mean(results, na.rm = T),
                  used.dists[h])
            }else{
              all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] <-
                c(median(results, na.rm = T),
                  used.dists[h])
            }
          }}else{
            results <- numeric(p)
            for(j in 1:p){
              test_data <- data %>%
                dplyr::select(j) %>%
                as.data.frame()

              training_data <- data %>%
                dplyr::select(-j) %>%
                as.data.frame()

              if(length(bool[bool != T]) == 0){
                d <- stats::dist(scale(training_data), method = dists[h])
              }else{
                d <- kmed::distmix(training_data, method = mixed_dist[h])
              }

              clust <- get(hc.func)(d)
              tree <- convert_to_phylo(clust)
              results[j] <- L_score(tree, test_data, ncores = 1)
            }
            if(median == F){
              all_results[[paste0(hc.func, ".", used.dists[h])]] <-
                c(mean(results, na.rm = T), used.dists[h])
            }else{
              all_results[[paste0(hc.func, ".", used.dists[h])]] <-
                c(median(results, na.rm = T), used.dists[h])
            }
          }
      }
    }
  }
  if(no_dists == F){
    all_cvs <- do.call(rbind, all_results)
    cvs_data.frame <- as.data.frame(all_cvs)
    cvs_data.frame$V1 <- as.numeric(cvs_data.frame$V1)
    colnames(cvs_data.frame) <- c("loss", "distance")
    return(cvs_data.frame)
  }else{
    all_cvs <- do.call(rbind, all_results)
    cvs_data.frame <- as.data.frame(all_cvs)
    return(cvs_data.frame)
  }
}
