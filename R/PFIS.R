#' @title Phylogenetical Feature Importance Score for Hierarchical clustering methods
#' @export
PFIS <- function(data,
                 cl.list,
                 dists = "euclidean",
                 mixed_dist = "gower",
                 tol = 1e-20,
                 seed = 99,
                 ncores = 2,
                 plot = TRUE){

  list.names <- names(cl.list)
  p <- ncol(data)
  l <- length(list.names)
  all_results <- list()
  types <- sapply(data, class)
  bool <- (types == "integer" | types == "numeric")
  row.names(data) <- 1:nrow(data)
  no_dists <- T

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

      # testing conditions
      if(length(bool[bool != T]) == 0){
        d <- stats::dist(scale(data), method = dists)
      }else{
        d <- kmed::distmix(data, method = mixed_dist)
      }

      if(anyNA(methods) == F){
        for(c in (1:m)){
          clust <- get(hc.func)(d, method = methods[c])
          tree <- convert_to_phylo(clust)
          results <- L_score(tree,
                             data,
                             score = TRUE,
                             ncores = 2)
          all_results[[paste0(hc.func, ".", methods[c])]] <- 1 - results
        }
      }else{
        clust = get(hc.func)(d)
        tree = convert_to_phylo(clust)
        results <- L_score(tree,
                           data,
                           score = TRUE,
                           ncores = 2)
        all_results[[hc.func]] = 1 - results
      }

    }else{
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        # testing conditions
        if(length(bool[bool != T]) == 0){
          d <- stats::dist(scale(data), method = dists[h])
        }else{
          d <- kmed::distmix(data, method = mixed_dist[h])
        }

        if(anyNA(methods) == F){
          for(c in (1:m)){
            clust = get(hc.func)(d, method = methods[c])
            tree = convert_to_phylo(clust)
            results <- L_score(tree,
                               data,
                               score = TRUE,
                               ncores = 2)
            all_results[[paste0(hc.func, ".",
                                methods[c], ".",
                                used.dists[h])]] <- c(1 - results,
                                                      used.dists[h])
          }
        }else{
          clust = get(hc.func)(d)
          tree = convert_to_phylo(clust)
          results <- L_score(tree,
                             data,
                             score = TRUE,
                             ncores = 2)
          all_results[[paste0(hc.func, ".",
                              used.dists[h])]] <- c(1 - results,
                                                    used.dists[h])
        }
      }
    }
  }
  if(no_dists == FALSE){
    all_cvs <- do.call(rbind, all_results)
    cvs_data.frame <- as.data.frame(all_cvs)
    cvs_data.frame[, 1:p] <- sapply(cvs_data.frame, as.numeric)
    colnames(cvs_data.frame)[1:p] <- colnames(data)
    if(plot){
      p1 <- PFIS_plot(cvs_data.frame, no_dists = FALSE)
      methods::show(p1)
      return(list(importance_data = cvs_data.frame,
                  plot = p1))
    }
    return(cvs_data.frame)
  }
  all_cvs <- do.call(rbind, all_results)
  cvs_data.frame <- as.data.frame(all_cvs)
  colnames(cvs_data.frame) <- colnames(data)
  if(plot){
    p1 <- PFIS_plot(cvs_data.frame, no_dists = TRUE)
    methods::show(p1)
    return(list(importance_data = cvs_data.frame,
                plot = p1))
  }
  return(cvs_data.frame)
}

#' @title Plot PFIS for each feature of the data
#' @export
PFIS_plot <- function(importance_data, no_dists){
  if(no_dists == TRUE){
    p1 <- importance_data %>%
      reshape2::melt() %>%
      dplyr::mutate(variable =
                      as.factor(
                        colnames(importance_data))) %>%
      dplyr::mutate(variable =
                      forcats::fct_reorder(variable,
                                           value,
                                           .fun = 'median',
                                           .desc = T)) %>%
      ggplot2::ggplot(ggplot2::aes(x = variable, y = value, group = 1)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Features",
                    y = "PFIS")
    return(p1)
  }
  p <- ncol(importance_data)
  n <- nrow(importance_data)
  dists <- importance_data %>%
    dplyr::pull(p)
  p1 <- importance_data %>%
    select(-p) %>%
    reshape2::melt() %>%
    dplyr::mutate(dist = rep(as.factor(dists), p - 1),
                  obs = rep(as.factor(1:n), p -1),
                  variable = as.factor(variable)) %>%
    dplyr::mutate(variable = forcats::fct_reorder(variable,
                                                  value,
                                                  .fun = 'median',
                                                  .desc = T)) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = variable,
      y  = value,
      group = obs)) +
    ggplot2::geom_line(ggplot2::aes(color = dist)) +
    ggplot2::geom_point(ggplot2::aes(color = dist)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Variables names",
                  y = "PFIS",
                  colour = "Distances")
  return(p1)

}
