# L_score main function
#' @importFrom foreach %dopar%
#' @export
L_score <- function(dend,
                   original_data,
                   tol  = 1e-20,
                   seed = 99,
                   ncores = 2,
                   score = FALSE){
  types <- sapply(original_data, class)
  p <- ncol(original_data)
  n <- nrow(original_data)
  total <- p*n
  score.matrix <- matrix(0, nrow = n, ncol = p)
  dend$edge.length[which(dend$edge.length %in% c(0))] <- 10^(-3)
  dend$edge.length[which(dend$edge.length < 0)] <- 10^(-3)
  names <- row.names(original_data)
  # paralellizing
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c("row_computing", "factorial.missing",
                        "onehotencoder"))
  doParallel::registerDoParallel(cl)
  loss.matrix <- foreach::foreach(j = 1:p, .combine = cbind,
                           .export = c("row_computing",
                                       "factorial.missing",
                                       "onehotencoder"),
                           .packages = c("magrittr")) %dopar% {
                             lines = row_computing(types,
                                                   dend,
                                                   original_data,
                                                   tol,
                                                   seed,
                                                   j)
                             lines
                                       }
  parallel::stopCluster(cl)
  if(p == 1){
      loss.matrix <- as.matrix(loss.matrix)
  }
  partial_loss <- colSums(loss.matrix)
  if(score == TRUE){
    loss <- partial_loss/n
  }else{
    loss <- sum(partial_loss)/total
  }
  return(loss)
}
