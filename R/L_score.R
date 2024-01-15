# L_score main function
#' @title Leave one feature out cross validation loss
#' @description Leave one feature out cross validation loss for a given dendrogram and dataset. This is an auxiliary function to
#' the main functions CVL and PFIS but can be used separately.
#' @param dend Dendrogram of interest (an object of class "phylo")
#' @param original_data Entire dataset (data.frame or matrix) or feature of interest (numeric or factor vector) to compute cross validate loss
#' @param seed Fixed random seed for SIMMAP algorithm. Default is 99.
#' @param ncores Number of cores to paralelize the computing process. Default is 2.
#' @param score Set whether a vector of leave one feature out loss should be computed (TRUE) or compute the overall mean of the loss
#' directly (FALSE)
#' @return A numeric vector of leave one feature out losses.
#' @importFrom foreach %dopar%
#' @export
L_score <- function(dend,
                   original_data,
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
