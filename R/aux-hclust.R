# auxilizary functions
# onehot encoder for categorical data
onehotencoder <- function(miss_data){
  nlvls <- nlevels(miss_data)
  onehot <- matrix(0, nrow = length(miss_data), ncol = nlvls)
  factors <- as.numeric(miss_data)
  for (i in 1:length(factors)){
    if(is.na(factors[i]) == TRUE){
      onehot[i, ] <- rep(1/nlvls, nlvls)
    }else{
      onehot[i, factors[i]] <- 1
    }
  }
  row.names(onehot) <- names(miss_data)
  colnames(onehot) <- levels(miss_data)
  return(onehot)
}

# computing pi for categorical loss
factorial.missing <- function(dend,
                              miss_data,
                              row,
                              nsim = 50){
  onehot <- onehotencoder(miss_data)
  fit_discrete <- phytools::make.simmap(phytools::as.multiPhylo(dend),
                                        onehot, model = "ER",
                                        message = FALSE, nsim = nsim)
  results <- phytools::describe.simmap(fit_discrete, plot = F)
  which_row <- which(rownames(results$tips) == row)
  pi <- results$tips[which_row, ]
  return(pi)
}

# computing distances for each type of feature
row_computing <- function(types,
                          dend,
                          original_data,
                          tol,
                          seed,
                          col){

  n <- nrow(original_data)
  lines <- numeric(n)
  names <- row.names(original_data)

  if (types[col] == "integer" | types[col] == "numeric"){
    # scaling each data
    current_data <- original_data %>%
      dplyr::select(col) %>%
      as.matrix() %>%
      scale()
    rownames(current_data) <- names
    # reordering according to dendrogram
    current_data <- as.matrix(current_data[dend$tip.label, ])
    fit <- mvMORPH::mvBM(dend,
                         current_data,
                         model = "BMM",
                         echo = T,
                         method = "pic")
    for (i in (1:n)){
      saved_value <- current_data[i]
      miss_data <- current_data
      miss_data[i, 1] <- NA
      imp <- mvMORPH::estim(dend, miss_data, fit)
      y.pred <- imp$estimates[i, 1]
      mse <- ((y.pred - saved_value)^2)
      lines[i] <-  mse
    }
  }else{
    current_data <- original_data %>%
      dplyr::pull(col)
    if (types[col] == "logical") current_data <- as.factor(current_data)
    set.seed(seed)
    for(i in (1:n)){
      saved_value <- numeric(nlevels(current_data))
      saved_value[as.numeric(current_data)[i]] <- 1
      miss_data <- current_data
      A <- nlevels(miss_data)
      names(miss_data) <- names
      miss_data[i] <- NA
      pi <- factorial.missing(dend, miss_data, i)
      print(pi)
      mse <- (sum((saved_value - pi)^2))
      lines[i] <- mse/A
    }
  }
  return(lines)
}


