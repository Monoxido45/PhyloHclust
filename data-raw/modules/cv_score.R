# code used to make the score function of interest
library(ape)
library(dendextend)
library(cluster)
library(tibble)
library(magrittr)
library(dplyr)
library(phytools)
library(mltools)
library(data.table)
library(factoextra)
setwd("~/estatistica_UFSCAR/cv_cluster/modules")
source("convert_to_parenthesis.R")
library(foreach)
library(doParallel)
library(tictoc)
library(mvMORPH)
library(kmed)
library(gtools)
library(MLmetrics)

# defining one hot encoder for simap
onehotencoder = function(miss_data){
  nlvls = nlevels(miss_data)
  onehot = matrix(0, nrow = length(miss_data), ncol = nlvls)
  factors = as.numeric(miss_data)
  for (i in 1:length(factors)){
    if(is.na(factors[i]) == TRUE){
      onehot[i, ] = rep(1/nlvls, nlvls)
    }else{
      onehot[i, factors[i]] = 1
    }
  }
  row.names(onehot) = names(miss_data)
  colnames(onehot) = levels(miss_data)
  return(onehot)
}

factorial.missing = function(dend, miss_data, row){
  onehot = onehotencoder(miss_data)
  fit_discrete = make.simmap(as.multiPhylo(dend), onehot, model = "ER",
                             message = FALSE)
  results = describe.simmap(fit_discrete, plot = F)
  which_row = which(rownames(results$tips) == row)
  pi = results$tips[which_row, ]
  return(pi)
}


# without variance
row_computing = function(types, dend, original_data, tol,seed, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  tree.order = dend$tip.label
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    current_data = as.matrix(original_data[, j])
    # not scaling
    colnames(current_data) = cname
    rownames(current_data) = names
    current_data = as.matrix(current_data[dend$tip.label, ])
    fit = mvBM(dend, current_data, model = "BMM", echo = T, method = "pic")
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = current_data
      names(miss_data) = names
      miss_data[i, 1] = NA
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[i, 1]
      mse = ((y.pred - saved_value)^2)
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    set.seed(seed)
    for (i in (1:n)){
      saved_value = numeric(nlevels(current_data))
      saved_value[as.numeric(current_data)[i]] = 1
      miss_data = current_data
      A = nlevels(miss_data)
      names(miss_data) = names
      miss_data[i] = NA
      pi = factorial.missing(dend, miss_data, i)
      print(pi)
      mse = (sum((saved_value - pi)^2))
      lines[i] = mse/A
    }
  }
  return(lines)
}

# scaled row_computing
scaled_row_computing = function(types, dend, original_data, tol,seed, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  tree.order = dend$tip.label
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    current_data = scale(as.matrix(original_data[, j]))
    # scaling
    colnames(current_data) = cname
    rownames(current_data) = names
    current_data = as.matrix(current_data[dend$tip.label, ])
    fit = mvBM(dend, current_data, model = "BMM", echo = T, method = "pic")
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = current_data
      names(miss_data) = names
      miss_data[i, 1] = NA
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[i, 1]
      mse = ((y.pred - saved_value)^2)
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    set.seed(seed)
    for (i in (1:n)){
      saved_value = numeric(nlevels(current_data))
      saved_value[as.numeric(current_data)[i]] = 1
      miss_data = current_data
      A = nlevels(miss_data)
      names(miss_data) = names
      miss_data[i] = NA
      pi = factorial.missing(dend, miss_data, i)
      print(pi)
      mse = (sum((saved_value - pi)^2))
      lines[i] = mse/A
    }
  }
  return(lines)
}
# Defining FOM function for simple computation of figure of merit score

FOM = function(data, nlvls, cl.list, dists = NA, mixed_dist = NA, scale = F){
  list.names = names(cl.list)
  n = nrow(data)
  p = ncol(data)
  l = length(list.names)
  all_results = list()
  dists.names = numeric(0)
  
  types = sapply(data, class)
  bool =  (types == "integer" | types == "numeric")
  if(length(bool[bool != T]) == 0){
    if(length(dists) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = dists
    }
  }else{
    if(is.na(mixed_dist) == T | length(mixed_dist) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = mixed_dist
    }
  }
  
  for(k in (1:l)){
    hc.func = list.names[k]
    methods = cl.list[[hc.func]]
    m = length(methods)
    if(no_dists == T){
      if(is.na(methods) == F){
        for(c in (1:m)){
          results = matrix(nrow = p, ncol = n)
          for(j in 1:p){
            training_data = data[, -j]
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if((is.na(mixed_dist) == F) & (length(mixed_dist) == 1)){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            clust = get(hc.func)(d, method = methods[c])
            lvls = factor(cutree(clust, k = nlvls))
            if(scale == F){
            test_data = data.frame(variable = data[, j],
                                   factors = lvls)
            }else{
              test_data = data.frame(variable = scale(data[, j]),
                                     factors = lvls)
            }
            rownames(test_data) = rownames(data)
            means = with(test_data, tapply(variable, 
                                           factors, mean))
            results[j, ] = (test_data[, 1] - means[as.numeric(test_data[, 2])])^2
          }
          sums = sqrt(1/n*rowSums(results))
          final_result = sum(sums)
          all_results[[paste0(hc.func, ".", methods[c])]] = final_result
        }}else{
          results = matrix(nrow = p, ncol = n)
          for(j in 1:p){
            training_data = data[, -j]
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if(length(mixed_dist) == 1){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            clust = get(hc.func)(d)
            lvls = factor(cutree(clust, k = nlvls))
            if(scale == F){
              test_data = data.frame(variable = data[, j],
                                     factors = lvls)
            }else{
              test_data = data.frame(variable = scale(data[, j]),
                                     factors = lvls)
            }
            rownames(test_data) = rownames(data)
            means = with(test_data, tapply(variable, 
                                           factors, mean))
            results[j, ] = (test_data[, 1] - means[as.numeric(test_data[, 2])])^2
          }
          sums = sqrt(1/n*rowSums(results))
          final_result = sum(sums)
          all_results[[hc.func]] = final_result
        }
    }else{
      if(length(bool[bool != T]) == 0){
        func = "dist"
      }else{
        func = "distmix"
      }
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        if(is.na(methods) == F){
          for(c in (1:m)){
            results = matrix(nrow = p, ncol = n)
            for(j in 1:p){
              training_data = data[, -j]
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d, method = methods[c])
              lvls = factor(cutree(clust, k = nlvls))
              if(scale == F){
                test_data = data.frame(variable = data[, j],
                                       factors = lvls)
              }else{
                test_data = data.frame(variable = scale(data[, j]),
                                       factors = lvls)
              }
              rownames(test_data) = rownames(data)
              means = with(test_data, tapply(variable, 
                                             factors, mean))
              results[j, ] = (test_data[, 1] - means[as.numeric(test_data[, 2])])^2
            }
            sums = sqrt(1/n*rowSums(results))
            final_result = sum(sums)
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] = c(final_result,
                                                                              used.dists[h])
          }}else{
            results = matrix(nrow = p, ncol = n)
            for(j in 1:p){
              training_data = data[, -j]
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d)
              lvls = factor(cutree(clust, k = nlvls))
              if(scale == F){
                test_data = data.frame(variable = data[, j],
                                       factors = lvls)
              }else{
                test_data = data.frame(variable = scale(data[, j]),
                                       factors = lvls)
              }
              rownames(test_data) = rownames(data)
              means = with(test_data, tapply(variable, 
                                             factors, mean))
              results[j, ] = (test_data[, 1] - means[as.numeric(test_data[, 2])])^2
            }
            sums = sqrt(1/n*rowSums(results))
            final_result = sum(sums)
            all_results[[paste0(hc.func, ".", used.dists[h])]] = c(final_result,
                                                                        used.dists[h])
          }
      }
    }
  }
  if(length(dists.names) > 0){
    all_foms = do.call(rbind, all_results)
    foms_data.frame = as.data.frame(all_foms)
    foms_data.frame[, 1] = as.numeric(foms_data.frame[, 1])
    return(foms_data.frame)
  }else{
    all_foms = as.data.frame(do.call(rbind, all_results))
    return(all_foms)
  }
}



# using mvMORPH
L_score = function(dend, original_data, scale = F, tol  = 1e-20, seed = 99){
  types = sapply(original_data, class)
  p = ncol(original_data)
  n = nrow(original_data)
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  dend$edge.length[which(dend$edge.length %in% c(0))] = 10^(-3)
  dend$edge.length[which(dend$edge.length < 0)] = 10^(-3)
  names = row.names(original_data)
  row.names(original_data)  = c(1:nrow(original_data))
  
  # paralellizing
  if(scale == F){
  cores = detectCores()
  cl = makeCluster(cores[1] - 3)
  clusterExport(cl, c("row_computing", "factorial.missing",
                      "onehotencoder"))
  registerDoParallel(cl)
  score.matrix = foreach(j = 1:p, .combine = cbind,
                         .export = c("row_computing", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                          lines = row_computing(types, dend, original_data, tol, seed, j)
                          lines
                         }
  stopCluster(cl)
  if(p == 1){
    score.matrix = as.matrix(score.matrix)
  }
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
  }else{
    cores = detectCores()
    cl = makeCluster(cores[1] - 1)
    clusterExport(cl, c("scaled_row_computing", "factorial.missing",
                        "onehotencoder"))
    registerDoParallel(cl)
    score.matrix = foreach(j = 1:p, .combine = cbind,
                           .export = c("scaled_row_computing", "factorial.missing","onehotencoder"),
                           .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                             lines = scaled_row_computing(types, dend, original_data, tol, seed, j)
                             lines
                           }
    stopCluster(cl)
    if(p == 1){
      score.matrix = as.matrix(score.matrix)
    }
    partial_score = colSums(score.matrix)
    score = sum(partial_score)/total
    return(score)
  }
}


L_cross_val = function(original_data, cl.list, dists = NA, mixed_dist = NA, tol = 1e-20, seed = 99, 
                       scale = F, median = F){
  list.names = names(cl.list)
  p = ncol(original_data)
  l = length(list.names)
  all_results = list()
  types = sapply(original_data, class)
  bool =  (types == "integer" | types == "numeric")
  if(length(bool[bool != T]) == 0){
    if(is.na(dists) == T | length(dists) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = dists
    }
  }else{
    if(is.na(mixed_dist) == T | length(mixed_dist) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = mixed_dist
    }
  }
  
  for(k in (1:l)){
    hc.func = list.names[k]
    methods = cl.list[[hc.func]]
    m = length(methods)
    if(no_dists == T){
      if(is.na(methods) == F){
        for(c in (1:m)){
          results = numeric(p)
          for(j in 1:p){
            test_data = as.data.frame(original_data[, j])
            colnames(test_data) = colnames(original_data)[j]
            rownames(test_data) = rownames(original_data)
            training_data = original_data[, -j]
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if((is.na(mixed_dist) == F) & (length(mixed_dist) == 1)){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            clust = get(hc.func)(d, method = methods[c])
            tree = convert_to_phylo(clust)
            results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                  silent = T, error=function(msg){
                                    return(NA)
                                  })
          }
          if(median == F){      
            all_results[[paste0(hc.func, ".", methods[c])]] = sum(results, na.rm = T)/(length(results) - 
                                                                                         sum(is.na(results)))
          }else{
            all_results[[paste0(hc.func, ".", methods[c])]] = median(results, na.rm = T)
          }
        }}else{
          results = numeric(p)
          for(j in 1:p){
            test_data = as.data.frame(original_data[, j])
            colnames(test_data) = colnames(original_data)[j]
            rownames(test_data) = rownames(original_data)
            training_data = original_data[, -j]
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if(length(mixed_dist) == 1){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            
            clust = get(hc.func)(d)
            tree = convert_to_phylo(clust)
            results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                  silent = T, error=function(msg){
                                    return(NA)
                                  })
          }
          if(median == F){      
            all_results[[hc.func]] = sum(results, na.rm = T)/(length(results) - 
                                                                sum(is.na(results)))
          }else{
            all_results[[hc.func]] = median(results, na.rm = T)
          }
        }
    }else{
      if(length(bool[bool != T]) == 0){
        func = "dist"
      }else{
        func = "distmix"
      }
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        if(is.na(methods) == F){
          for(c in (1:m)){
            results = numeric(p)
            for(j in 1:p){
              test_data = as.data.frame(original_data[, j])
              colnames(test_data) = colnames(original_data)[j]
              rownames(test_data) = rownames(original_data)
              training_data = original_data[, -j]
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d, method = methods[c])
              tree = convert_to_phylo(clust)
              results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                    silent = T, error=function(msg){
                                      return(NA)
                                    })
            }
            if(median == F){      
              all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] = 
                c(sum(results, na.rm = T)/(length(results) - sum(is.na(results))),
                  used.dists[h])
            }else{
              all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] = 
                c(median(results, na.rm = T), used.dists[h])
            }
          }}else{
            results = numeric(p)
            for(j in 1:p){
              test_data = as.data.frame(original_data[, j])
              colnames(test_data) = colnames(original_data)[j]
              rownames(test_data) = rownames(original_data)
              training_data = original_data[, -j]
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d)
              tree = convert_to_phylo(clust)
              results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                    silent = T, error=function(msg){
                                      return(NA)
                                    })
            }
            if(median == F){      
              all_results[[paste0(hc.func, ".", used.dists[h])]] = 
                c(sum(results, na.rm = T)/(length(results) - sum(is.na(results))),  
                  used.dists[h])
            }else{
              all_results[[paste0(hc.func, ".", used.dists[h])]] = 
                c(median(results, na.rm = T), used.dists[h])
            }
          }
      }
    }
  }
  if(no_dists == F){
    all_cvs = do.call(rbind, all_results)
    cvs_data.frame = as.data.frame(all_cvs)
    cvs_data.frame$V1 = as.numeric(cvs_data.frame$V1)
    return(cvs_data.frame)
  }else{
    all_cvs = do.call(rbind, all_results)
    cvs_data.frame = as.data.frame(all_cvs)
    return(cvs_data.frame)
  }
}



all_supervised_comparing = function(data, clust_list, test_index, dist = NA, scale = F, 
                                    mixed_dists = NA){
  list.names = names(clust_list)
  training_data = data[, -test_index]
  types = sapply(training_data, class)
  bool = (types == "integer" | types == "numeric")
  size = 1
  comp_list = list()
  if(length(bool[bool != T]) == 0){
    if(is.na(dist) == F){size = length(dist)
    mixed_dists = NA
    dist.names = dist}
  }else{
    if(is.na(mixed_dists) == F){size = length(mixed_dists)
    dist = NA
    dist.names = mixed_dists}
  }
  for(i in (1:length(list.names))){
    clust.method = list.names[i]
    all.methods = test.list[[i]]
    if(is.na(dist) & is.na(mixed_dists)){mat = matrix(nrow = length(all.methods),
                                                      ncol = 3)}else{
                                        mat = matrix(nrow = length(all.methods)*size,
                                                     ncol = 4)
                                                      }
    row.names(mat) = numeric(length = size * length(all.methods))
    for(j in (1:length(all.methods))){
      if(is.na(dist) & is.na(mixed_dists)){
        mat[j, ] = supervised_comparing(data, clust.method, all.methods[j], test_index, 
                                        scale = scale)
        row.names(mat)[j] = paste0(clust.method,".", all.methods[j])
      }else{
        for(k in (1:size)){
          mat[k + ((j - 1)* size), c(1,2,3)] = supervised_comparing(data, clust.method, 
                                          all.methods[j], test_index, dist = dist[k], 
                                          mixed_dists = mixed_dists[k],
                                          scale = scale)
          mat[k + ((j - 1)* size),4] = dist.names[k]
          row.names(mat)[k + ((j - 1)* size)] = paste0(clust.method,".", all.methods[j])
        }
      }
    }
    comp_list[[clust.method]] = mat
  }
comp_data.frame = as.data.frame(do.call(rbind, comp_list))
  if(is.na(dist) == F | is.na(mixed_dists) == F){
    comp_data.frame = comp_data.frame %>%
      mutate(V1 = as.numeric(V1),
             V2 = as.numeric(V2),
             V3 = as.numeric(V3),
             V4 = as.factor(V4))
  }
return(comp_data.frame)
}


supervised_comparing = function(data, clust, clust_method, test_index, dist = NA, 
                                mixed_dists = NA, p_mink = 3.5, scale = F){
  n = nrow(data)
  p = ncol(data) - 1
  training_data = data[, - test_index]
  test_data = data[, test_index]
  lvls = levels(test_data)
  nlvls= length(lvls)
  types = sapply(training_data, class)
  relab_y = mapvalues(test_data, from = lvls, to = 1:nlvls)
  new_lvls = levels(relab_y)
  bool =  (types == "integer" | types == "numeric")
  if (length(bool[bool != T]) == 0){
    if(is.na(dist) == TRUE){ d = dist(scale(training_data))}else{
      d = dist(scale(training_data), method = dist, p = p_mink)}
  }else{
    if(is.na(mixed_dists) == T){d = daisy(training_data, metric = "gower")}else{ 
      d = distmix(training_data, method = mixed_dists)}
  }
  if(clust == "hclust"){
      clust = hclust(d, method = clust_method)
      clusters = factor(cutree(clust, k =  nlvls))
    }else{
      if(clust == "agnes"){
        clust = agnes(d, diss = T, method = clust_method)
        clusters = factor(cutree(clust, k = nlvls))
      }else{
        clust = diana(d, diss = T)
        clusters = factor(cutree(clust, k =  nlvls))
      }
    }
  perms = permutations(n = nlvls, r = nlvls,v = new_lvls)
  hits = numeric(nrow(perms))
  clusts = matrix(0, nrow = nrow(perms), ncol = length(clusters))
  for(i in (1:nrow(perms))){
    new_clust = mapvalues(clusters, from = levels(clusters), to = perms[i, ])
    clusts[i, ] = new_clust
    temp_hit = (sum(new_clust == relab_y)/length(relab_y))
    hits[i] = temp_hit
  }
  max.index = which.max(hits)
  F1_s = F1_Score(relab_y, clusts[max.index, ])
  maximum_index = hits[max.index]
  dend = convert_to_phylo(clust)
  score = tryCatch(L_score(dend, training_data, scale = scale),  silent = T, error=function(e) e)
  if(inherits(score, "error")){
    score = NA
  }
  names = c("Hits proportion", "F1", "Score")
  comparisson = setNames(c(maximum_index, F1_s, score), names)
  return(comparisson)
}

supervised_comparing_L_cross_val = function(data, clust_list, test_index, dists = NA, mixed_dist = NA,
                                            tol = 1e-20, scale = F, median = F){
  list.names = names(clust_list)
  training_data = data[, -test_index]
  test_data = data[, test_index]
  lvls = levels(test_data)
  nlvls= length(lvls)
  types = sapply(training_data, class)
  bool = (types == "integer" | types == "numeric")
  relab_y = mapvalues(test_data, from = lvls, to = 1:nlvls)
  new_lvls = levels(relab_y)
  cross_val_results = L_cross_val(training_data, clust_list, dists = dists, mixed_dist = mixed_dist, 
              tol = tol, seed = 99, scale =  scale, median = median)
  all_hits = numeric(0)
  all_F1 = numeric(0)
  if(length(bool[bool != T]) == 0){
    if(is.na(dists) == T | length(dists) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = dists
    }
  }else{
    if(is.na(mixed_dist) == T | length(mixed_dist) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = mixed_dist
    }
  }
  for(k in 1:length(list.names)){
    hc.func = list.names[k]
    methods = clust_list[[hc.func]]
    m = length(methods)
    if(no_dists == T){
      if(is.na(methods) == F){
        for(c in (1:m)){
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if((is.na(mixed_dist) == F) & (length(mixed_dist) == 1)){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            
            clust = get(hc.func)(d, method = methods[c])
            clusters = factor(cutree(clust, k = nlvls))
            perms = permutations(n = nlvls, r = nlvls,v = new_lvls)
            hits = numeric(nrow(perms))
            clusts = matrix(0, nrow = nrow(perms), ncol = length(clusters))
            for(i in (1:nrow(perms))){
              new_clust = mapvalues(clusters, from = levels(clusters), to = perms[i, ])
              clusts[i, ] = new_clust
              temp_hit = (sum(new_clust == relab_y)/length(relab_y))
              hits[i] = temp_hit
            }
            max.index = which.max(hits)
            F1_s = F1_Score(relab_y, clusts[max.index, ])
            maximum_index = hits[max.index]
            all_hits = c(all_hits, maximum_index)
            all_F1 = c(all_F1, F1_s)
        }
        }else{
          # testando condi??es para ver se as vari?veis sao mistas ou nao
          if(length(bool[bool != T]) == 0){
            if((is.na(dists) == F) & (length(dists) == 1)){
              d = dist(scale(training_data), method = dists)
            }else{
              d = dist(scale(training_data))
            }
          }else{
            if(length(mixed_dist) == 1){
              d = distmix(scale(training_data), method = mixed_dist)
            }else{
              d = distmix(scale(training_data))
            }
          }
          clust = get(hc.func)(d)
          clusters = factor(cutree(clust, k = nlvls))
          perms = permutations(n = nlvls, r = nlvls,v = new_lvls)
          hits = numeric(nrow(perms))
          clusts = matrix(0, nrow = nrow(perms), ncol = length(clusters))
          for(i in (1:nrow(perms))){
            new_clust = mapvalues(clusters, from = levels(clusters), to = perms[i, ])
            clusts[i, ] = new_clust
            temp_hit = (sum(new_clust == relab_y)/length(relab_y))
            hits[i] = temp_hit
          }
          max.index = which.max(hits)
          F1_s = F1_Score(relab_y, clusts[max.index, ])
          maximum_index = hits[max.index]
          all_hits = c(all_hits, maximum_index)
          all_F1 = c(all_F1, F1_s)
        }
    }else{
      if(length(bool[bool != T]) == 0){
        func = "dist"
      }else{
        func = "distmix"
      }
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        if(is.na(methods) == F){
          for(c in (1:m)){
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d, method = methods[c])
              clusters = factor(cutree(clust, k = nlvls))
              perms = permutations(n = nlvls, r = nlvls,v = new_lvls)
              hits = numeric(nrow(perms))
              clusts = matrix(0, nrow = nrow(perms), ncol = length(clusters))
              for(i in (1:nrow(perms))){
                new_clust = mapvalues(clusters, from = levels(clusters), to = perms[i, ])
                clusts[i, ] = new_clust
                temp_hit = (sum(new_clust == relab_y)/length(relab_y))
                hits[i] = temp_hit
              }
              max.index = which.max(hits)
              F1_s = F1_Score(relab_y, clusts[max.index, ])
              maximum_index = hits[max.index]
              all_hits = c(all_hits, maximum_index)
              all_F1 = c(all_F1, F1_s)
            }
          }else{
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d)
              clusters = factor(cutree(clust, k = nlvls))
              perms = permutations(n = nlvls, r = nlvls,v = new_lvls)
              hits = numeric(nrow(perms))
              clusts = matrix(0, nrow = nrow(perms), ncol = length(clusters))
              for(i in (1:nrow(perms))){
                new_clust = mapvalues(clusters, from = levels(clusters), to = perms[i, ])
                clusts[i, ] = new_clust
                temp_hit = (sum(new_clust == relab_y)/length(relab_y))
                hits[i] = temp_hit
              }
              max.index = which.max(hits)
              F1_s = F1_Score(relab_y, clusts[max.index, ])
              maximum_index = hits[max.index]
              all_hits = c(all_hits, maximum_index)
              all_F1 = c(all_F1, F1_s)
            }
          }
  }
  }
  cross_val_results$hits = all_hits
  cross_val_results$F1 = all_F1
  return(cross_val_results)
}


L_cross_val_per_var = function(original_data, cl.list, dists = NA, mixed_dist = NA, tol = 1e-20, 
                               seed = 99, scale = F){
  list.names = names(cl.list)
  p = ncol(original_data)
  l = length(list.names)
  all_results = list()
  types = sapply(original_data, class)
  bool =  (types == "integer" | types == "numeric")
  if(length(bool[bool != T]) == 0){
    if(is.na(dists) == T | length(dists) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = dists
    }
  }else{
    if(is.na(mixed_dist) == T | length(mixed_dist) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = mixed_dist
    }
  }
  
  for(k in (1:l)){
    hc.func = list.names[k]
    methods = cl.list[[hc.func]]
    m = length(methods)
    if(no_dists == T){
      if(is.na(methods) == F){
        for(c in (1:m)){
          results = numeric(p)
          for(j in 1:p){
            test_data = as.data.frame(original_data[, j])
            colnames(test_data) = colnames(original_data)[j]
            rownames(test_data) = rownames(original_data)
            training_data = original_data[, -j]
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if((is.na(mixed_dist) == F) & (length(mixed_dist) == 1)){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            clust = get(hc.func)(d, method = methods[c])
            tree = convert_to_phylo(clust)
            results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                  silent = T, error=function(msg){
                                    return(NA)
                                  })
          }
          all_results[[paste0(hc.func, ".", methods[c])]] = results
        }}else{
          results = numeric(p)
          for(j in 1:p){
            test_data = as.data.frame(original_data[, j])
            colnames(test_data) = colnames(original_data)[j]
            rownames(test_data) = rownames(original_data)
            training_data = original_data[, -j]
            types = sapply(training_data, class)
            bool =  (types == "integer" | types == "numeric")
            
            
            # testando condi??es para ver se as vari?veis sao mistas ou nao
            if(length(bool[bool != T]) == 0){
              if((is.na(dists) == F) & (length(dists) == 1)){
                d = dist(scale(training_data), method = dists)
              }else{
                d = dist(scale(training_data))
              }
            }else{
              if(length(mixed_dist) == 1){
                d = distmix(scale(training_data), method = mixed_dist)
              }else{
                d = distmix(scale(training_data))
              }
            }
            
            clust = get(hc.func)(d)
            tree = convert_to_phylo(clust)
            results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                  silent = T, error=function(msg){
                                    return(NA)
                                  })
          }
          all_results[[hc.func]] = results
        }
    }else{
      if(length(bool[bool != T]) == 0){
        func = "dist"
      }else{
        func = "distmix"
      }
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        if(is.na(methods) == F){
          for(c in (1:m)){
            results = numeric(p)
            for(j in 1:p){
              test_data = as.data.frame(original_data[, j])
              colnames(test_data) = colnames(original_data)[j]
              rownames(test_data) = rownames(original_data)
              training_data = original_data[, -j]
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d, method = methods[c])
              tree = convert_to_phylo(clust)
              results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                    silent = T, error=function(msg){
                                      return(NA)
                                    })
            }
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] = c(results, 
                                                                                    used.dists[h])
          }}else{
            results = numeric(p)
            print(hc.func)
            for(j in 1:p){
              test_data = as.data.frame(original_data[, j])
              colnames(test_data) = colnames(original_data)[j]
              rownames(test_data) = rownames(original_data)
              training_data = original_data[, -j]
              d = get(func)(scale(training_data), method = used.dists[h])
              clust = get(hc.func)(d)
              tree = convert_to_phylo(clust)
              results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                    silent = T, error=function(msg){
                                      return(NA)
                                    })
            }
            all_results[[paste0(hc.func, ".", used.dists[h])]] = c(results, used.dists[h])
          }
      }
    }
  }
  if(no_dists == F){
    all_cvs = do.call(rbind, all_results)
    cvs_data.frame = as.data.frame(all_cvs)
    cvs_data.frame[, 1:p] <- sapply(cvs_data.frame, as.numeric)
    return(cvs_data.frame)
  }else{
    all_cvs = do.call(rbind, all_results)
    cvs_data.frame = as.data.frame(all_cvs)
    return(cvs_data.frame)
  }
}

# alternative version of L_cross_val_per_var

L_cross_val_per_var_alt = function(original_data, cl.list, dists = NA, mixed_dist = NA, tol = 1e-20, 
                               seed = 99, scale = F){
  list.names = names(cl.list)
  p = ncol(original_data)
  l = length(list.names)
  all_results = list()
  types = sapply(original_data, class)
  bool =  (types == "integer" | types == "numeric")
  if(length(bool[bool != T]) == 0){
    if(is.na(dists) == T | length(dists) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = dists
    }
  }else{
    if(is.na(mixed_dist) == T | length(mixed_dist) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = mixed_dist
    }
  }
  
  for(k in (1:l)){
    hc.func = list.names[k]
    methods = cl.list[[hc.func]]
    m = length(methods)
    if(no_dists == T){
      if(is.na(methods) == F){
        for(c in (1:m)){
          results = numeric(p)
          types = sapply(original_data, class)
          bool =  (types == "integer" | types == "numeric")
          
          # testando condi??es para ver se as vari?veis sao mistas ou nao
          if(length(bool[bool != T]) == 0){
            if((is.na(dists) == F) & (length(dists) == 1)){
              d = dist(scale(original_data), method = dists)
            }else{
              d = dist(scale(original_data))
            }
          }else{
            if((is.na(mixed_dist) == F) & (length(mixed_dist) == 1)){
              d = distmix(scale(original_data), method = mixed_dist)
            }else{
              d = distmix(scale(original_data))
            }
          }
          clust = get(hc.func)(d, method = methods[c])
          tree = convert_to_phylo(clust)
          for(j in 1:p){
            test_data = as.data.frame(original_data[, j])
            colnames(test_data) = colnames(original_data)[j]
            rownames(test_data) = rownames(original_data)
            results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                  silent = T, error=function(msg){
                                    return(NA)
                                  })
          }
          all_results[[paste0(hc.func, ".", methods[c])]] = 1 - (results)
        }}else{
          results = numeric(p)
          types = sapply(original_data, class)
          bool =  (types == "integer" | types == "numeric")
          
          
          # testando condi??es para ver se as vari?veis sao mistas ou nao
          if(length(bool[bool != T]) == 0){
            if((is.na(dists) == F) & (length(dists) == 1)){
              d = dist(scale(original_data), method = dists)
            }else{
              d = dist(scale(original_data))
            }
          }else{
            if(length(mixed_dist) == 1){
              d = distmix(scale(original_data), method = mixed_dist)
            }else{
              d = distmix(scale(original_data))
            }
          }
          
          clust = get(hc.func)(d)
          tree = convert_to_phylo(clust)
          for(j in 1:p){
            test_data = as.data.frame(original_data[, j])
            colnames(test_data) = colnames(original_data)[j]
            rownames(test_data) = rownames(original_data)
            results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                  silent = T, error=function(msg){
                                    return(NA)
                                  })
          }
          all_results[[hc.func]] = 1 - results
        }
    }else{
      if(length(bool[bool != T]) == 0){
        func = "dist"
      }else{
        func = "distmix"
      }
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        if(is.na(methods) == F){
          for(c in (1:m)){
            results = numeric(p)
            d = get(func)(scale(original_data), method = used.dists[h])
            clust = get(hc.func)(d, method = methods[c])
            tree = convert_to_phylo(clust)
            for(j in 1:p){
              test_data = as.data.frame(original_data[, j])
              colnames(test_data) = colnames(original_data)[j]
              rownames(test_data) = rownames(original_data)
              results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                    silent = T, error=function(msg){
                                      return(NA)
                                    })
            }
            all_results[[paste0(hc.func, ".", methods[c], ".", used.dists[h])]] = c(1 - results, 
                                                                                    used.dists[h])
          }}else{
            results = numeric(p)
            d = get(func)(scale(original_data), method = used.dists[h])
            clust = get(hc.func)(d)
            tree = convert_to_phylo(clust)
            for(j in 1:p){
              test_data = as.data.frame(original_data[, j])
              colnames(test_data) = colnames(original_data)[j]
              rownames(test_data) = rownames(original_data)
              results[j] = tryCatch(L_score(tree, test_data, scale = scale),  
                                    silent = T, error=function(msg){
                                      return(NA)
                                    })
            }
            all_results[[paste0(hc.func, ".", used.dists[h])]] = c(1 - results, 
                                                                used.dists[h])
          }
      }
    }
  }
  if(no_dists == F){
    all_cvs = do.call(rbind, all_results)
    cvs_data.frame = as.data.frame(all_cvs)
    cvs_data.frame[, 1:p] <- sapply(cvs_data.frame, as.numeric)
    return(cvs_data.frame)
  }else{
    all_cvs = do.call(rbind, all_results)
    cvs_data.frame = as.data.frame(all_cvs)
    return(cvs_data.frame)
  }
}

