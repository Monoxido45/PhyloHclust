# testing the silhoute score and cophenetic score for all data
# loading some data
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
library(mvMORPH)
library(kmed)
library(gtools)
library(MLmetrics)
library(tictoc)
library(RColorBrewer)
library(MASS)
library(plyr)
library(fossil)
source("modules/convert_to_parenthesis.R")
source("modules/cv_score.R")


data(ruspini)
pr4 <- pam(ruspini, 4)
si <- silhouette(pr4)

# Loading all  data to test new losses ------------------------------------
# simulated data
data.generator = function(n, mu_0, mu_1, p = 2, S = NULL, seed = 500){
  set.seed(seed)
  Y = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  if(is.null(S) == T) S = diag(nrow = p, ncol = p)
  X = matrix(nrow = n, ncol = p)
  X[Y == 0, ] = round(mvrnorm(sum(Y == 0), mu_0, S), 4)
  X[Y == 1, ] = round(mvrnorm(sum(Y == 1), mu_1, S), 4)
  
  sim.data = as.data.frame(cbind(X , Y))
  sim.data$Y = as.factor(Y)
  if(p == 2){
    colnames(sim.data) = c("X1", "X2", "Y")}
  row.names(sim.data) = c(1:nrow(sim.data))
  return(sim.data)
}

n = 250
p = 2
mu_0 = c(0, 0)
mu_1 = c(2, 2)
sim.data = data.generator(n, mu_0, mu_1, seed = 150)
head(sim.data, 4)

# iris dataset
flowers = iris

# pima indian diabetes dataset
n = 350
library(data.table)
prima_data <- as.data.frame(fread('https://raw.githubusercontent.com/jbrownlee/Datasets/master/pima-indians-diabetes.csv'))
selected_rows = sample(1:nrow(prima_data), n, replace = F)
# selecting rows (too many samples)
prima_data$V9 = as.factor(prima_data$V9)
prima_data = prima_data[selected_rows, ]
row.names(prima_data) = 1:nrow(prima_data)
head(prima_data)



# wheat seeds dataset
wheat_data = read.delim("../data/seeds_dataset.txt")
wheat_data$X1 = as.factor(wheat_data$X1)
head(wheat_data)

# ionosphere dataset
ionosphere_data = as.data.frame(read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data", sep = ","))
ionosphere_data$V35 = as.factor(ionosphere_data$V35)

# droping second variable
ionosphere_data = ionosphere_data[, -2]
pca_ionosphere = prcomp(ionosphere_data[, -34], center = TRUE,scale. = TRUE)
summary(pca_ionosphere)

# selecting first 12 componentes
ionosphere_components = as.data.frame((pca_ionosphere$x)[, 1:12])
ionosphere_components$label = ionosphere_data[, 34]
head(ionosphere_components)

# glass data
glass.data = as.data.frame(read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/glass/glass.data", sep = ","))
glass.data = glass.data[, -1]
glass.data$V11 = as.factor(glass.data$V11)
head(glass.data)

# haberman data
haberman.data = as.data.frame(read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/haberman/haberman.data", sep = ","))
haberman.data$V4 = as.factor(haberman.data$V4)
head(haberman.data)

# wine data
wine.data = as.data.frame(read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data", sep = ","))
wine.data$V1 = as.factor(wine.data$V1)
head(wine.data)



# Creating all scores for all possible combination of methods --------
# all combinations
test.list = list(hclust = c("ward.D", "single","ward.D2", "average", "complete", "mcquitty", "median"), 
                 agnes = c("weighted", "average", "ward"), diana = NA)
dists = c("euclidean", "manhattan", "canberra")

# function to display all correlations and results
display_silhouette_coph = function(data, clust.list, test_index, 
                                  dist = NA, scale = F){
  # adjusting data
  training_data = data[, -test_index]
  test_data = data[, test_index]
  
  coph_sil_res <- silhouette_coph_val(training_data, test_data, clust.list, dists = dist, 
                                      mixed_dist = mixed_dist)

  cat("Correlations between silhouette score and F1: \n")
  cat("silhouette vs F1: ", cor(coph_sil_res$silhouette, coph_sil_res$f1, use = "complete.obs",
                                method = "spearman"), "\n")
  
  cat("Correlations between Dunn index and F1: \n")
  cat("cophenetic vs F1: ", cor(coph_sil_res$dun, coph_sil_res$f1,  use = "complete.obs",
                                method = "spearman"), "\n")
  
  cat("Correlations between connectivity and F1: \n")
  cat("cophenetic vs F1: ", cor(coph_sil_res$con, coph_sil_res$f1,  use = "complete.obs",
                                method = "spearman"), "\n")
  
  cat("Correlations between cophenetic score and F1: \n")
  cat("cophenetic vs F1: ", cor(coph_sil_res$coph, coph_sil_res$f1, use = "complete.obs",
                                method = "spearman"), "\n")
  
  cat("Correlations between Rand index and F1: \n")
  cat("cophenetic vs F1: ", cor(coph_sil_res$rand, coph_sil_res$f1,  use = "complete.obs",
                                method = "spearman"), "\n")
  return(coph_sil_res)
}

silhouette_coph_val <- function(training_data, test_data, clust_list, 
                           dists, mixed_dist = mixed_dist){
  # levels
  lvls = levels(test_data)
  nlvls= length(lvls)
  types = sapply(training_data, class)
  relab_y = mapvalues(test_data, from = lvls, to = 1:nlvls)
  new_lvls = levels(relab_y)
  
  list.names = names(clust_list)
  all_F1 = numeric(0)
  all_sil = numeric(0)
  all_coph = numeric(0)
  all_rand = numeric(0)
  all_dun = numeric(0)
  all_con = numeric(0)
  
  types = sapply(training_data, class)
  bool =  (types == "integer" | types == "numeric")
  
  if(length(bool[bool != T]) == 0){
    if(any(is.na(dists) == T) | length(dists) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = dists
    }
  }else{
    if(any(is.na(mixed_dist) == T) | length(mixed_dist) == 1){
      no_dists = T
    }else{
      no_dists = F
      used.dists = mixed_dist
    }
  }
  
  # creating the clusters
  for(k in 1:length(list.names)){
    hc.func = list.names[k]
    methods = clust_list[[hc.func]]
    m = length(methods)
    if(no_dists == T){
      if(is.na(methods) == F){
        for(c in (1:m)){
          # testing if variables are mixed or not
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
          clusts = matrix(0, nrow = nrow(perms), ncol = length(clusters))
          
          for(i in (1:nrow(perms))){
            new_clust = mapvalues(clusters, from = levels(clusters), to = perms[i, ])
            clusts[i, ] = new_clust
            temp_hit = (sum(new_clust == relab_y)/length(relab_y))
            hits[i] = temp_hit
          }
          max.index = which.max(hits)
          # computing F1 in the indexes that maximizes proportion of hits
          F1_s = F1_Score(relab_y, clusts[max.index, ])
          
          # silhouette score
          sil_obj <- silhouette(clusts[max.index, ], d) |> summary()
          ave_sil <- -sil_obj$avg.width
          
          # dunn index
          dun <- -clValid::dunn(d, clusts[max.index, ])
          
          # connectiviy
          con <- clValid::connectivity(d, clusts[max.index, ])
          
          # rand index
          y <- relab_y |> as.character() |> as.numeric()
          rand_val <- -adj.rand.index(y, clusts[max.index, ])
          
          # cophenetic
          coph_obj = cophenetic(as.hclust(clust))
          coph_val = -cor(d, coph_obj)
          
          all_F1 = c(all_F1, F1_s)
          all_dun = c(all_dun, dun)
          all_sil = c(all_sil, ave_sil)
          all_rand = c(all_rand, rand_val)
          all_coph = c(all_coph, coph_val)
          all_con = c(con, all_con)
        }
      }else{
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
        # computing F1 in the indexes that maximizes proportion of hits
        F1_s = F1_Score(relab_y, clusts[max.index, ])
        
        # silhouette score
        sil_obj <- silhouette(clusts[max.index, ], d) |> summary()
        ave_sil <- -sil_obj$avg.width
        
        # dunn index
        dun <- -clValid::dunn(d, clusts[max.index, ])
        
        # connectiviy
        con <- clValid::connectivity(d, clusts[max.index, ])
        
        # rand index
        y <- relab_y |> as.character() |> as.numeric()
        rand_val <- -adj.rand.index(y, clusts[max.index, ])
        
        # cophenetic
        coph_obj = cophenetic(as.hclust(clust))
        coph_val = -cor(d, coph_obj)
        
        all_F1 = c(all_F1, F1_s)
        all_sil = c(all_sil, ave_sil)
        all_coph = c(all_coph, coph_val)
        all_dun = c(all_dun, dun)
        all_rand = c(all_rand, rand_val)
        all_con = c(con, all_con)
      }
    
      }else{
      if(length(bool[bool != T]) == 0){
        func = "dist"
      }else{
        func = "distmix"
      }
      dists.length = length(used.dists)
      for(h in 1:dists.length){
        if(any(is.na(methods) == F)){
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
            
            # computing F1 in the indexes that maximizes proportion of hits
            F1_s = F1_Score(relab_y, clusts[max.index, ])
            
            # silhouette score
            sil_obj <- silhouette(clusts[max.index, ], d) |> summary()
            ave_sil <- -sil_obj$avg.width
            
            # cophenetic
            coph_obj = cophenetic(as.hclust(clust))
            coph_val = -cor(d, coph_obj)
            
            # dunn index
            dun <- -clValid::dunn(d, clusts[max.index, ])
            
            # connectiviy
            con <- clValid::connectivity(d, clusts[max.index, ])
            
            # rand index
            y <- relab_y |> as.character() |> as.numeric()
            rand_val <- -adj.rand.index(y, clusts[max.index, ])
            
            all_F1 = c(all_F1, F1_s)
            all_sil = c(all_sil, ave_sil)
            all_coph = c(all_coph, coph_val)
            all_rand = c(all_rand, rand_val)
            all_dun = c(all_dun, dun)
            all_con = c(all_con, con)
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
          # computing F1 in the indexes that maximizes proportion of hits
          F1_s = F1_Score(relab_y, clusts[max.index, ])
          
          # silhouette score
          sil_obj <- silhouette(clusts[max.index, ], d) |> summary()
          ave_sil <- -sil_obj$avg.width
          
          # dunn index
          dun <- -clValid::dunn(d, clusts[max.index, ])
          
          # connectiviy
          con <- clValid::connectivity(d, clusts[max.index, ])
          
          # cophenetic
          coph_obj = cophenetic(as.hclust(clust))
          coph_val = -cor(d, coph_obj)
          
          # rand index
          y <- relab_y |> as.character() |> as.numeric()
          rand_val <- -adj.rand.index(y, clusts[max.index, ])
          
          all_F1 = c(all_F1, F1_s)
          all_sil = c(all_sil, ave_sil)
          all_coph = c(all_coph, coph_val)
          all_dun = c(all_dun, dun)
          all_rand = c(all_rand, rand_val)
          all_con = c(all_con, con)
        }
      }
    }
  }
  results_data = data.frame(silhouette = all_sil,
                            f1 = all_F1,
                            dun = all_dun,
                            con = all_con,
                            coph = all_coph,
                            rand = all_rand)
  return(results_data)
}




# testing for all data ----------------------------------------------------
# simulated data
coph_data <- display_silhouette_coph(sim.data, test.list, test_index = 3, dist = dists)

iris_exp <- display_silhouette_coph(flowers, test.list, test_index = 5, dist = dists)

pima_exp <- display_silhouette_coph(prima_data, test.list, test_index = 9, dist = dists)

wheat_exp = display_silhouette_coph(wheat_data, test.list, test_index = 8, dist = dists)

ionosphere_exp = display_silhouette_coph(ionosphere_components, test.list, test_index = 13, 
                                           dist = dists)

glass_exp = display_silhouette_coph(glass.data, test.list, test_index = 10, 
                                     dist = dists)

haberman_exp = display_silhouette_coph(haberman.data, test.list, test_index = 4, 
                                        dist = dists)

wine_exp = display_silhouette_coph(wine.data, test.list, test_index = 1, 
                                    dist = dists)
