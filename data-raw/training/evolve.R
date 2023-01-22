# simulating data from bayesian network model
set.seed(2075)
d <- 6
k <- 8
cov_stab <- seq(from = 0.5, to = 1.5, length.out = 6)
evolve <- function(obs, time) 
  obs + rnorm(d, 0, cov_stab^time)

param <- list()
param[[1]] <- matrix(0, nrow = 1, ncol = d)
for(time in 1:(k-1))
{
  param[[time + 1]] <- matrix(NA, nrow = 2^time, ncol = d)
  n_param_ant = 2^(time-1)
  for(jj in 1:n_param_ant)
  {
    mut_1 = evolve(param[[time]][jj,], time)
    mut_2 = evolve(param[[time]][jj,], time)
    param[[time + 1]][jj,] = mut_1
    param[[time + 1]][jj + n_param_ant,] = mut_2
  }
}


sim_data <- param[[k]] |>
  as.data.frame()


# importing packages
library(ape)
library(dendextend)
library(cluster)
library(tidyverse)
library(phytools)
library(mltools)
library(data.table)
library(factoextra)
setwd("~/estatistica_UFSCAR/cv_cluster/modules")
source("convert_to_parenthesis.R")
source("cv_score.R")
library(tictoc)
library(mvMORPH)
library(RColorBrewer)
library(MASS)
library(plyr)
library(clValid)

# implementing variable importance score from badih also called lovo score
lovo <- function(data, hclust_type = "ward.D2", dist_type = "euclidean", 
                 k = 3){
  p <- ncol(data)
  
  1:p |>
    map_dbl(function(x){
      new_data <- data[, -x]
      n <- nrow(data)
      
      # obtaining partition for leave one out
      part <- new_data |> 
        scale() |>
        dist(method = dist_type) |>
        hclust(method = hclust_type) |>
        cutree(k = k)
      
      new_data <- new_data |>
        mutate(part = part)
      
      # within cluster heterogeneity for variable p
      lovo <- 1:k |>
        map_dbl(function(t){
          trace_t <- new_data |> 
            filter(part == t) |>
            var() |>
            diag() |>
            sum()
          
          nt <- new_data |> 
            filter(part == t) |>
            nrow()
          
          nt/n * trace_t
        }) |>
        mean()
      
      return(lovo)
    })
}


test.list = list(hclust = c("ward.D2"))
dists = c("euclidean")
# variable importance by our method
import_sim = L_cross_val_per_var_alt(sim_data, test.list, 
                                   dists, scale = T)

# variable importance by lovo
lovo_sim <- lovo(sim_data)

# variable importance by XGboost
# current dendrogram
dend_ward = sim_data %>%
  scale() %>%
  dist() %>%
  hclust(method = "ward.D2")

sim_data %<>%
  mutate(clust3 = as.factor(cutree(dend_ward, k = 3)))

data_matrix <- model.matrix(clust3 ~ ., data = sim_data)[,-1]

xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = 3)

xgb <- xgboost::xgboost(
  params = xgb_params,
  data = data_matrix, 
  label = as.numeric(sim_data$clust3) - 1,
  max.depth = 2,
  nrounds = 50)

# obtaining importance score
importance_matrix <- xgboost::xgb.importance(colnames(data_matrix),
                                   model = xgb) |>
  arrange(Feature)




# simulated data with dendrogram
sim_data %<>%
  mutate(clust3 = as.factor(cutree(dend_ward, k = 3)))




melted_sim = reshape2::melt(import_sim)
melted_sim$variable = as.factor(colnames(import_sim))
p1 = melted_sim %>%
  mutate(variable = fct_reorder(variable, value, .desc = T)) %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Features",
       y = "Importance score") +
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(min(melted_sim$value) - 0.05, max(melted_sim$value) + 0.05)

p1

# RF
nb = NbClust::NbClust(sim_data, distance = "euclidean",
                      min.nc = 2, max.nc = 10, method = "ward.D2")

factoextra::fviz_nbclust(nb) +
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))


rf_clust_3 = ranger::ranger(formula = clust3 ~ .,
                            data = sim_data,
                            num.trees = 500,
                            importance = "impurity",
                            write.forest = TRUE,
                            verbose = FALSE,
                            probability = T)

import = tibble(variable = c(names(ranger::importance(rf_clust_3)), 
                             names(ranger::importance(rf_clust_3)),
                             names(ranger::importance(rf_clust_3)),
                             names(ranger::importance(rf_clust_3)),
                             names(ranger::importance(rf_clust_3))),
                importance = c(ranger::importance(rf_clust_3),
                               importance_matrix$Gain,
                               lovo_sim,
                               melted_sim$value,
                               1/cov_stab)) %>%
  mutate(cluster = as.factor(rep(c("RF", "XGboost", "LOVO" ,"PFIS", "Truth"), 
                                 each = dim(sim_data)[2] - 1)))

import$cluster_f = factor(import$cluster, levels=c('RF', "XGboost",
                                                   'LOVO', 'PFIS','Truth'))

facet_names <- c(
  `RF` = "RF",
  `XGboost` = "XGboost",
  `LOVO` = "LOVO",
  `PFIS` = expression("PFIS"~ "(Our approach)"),
  `Truth` = expression("Ground truth (" ~ sigma[j]^{-1} ~ ")")
)

import = mutate_at(import, .vars = "cluster_f",
                   .funs = factor, labels = facet_names)

import %>%
  ggplot(aes(x = variable, y  = importance, group = 1)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Features",
       y = "Importance score") +
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~cluster_f, scales = "free_y", ncol = 5,
             labeller = label_parsed)

dend_ward %<>% convert_to_phylo()

par(mfrow= c(2, 3))
x = setNames(round(sim_data[, 1], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.75), outline=FALSE, lwd = c(2,5), leg.txt="X1", ftype = "off")


x = setNames(round(sim_data[, 2], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.75), outline=FALSE, lwd = c(2,5), leg.txt="X2", ftype = "off")


x = setNames(round(sim_data[, 3], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.75), outline=FALSE, lwd = c(2,5), leg.txt="X3", ftype = "off")


x = setNames(round(sim_data[, 4], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.75), outline=FALSE, lwd = c(2,5), leg.txt="X4", ftype = "off")

x = setNames(round(sim_data[, 5], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.75), outline=FALSE, lwd = c(2,5), leg.txt="X5", ftype = "off")


x = setNames(round(sim_data[, 6], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.75), outline=FALSE, lwd = c(2,5), leg.txt="X6", ftype = "off")

