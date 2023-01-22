# code to make graphical abstract and boxplot/evolutionary dendrogram comparisson
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

# making simulated dataset for graphical abstract
# simulating dataset from multivariate normal
# according to a latent variable with 3 levels
set.seed(1275)
n = 15
y = sample(c(1, 2, 3), size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3))
mu_1 = c(1, 1, 1, 3)
mu_2 = c(4, -2, 1, 3)
mu_3 = c(-2, -2, 1, 3)
S = rbind(c(1, 0.5, 0.15, 0.15), c(0.5, 1, 0.15, 0.15), 
          c(0.15, 0.15, 1, 0.5), c(0.15, 0.15, 0.5, 1))
X = ((y == 1)*(mvrnorm(15, mu_1, S)) + (y == 2)*(mvrnorm(15, mu_2, S)) +
  (y == 3)*(mvrnorm(15, mu_2, S))) %>% 
  as.data.frame()
colnames(X) = c("X1", "X2", "X3", "X4")

# generating 3 different dendrograms with euclidean distance
ward_tree = hclust(dist(scale(X)), method = "ward.D2") %>%
  convert_to_phylo()
comp_tree = hclust(dist(scale(X)), method = "complete") %>%
  convert_to_phylo()
mcquitty_tree = hclust(dist(scale(X)), method = "mcquitty") %>%
  convert_to_phylo()


par(mfrow = c(1, 3))
plot(ward_tree,cex = 1)

plot(comp_tree,cex = 1)

plot(mcquitty_tree, cex = 1)

test.list = list(hclust = c("ward.D2", "complete","mcquitty"))
dists = c("euclidean")

tol = 1e-20

scores_arr = L_cross_val(X, cl.list = test.list, dists = dists,
                         tol = tol, seed = 99, scale =  T)

# importance score
test.list = list(hclust = c("mcquitty"))
dists = c("euclidean")
import_x = L_cross_val_per_var_alt(X, test.list, 
                        dists, scale = T)

melted_sim = reshape2::melt(import_x)
melted_sim$variable = as.factor(colnames(X))
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



# evolutionary dendrograms for these two cases:
par(mfrow = c(1, 2))

# ploting X2 first
x = setNames(X$X2,gsub(" ", "", rownames(X)))
reordered_x = x[mcquitty_tree$tip.label]

obj = contMap(mcquitty_tree, reordered_x, plot=FALSE, cols =  viridis::viridis(15))
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize= 0.5,outline=FALSE, lwd = c(2,5), leg.txt = "X2")

# X1 now
x = setNames(X$X1,gsub(" ", "", row.names(X)))
reordered_x = x[mcquitty_tree$tip.label]

obj = contMap(mcquitty_tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize= 0.5,outline=FALSE, lwd = c(2,5), leg.txt = "X1")

# generating boxplots for USArrest according to dendrogram partitions
arr_dend = hclust(dist(scale(USArrests), method = "manhattan"), method = "mcquitty")
# scaling and secting several partitions
arr_data = USArrests
arr_data %<>%
  scale() %<>%
  as.data.frame() %<>%
  mutate("k = 2" = as.factor(cutree(arr_dend, k = 2)),
         "k = 3" = as.factor(cutree(arr_dend, k = 3)),
         "k = 4" = as.factor(cutree(arr_dend, k = 4)))


melted_data = arr_data %>%
  pivot_longer(Murder, names_to = "variable", values_to = "values") %>%
  pivot_longer("k = 2":"k = 4", names_to = "k", values_to = "cluster")


melted_data %>%
  ggplot(aes(y = variable, x = values, fill = cluster)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Standardized values",
       y = "",
       fill = "Cluster") +
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~k)

mcquitty.tree = convert_to_phylo(arr_dend)
par(mfrow = c(1,1))
x = setNames(scale(arr_data[,1]),gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]
obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd=c(3,7),leg.txt="Murder", legend = FALSE)
nodelabels(text= c("A", "B", "C", "D", "E", "F"),node= c(52, 82, 53, 59, 83, 89))

# comparing our importance score to random forest importance
# using wheat seed dataset
# in two different situations concerning clustering
# importing wheat seeds again
setwd("~/estatistica_UFSCAR/cv_cluster")
wheat_data = read.delim("data/seeds_dataset.txt",
                        header = F)
wheat_data$V8 = as.factor(wheat_data$V8)
head(wheat_data)
dados_quant = wheat_data[, -8]
# obtaining dendrogram
wheat_dend = wheat_data %>%
  dplyr::select(-V8) %>%
  scale() %>%
  dist() %>%
  hclust(method = "ward.D2")

# good partition and bad partition
wheat_data %<>%
  mutate(clust3 = as.factor(cutree(wheat_dend, k = 3)),
         clust6 = as.factor(cutree(wheat_dend, k = 6)))




# applying random forest
rf_clust_3 = ranger::ranger(formula = clust3 ~ . - clust6 - V8,
                            data = wheat_data,
                            num.trees = 500,
                            importance = "impurity",
                            write.forest = TRUE,
                            verbose = FALSE,
                            probability = T)

rf_clust_6 = ranger::ranger(formula = clust6 ~ . - clust3 - V8,
                            data = wheat_data,
                            num.trees = 500,
                            importance = "impurity",
                            write.forest = TRUE,
                            verbose = FALSE,
                            probability = T)


test.list = list(hclust = c("ward.D2"))
dists = c("euclidean")

import_x = L_cross_val_per_var_alt(dados_quant, test.list, 
                                   dists, scale = T)

hcpl_import = import_x %>% reshape2::melt()

hcpl_import = setNames(hcpl_import$value, hcpl_import$variable)

ggplot_data = reshape2::melt(import_x)
ggplot_data$variable = as.factor(colnames(X))


import = tibble(variable = c(names(ranger::importance(rf_clust_3)), 
                             names(ranger::importance(rf_clust_6)),
                             names(hcpl_import)),
importance = c(ranger::importance(rf_clust_3),
               ranger::importance(rf_clust_6), 
               hcpl_import)) %>%
  mutate(cluster = as.factor(rep(c("RF k = 3", "RF k = 6", "Our score"),
                                 each = dim(wheat_data)[2] - 3)))


import %>%
  ggplot(aes(x = variable,
             y = importance))+
  geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") + 
  coord_flip()+
  labs(y = "Importance score",
       x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15,
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~cluster, scales = "free_x")



# plotting evolutionary dendrogram only for V6 and V4
# dendrogram for V6
ward.tree = wheat_dend %>% convert_to_phylo()
par(mfrow = c(1,2))
x = round(setNames(dados_quant$V4,gsub(" ", "", rownames(dados_quant))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.4,0.6),outline=FALSE,lwd = c(2,5), leg.txt="V4", ftype = "off")

# dendrogram for V4
x = round(setNames(dados_quant$V6, gsub(" ", "", rownames(dados_quant))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.4,0.6),outline=FALSE,lwd = c(2,5), leg.txt="V6", ftype = "off")


test.list = list(hclust = c("ward.D2"))
dists = c("euclidean")
import_x = L_cross_val_per_var_alt(dados_quant, test.list, 
                                   dists, scale = T)

melted_sim = reshape2::melt(import_x)
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

ggplot_data = reshape2::melt(import_x)
ggplot_data$variable = as.factor(colnames(X))
p1 = ggplot_data %>%
  mutate(variable = fct_reorder(variable, value, .desc = T)) %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3")+
  theme_minimal() +
  labs(x = "Features",
       y = "Importance score") +
  coord_flip()+
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))
p1

