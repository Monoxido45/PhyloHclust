# code exploring USArrests using our score
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

# importing USArrests
data("USArrests")
arr.data = USArrests
row.names(arr.data) = 1:nrow(arr.data)

# list of possible combinations:
# removing median and centroid linkage
test.list = list(hclust = c("ward.D", "single","ward.D2", "complete", "mcquitty", "average"), diana = NA)
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20

scores_arr = L_cross_val(arr.data, cl.list = test.list, dists = dists,
                         tol = tol, seed = 99, scale =  T)

id = which.min(scores_arr$V1)
min_id = scores_arr[id,]

# ordenando os scores
scores_arr %<>% mutate(names = row.names(scores_arr))
scores_arr %>% arrange(V1) %>% head(8)

arr_dend = hclust(dist(scale(USArrests), method = "manhattan"), method = "mcquitty")

# selecting the best 
mcquitty.tree = convert_to_phylo(arr_dend)

# dendrogram for murder
par(mfrow = c(1,2))
x = setNames(scale(arr.data[,1]),gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd=c(3,7),leg.txt="Murder")

# dendrogram for assault
x = setNames(scale(arr.data[,2]),gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd=c(3,7),leg.txt="Assault")

# dendrogram for Urbanpop
x = setNames(scale(arr.data[,3]),gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd=c(3,7),leg.txt="Urbanpop")

# dendrogram for rape
x = setNames(scale(arr.data[,4]),gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6, 0.8),outline=FALSE,lwd=c(3,7),leg.txt="Rape")

# importance by variable
test.list = list(hclust = c("mcquitty"))
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20
L_per_var = L_cross_val_per_var(arr.data, test.list, dists, scale = T)

ggplot_data = reshape2::melt(L_per_var[3, -5])
ggplot_data$variable = as.factor(colnames(arr.data))
p1 = ggplot_data %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Cross validated score values",
       title = "Scatterplot of cross validated score values versus 
       variables from USArrests dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0, max(ggplot_data$value) + 0.05)
p1


# correlation matrix
library(ggcorrplot)
ggcorrplot(cor(arr.data), hc.order = TRUE, type = "lower",
           outline.col = "white", lab = T) +
  ggsave(filename  = "corr_usa_arrests.pdf",
path = "C:/Users/lucru/Estat?stica_UFSCar/cv_cluster/figures",
width = 20.75, height = 12.5, units = "cm")


# using an alternative version of importance by variable
test.list = list(hclust = c("mcquitty"))
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20
L_per_var = L_cross_val_per_var_alt(arr.data, test.list, dists, scale = T)

ggplot_data = reshape2::melt(L_per_var[3, -5])
ggplot_data$variable = as.factor(colnames(arr.data))
p1 = ggplot_data %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Importance value",
       title = "Scatterplot of variable importance for USArrests dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0, max(ggplot_data$value) + 0.05)
p1


# adding and testing wheat seeds dataset
# wheat seeds dataset
wheat_data = read.delim("data/seeds_dataset.txt",
                        header = F)
wheat_data$V8 = as.factor(wheat_data$V8)
head(wheat_data)


# using an alternative version of importance by variable
test.list = list(hclust = c("mcquitty"))
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20
L_per_var_wheat = L_cross_val_per_var_alt(wheat_data[, -8], test.list, dists, scale = T)

ggplot_data = reshape2::melt(L_per_var_wheat[3, -8])
ggplot_data$variable = as.factor(colnames(wheat_data[,-8]))
p1 = ggplot_data %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
labs(x = "Variables names",
       y = "Cross validated score values",
       title = "Wheat seeds dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0, max(ggplot_data$value) + 0.05)
p1


# variable alternative importance for simulated dataset and wheat seeds
# new test.list with more methods
test.list = list(hclust = c("ward.D", "single","ward.D2", "average", "complete", "mcquitty", "median"), 
                 agnes = c("weighted", "average", "ward"), diana = NA)

# wheat seeds dataset
wheat.l_cross_per_var = L_cross_val_per_var_alt(wheat_data[, -8], test.list, dists, scale = T)
  

# simulated dataset
# simulating in a different way
mu_0 = c(2, 1)
mu_1 = c(8, 1)
S = rbind(c(1, 0.5), c(0.5,1))
n = 225
set.seed(600)

Y = sample(c(0,1), n, replace = T, prob = c(0.5, 0.5))
X = matrix(nrow = n, ncol = 2)
X[Y == 0,] = round(mvrnorm(sum(Y == 0), mu_0, S), 4)
X[Y == 1,] = round(mvrnorm(sum(Y == 1), mu_1, S), 4)
sim.data = as.data.frame(cbind(X, Y))
sim.data$Y = as.factor(Y)
colnames(sim.data) = c("X1", "X2", "Y")


sim.data.l_cross_per_var_scaled = L_cross_val_per_var_alt(sim.data[, -3], test.list, dists, scale = T)

sim.data.l_cross_per_var_scaled2 = L_cross_val_per_var(sim.data[, -3], test.list, dists, scale = T)



nb.cols = 33
mycolors = colorRampPalette(brewer.pal(33, "Paired"))(nb.cols)

ggplot_data = reshape2::melt(wheat.l_cross_per_var[, -8])
ggplot_data$dists = rep(wheat.l_cross_per_var$V8, 7)
ggplot_data$obs = rep(as.factor(1:33), 7)

library(forcats)

p1 = ggplot_data %>%
  mutate(variable = as.factor(variable),
         dists = as.factor(dists)) %>%
  mutate(variable = fct_reorder(variable, value, .fun = 'median', .desc = T)) %>%
  ggplot(aes(x = variable, y  = value, group = obs)) +
  geom_line(aes(color = obs)) +
  geom_point(aes(color = obs)) +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Importance score values",
       colour = "Combinations",
       title = "Wheat seeds dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
       plot.title = element_text(hjust = 0.5))
p1


ggplot_data = reshape2::melt(sim.data.l_cross_per_var_scaled[, -3])
library(plyr)
ggplot_data$dists = rep(sim.data.l_cross_per_var_scaled$V3, 2)
ggplot_data$obs = rep(as.factor(1:33), 2)

p2 = ggplot_data %>%
  mutate(variable = revalue(as.factor(variable), c("V1" = "X1", "V2" = "X2")),
         dists = as.factor(dists)) %>%
  mutate(variable = fct_reorder(variable, value, .fun = 'median', .desc = T)) %>%
  ggplot(aes(x = variable, y  = value, group = obs)) +
  geom_line(aes(color = obs)) +
  geom_point(aes(color = obs)) +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Importance score values",
       colour = "Combinations",
       title = "Simulated dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))
p2 + 
  ggsave(filename  = "simulated_dataset_importances.pdf",
         path = "C:/Users/lucru/Estat?stica_UFSCar/cv_cluster/figures",
         width = 20.75, height = 14.5, units = "cm")

library(ggpubr)
ggarrange(p1, p2, common.legend = T, nrow = 2, legend = "right")+
  ggsave(filename  = "wheat_seeds_simulated_importances_v2.pdf",
         path = "C:/Users/lucru/Estat?stica_UFSCar/cv_cluster/figures",
         width = 22.75, height = 14.5, units = "cm")

# clustering by ward.D2 in wheat seeds
wheat_data_temp = wheat_data[,-8]
wheat_dend = hclust(dist(scale(wheat_data_temp), method = "euclidean"),
                    method = "ward.D2")

wheat_data_temp %<>%
  scale() %<>%
  as.data.frame() %<>%
  mutate(clust = as.factor(cutree(wheat_dend, k = 3)))

# scaling for boxplot
dados_quant = wheat_data_temp[,-8]
ggplot_data = melt(dados_quant)

# boxplot separating by ward.d2
ggplot_data$factor = rep(wheat_data_temp$clust, 7)
ggplot_data %>%
  ggplot(aes(y = variable, x = value, fill = factor)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Values",
       y = "Variables",
       fill = "Cluster",
       title = "Separation for wheat seeds") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1") +
  ggsave(filename = "boxplots_vars_clust_wheat_seeds.pdf",
         path = "C:/Users/lucru/Estat?stica_UFSCar/cv_cluster/figures", 
         width = 20.75, height = 12.5, units = "cm")


# making dendrograms for V1, V2, V3 and V6

# selecting the best 
ward.tree = convert_to_phylo(wheat_dend)

# dendrogram for V1
par(mfrow = c(2,4))
x = round(setNames(wheat_data_temp[,1],gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V1", ftype = "off")

# dendrogram for V2
x = round(setNames(wheat_data_temp[,2], gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V2", ftype = "off")

# dendrogram for V3
x = round(setNames(wheat_data_temp[,3], gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V3", ftype = "off")

# dendrogram for V4
x = round(setNames(wheat_data_temp[,4], gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V4", ftype = "off")

# dendrogram for V5
x = round(setNames(wheat_data_temp[,5], gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V5", ftype = "off")


# dendrogram for V6
x = round(setNames(wheat_data_temp[,6], gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V6", ftype = "off")

# dendrogram for V7
x = round(setNames(wheat_data_temp[,7], gsub(" ", "", rownames(wheat_data_temp))), 2)
reordered_x = x[ward.tree$tip.label]

obj = contMap(ward.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.6,0.8),outline=FALSE,lwd = c(2,5), leg.txt="V7", ftype = "off")


# selecting two combinations and analysing colored dendrogram for simulated dataset
# first combination with wide differences
# single linkage
sim_single_dend = hclust(dist(scale(sim.data[, -3]), method = "euclidean"), method = "single")

# converting to tree
single_sim.tree = convert_to_phylo(sim_single_dend)

# dendrogram for V1
par(mfrow = c(1,2))
x = setNames(sim.data[,1],gsub(" ", "", rownames(sim.data)))
reordered_x = x[single_sim.tree$tip.label]

obj = contMap(single_sim.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.5,0.9),outline=FALSE, lwd = c(2,5), leg.txt="X1", ftype = "off")

# dendrogram for V2
x = setNames(sim.data[,2], gsub(" ", "", rownames(sim.data)))
reordered_x = x[single_sim.tree$tip.label]

obj = contMap(single_sim.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.5,0.9),outline=FALSE, lwd = c(2,5), leg.txt="X2", ftype = "off")


# ward.D2 linkage
sim_ward_dend = hclust(dist(scale(sim.data[, -3]), method = "euclidean"), method = "ward.D2")

# converting to tree
ward_sim.tree = convert_to_phylo(sim_ward_dend)

# dendrogram for V1
par(mfrow = c(1,2))
x = setNames(sim.data[,1],gsub(" ", "", rownames(sim.data)))
reordered_x = x[ward_sim.tree$tip.label]

obj = contMap(ward_sim.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.5,0.9),outline=FALSE, lwd = c(2,5), leg.txt="X1", ftype = "off")

# dendrogram for V2
x = setNames(sim.data[,2], gsub(" ", "", rownames(sim.data)))
reordered_x = x[ward_sim.tree$tip.label]

obj = contMap(ward_sim.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.5,0.9),outline=FALSE, lwd = c(2,5), leg.txt="X2", ftype = "off")


# toy example for ceramic samples data set
ceramic_data = read.csv("/home/kuben/estatistica_UFSCAR/cv_cluster/data/Chemical Composion of Ceramic.csv")
head(ceramic_data)

# ignoring ceramic name and only using the binary categorical variable
ceramic_data %<>%
  dplyr::select(-(Ceramic.Name)) %>%
  mutate(Part = as.factor(Part))

# test list with only euclidean
test.list = list(hclust = c("single","ward.D2", "ward.D", "complete", "mcquitty", 
                            "average"), 
                 diana = NA)
dists = c("euclidean")
tol = 1e-20


# computing scores only for ward.D2 linkage
test.list = list(hclust = c("mcquitty"))
dists = c("euclidean")

# computing scores:
ceramic.l_cross_per_var = L_cross_val_per_var_alt(ceramic_data[, -1], test.list, 
                                                dists, scale = T)

ggplot_data = reshape2::melt(ceramic.l_cross_per_var)
ggplot_data$variable = as.factor(colnames(ceramic_data[, -1]))
p1 = ggplot_data %>%
  mutate(variable = fct_reorder(variable, value, .desc = T)) %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Features",
       y = "Phylogenetic Feature Importance Score (PFIS)") +
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(min(ggplot_data$value) - 0.05, max(ggplot_data$value) + 0.05)

ggsave(p1, filename = "importances_ceramic_samples_data.pdf",
         path = "/home/kuben/estatistica_UFSCAR/cv_cluster/figures", 
         width = 20.75, height = 10.5, units = "cm")

# plotting the best and worst variables according to importance
# ward.D2 linkage
ceramic_mcquitty_dend = hclust(dist(scale(ceramic_data[, -1]), 
                            method = "euclidean"), method = "mcquitty")

# converting to tree
mcquitty_ceramic.tree = convert_to_phylo(ceramic_mcquitty_dend)

# dendrogram for V1
par(mfrow = c(1,2))
x = setNames(ceramic_data$Al2O3,gsub(" ", "", rownames(ceramic_data)))
reordered_x = x[mcquitty_ceramic.tree$tip.label]

obj = contMap(mcquitty_ceramic.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.9),outline=FALSE, lwd = c(2,5), leg.txt = "Al203", ftype = "off")

# dendrogram for V2
x = setNames(ceramic_data$MgO,gsub(" ", "", rownames(ceramic_data)))
reordered_x = x[mcquitty_ceramic.tree$tip.label]

obj = contMap(mcquitty_ceramic.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,fsize=c(0.5,0.9),outline=FALSE, lwd = c(2,5), leg.txt = "MgO", ftype = "off")


set.seed(12500)
sim_tree = rtree(n = 4, tip.label = c("1", "0", "0", "?"))
plotTree(sim_tree, fsize= 1.5)
# adding the phylogenetic simulated example
# simulating tree
set.seed(1235)
sim_tree = rtree(n = 24, tip.label = letters[1:24])
plotTree(sim_tree, lwd = c(2, 5))
# generating discrete states:
disc_states = rTraitDisc(sim_tree, k = 2, freq = 0.45, rate = 0.65)

# plotting:
fitER<-ace(disc_states, sim_tree, model="ARD", type="discrete")
fitER
cols<-setNames(c("red","blue"),levels(disc_states))


plotTree(sim_tree, type="fan",fsize=0.9,ftype="i",lwd=1)
nodelabels(node=1:sim_tree$Nnode+Ntip(sim_tree),
           pie=fitER$lik.anc,
           piecol=cols,
           cex=0.4)
tiplabels(pie=to.matrix(disc_states[sim_tree$tip.label],
                        levels(disc_states)),piecol=cols,cex=0.25)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

# generating continuous
cont_states = rTraitCont(sim_tree, sigma = 2)


obj<-contMap(sim_tree, cont_states, plot=FALSE)
n = length(obj$cols)
obj$cols[1:n] = viridis::viridis(n)
plot(obj,lwd=7)

# plotting iris clustering:
# for categorical ancestral states
flowers = iris
d = daisy(flowers, metric = "gower")
complete.hc = hclust(d, method = "complete")
complete.tree = convert_to_phylo(complete.hc)
complete.tree$edge.length[which(complete.tree$edge.length %in% c(0))] = 10^(-3)
x = setNames(flowers$Species, row.names(flowers))
fitER<-ace(x, complete.tree, model="ER", type="discrete")
cols<-setNames(c("red","blue", "green"),levels(x))

plotTree(complete.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:complete.tree$Nnode+Ntip(complete.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(x[complete.tree$tip.label],
                        levels(x)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

# for continuous traits
par(mfrow= c(2, 2))
flowers_scale = scale(flowers[,-5])
x = setNames(round(flowers_scale[, 1], 2), gsub(" ", "", rownames(flowers)))
reordered_x = x[complete.tree$tip.label]

obj = contMap(complete.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="Sepal.length")


x = setNames(round(flowers_scale[, 2], 2), gsub(" ", "", rownames(flowers)))
reordered_x = x[complete.tree$tip.label]

obj = contMap(complete.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="Sepal.width")


x = setNames(round(flowers_scale[, 3], 2), gsub(" ", "", rownames(flowers)))
reordered_x = x[complete.tree$tip.label]

obj = contMap(complete.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="Petal.length")


x = setNames(round(flowers_scale[, 4], 2), gsub(" ", "", rownames(flowers)))
reordered_x = x[complete.tree$tip.label]

obj = contMap(complete.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="Petal.width")



