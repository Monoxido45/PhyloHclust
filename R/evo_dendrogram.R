#' @title Plot evolutionary dendogram
#' @description
#' Plot evolutionary dendrogram of a selected variable in a hierarchical clustering of a given dataset of interest.
#' @param data Data frame object of interest.
#' @param var_col Columns from the dataset to plot the evolutionary dendrogram.
#' @param hclust_obj Hierarchical clustering object obtained by applying a hierarchical clustering method in the
#' dataset of interest.
#' @param tip_names Set whether to plot row names in the evolutionary dendrogram. Default is TRUE.
#' @param scale_data Set whether to scale all continuous variable in the dataset when plotting.
#' @param fsize Figure size of each evolutionary dendrogram relative to the default dimensions from base plot. Default
#' is c(0.9, 0.8).
#' @param outline Set whether to draw border around each evolutionary dendrogram plot. Default is FALSE.
#' @param lwd Line width relative to the default. Default is 1.
#' @param leg_txt Legend text to place in legend corner. Default is NA.
#' @param palette_discrete Palette used to Ancestral State reconstruction of categorical variables. Default is "Set1".
#' @param cex  A numerical value giving the amount by which plotting text and symbols should be magnified relative to the
#' default. Default is c(0.4, 0.25).
#' @param ... Additional R base plot arguments to be passed to graphical parameters.
#' @return Evolutionary dendrogram plot.
#' @export
evo_dendrogram <- function(data,
                           var_col,
                           hclust_obj,
                           tip_names = TRUE,
                           scale_data = TRUE,
                           fsize = c(0.9, 0.8),
                           outline = FALSE,
                           lwd = 1,
                           leg_txt = NA,
                           palette_discrete = "Set1",
                           cex = c(0.4, 0.25),
                           ...){
  hclust_tree <- convert_to_phylo(hclust_obj)
  # plotting evo dendrogram for numeric variable
  if(class(data[,var_col]) == "numeric"){
    if(tip_names){
      if(scale_data == TRUE){
        scaled_data = scale(data[,var_col])
      }
      x = setNames(scaled_data, gsub(" ", "", rownames(data)))
      reordered_x = x[hclust_tree$tip.label]
    }
    obj <- phytools::contMap(hclust_tree, reordered_x, plot=FALSE)
    obj <- phytools::setMap(obj, invert=TRUE)

    n = length(obj$cols)
    obj$cols[1:n] = viridis::viridis(n)
    plot(obj,fsize = fsize,
         outline = outline,
         lwd = lwd,
         leg.txt= ifelse(is.na(leg_txt), paste0("Var_", var_col), leg_txt),
         ...)
  }else{ # plotting evo dendgram for categorical variable
    if(tip_names){
      x = setNames(as.factor(data[,var_col]),
                   gsub(" ", "", rownames(data)))
      reordered_x = x[hclust_tree$tip.label]
    }
    fitER <- ape::ace(x, hclust_tree, model="ER", type="discrete")
    lvls <- levels(x)
    cols <- setNames(palette.colors(n = length(lvls), palette_discrete), levels(x))

    phytools::plotTree(hclust_tree, type="fan", fsize = fsize[1],
                       ftype="i", lwd = lwd, ...)
    ape::nodelabels(node=1:hclust_tree$Nnode + Ntip(hclust_tree),
                    pie=fitER$lik.anc,
                    piecol=cols,
                    cex = cex[1])
    ape::tiplabels(pie = phytools::to.matrix(x[hclust_tree$tip.label],
                                             lvls),
                   piecol=cols,
                   cex = cex[2])
    phytools::add.simmap.legend(colors = cols,
                                prompt = FALSE,
                                x = 0.9*par()$usr[1],
                                y = 0.8*par()$usr[3],
                                fsize = fsize[2])
  }
}
