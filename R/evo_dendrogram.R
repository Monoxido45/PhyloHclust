#' @title Plot evolution dendogram for given data and selected variables
#' @export
evo_dendrogram <- function(data,
                           var_col,
                           hclust_obj,
                           tip_names = TRUE,
                           scale_data = TRUE,
                           fsize = c(0.9, 0.8),
                           outline = FALSE,
                           lwd = c(3,7),
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
