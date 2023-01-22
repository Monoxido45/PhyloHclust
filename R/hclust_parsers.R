# functions to convert clustering objects into dendrogram object
# and parsing dendrogram into phylo object
to.dend = function(cl.obj){
  # converting generical clustering object into hclust object first
  cluster.obj = cl.obj %>% stats::as.hclust()
  dend = cluster.obj %>% stats::as.dendrogram(hang = -1, check = TRUE)
  return(dend)
}

# parsing dendrogram into parentetic format
convert_to_par <- function(dend, first_it = TRUE)
{
  if (first_it == TRUE){
    dist = as.double(attr(dend, "height"))
    dend.object1 = dend[[1]]
    dist1 = as.double(attr(dend.object1, "height"))
    dend.object2 = dend[[2]]
    dist2 = as.double(attr(dend.object2, "height"))
    first_it = FALSE
    if (is.list(dend.object1) == TRUE & is.list(dend.object2) == TRUE){
      return(paste0("((",
                    convert_to_par(dend.object1, first_it),
                    "):",
                    dist - dist1,
                    ",",
                    "(",
                    convert_to_par(dend.object2, first_it),
                    "):",
                    dist - dist2,
                    ");"))}else{
                      if(is.list(dend.object1) == TRUE & is.list(dend.object2) == FALSE){
                        label = attr(dend.object2, "label")
                        return(paste0("((",
                                      convert_to_par(dend.object1, first_it),
                                      "):",
                                      dist - dist1,
                                      ",",
                                      label,
                                      ":",
                                      dist - dist2,
                                      ");"))
                      }else{
                        label = attr(dend.object1, "label")
                        return(paste0("(",
                                      label,
                                      ":",
                                      dist - dist1,
                                      ",",
                                      "(",
                                      convert_to_par(dend.object2, first_it),
                                      "):",
                                      dist - dist2,
                                      ");"))
                      }
                    }
  }else{
    dist = as.double(attr(dend, "height"))
    dend.object1 =  dend[[1]]
    dist1 = as.double(attr(dend.object1, "height"))
    dend.object2 = dend[[2]]
    dist2 = as.double(attr(dend.object2, "height"))
    if(is.list(dend.object1) == TRUE & is.list(dend.object2) == TRUE){
      return(paste0("(",
                    convert_to_par(dend.object1, first_it),
                    "):", dist - dist1,
                    ",",
                    "(",
                    convert_to_par(dend.object2, first_it),
                    "):", dist - dist2))
    }else{
      if(is.list(dend.object1) == FALSE & is.list(dend.object2) == TRUE){
        label1 = attr(dend.object1,"label")
        return(paste0(label1, ":",
                      dist - dist1,
                      ",",
                      "(",
                      convert_to_par(dend.object2, first_it),
                      "):", dist - dist2))
      }else{
        if(is.list(dend.object1) == TRUE & is.list(dend.object2) == FALSE){
          label2 = attr(dend.object2, "label")
          return(paste0("(",
                        convert_to_par(dend.object1, first_it),
                        "):", dist - dist1, ",",
                        label2, ":",
                        dist - dist2))
        }else{
          label1 = attr(dend.object1, "label")
          label2 = attr(dend.object2, "label")
          return(paste0(label1, ":", dist - dist1,
                        ",",
                        label2, ":", dist - dist2))
        }
      }
    }
  }
}

# finally, converting hclust object into phylo object
#' Converting clustering object into phylo object
#' @param cl.obj Clustering object
#' @return Returns a phylo object
convert_to_phylo = function(cl.obj){
  dend = to.dend(cl.obj)
  dend.str = convert_to_par(dend)
  tree = ape::read.tree(text = dend.str)
  return(tree)
}

