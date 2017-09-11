#' Compute an abundance weighted FPD
#' 
#' @description This function is designed to be called internally
#' by \code{rarefy_FPD} and weights each species distance based on
#' their local abundance. 
#' 
#' @inheritParams make_local_phylo_dist
#' @param fpd.mat A functional-phylogenetic distance matrix
#' @param metric The metric of novelty you'd like to calculate.
#' Currently supports mean pairwise distance and nearest neighbor
#' distance
#' @param na.rm A logical indicating whether to remove missing
#' values
#' 
#' @return A numeric
#' 
#' @author Sam Levin
#' 
#' @importFrom dplyr filter
#' @importFrom stats weighted.mean
#' 
#' @export

AW_calc <- function(focal.species, community.data, fpd.mat, 
                    metric = c('NND', 'MPD'), na.rm = TRUE){
  
  focal.pos <- which(rownames(fpd.mat) == focal.species)
  fpd.mat[focal.pos] <- 0
  
  if(!identical(rownames(fpd.mat), community.data$community)){
    stop('The developer has made a mistake in AW_calc or
         rarefy_FPD (lines 98-113)')
  }
  if(na.rm){
    
    idx <- which(!is.na(community.data[ ,2]))
    community.data <- community.data[idx, ]
    fpd.mat <- data.frame(fpd.mat[idx, ])
    
    rownames(fpd.mat) <- community.data$community
  }
  weights <- as.matrix(community.data$percentcover) %*%
             t(as.matrix(community.data$percentcover))
  colnames(weights) <- community.data$community
  
  fpd.mat[which(rownames(fpd.mat) == focal.species),
          1] <- NA
  
  if(metric == 'MPD'){
    out <- weighted.mean(as.vector(fpd.mat[ ,1]),
                         as.vector(weights[ ,focal.species]),
                         na.rm = na.rm)
  }
  
  if(metric == 'NND'){
    out <- weighted.min(as.vector(fpd.mat[ ,1]),
                         as.vector(weights[ ,focal.species]),
                        na.rm = na.rm)
  }
  
  return(out)
}
