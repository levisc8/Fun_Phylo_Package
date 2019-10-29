#' Make a species pool-level phylogenetic distance matrix
#' 
#' @description Create a phylogeny from a species pool and then calculate
#' patristic distances between all species
#' 
#' @param spp.vector A character vector of species names. Species names
#' must match those in phylogeny.
#' @param phylo A larger phylogeny that contains all of the 
#' species in the \code{spp.vector}. If this condition is not met,
#' the function will try to merge in the missing species using
#' \code{pez::congeneric.merge} while warning the user. 
#' @param square_root Square root transform branch lengths in distance matrix?
#' 
#' @return A phylogenetic distance matrix in the form of a
#' \code{data.frame}
#' 
#' @author Sam Levin
#' 
#' @importFrom magrittr %>%
#' @importFrom ape cophenetic.phylo drop.tip
#' @importFrom pez congeneric.merge
#' @export

make_regional_phylo_dist <- function(spp.vector, phylo, square_root = TRUE) {
  
  spp.vector <- gsub("-","\\.",spp.vector)
  
  if(!all(spp.vector %in% phylo$tip.label)){
    warning('Attempting merge missing species into tree.\n',
            'Make sure to check output for accuracy.')
    phylo <- pez::congeneric.merge(phylo, spp.vector)
  }
  
  phylo <- ape::drop.tip(phylo, setdiff(phylo$tip.label, spp.vector))
  bigDist <- ape::cophenetic.phylo(phylo) 
  
  if(square_root) {
    bigDist <- sqrt(bigDist)
  }

  if(dim(bigDist)[1] != dim(bigDist)[2]){
    stop('resulting matrix is not square')
  }
  
  return(bigDist)
}
