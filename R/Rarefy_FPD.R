#' @title Rarefied Function-Phylogenetic Distances
#' @description 
#' \code{rare_FPD} resamples a phylogenetic and functional distance matrix
#' a set number of times so that communities of different sizes can be 
#' compared to each other.
#' 
#' @param focal.species The name of the focal species in the 
#' community to be rarefied.
#' @param phylo.mat,fun.mat A phylogenetic or functional distance
#' matrix with species names as row names and column names. 
#' Can either be a \code{data frame}, \code{matrix}, or \code{dist}.
#' @param metric Either MPD (\code{MPD}) or nearest neighbor distance
#' (\code{NND}). The default is to calculate both metrics.
#' @param n.resamp The number of rarefied resamplings to perform. 
#' Default is 1000. 
#' @param n.rare The number of species to rarefy the community down to.
#' Must be less than the total in the community.
#' @param a The phylogenetic scaling factor for the calculation of
#' functional-phylogenetic distances
#' @param p The power to raise each distance value to in the functional-
#' phylogenetic distance calculation
#' 
#' @author Sam Levin
#' 
#' @importFrom magrittr %>%
#' @export
#' 
rarefy_FPD <- function(focal.species, phylo.mat, fun.mat,
                       metric = c("MPD", "NND"),
                       n.resamp = 1000, n.rare, 
                       a, p) {
  
  mpd.tf <- "MPD" %in% metric
  nnd.tf <- "NND" %in% metric
  
  if(!inherits(phylo.mat, c('dist', 'matrix', 'data.frame')) |
     !inherits(fun.mat,  c('dist', 'matrix', 'data.frame'))) {
    stop('coerce phylogenetic and/or functional distance matrices to\n',
         'one of the following classes: "dist", "matrix", or "data.frame"')
  }
  
  if(inherits(phylo.mat, 'dist') | inherits(phylo.mat,'matrix')){
    phylo.mat <- as.data.frame(as.matrix(phylo.mat))
  }

  if(inherits(fun.mat, 'dist') | inherits(fun.mat,'matrix')){
    fun.mat <- as.data.frame(as.matrix(fun.mat))
  }

  if(dim(phylo.mat)[1] < dim(fun.mat)[1]){
    stop('More species with trait data than are in phylogeny')
  }
  
  phylo.mat <- phylo.mat[rownames(phylo.mat) %in% sort(rownames(fun.mat)),
                         names(phylo.mat) %in% sort(names(fun.mat))] %>%
    .[sort(rownames(.)), sort(names(.))]
  fun.mat <- fun.mat[sort(rownames(fun.mat)), sort(names(fun.mat))]
  
  if(!identical(names(phylo.mat), names(fun.mat)) & 
     !identical(rownames(phylo.mat), rownames(fun.mat))){ 
    stop('The developer is dumb and has made a mistake.\n',
         'Email a reproducible example to levisc8@gmail.com')
  }
  fpd <- func_phy_dist(FDist = as.matrix(fun.mat),
                       PDist = as.matrix(phylo.mat),
                       phyloWeight = a, p = p)
  
  if(mpd.tf){
    rare.mpd <- rep(NA, n.resamp)
  }
  if(nnd.tf){
    rare.bl <- rep(NA, n.resamp)
  }

  diag(fpd) <- NA
  focal.column <- fpd[ ,focal.species]

  for(i in 1:n.resamp){
    resamp.x <- base::sample(1:length(focal.column),
                             size = n.rare,
                             replace = FALSE)

    if(mpd.tf) {
      rare.mpd[i] <- mean(focal.column[resamp.x], na.rm = TRUE)
    }
    if(nnd.tf){
      rare.bl[i] <- min(focal.column[resamp.x], na.rm = TRUE)
    }
  }

  out <- list()
  if(mpd.tf){
    out$rare.mpd <- mean(rare.mpd)
    out$sample.mpds <- rare.mpd
  }
  if(nnd.tf){
    out$rare.nnd <- mean(rare.bl)
    out$sample.nnds <- rare.bl
  }

  return(out)

}

