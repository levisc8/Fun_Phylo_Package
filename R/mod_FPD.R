#' @title Functional-Phylogenetic Distances (FPD)
#' @description Calculate hybrid functional phylogenetic distance matrices
#' based on the method of Cadotte et al 2013
#' 
#' @param FDist A matrix of functional trait distances
#' @param PDist A matrix of phylogenetic distances
#' @param phyloWeight The phylogenetic scaling parameter \emph{a} as
#' described in Cadotte et al 2013. A value of 1 denotes only phylogenetic
#' information while a value of 0 indicates only functional trait 
#' information. The value must be between 0 and 1. 
#' @param p the power to raise each matrix to when combining the two 
#' distance matrices. 
#' @param ... Ignored here
#' 
#' @return A matrix with weighted interspecific distances
#' 
#' @note The function does not check to make sure names in 
#' \code{FDist} and \code{PDist match}.
#' 
#' @seealso \code{\link{rarefy_FPD}}
#' @author Sam Levin
#' @export


func_phy_dist <- function (FDist, PDist, phyloWeight, p, ...)
{
  if (phyloWeight < 0 | phyloWeight > 1)
    stop("'phyloWeight' must be between 0 and 1")
  if (!is.numeric(p))
    stop("'p' must be a numeric")

  FDist <- FDist/max(FDist)
  PDist <- PDist/max(PDist)
  (phyloWeight * PDist^p + (1 - phyloWeight) * FDist^p)^(1/p)
}


