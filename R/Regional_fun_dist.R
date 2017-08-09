#' Make a regional scale functional distance matrix
#' 
#' @description This function is very similar to 
#' \code{make_local_traits_distance}, but will not subset to plot/focal species
#' level by default.
#' 
#' @inheritParams make_local_trait_dist
#' 
#' @return An object of class \code{dist}
#' 
#' @note These
#' can easily be coerced to other data structures, but must
#' be coerced to a \code{matrix} first.
#' 
#' @author Sam Levin
#' 
#' @seealso \code{\link{make_local_trait_dist}}, \code{\link{make_regional_trait_ktab}}
#' 
#' @importFrom ade4 dist.ktab
#' @export 


make_regional_trait_dist <- function(trait.data, traits){
  
  ktab <- make_regional_trait_ktab(trait.data, traits)
  
  out <- ade4::dist.ktab(ktab$KTab, type = ktab$VarTypes,
                         option = 'scaledBYrange')
  
  return(out)
}