#' Make a functional trait distance matrix
#' 
#' @description Creates a Gower distance matrix using the method
#' of Pavoine et al 2009 to handle many types of data
#' 
#' @inheritParams make_local_trait_ktab
#' @param scale The spatial scale at which to analyze the data
#' 
#' @return An object of class \code{dist}.
#' 
#' @note These
#' can easily be coerced to other data structures, but must
#' be coerced to a \code{matrix} first.
#' 
#' @author Sam Levin
#' 
#' @seealso \code{\link{make_regional_trait_dist}}, \code{\link{make_local_trait_ktab}}
#' 
#' @importFrom ade4 dist.ktab
#' @export 
#' 

make_local_trait_dist <- function(focal.species, community.data, trait.data, traits, scale){
  

  ktab <- make_local_trait_ktab(focal.species, community.data, trait.data, 
                                 traits)
  
  out <- ade4::dist.ktab(ktab$KTab, type = ktab$VarTypes,
                         option = scale)
  
  return(out)
  
  
}
