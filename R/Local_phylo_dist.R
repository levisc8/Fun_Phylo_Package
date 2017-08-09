#' Make Local Scale Phylogenetic Distance Matrix
#' 
#' @description Function for creating a phylogenetic distance matrix
#' from data that does not require you to subset first. 
#' 
#' @param focal.species The focal species or plot identifier that
#' you want to use to create a the phylogeny
#' @param community.data The community data file. Each observation
#' of a species in a community should be its own row. For now,
#' the column specifiy which community is which is called 
#' \code{exotic_species}, but this will be generalized later
#' using an argument to specifiy which column marks community 
#' separators and which column holds the actual species names
#' @param phylo A larger phylogeny that contains all of the 
#' species in the \code{community.data} file. 
#' 
#' @return A phylogenetic distance matrix in the form of a
#' \code{data.frame}
#' 
#' @author Sam Levin
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom ape drop.tip
#' @importFrom ape cophenetic.phylo
#' @export

make_local_phylo_dist <- function(focal.species, community.data,
                                  phylo){
  local.com <- dplyr::filter(community.data,
                             exotic_species == focal.species) %>%
    .$community %>% as.character()
  
  out <- ape::drop.tip(phylo, setdiff(phylo$tip.label, local.com)) %>%
    ape::cophenetic.phylo() %>% data.frame()
  
  return(out)
  
}
