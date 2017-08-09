#' Make a KTab for a whole data set
#' 
#' @description A wrapper around functions from \code{ade4} to create
#' the specialized data structure known as a \code{ktab}. These
#' are used in subsequent functions to calculate a generalized
#' Gower distance based on the methods of Pavoine et al 2009.
#' 
#' @inheritParams make_local_trait_ktab
#' 
#' @details \code{ade4} provides many useful functions for calculating
#' trait based distances. The aim of this particular function is to create
#' an easy-to-automate version of these for use in stepwise model selection
#' procedures which can handle odd boundary cases not necessarily forseen
#' in previous implementations. This particular function makes it easy
#' to loop across a variety of traits without having to stop and figure 
#' out which type of variable each particular portion of the \code{ktab}
#' is. The \code{VarType} output is designed to seamlessly integrate
#' with \code{dist.ktab}'s \code{type} argument to allow for effortless
#' integration.
#' 
#' @return A list consisting of 
#' \itemize{
#'    \item{\code{KTab}}{The \code{ktab} which can be passed to 
#'    \code{dist.ktab}}
#'    \item{\code{VarTypes}}{A character vector that can be passed to 
#'    \code{dist.ktab}}
#' }
#' 
#' @note The output from this function is designed to be passed directly to
#' \code{make_regional_trait_dist}. 
#' 
#' @seealso \code{\link{make_regional_trait_dist}}, \code{\link{make_local_trait_ktab}}
#' @author Sam Levin
#' @importFrom ade4 prep.binary prep.circular ktab.list.df
#' @export
#' 
make_regional_trait_ktab <- function(trait.data, traits){

  ContTraitNames <- c('SLA','Tough','WoodDens','Height')

  GrowthForm <- c("Woody", "Stemmed_Herb",
                  "Tree", "Rosette",
                  "Vine", "SubShrub",
                  "Shrub",
                  "Elongated_Leafy_Rhizomatous",
                  "N_Fixer")

  DispersalTraitNames <- c("Subterranean",
                           "Unassisted",
                           "Wind",
                           "ExoZoochory",
                           "EndoZoochory",
                           "Ballistic",
                           "Hoarding",
                           "Myrmecochory",
                           "Water",
                           "Clonal")

  spp.names <- trait.data$Species.Name

  ContNames <- ContTraitNames[ContTraitNames %in% traits]
  GFNames <- GrowthForm[GrowthForm %in% traits]
  DispNames <- DispersalTraitNames[DispersalTraitNames %in% traits]

  if(length(ContNames) > 0){
    ContTraits <- data.frame(trait.data[ ,ContNames])
  }
  if(length(GFNames) > 0){
    GFTraits <- data.frame(trait.data[ ,GFNames])
  }
  if(length(DispNames) > 0){
    DispTraits <- data.frame(trait.data[ ,DispNames])
  }
  if("Flower.Period" %in% traits){
    Flower.Period <- data.frame(trait.data[ ,'Flower.Period'])
  }

  VarTypes <- NULL

  if(exists("ContTraits")){
    VarTypes <- c(VarTypes, "Q")
  } else {
    ContTraits <- NULL
  }

  if(exists("GFTraits")){
    GFKTab <- ade4::prep.binary(GFTraits,
                                col.blocks = dim(GFTraits)[2],
                                label = "Growth_Form")
    VarTypes <- c(VarTypes, "B")
  } else {
    GFKTab <- NULL
  }

  if(exists("DispTraits")){
    DispKTab <- ade4::prep.binary(DispTraits,
                                  col.blocks = dim(DispTraits)[2],
                                  label="Dispersal_Traits")
    VarTypes <- c(VarTypes, "B")

  } else {
    DispKTab <- NULL
  }

  if(exists("Flower.Period")){
    FlowerKTab <- ade4::prep.circular(Flower.Period,
                                      rangemin=1, rangemax=12)
    VarTypes <- c(VarTypes, "C")

  } else {
    FlowerKTab <- NULL
  }


  forKTab <- check_empty_elements(ContTraits, GFKTab,
                                  DispKTab, FlowerKTab)

  KTab <- ade4::ktab.list.df(forKTab,
                             rownames = spp.names)

  out <- list(KTab = KTab, VarTypes = VarTypes)

  return(out)

}
