#' Create a subsetted KTab
#' @description A wrapper around functions from \code{ade4} to create
#' the specialized data structure known as a \code{ktab}. These
#' are used in subsequent functions to calculate a generalized
#' Gower distance based on the methods of Pavoine et al 2009.
#'  
#' @inheritParams make_local_phylo_dist
#' @param trait.data A data frame that contains trait data information.
#' Currently, this only supports names from \emph{Levin et al 2017}, but
#' will hopefully soon be generalized to include all trait names listed
#' in \emph{TRY (Kattge et al 2017)}.
#' @param traits A character vector of trait names that should match the names
#' in \code{trait.data}. An effect method of subsetting this is to simply call
#' \code{names(trait.data)[first.trait:last.trait]}
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
#' @author Sam Levin
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom ade4 prep.circular prep.binary ktab.list.df
#' @export
#' 

make_local_traits_ktab <- function(focal.species,community.data,
                                   trait.data, traits){

  WoodyTraitNames <- c('SLA','Tough','WoodDens')
  HerbTraitNames <- c("Height","SLA","Tough")

  GrowthForm <- c("Stemmed_Herb",
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

  wood <- trait.data[trait.data$Species.Name == focal.species, 'Woody']
  whole.com <- dplyr::filter(community.data,
                             exotic_species == focal.species) %>%
    .$community

  local.com <- trait.data[trait.data$Woody == wood &
                            trait.data$Species.Name %in% whole.com, ]

  if(wood == 0){
    ContNames <- HerbTraitNames[HerbTraitNames %in% traits]
  } else {
    ContNames <- WoodyTraitNames[WoodyTraitNames %in% traits]
  }
  GFNames <- GrowthForm[GrowthForm %in% traits]
  DispNames <- DispersalTraitNames[DispersalTraitNames %in% traits]

  if(length(ContNames) > 0){
    ContTraits <- data.frame(local.com[ ,ContNames])
  }
  if(length(GFNames) > 0){
    GFTraits <- data.frame(local.com[ ,GFNames])
  }
  if(length(DispNames) > 0){
    DispTraits <- data.frame(local.com[ ,DispNames])
  }
  if("Flower.Period" %in% traits){
    Flower.Period <- data.frame(local.com[ ,'Flower.Period'])
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
                             rownames = local.com$Species.Name)

  out <- list(KTab = KTab, VarTypes = VarTypes)

  return(out)

}
