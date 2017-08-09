# tests for server

rm(list = ls(all = T))
#library(shiny)
library(FunPhylo)
#source('helpers.R')

data('tyson')

communities <- tyson$communities
phylo <- tyson$phylo
tyson <- tyson$spp.list
demog <- tyson$demo.data
trait.tyson <- tyson$traits

traits <- names(trait.data[-1])
# test trait ktab functions
for(i in 1:10){
  trait.test <- traits[base::sample(1:24,9)]
  
  for(x in unique(demog$Species)){
    cat(x,"    ", i,'\n')
    test.local <- make_local_traits_ktab(x, community.data = communities,
                                         trait.data = trait.data, traits = traits)
  }
  
  test.regional <- make_regional_traits_ktab(trait.data, traits)
  
  
}


# test rarefying functions

for(x in unique(demog$Species)){
  phylo.mat <- make_local_phylo_dist(x, communities, phylo)
  fun.mat <- FunPhylo:::make_local_trait_dist(x, communities, trait.data,
                                              traits = traits,
                                              scale = 'scaledBYrange')
  
  FPD <- rarefy_FPD(x, phylo.mat = phylo.mat,
                    fun.mat = fun.mat,
                    n.rare = 11, a = .5, p = 2)
  
  demog[demog$Species == x, 'out'] <- as.numeric(FPD$rare.nnd)
  
}



