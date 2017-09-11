rm(list = ls(all = T))
library(FunPhylo)
data(tyson)
demog <- tyson$demo.data
communities <- tyson$communities
trait.data <- tyson$traits
traits <- names(trait.data)[-1]
phylo = tyson$phylo

out <- numeric()
for(x in unique(demog$Species)){
  y <- make_local_phylo_dist(x, communities, phylo)
  z <- make_local_trait_dist(x, communities, trait.data,
                             traits, 'scaledBYrange')
  rare <- rarefy_FPD(x, y, z, metric = 'NND',
                     n.rare = 11, a = .5, p = 2, abundance.weighted = T,
                     community.data = communities)
  
  out <- c(out, rare$rare.nnd)
}
