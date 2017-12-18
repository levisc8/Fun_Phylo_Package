context('Rarefying and abundance weighting')

rm(list = ls(all = T))
library(FunPhylo)
library(data.table)
data(tyson)
demog <- tyson$demo.data
communities <- data.table(tyson$communities, key = 'exotic_species')
trait.data <- tyson$traits
traits <- names(trait.data)[-1]
phylo = tyson$phylo


test_that('Rarefying and abundance weighting work', {
  y <- make_local_phylo_dist('Ailanthus_altissima', communities, phylo)
  z <- make_local_trait_dist('Ailanthus_altissima', communities, trait.data,
                             traits, 'scaledBYrange')
  rare <- rarefy_FPD('Ailanthus_altissima', y, z,
                     n.rare = 11, a = .5, p = 2, abundance.weighted = T,
                     community.data = communities)
  
  expect_equal(class(rare), 'list')
  expect_equal(length(rare$sample.mpds), 1000)
  expect_true(is.numeric(rare$rare.mpd))
  expect_true(all(is.numeric(rare$rare.mpd)))
  expect_true(all(is.numeric(rare$rare.nnd)))
  
  UW <- rarefy_FPD('Ailanthus_altissima', y, z,
                   n.rare = 11, a = .5, p = 2)
  expect_equal(class(UW), 'list')
  expect_equal(length(UW$sample.mpds), 1000)
  expect_true(is.numeric(UW$rare.mpd))
  expect_true(all(is.numeric(UW$rare.mpd)))
  expect_true(all(is.numeric(UW$rare.nnd)))
  
})



