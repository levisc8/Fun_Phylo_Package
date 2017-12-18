
context('FunPhylo tests')
#library(shiny)
library(FunPhylo)
#source('helpers.R')

data('tyson')

communities <- tyson$communities
phylo <- tyson$phylo
spp.list <- tyson$spp.list
demog <- tyson$demo.data
trait.data <- tyson$traits

traits <- names(trait.data[-1])
# test trait ktab functions
test_that('KTabs are constructed correclty', {
  trait.test <- traits[base::sample(1:24,9)]
  
  for(x in unique(demog$Species)){
    test.local <- make_local_trait_ktab(x, community.data = communities,
                                        trait.data = trait.data, traits = traits)
    expect_true(any(grepl('Q|B|C', test.local$VarTypes)))
    
  }
  
  test.regional <- make_regional_trait_ktab(trait.data, traits)
  
  expect_equal(class(test.regional), 'list')
  expect_true(all(grepl('Q|B|C', test.regional$VarTypes)))

})



