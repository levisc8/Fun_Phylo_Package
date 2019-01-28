# Script to build the summarized trait data, phylogeny, and community data
# data sets. These then get stored in tyson.rda in FunPhylo

library(dplyr)
library(ape)
library(pez)
library(stringr)
library(BIOMASS)
library(FunPhylo)
data(tyson)

Marko_Data <- read.csv('C:/Users/sl13sise/Desktop/Tyson/Writing/Tyson_Traits_Analysis/Data/TRCP_Traits_Marko.csv')
Marko_Data$Species <- gsub(" ", "_", Marko_Data$Species)


# The data are scattered between a couple traits files. Rather
# than do a merge by hand, I can check to make sure all records are matched using R 
# and then join them in here. 
RawTraits1 <- read.csv('../../Data/Full Traits w_Sites.csv',
                       stringsAsFactors = FALSE)
RawTraits2 <- read.csv('C:/Users/sl13sise/Desktop/Tyson/Writing/Traits with Sites incomplete.csv',
                       stringsAsFactors = FALSE)

# # check to make sure columns are identical
identical(RawTraits1$Species.Name,RawTraits2$Species.Name)
# 
# # crap!
setdiff(RawTraits1$Species.Name,RawTraits2$Species.Name)
# 
# # A ha! found the culprit!
RawTraits2$Species.Name<-gsub("Plantanus occidentalis","Platanus occidentalis",RawTraits2$Species.Name)
# 
# # check again
identical(RawTraits1$Species.Name,RawTraits2$Species.Name)
# 
# # time to join.
RawTraits <- RawTraits2 %>% 
  dplyr::select(Tough, Chlor, Legume:Flower.Period) %>% 
  cbind(RawTraits1, .) %>%
  select(-c(Growth.form, Leaf.Type, Source.LT.GF))



# we need to summarize the continuous data and
# create dummy variables for all of our categorical data. 
SummaryData <- SummaryData1 <- RawTraits %>% 
  group_by(Species.Name) %>%
  summarise(Height = mean(Height, na.rm = TRUE),
            SLA = mean(SLA.Ratio, na.rm = TRUE),
            Tough = mean(Tough, na.rm = TRUE),
            Chlor = mean(Chlor, na.rm = TRUE),
            #Wet.Dry = mean(Wet.Dry.Ratio, na.rm = TRUE),
            Woody = first(GF1),
            GF2 = first(GF2),
            LifeSpan = first(LifeSpan),
            LeafType = first(LeafType),
            N_Fixer = first(Legume),
            Flower.Period = mean(
              match(
                str_sub(Flower.Period,
                        1, 
                        3),
                month.abb
              )
            ),
            D1 = first(Dispersal1),
            D2 = first(Dispersal2),
            D3 = first(Dispersal3),
            Clonal = mean(Clonal, na.rm = TRUE)
  )



# Read in data that's already prepared. Alternatively, comment out this part
# and run commented sections before and after. Removing Euphorbia maculata
# and stenaria nigracans because we don't have data for any of their
# continuous traits, which causes problems in the distance calculations that we
# end up using the data for


names(SummaryData)[names(SummaryData) == "LeafType"] <- "Evergreen"
names(SummaryData)[names(SummaryData) == "GF1"] <- "Woody"

TankTree <- read.tree('C:/Users/sl13sise/Dropbox/MLU/Frequent_Data/PhytoPhylo.tre')


# Some species had height info collected incorrectly. For example, Carduus nutans
# height info was collected as height of the flowering stem, so it averaged ~120cm.
# It's a monocarpic perennial rosette. In the absence of height data for the 
# rosette leaves, we've decided to exclude it from calculations involving plant
# height

HtProblems<-read.csv('../../Data/Species_with_Ht_Problems.csv', 
                     stringsAsFactors = FALSE) %>%
  as.vector()

SummaryData[SummaryData$Species.Name %in% HtProblems$Species, 'Height'] <- NaN

SummaryData$Species.Name <- trimws(SummaryData$Species.Name)

SummaryData$Species.Name <- gsub(" ", "_", SummaryData$Species.Name)

# Now, create some additional columns for taxonomy and add wood density information
# in
SummaryData$Genus <- lapply(unique(SummaryData$Species.Name),
                            FUN = function(x) strsplit(x,"_")[[1]][1]) %>%
  as.character()

SummaryData$Species <- lapply(unique(SummaryData$Species.Name),
                              FUN = function(x) strsplit(x,"_")[[1]][2]) %>%
  as.character()

wood_dens <- getWoodDensity(genus=SummaryData$Genus, species=SummaryData$Species,
                            stand = NULL, family = NULL, region = "World",
                            addWoodDensityData = NULL) %>%
  filter(levelWD!='dataset')

wood_dens$FullName <- paste(wood_dens$genus, wood_dens$species, sep =  "_")

SummaryData$WoodDens <- NA_real_

for(x in wood_dens$FullName) {
  SummaryData[SummaryData$Species.Name == x,'WoodDens'] <- wood_dens[wood_dens$FullName == x, 'meanWD']
}

# create dummy variables for each categorical variable
SummaryData$Woody <- ifelse(SummaryData$Woody == "Woody",
                            1,
                            0)

for(i in unique(SummaryData$GF2)) {
  SummaryData[ ,i] <- ifelse(SummaryData$GF2 == i,
                             1,
                             0)
}

for(i in unique(SummaryData$LifeSpan)) {
  SummaryData[ ,i] <- ifelse(SummaryData$LifeSpan == i,
                             1,0)
}

SummaryData$EverGreen <- ifelse(SummaryData$Evergreen == "Evergreen", 1, 0)
SummaryData$N_Fixer <- ifelse(SummaryData$N_Fixer == "N_Fixer", 1, 0)

for(i in c("Subterranean", unique(SummaryData$D1))) {
  for(j in 1:length(SummaryData$D1)) {
    SummaryData[j, i] <- ifelse(i %in% c(SummaryData$D1[j],
                                         SummaryData$D2[j],
                                         SummaryData$D3[j]),
                                1,
                                0)
  }
}

# add in Marko Spasojevic's trait data from tyson. 
# If you use this data, you should cite Spasojevic et al. 2016 When does intraspecific
# trait variation contribute to beta diversity

for(x in unique(Marko_Data$Species)) {
  
  # If they have data on a specific trait and we do not, then we use their values.
  # Doing this for SLA, leaf toughness, and wood density
  if(x %in% SummaryData$Species.Name) {
    
    if(is.na(SummaryData$SLA[SummaryData$Species.Name == x])) {
      
      SummaryData$SLA[SummaryData$Species.Name == x] <- Marko_Data$SLA[Marko_Data$Species == x]
    }
    
    if(is.na(SummaryData$Tough[SummaryData$Species.Name == x])) {
      
      SummaryData$Tough[SummaryData$Species.Name == x] <- Marko_Data$LeafToughness[Marko_Data$Species == x] / 1000
    }
    
    if(is.na(SummaryData$WoodDens[SummaryData$Species.Name == x])) {
      
      SummaryData$WoodDens[SummaryData$Species.Name == x] <- Marko_Data$WoodDensity[Marko_Data$Species == x]
    }
  }
}

# Finally, our analyses excluded a few traits. Namely, evergreen, lifespan,
# and chlorophyll content. This will remove those from the final data set.
# Note that some of the variable names are for columns that were used to create
# dummy binary variables (e.g. GF2, GF1)

RmNames <- c('Chlor', 'GF2', 'GF1', 'LifeSpan', 'EverGreen', 'D1', 'D2', 'D3',
             'Genus', 'Species', 'Annual', 'Perennial', 'Monocarpic_Perennial',
             'Evergreen', 'Unassisted ')

FinalData <- SummaryData[ , !names(SummaryData) %in% RmNames] %>%
  select(c(Species.Name:Woody, Flower.Period:WoodDens, N_Fixer, Stemmed_Herb:Water))

# Removing Euphorbia maculata and stenaria nigracans because we don't have data
# for any of their continuous traits, which causes problems in the distance
# calculation (I think). Additionally, Allium vineale is a monocot, so I am removing
# that as well

RmSpp <- c('Allium_vineale', 'Euphorbia_maculata', 'Stenaria_nigricans')

FinalTraitData <- FinalData[!FinalData$Species.Name %in% RmSpp, ]
FinalTraitData$N_Fixer <- ifelse(FinalTraitData$N_Fixer == 'N_Fixer', 1, 0)

# Next, TRC species list. This will get used for making the phylogeny

spp_list <- read.csv("C:/Users/sl13sise/Dropbox/Thesis_SL/Data/Complete Tyson Flora List w_o Grasses_Sedges.csv",
                     stringsAsFactors = FALSE) %>%
  filter(Monocot != 1)

# trim whitespaces from ends and make sure anything w/ 2 spaces only has one.
# then replace spaces w/ underscores. Finally, replace dashes w/ periods
spp_list$Species <- trimws(spp_list$Species)
spp_list$Species <- gsub('  ', " ", spp_list$Species)
spp_list$Species <- gsub(' ', '_', spp_list$Species)
spp_list$Species <- gsub('-', '\\.', spp_list$Species)

setdiff(spp_list$Species, tyson$spp.list$Species)

# Tyson species list corrections
spp_list[spp_list$Species == 'Achillea_millefolium', c('Exotic', 'Invasive')] <- 1
spp_list$Species[spp_list$Species == 'Diodia_teres'] <- 'Diodella_teres'
spp_list$Species[spp_list$Species == 'Psoralidium_tenuiflorum'] <- 'Psoralea_tenuiflora'

# Next, add in species in community data that are missing from the complete species
# list.

communities <- read.csv("C:/Users/sl13sise/Dropbox/Thesis_SL/Data/Branch.length.data.csv",
                        stringsAsFactors = FALSE)

# convert -'s to periods, and make sure there is no leading or trailing
# whitespace in the species names
communities$exotic_species <- gsub(" ", "_", communities$exotic_species)
communities$community <- gsub('\\.', ' ', communities$community)
communities$community <- trimws(communities$community)
communities$community <- gsub('-', '\\.', communities$community)
communities$community <- gsub(' ', '_', communities$community)

# finally, correct the _sp to spp for consistency's sake
communities$community <- gsub('_sp$', '_spp', communities$community)

to_add <- communities[!communities$community %in% spp_list$Species, ] %>%
  .[!duplicated(.$community), ] %>%
  select(community, alien, Invasive) %>%
  setNames(
    c(
      'Species',
      'Exotic',
      'Invasive'
    )
  ) %>%
  mutate(Monocot = 0) %>%
  filter(Species != 'Allium_vineale')

# Finally, add Cerastium spp. so it'll match the trait data.
cer_spp <- data.frame(Species = 'Cerastium_spp.', Exotic = NA, Invasive = NA, Monocot = 0)

spp_list <- rbind(spp_list, cer_spp, to_add)

# Done w/ species list! Now time for the phylogeny

FullTree <- congeneric.merge(TankTree, spp_list$Species)

TyPhy <- drop.tip(FullTree, setdiff(FullTree$tip.label, spp_list$Species))

# rm(FullTree, TankTree)

# Done with the phylogeny!

# Now, fix the names in each community data set

# Check to see how many are missing from TyPhy
unique(communities$community)[! unique(communities$community) %in% TyPhy$tip.label]

# Allium is the only one missing :) This should be true, as we removed it 
# due to it being a monocot

tyson$communities <- communities
tyson$traits <- FinalTraitData
tyson$phylo <- TyPhy
tyson$spp.list <- spp_list

usethis::use_data(tyson,
                  overwrite = TRUE)
