# Top ---------------------------------------------------------------------
# cleaning occurence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 16, 2024

library(tidyverse)
library(CoordinateCleaner)

# load CSV occurrence data files
gbif_cor <- read.csv(file = "occ_coronaria.csv") # load coronaria data
gbif_fusca <- read.csv(file = "occ_fusca.csv") # load fusca data


occ_cor <- gbif_cor %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot cardens 
  cc_sea() %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T)
  
