# Top ---------------------------------------------------------------------
# cleaning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 16, 2024

library(tidyverse) #grammar, data management
library(CoordinateCleaner) #helpful functions to clean data
library(terra) #working with vector/raster data
library(geodata) #download basemaps
library(scales) #alpha adjust colours

# set wd
setwd("./occ_data/") 


# load occurrence csv files ------------------------------------------------
gbif_cor <- read.csv(file = "occ_coronaria.csv") # load coronaria data
gbif_fusca <- read.csv(file = "occ_fusca.csv") # load fusca data


# filter occurrence data in df ---------------------------------------------

# filter M. corornia data
occ_cor <- gbif_cor %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot gardens 
  cc_sea() %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>%
  filter(decimalLongitude !=-123.10000) %>% #remove one record on west coast
  dplyr::select(species, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
         )

# filter M. fusca data in df
occ_fus <- gbif_fusca %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot cardens 
  cc_sea() %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>% 
  select(species, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )

# map occurrence points (take a peak) -------------------------------------

# download/load maps
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = '../occ_data/base_maps') #Canada basemap w. Provinces

canUS_map <- rbind(us_map, ca_map) #combine US and Canada vector map

# plot basemap
plot(canUS_map, xlim = c(-180, -50))

# plot Malus coronia occurrences
points(occ_cor, pch = 16,
       col = alpha("red", 0.2))

# plot Malus fusca occurrences
points(occ_fus$decimalLongitude, occ_fus$decimalLatitude, pch = 16,
       col = alpha('blue', 0.2))

# add legend
legend(x= -170,
       y = 40,
       legend = c('M. coronaria', 'M. fusca'),
       fill = c('red', 'blue'))

# clear graphics plot
dev.off()


# Zoom into each species --------------------------------------------------

# M coronaria
plot(canUS_map, xlim = c(-100, -70), ylim = c(30, 50))
# plot Malus coronia occurrences
points(occ_cor$decimalLongitude, occ_cor$decimalLatitude, pch = 16,
       col = alpha("red", 0.2))

dev.off()

# M fusca
plot(canUS_map, xlim = c(-170, -110), ylim = c(30, 60))
# plot Malus fusca occurrences
points(occ_fus$decimalLongitude, occ_fus$decimalLatitude, pch = 16,
       col = alpha("blue", 0.2))

# Zoom into PNW
plot(canUS_map, xlim = c(-140, -110), ylim = c(43, 55))
# plot Malus fusca occurrences
points(occ_fus$decimalLongitude, occ_fus$decimalLatitude, pch = 16,
       col = alpha("blue", 0.2))


# save cleaned plot data --------------------------------------------------
# if happy with filter save df as .rdata file
# load in following analysis
saveRDS(occ_cor, file = "occ_cor.Rdata")
saveRDS(occ_fus, file = "occ_fus.Rdata")
