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
getwd()
setwd("../occ_data/") 


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
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, botanical gardens, and other institutions 
  cc_sea() %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>%
  filter(decimalLongitude >= -100) %>% # remove some records on west coast, two from bot gardens
  filter(!(decimalLatitude < 35 & decimalLongitude < -86)) %>% # remove records from Texas, Oklahoma, Louisiana
  filter(!(decimalLatitude > 45)) %>% # remove record from New Brunswick
  filter(!(decimalLongitude < -98)) %>% # remove iNat record from Kansas, northern Kansas record verified by taxonomist
  dplyr::select(species, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
         )


# Note it is helpful to plot the occurrences bellow, and then add more conditions to clean inaccurate points
# Pay special attention to points at the edge of the range of occurrences, as these are most likely to be suspicious
# as well as influence the model in strange ways, in comparison to inaccurate points well within the other occurrences

# Save M. coronaria occurrence dataframe ----------------------------------
saveRDS(occ_cor, file = "occ_cor.Rdata") 

# compare pre/post 1970 occurrence data
# if no points between pre/post are drastically different then keep all data
# consider removing if different from post-1970

occ_cor_pre <- occ_cor %>% # pre 1970
  filter(year < 1970)

occ_cor_post <- occ_cor %>% # post 1970
  filter(year >= 1970)

occ_cor_nd <- occ_cor %>% 
  filter(year %in% NA)

# filter M. fusca data in df
occ_fus <- gbif_fusca %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot cardens 
  cc_sea() %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>% 
  filter(!(decimalLongitude > -113)) %>%  # remove record from Idaho herbarium 
  filter(!(decimalLatitude > 64)) %>% # remove record from Interior Alaska
  filter(!(stateProvince == 	'Idaho')) %>% # remove phishy record from Idaho
  dplyr::select(species, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )


# Save M. fusca occurrence dataframe --------------------------------------
saveRDS(occ_fus, file = "occ_fus.Rdata")

# compare pre/post 1970 occurrences
occ_fus_pre <- occ_fus %>% 
  filter(year < 1970)

occ_fus_post <- occ_fus %>% 
  filter(year >= 1970)

occ_fus_nd <- occ_fus %>% 
  filter(year %in% NA)

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
points(occ_cor$decimalLongitude, occ_cor$decimalLatitude, pch = 16,
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
plot(canUS_map, xlim = c(-100, -60), ylim = c(25, 50))
# plot Malus coronia occurrences
# pre-1970
points(occ_cor_pre$decimalLongitude, occ_cor_pre$decimalLatitude, pch = 16,
       col = alpha("red", 0.2))
# post-1970
points(occ_cor_post$decimalLongitude, occ_cor_post$decimalLatitude, pch = 16,
       col = alpha("blue", 0.2))
# No date
points(occ_cor_nd$decimalLongitude, occ_cor_nd$decimalLatitude, pch = 16,
       col = alpha('black', 0.2))

legend(x= -75,
       y = 33,
       title = 'M. coronaria',
       legend = c('Pre-1970 (n=422)', 'Post-1970 (n=558)', 'No Date (n=52)'),
       fill = c('red', 'blue', 'black'))

dev.off() # close graphics plot window

# M fusca
plot(canUS_map, xlim = c(-170, -110), ylim = c(30, 65))
# plot Malus fusca occurrences
# post-1970
points(occ_fus_post$decimalLongitude, occ_fus_post$decimalLatitude, pch = 16,
       col = alpha("blue", 0.2))
# pre-1970
points(occ_fus_pre$decimalLongitude, occ_fus_pre$decimalLatitude, pch = 16,
       col = alpha("red", 0.2))
# No date
points(occ_fus_nd$decimalLongitude, occ_fus_nd$decimalLatitude, pch = 16,
       col = alpha('black', 0.2))


legend(x= -165,
       y = 40,
       title = 'M. fusca',
       legend = c('Pre-1970 (n=191)', 'Post-1970 (n=1230)', 'No Date (n=45)'),
       fill = c('red', 'blue', 'black'))

dev.off()

# Zoom into PNW
plot(canUS_map, xlim = c(-140, -110), ylim = c(43, 55))
# plot Malus fusca occurrences
points(occ_fus$decimalLongitude, occ_fus$decimalLatitude, pch = 16,
       col = alpha("blue", 0.2))


