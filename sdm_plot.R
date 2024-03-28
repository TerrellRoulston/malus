# Top ---------------------------------------------------------------------
# SDM plotting and map making
# Started Mar 27, 2024

library(tidyverse)
library(terra)




# Load basemaps and sdm raster layers -------------------------------------
# Occurrence Points in SpatVectors
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = 'occThin_fus.Rdata') # M. fusca



# Download/load basemaps
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

us_map_0 <- gadm(country = 'USA', level = 0, resolution = 2,
                 path = "../occ_data/base_maps") #USA basemap without States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = '../occ_data/base_maps') # Canada basemap w. Provinces

ca_map_0 <- gadm(country = 'CA', level = 0, resolution = 2,
                 path = '../occ_data/base_maps') # Canada basemap without Provinces

mex_map <-gadm(country = 'MX', level = 1, resolution = 2,
               path = '../occ_data/base_maps') # Mexico basemap w. States

mex_map_0 <-gadm(country = 'MX', level = 0, resolution = 2,
               path = '../occ_data/base_maps') # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map
canUSMex_map_0 <- rbind(us_map_0, ca_map_0, mex_map_0) # Combine Country boundaries only


NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in Western Hemisphere

canUSMex_map <- crop(canUSMex_map, NA_ext) # crop to Western Hemisphere
canUSMex_map_0 <- crop(canUSMex_map_0, NA_ext) # crop to Western Hemisphere

# Great Lakes shapefiles for making pretty maps
great_lakes <- vect('C:/Users/terre/Documents/UBC/Botanical Garden/Malus Project/maps/great lakes/combined great lakes/')
great_lakes <- crop(great_lakes, NA_ext)


plot(canUSMex_map_0)
rect()
lines(canUSMex_map, col = 'grey') # plot basemap
lines(canUSMex_map_0, lwd = 1.25, col = 'black')
polys(great_lakes, col = 'lightblue')

dev.off()


