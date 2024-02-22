# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 21, 2024

library(tidyverse) # grammar, data management
library(ENMeval) # ecological niche models, aka SDMs, MaxEnt
library(raster) # working with rasters and vectors
library(sf) # converting raster objects
library(dismo)
library(geodata)
library(predicts)


# Load thinned occurrence records -----------------------------------------
getwd() #check wd and return or forward (../ or ./) to occ_data
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata')
occThin_fus <- readRDS(file = 'occThin_fus.Rdata')
# Load WorldClim data
setwd("../wclim_data/")
# Note DO NOT PUSH wclim data**
wclim <- geodata::worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "../wclim_data/")


# Download/load maps ------------------------------------------------------
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = '../occ_data/base_maps') #Canada basemap w. Provinces

canUS_map <- rbind(us_map, ca_map) # combine US and Canada vector map


# Prepare maps anddata for MESS -------------------------------------
e <- ext(-180, 0, 10, 90) # limit extent to western hemisphere
canUS_map <- crop(canUS_map, e) # crop raster to extent

wclimCanUS <- crop(wclim, canUS_map) #crop worldclim data to Canada/US extent

envs <- raster::stack(wclimCanUS) #convert worldclim data to a raster stack

plot(envs[[1]]) #plot bioclimatic variable 1 to make sure raster stack worked correctly

occThin_cor.x <- st_as_sf(occThin_cor) # convert spat vector to sf object
occThin_fus.x <- st_as_sf(occThin_fus)

occThin_cor.sp <- sp::SpatialPoints(occThin_cor.z)

# MESS --------------------------------------------------------------------
# Extract climatic data from cells where occurrence points are present
# These values are a 'reference' for the model
occThin_cor.z <- raster::extract(envs, occThin_cor.x) # M coronaria
occThin_fus.z <- raster::extract(envs, occThin_fus.x) # M fusca


# Compare environmental conditions of reference points to background points 
# similarity() calculates the environmental similarity of background environments
# to the the environments where species occurrences exist
# positive values at more similar and negative more dissimilar 

# occSim_cor <- ENMeval::similarity(envs, occThin_cor.z) 
# occSim_fus <- ENMeval::similarity(envs, occThin_fus.z)

# memory intensive so save and load data
setwd("../mess_data/")
# saveRDS(occSim_cor, file = 'occSim_cor.Rdata')
# saveRDS(occSim_fus, file = 'occSim_fus.Rdata')

occSim_cor <- readRDS(file = 'occSim_cor.Rdata')
occSim_fus <- readRDS(file = 'occSim_fus.Rdata')

# Note similarity objects have three variables
# $similarity_min is the minimum similarity values, i.e. MESS
# $mod = which variables are most dissimilar to reference
# $mos = which variables are most similar 


# MESS
occMESS_fus <- occSim_fus$similarity_min
occMESS_cor <- occSim_cor$similarity_min


plot(occMESS_cor) 
points(occThin_cor)

dev.off() # close graphics device

# prepare colours for the 19 climatic variables
# mode
category_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                     "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3")

# which variable 
terra::plot(occSim_cor$mod, col = category_colors,
            xlim = c(-100, -70), ylim = c(30, 50), cex = 2)
terra::plot(canUS_map, xlim = c(-100, -70), ylim = c(30, 50), add = T)
terra::points(occThin_cor, pch = 3, add = T, col = 'red')


plot(occThin_cor)
dev.off()
