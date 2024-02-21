# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 21, 2024

library(tidyverse) # grammar, data management
library(ENMeval) # ecological niche models, aka SDMs, MaxEnt
library(raster) # working with rasters and vectors
library(sf) # converting raster objects
library(dismo)
# Load thinned occurrences points
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata')
occThin_fus <- readRDS(file = 'occThin_fus.Rdata')

# Load WorldClim data
setwd("../wclim_data/")
# Note DO NOT PUSH wclim data**
wclim <- worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "../wclim_data/")

envs <- raster::stack(wclim) #convert bioclimatic data to a raster stack

# Extract climatic variables at the occurrence points 
# These values are a 'reference' for the model

occThin_cor.z <- raster::extract(wclim, occThin_cor) # M coronaria
occThin_fus.z <- raster::extract(wclim, occThin_fus) # M fusca

occSim_cor <- similarity(wclim, occThin_cor.z) # compare enviromental conditions of reference points to background points 
plot(occSim_cor$similarity_min)
