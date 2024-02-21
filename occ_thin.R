# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 20, 2024

library(tidyverse)
library(dismo)
library(terra) 
library(raster)
library(sp)


# Load cleaned occurrence data ---------------------------------------------
setwd("../occ_data/")
occ_cor <- readRDS(file = "occ_cor.Rdata") # n = 518
occ_fus <- readRDS(file = 'occ_fus.Rdata') # n = 985

# vectorize occurrence df to coordinates for bellow
occ_cor <- vect(occ_cor, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_fus <- vect(occ_fus, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

# Download WorldClim Bioclimatic raster -----------------------------------
setwd("../wclim_data/")
# Note DO NOT PUSH wclim data**
wclim <- worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "../wclim_data/")

# plot a raster to check it downloaded properly
plot(wclim$wc2.1_2.5m_bio_1)

occ_corALL <- occ_cor # keep for later
occ_fusALL <- occ_fus # keep for later


# Thin data using sampler -------------------------------------------------

set.seed(1337) # set random generator seed to get consistent results

# M. coronaria thinning
occ_cor <- spatSample(occ_cor, size = 1, 
                      strata = wclim) #sample one occurrence from each climatic cell

# M. fusca thinning
occ_fus <- spatSample(occ_fus, size = 1,
                      strata = wclim) #sample one occurrence from each climatic cell

# take a look at the thinned points

# plot(canUS_map, xlim = c(-85, -80), ylim = c(40.75, 41.25))
# plot(wclim$wc2.1_2.5m_bio_1, legend = F, add = T)
# plot(canUS_map, add = T)
# points(occ_corALL, pch = 16,col = "black")
# points(occ_cor, col = 'blue', pch = 1)

# dev.off()


# Extract environmental data from wclim ------------------------------------
# Extract and add wclim raster data to spatial occurrence data

occ_cor <- cbind(occ_cor, extract(wclim, occ_cor)) # M. coronaria
occ_fus <- cbind(occ_fus, extract(wclim, occ_fus)) # M. fusca

