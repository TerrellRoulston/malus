# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 20, 2024

library(tidyverse) # grammar and data management 
library(terra) # working with spatial data
library(geodata) # basemaps and climate data


# Load cleaned occurrence data ---------------------------------------------

occ_cor <- readRDS(file = "./occ_data/occ_cor.Rdata")
occ_fus <- readRDS(file = './occ_data/occ_fus.Rdata') 

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
plot(wclim$wc2.1_2.5m_bio_1, main = 'Annual Mean Temperature')


# Thin data using sampler -------------------------------------------------

set.seed(1337) # set random generator seed to get reproducible results

# M. coronaria thinning
occ_cor <- spatSample(occ_cor, size = 1, 
                      strata = wclim) #sample one occurrence from each climatic cell

# M. fusca thinning
occ_fus <- spatSample(occ_fus, size = 1,
                      strata = wclim) #sample one occurrence from each climatic cell


# Save thinned occurrence points for further analysis ----------------------
setwd("../occ_data/")
saveRDS(occ_cor, file = 'occThin_cor.Rdata')
saveRDS(occ_fus, file = 'occThin_fus.Rdata')

# Extract environmental data from wclim ------------------------------------
# Extract and add wclim raster data to spatial occurrence data

#occ_cor <- cbind(occ_cor, extract(wclim, occ_cor)) # M. coronaria
#occ_fus <- cbind(occ_fus, extract(wclim, occ_fus)) # M. fusca

