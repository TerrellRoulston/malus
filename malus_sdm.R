# Top ---------------------------------------------------------------------
# MaxEnt Species Distribution Modeling (SDM) for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024

library(tidyverse)
library(terra) 
library(predicts)
library(geodata)
library(rJava)



# Load data for SDM -------------------------------------------------------

getwd()
setwd("./occ_data/")

# Predictors values in dataframes
cor_sdmData <- readRDS(file = 'cor_sdmData.Rdata') # M. coronaria
cor_pb <- cor_sdmData[1] %>% # save object with 0/1 background/presence vector
  unlist()
fus_sdmData <- readRDS(file = 'fus_sdmData.Rdata') # M. fusca
fus_pb <- fus_sdmData[1] %>% # save object with 0/1 background/presence vector
  unlist()

# Background points in SpatVectors
setwd("../occ_data/")
cor_bg_vec <- readRDS(file = 'cor_bg_vec.Rdata')
fus_bg_vec <- readRDS(file = 'fus_bg_vec.Rdata')

# Presence Points in SpatVectors
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata')
occThin_fus <- readRDS(file = 'occThin_fus.Rdata')

# Download/load basemaps
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = '../occ_data/base_maps') # Canada basemap w. Provinces

canUS_map <- rbind(us_map, ca_map) # Combine US and Canada vector map


NA_ext <- ext(-180, -40, 18, 85)

canUS_map <- crop(canUS_map, NA_ext) #crop to western hemisphere

plot(canUS_map)# plot basemap


# Worldclim Predictors
setwd("../wclim_data/")
# Note DO NOT PUSH wclim data**
wclim <- geodata::worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "../wclim_data/")

# Crop wclim data to North America
wclimNA <- crop(wclim, canUS_map, mask = T)
# Sub-sample Test and Training Occurrences
occTrain_fus
occTest_fus

# Model from Predictors in Raster -------------------------------------------------------
# M. coronaria

cor_model <- predicts::MaxEnt(x = wclimNA, p = occThin_cor, a = cor_bg_vec)

# M. fusca

fus_model <- predicts::MaxEnt(x = wclimNA, p = occThin_fus, a = fus_bg_vec)

# Run prediction in parrallel to help speed up computation
# Terra implements parallel functions natively
# But load library for additional functions like <decectCores()>
library(parallel)
cn <- detectCores(logical = F) # number of physical RAM cores in your computure

# M. coronaria
cor_pred <- terra::predict(wclimNA, cor_model, cores = cn - 1) # leave one core for other processes, browser etc.
#plot
plot(cor_pred, main = 'Malus coronaria (Historical)')
points(occThin_cor, cex = 0.5) # Overlay thinned occurrences


# M. fusca
fus_pred <- terra::predict(wclimNA, fus_model, cores = cn - 1) # leave one core for other processes, browser etc.
#plot
plot(fus_pred, main = 'Malus fusca (Historical)')
points(occThin_fus, cex = 0.5) # Overlay thinned occurrences


# Model from predictors values in Dataframes ---------------------------------------------------
# ***NOTE***
# See dendograms (in malus_bg) for reference
# Select variables that are not colinear together, or one var from each colinear group
# See https://www.worldclim.org/data/bioclim.html for definition of variables from codes


# M. coronaria
names(cor_sdmData) # peak at predictor variables for assistance
cor_sdmData.x <- cor_sdmData %>% 
  dplyr::select('wc2.1_2.5m_bio_8', 'wc2.1_2.5m_bio_2', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
#BIO8 = Mean Temperature of Wettest Quarter; BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO12 =  Annual Mean Precipitation; BIO1 = Annual Mean Temperature

# M. fusca
names(fus_sdmData) # peak at predictor variables for assistance
fus_sdmData.x <- fus_sdmData %>% 
  dplyr::select('wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_8', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
#BIO15 = Precipitation Seasonality (Coefficient of Variation); BIO8 = Mean Temperature of Wettest Quarter
#BIO12 =  Annual Mean Precipitation; BIO1 = Annual Mean Temperature


# model using dataframes
cor_max <- predicts::MaxEnt(x = cor_sdmData.x, p = cor_pb)


