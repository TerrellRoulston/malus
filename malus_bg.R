# Top ---------------------------------------------------------------------
# Setting Background area and predictor vars for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024

library(tidyverse)
library(dismo)
library(terra) 
library(raster)
library(sp)
library(predicts)
library(geodata)
library(ENMTools)


# Ecoregion prep ----------------------------------------------------------
# Download NA Ecoregion shapefile from: https://www.epa.gov/eco-re 
# Load shapefile from local files
ecoNA <- vect(x = "C:/Users/terre/Documents/UBC/Botanical Garden/Malus Project/maps/eco regions/na_cec_eco_l2/", layer = 'NA_CEC_Eco_Level2')
ecoNA <- project(ecoNA, 'WGS84') # project ecoregion vector to same coords ref as basemap


# download/load maps
getwd()
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = '../occ_data/base_maps') #Canada basemap w. Provinces

canUS_map <- rbind(us_map, ca_map) #combine US and Canada vector map

# plot basemap
plot(canUS_map, xlim = c(-180, -50))
# plot ecoregions 
lines(ecoNA, col = 'red')

dev.off()

# load occurrence data
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata')
occThin_fus <- readRDS(file = 'occThin_fus.Rdata')

# Load WorldClim data
setwd("../wclim_data/")
# Note DO NOT PUSH wclim data**
wclim <- geodata::worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "../wclim_data/")


# M. coronaria eco regions ------------------------------------------------
# extract eco region polygon that contain M. coronaria occurrence points
eco_cor <- extract(ecoNA, occThin_cor) # extract what polygons contained points 

# return vector of eco region codes of the polygons that contain occurrences
eco_cor_code <- eco_cor$NA_L2CODE %>% unique() 
eco_cor_code <- eco_cor_code[eco_cor_code != '0.0']  #remove the 'water' '0.0' ecoregion

#CODES: "8.1" "8.2" "5.3" "8.4" "8.3" "9.4" "8.5" "5.2"

ecoNA_cor <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_cor_code) # subset eco region spat vector by the codes

plot(ecoNA_cor) # plot the subseted M. coronaria eco regions
points(occThin_cor, pch = 3, col = 'red') # plot M. coronaria points

# M. fusca eco regions ----------------------------------------------------
# # extract eco region polygon that contain M. fusca occurrence points
eco_fus <- extract(ecoNA, occThin_fus) # extract what polygons contained points

# return vector of eco region codes of the polygons that contain occurrences
eco_fus_code <- eco_fus$NA_L2CODE %>% unique() 
eco_fus_code <- eco_fus_code[eco_fus_code != '0.0'] # remove NA value

#CODES "7.1""6.2"  "10.1" "11.1" "10.2" "3.1" 

ecoNA_fus <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_fus_code)

plot(ecoNA_fus, col = 'red')
points(occThin_fus, pch = 3, col = 'red') # plot M. coronaria points

setwd('../occ_data/eco_regions')
saveRDS(ecoNA_cor, file = 'ecoNA_cor.Rdata')
saveRDS(ecoNA_fus, file = 'ecoNA_fus.Rdata')


# Crop WorldClim to Ecoregions and Create Background ----------------------

# crop+mask extent of WorldClim data to the Malus ecoregions
wclim_cor <- terra::crop(wclim, ecoNA_cor, mask = T)
wclim_fus <- terra::crop(wclim, ecoNA_fus, mask = T)

# Save cropped wclim data for downsteam SDM workflow
saveRDS(wclim_cor, file = 'wclim_cor.Rdata')
saveRDS(wclim_fus, file = 'wclim_fus.Rdata')


set.seed(1337) # set a seed to ensure consistent results

# NOTE: Set arguments as.raster = T to return raster
# OR as.points to return spatvector = T to return spatvector

# M. coronaria background
# SpatVector

# Note upped bg points from 5000 to 20000 to be more suitable to better reflect a mean probability of presence 1 - a/2
cor_bg_vec <- spatSample(wclim_cor, 20000, 'random', na.rm = T, as.points = T) #ignore NA values
plot(wclim_cor[[1]])
points(cor_bg_vec, cex = 0.01)

expanse(wclim_cor[[1]], unit = 'km') # total area of raster in km^2
# 5683684 km^2
# 5000/5683684 = 0.000879 samples/km

# M. fusca background
# SpatVector
fus_bg_vec <- spatSample(wclim_fus, 20000, 'random', na.rm = T, as.points = T) #ignore NA values
plot(wclim_fus[[1]])
points(fus_bg_vec, cex = 0.01)

expanse(wclim_fus[[1]], unit = 'km') # total area of raster in km^2
# 4659175 km^2
# 5000/4659175 = 0.00107 samples/km

# Save background SpatVectors
setwd("../occ_data/")
saveRDS(cor_bg_vec, file = 'cor_bg_vec.Rdata')
saveRDS(fus_bg_vec, file = 'fus_bg_vec.Rdata')


# Load background SpatVectors
setwd("../occ_data/")
cor_bg_vec <- readRDS(file = 'cor_bg_vec.Rdata')
fus_bg_vec <- readRDS(file = 'fus_bg_vec.Rdata')



# Extracting presence-background raster values ----------------------------

cor_predvals <- extract(wclim_cor, occThin_cor) # M. coronaria
cor_predvals <- cor_predvals[-1] # drop ID column

fus_predvals <- extract(wclim_fus, occThin_fus) # M. fusca
fus_predvals <- fus_predvals[-1] # drop ID column

cor_bgvals <- values(cor_bg) # Extract raster values for bg points
fus_bgvals <- values(fus_bg) # Extract raster values for bg points


# Create a df for presence-background raster values for SDM ---------------

cor_pb <- c(rep(1, nrow(cor_predvals)), rep(0, nrow(cor_bgvals))) #T/F presence or background string
fus_pb <- c(rep(1, nrow(fus_predvals)), rep(0, nrow(fus_bgvals))) #T/F presence or background string

# combine presence and background dataframes for SDM

cor_sdmData <- data.frame(cbind(cor_pb, rbind(cor_predvals, cor_bgvals)))
fus_sdmData <- data.frame(cbind(fus_pb, rbind(fus_predvals, cor_bgvals)))

#Save df for downstream SDM work
setwd("../occ_data/")
saveRDS(cor_sdmData, file = 'cor_sdmData.Rdata')
saveRDS(fus_sdmData, file = 'fus_sdmData.Rdata')

# Check for colinearity of predictor varirables for presence-bg -----------
# Want to select variables that are not colinear to avoid issues with model fitting
# Dendograms useful for indentifying groupings of variables
# Select one variable from each group to use in modeling


# M. coronaria
pairs(cor_sdmData[,2:5])
pairs(cor_sdmData[,6:9])
pairs(cor_sdmData[,10:13])
pairs(cor_sdmData[,14:17])
pairs(cor_sdmData[,17:20])

pairs(cor_sdmData[,-1])

# Dendogram cluster of predictor colinearity
threshold <- 0.7 # set the threshold for colinearity 
cor_cors <- raster.cor.matrix(wclim_cor)
cor_dist <- as.dist(1 - abs(cor_cors)) # calculate distance of predictors

cor_clust <- hclust(cor_dist, method = 'single') # calculate cluster dendogram
cor_groups <- cutree(cor_clust, h = 1 - threshold) #calculate groupings of variables

plot(cor_clust, hang = -1, main = "M. coronaria Predictors")
rect.hclust(cor_clust, h = 1 - threshold)


#M. fusca
pairs(fus_sdmData[,2:5])
pairs(fus_sdmData[,6:9])
pairs(fus_sdmData[,10:13])
pairs(fus_sdmData[,14:17])
pairs(fus_sdmData[,17:20])

pairs(fus_sdmData[,-1])

# Dendogram cluster of predictor colinearity
threshold <- 0.7 # set the threshold for colinearity 
fus_cors <- raster.cor.matrix(wclim_fus)
fus_dist <- as.dist(1 - abs(fus_cors)) # calculate distance of predictors

fus_clust <- hclust(fus_dist, method = 'single') # calculate cluster dendogram
fus_groups <- cutree(fus_clust, h = 1 - threshold) #calculate groupings of variables

plot(fus_clust, hang = -1, main = "M. fusca Predictors")
rect.hclust(fus_clust, h = 1 - threshold)



# Predictor Kernel Density Plots -------------------------------------------------
# It is helpful to visualize two predictors pairs of
# Presence points and background points


cor_occ.temp <- cor_sdmData %>% filter(cor_pb == 1) %>% # Presence points
  dplyr::select(wc2.1_2.5m_bio_1) %>% # Mean annual temp
  drop_na() %>% 
  unlist()
cor_bg.temp <- cor_sdmData %>% filter(cor_pb == 0) %>% # Background points
  dplyr::select(wc2.1_2.5m_bio_1) %>% # Mean annual temp
  drop_na() %>% 
  unlist()

cor_occ.perc <- cor_sdmData %>% filter(cor_pb == 1) %>% # Presence points
  dplyr::select(wc2.1_2.5m_bio_12) %>% # Annual precipitation
  drop_na() %>% 
  unlist()

cor_bg.perc <- cor_sdmData %>% filter(cor_pb == 0) %>% # Background points
  dplyr::select(wc2.1_2.5m_bio_12) %>% # Annual precipitation
  drop_na() %>% 
  unlist()


library(plotly) # 3D surface Kernel bivariate plots
library(MASS)

cor_occ.3d <- kde2d(cor_occ.temp, cor_occ.perc)
cor_bg.3d <- kde2d(cor_bg.temp, cor_bg.perc)

#Plot 3D surface Kernel density estimation

plot_cor.occ_3d <- plot_ly(x=cor_occ.3d$x, y=cor_occ.3d$y, z=cor_occ.3d$z) %>% 
  add_surface() %>% 
  layout(scene = list(xaxis = list(title = 'Mean Annual Temp (C)', autotick = F, nticks = 5, tickvals = list(0,5,10,15,20)), 
                      yaxis = list(title = 'Mean Annual Percip. (mm)', tick0=0, tick1=2000, dtick=200), 
                      zaxis = list(title = 'Kernel Density', tick0=0, tick1=0.001, dtick=0.0002)),
         title = list(text = "<i>M. coronaria<i> Occurrence Points", 
                      y = 0.95, x = 0.5, 
                      xanchor = 'center', 
                      yanchor = 'top'))

plot_cor.bg_3d <- plot_ly(x=cor_bg.3d$x, y=cor_bg.3d$y, z=cor_bg.3d$z) %>% 
  add_surface() %>% 
  layout(scene = list(xaxis = list(title = 'Mean Annual Temp (C)', tick0=0, tick1=20, dtick=5), 
                      yaxis = list(title = 'Mean Annual Percip. (mm)', tick0=0, tick1=2000, dtick=200), 
                      zaxis = list(title = 'Kernel Density')),
         title = list(text = "<i>M. coronaria<i> Background Points", 
                      y = 0.95, x = 0.5, 
                      xanchor = 'center', 
                      yanchor = 'top'))



