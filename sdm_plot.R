# Top ---------------------------------------------------------------------
# SDM plotting and map making
# Started May 7th, 2024

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package
library(geodata) # Basemaps


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

# Predicted habitat suitability rasters
# M. coronaria
setwd('../sdm_output')
cor_pred_hist <- readRDS(file = 'cor_pred_hist.Rdata')

cor_pred_ssp245_30 <- readRDS(file = 'cor_pred_ssp245_30.Rdata')
cor_pred_ssp245_50 <- readRDS(file = 'cor_pred_ssp245_50.Rdata')
cor_pred_ssp245_70 <- readRDS(file = 'cor_pred_ssp245_70.Rdata')

cor_pred_ssp585_30 <- readRDS(file = 'cor_pred_ssp585_30.Rdata')
cor_pred_ssp585_50 <- readRDS(file = 'cor_pred_ssp585_50.Rdata')
cor_pred_ssp585_70 <- readRDS(file = 'cor_pred_ssp585_70.Rdata')

#M. fusca
fus_pred_hist <- readRDS(file = 'fus_pred_hist.Rdata')

fus_pred_ssp245_30 <- readRDS(file = 'fus_pred_ssp245_30.Rdata')
fus_pred_ssp245_50 <- readRDS(file = 'fus_pred_ssp245_50.Rdata')
fus_pred_ssp245_70 <- readRDS(file = 'fus_pred_ssp245_70.Rdata')

fus_pred_ssp585_30 <- readRDS(file = 'fus_pred_ssp585_30.Rdata')
fus_pred_ssp585_50 <- readRDS(file = 'fus_pred_ssp585_50.Rdata')
fus_pred_ssp585_70 <- readRDS(file = 'fus_pred_ssp585_70.Rdata')

# Thresholds
# M. coronaria
setwd('../sdm_output/thresholds')
corPred_threshold_1 <- readRDS(file = 'corPred_threshold_1.Rdata')
corPred_threshold_10 <- readRDS(file = 'corPred_threshold_10.Rdata')
corPred_threshold_50 <- readRDS(file = 'corPred_threshold_50.Rdata')

#M. fusca
fusPred_threshold_1 <- readRDS(file = 'fusPred_threshold_1.Rdata')
fusPred_threshold_10 <- readRDS(file = 'fusPred_threshold_10.Rdata')
fusPred_threshold_50 <- readRDS(file = 'fusPred_threshold_50.Rdata')



# Project to Lambert Conformal Conic --------------------------------------

projLam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Basemaps
canUSMex_map.lcc <- project(canUSMex_map, projLam)
canUSMex_map_0.lcc <- project(canUSMex_map_0, projLam)

# Great lakes
great_lakes.lcc <- project(great_lakes, projLam)

# Occurrences
occThin_cor.lcc <- project(occThin_cor, projLam)
occThin_fus.lcc <- project(occThin_fus, projLam)

# Habitat suitability
cor_pred_hist.lcc <- project(cor_pred_hist, projLam)

cor_pred_ssp245_30.lcc <- project(cor_pred_ssp245_30, projLam)
cor_pred_ssp245_50.lcc <- project(cor_pred_ssp245_50, projLam)
cor_pred_ssp245_70.lcc <- project(cor_pred_ssp245_70, projLam)

cor_pred_ssp585_30.lcc <- project(cor_pred_ssp585_30, projLam)
cor_pred_ssp585_50.lcc <- project(cor_pred_ssp585_50, projLam)
cor_pred_ssp585_70.lcc <- project(cor_pred_ssp585_70, projLam)

#M. fusca
fus_pred_hist.lcc <- project(fus_pred_hist, projLam)

fus_pred_ssp245_30.lcc <- project(fus_pred_ssp245_30, projLam)
fus_pred_ssp245_50.lcc <- project(fus_pred_ssp245_50, projLam)
fus_pred_ssp245_70.lcc <- project(fus_pred_ssp245_70, projLam)

fus_pred_ssp585_30.lcc <- project(fus_pred_ssp585_30, projLam)
fus_pred_ssp585_50.lcc <- project(fus_pred_ssp585_50, projLam)
fus_pred_ssp585_70.lcc <- project(fus_pred_ssp585_70, projLam)


# M. coronaria future habitat plot ----------------------------------------
# Predicted historical distribtuion
legend_labs <- c('Low Suitability (1st percentile)', 'Moderate Suitability (10th percentile)', 'High Suitability (50th percentile)')
fill_cols <- c('#D81B60', '#1E88E5', '#FFC107')

cor.xlim <- c(-5*10^5, 3.1*10^6)
cor.ylim <- c(-2*10^6, 2*10^6)

# Plot a legend that can be added on its own
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', xpd = NA, title = c(as.expression(bquote(bold('Habitat Suitability')))), legend = legend_labs, fill = fill_cols, cex = 2)



terra::plot(cor_pred_hist.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = 'Historical (1970-2000)',
            cex.main = 1.5,
            axes = F,
            box = T,
            mar = c(1, 0, 2, 0))
terra::plot(cor_pred_hist.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_hist.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)
#points(occThin_cor.lcc, col = 'black', cex = 0.3, pch = 4)


# Make a figure with three plots next to one another
par(mfrow = c(1, 3))

# SSP245
# 2030

terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2021-2040',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0.1, 0, 0))
terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2050

terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2041-2060',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 0, 0))
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2070

terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2061-2080',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(0, 0, 0,  0.1))
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# SSP585
par(mfrow = c(1, 3))
# 2030

terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2021-2040',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 0, 0))
terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2050

terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2041-2060',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 0, 0))
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2070

terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2061-2080',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(0, 0, 0,  0))
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)


# M. fusca future habitat plot --------------------------------------------

# Predicted historical distribtuion
legend_labs <- c('Low Suitability (1st percentile)', 'Moderate Suitability (10th percentile)', 'High Suitability (50th percentile)')
fill_cols <- c('#D81B60', '#1E88E5', '#FFC107')

fus.xlim <- c(-4*10^6, -1*10^6)
fus.ylim <- c(-2*10^6, 3*10^6)

fus.bc.xlim <- c(-2.5*10^6, -1.5*10^6)
fus.bc.ylim <- c(3.5*10^5, 2*10^6)

# Plot a legend that can be added on its own
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', xpd = NA, title = c(as.expression(bquote(bold('Habitat Suitability')))), legend = legend_labs, fill = fill_cols, cex = 2)



terra::plot(fus_pred_hist.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = 'Historical (1970-2000)',
            cex.main = 1.5,
            axes = F,
            box = T,
            mar = c(1, 5, 2, 5))
terra::plot(fus_pred_hist.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_hist.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# PNW PLOT
terra::plot(fus_pred_hist.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.bc.xlim, ylim = fus.bc.ylim, 
            main = 'PNW Historical (1970-2000)',
            cex.main = 1.5,
            axes = F,
            box = T,
            mar = c(1, 5, 2, 5))
terra::plot(fus_pred_hist.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_hist.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)
#


# Make a figure with three plots next to one another
par(mfrow = c(1, 3))

# SSP245
# 2030

terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '2021-2040',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 2, 0))
terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2050

terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '2041-2060',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 2, 0))
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2070

terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '2061-2080',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(0, 0, 2, 0))
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# SSP585
par(mfrow = c(1, 3))
# 2030

terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '2021-2040',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 0, 0))
terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2050

terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '2041-2060',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(0, 0, 0, 0))
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)

# 2070

terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '2061-2080',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(0, 0, 0,  0))
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map.lcc, add = T)


