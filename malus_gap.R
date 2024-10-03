# Top ---------------------------------------------------------------------
# Conservation Gap Analysis
# Started September 26th, 2024

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package
library(geodata) # basemaps
library(ggpubr) # arrange ggplot figures


# The following gap analysis is referencing Carver et al. (2021)
# https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.05430

# For the purposes of the analysis we are restricting the study area to just Canada


# Protected areas ---------------------------------------------------------
# Protected Area and OECM data from https://www.canada.ca/en/environment-climate-change/services/national-wildlife-areas/protected-conserved-areas-database.html
# Download .kmz file and vector using terra

pro_area <- terra::vect("C:/Users/terre/Documents/Acadia/Malus Project/maps/canada_pa/ProtectedConservedArea_2023/ProtectedConservedArea_2023/ProtectedConservedArea_2023.gdb")
pro_area <- pro_area[pro_area$BIOME == 'T'] # filter only the terrestrial protected areas and OECMs
pro_area <- pro_area[pro_area$PA_OECM_DF == '1'] # filter only protected areas (remove OECMs)
pro_area <- project(pro_area, 'WGS84') # project the PAs to be the same CRS as the SDM layers

# Take a peak
dev.new()
terra::plot(pro_area)

# Load suitable habitat data and occurrences ------------------------------
# I am using the complete suitability rasters and thresholds separately
# Could also load pre-categorized layers as well

# M. coronaria
getwd()
setwd('./sdm_output/')

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

# Load occurrence points
getwd()
setwd("./occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = 'occThin_fus.Rdata') # M. fusca


# Crop suitable habitat data to Canada ONLY -------------------------------
# Load a Canada Admin boundary
getwd()
setwd('../occ_data')
ca_bound <- gadm(country = 'CA', level = 0, resolution = 1,
               path = '../occ_data/base_maps')  

# Mask and crop SDM layer to Canada

# M. coronaria
# Historical
can_cor_pred_hist <- crop(cor_pred_hist, ca_bound, mask = T) 

# SSP245
can_cor_pred_ssp245_30 <- crop(cor_pred_ssp245_30, ca_bound, mask = T)
can_cor_pred_ssp245_50 <- crop(cor_pred_ssp245_50, ca_bound, mask = T)
can_cor_pred_ssp245_70 <- crop(cor_pred_ssp245_70, ca_bound, mask = T)

# SSP585
can_cor_pred_ssp585_30 <- crop(cor_pred_ssp585_30, ca_bound, mask = T)
can_cor_pred_ssp585_50 <- crop(cor_pred_ssp585_50, ca_bound, mask = T)
can_cor_pred_ssp585_70 <- crop(cor_pred_ssp585_70, ca_bound, mask = T)


# M. fusca
# Historical
can_fus_pred_hist <- crop(fus_pred_hist, ca_bound, mask = T) 

# SSP245
can_fus_pred_ssp245_30 <- crop(fus_pred_ssp245_30, ca_bound, mask = T)
can_fus_pred_ssp245_50 <- crop(fus_pred_ssp245_50, ca_bound, mask = T)
can_fus_pred_ssp245_70 <- crop(fus_pred_ssp245_70, ca_bound, mask = T)

# SSP585
can_fus_pred_ssp585_30 <- crop(fus_pred_ssp585_30, ca_bound, mask = T)
can_fus_pred_ssp585_50 <- crop(fus_pred_ssp585_50, ca_bound, mask = T)
can_fus_pred_ssp585_70 <- crop(fus_pred_ssp585_70, ca_bound, mask = T)

# Mask suitable habitat area to PA area -----------------------------------
# Mask Suitable Habitat rasters to the PA polygonss

#M. coronaria
#Historical 
pa_cor_pred_hist <- crop(can_cor_pred_hist, pro_area, mask = T) 

# SSP245
pa_cor_pred_ssp245_30 <- crop(cor_pred_ssp245_30, pro_area, mask = T)
pa_cor_pred_ssp245_50 <- crop(cor_pred_ssp245_50, pro_area, mask = T) 
pa_cor_pred_ssp245_70 <- crop(cor_pred_ssp245_70, pro_area, mask = T) 

# SSP585
pa_cor_pred_ssp585_30 <- crop(cor_pred_ssp585_30, pro_area, mask = T) 
pa_cor_pred_ssp585_50 <- crop(cor_pred_ssp585_50, pro_area, mask = T) 
pa_cor_pred_ssp585_70 <- crop(cor_pred_ssp585_70, pro_area, mask = T) 


# M. fusca
# Historical
pa_fus_pred_hist <- crop(can_fus_pred_hist, pro_area, mask = T) 

# SSP245
pa_fus_pred_ssp245_30 <- crop(fus_pred_ssp245_30, pro_area, mask = T)
pa_fus_pred_ssp245_50 <- crop(fus_pred_ssp245_50, pro_area, mask = T) 
pa_fus_pred_ssp245_70 <- crop(fus_pred_ssp245_70, pro_area, mask = T) 

# SSP585
pa_fus_pred_ssp585_30 <- crop(fus_pred_ssp585_30, pro_area, mask = T) 
pa_fus_pred_ssp585_50 <- crop(fus_pred_ssp585_50, pro_area, mask = T) 
pa_fus_pred_ssp585_70 <- crop(fus_pred_ssp585_70, pro_area, mask = T) 


# SRSin score -------------------------------------------------------------
# Sampling representativeness score in situ

# Extract occurrence points that fall within the SDM prediction

# M. coronaria
# Historical
index_cor_pos_hist_low <- extract(can_cor_pred_hist > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_hist_mod <- extract(can_cor_pred_hist > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_hist_high <- extract(can_cor_pred_hist > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

# SSP245
index_cor_pos_ssp245_30_low <- extract(can_cor_pred_ssp245_30 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_ssp245_30_mod <- extract(can_cor_pred_ssp245_30 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_ssp245_30_high <- extract(can_cor_pred_ssp245_30 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_cor_pos_ssp245_50_low <- extract(can_cor_pred_ssp245_50 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_ssp245_50_mod <- extract(can_cor_pred_ssp245_50 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_ssp245_50_high <- extract(can_cor_pred_ssp245_50 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_cor_pos_ssp245_70_low <- extract(can_cor_pred_ssp245_70 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_ssp245_70_mod <- extract(can_cor_pred_ssp245_70 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_ssp245_70_high <- extract(can_cor_pred_ssp245_70 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

# SSP585
index_cor_pos_ssp585_30_low <- extract(can_cor_pred_ssp585_30 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_ssp585_30_mod <- extract(can_cor_pred_ssp585_30 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_ssp585_30_high <- extract(can_cor_pred_ssp585_30 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_cor_pos_ssp585_50_low <- extract(can_cor_pred_ssp585_50 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_ssp585_50_mod <- extract(can_cor_pred_ssp585_50 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_ssp585_50_high <- extract(can_cor_pred_ssp585_50 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_cor_pos_ssp585_70_low <- extract(can_cor_pred_ssp585_70 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_cor_pos_ssp585_70_mod <- extract(can_cor_pred_ssp585_70 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_cor_pos_ssp585_70_high <- extract(can_cor_pred_ssp585_70 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

# M. fusca
# Historical
index_fus_pos_hist_low <- extract(can_fus_pred_hist > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_hist_mod <- extract(can_fus_pred_hist > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_hist_high <- extract(can_fus_pred_hist > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

# SSP245
index_fus_pos_ssp245_30_low <- extract(can_fus_pred_ssp245_30 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_ssp245_30_mod <- extract(can_fus_pred_ssp245_30 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_ssp245_30_high <- extract(can_fus_pred_ssp245_30 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_fus_pos_ssp245_50_low <- extract(can_fus_pred_ssp245_50 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_ssp245_50_mod <- extract(can_fus_pred_ssp245_50 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_ssp245_50_high <- extract(can_fus_pred_ssp245_50 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_fus_pos_ssp245_70_low <- extract(can_fus_pred_ssp245_70 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_ssp245_70_mod <- extract(can_fus_pred_ssp245_70 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_ssp245_70_high <- extract(can_fus_pred_ssp245_70 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

# SSP585
index_fus_pos_ssp585_30_low <- extract(can_fus_pred_ssp585_30 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_ssp585_30_mod <- extract(can_fus_pred_ssp585_30 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_ssp585_30_high <- extract(can_fus_pred_ssp585_30 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_fus_pos_ssp585_50_low <- extract(can_fus_pred_ssp585_50 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_ssp585_50_mod <- extract(can_fus_pred_ssp585_50 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_ssp585_50_high <- extract(can_fus_pred_ssp585_50 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_fus_pos_ssp585_70_low <- extract(can_fus_pred_ssp585_70 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_fus_pos_ssp585_70_mod <- extract(can_fus_pred_ssp585_70 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_fus_pos_ssp585_70_high <- extract(can_fus_pred_ssp585_70 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)


# Extract that occurrences that fall within the PAs

# M. coronaria
# Historical
index_pa_cor_pos_hist_low <- extract(pa_cor_pred_hist > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_hist_mod <- extract(pa_cor_pred_hist > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_hist_high <- extract(pa_cor_pred_hist > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

# SSP245
index_pa_cor_pos_ssp245_30_low <- extract(pa_cor_pred_ssp245_30 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_ssp245_30_mod <- extract(pa_cor_pred_ssp245_30 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_ssp245_30_high <- extract(pa_cor_pred_ssp245_30 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_pa_cor_pos_ssp245_50_low <- extract(pa_cor_pred_ssp245_50 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_ssp245_50_mod <- extract(pa_cor_pred_ssp245_50 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_ssp245_50_high <- extract(pa_cor_pred_ssp245_50 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_pa_cor_pos_ssp245_70_low <- extract(pa_cor_pred_ssp245_70 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_ssp245_70_mod <- extract(pa_cor_pred_ssp245_70 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_ssp245_70_high <- extract(pa_cor_pred_ssp245_70 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

# SSP585
index_pa_cor_pos_ssp585_30_low <- extract(pa_cor_pred_ssp585_30 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_ssp585_30_mod <- extract(pa_cor_pred_ssp585_30 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_ssp585_30_high <- extract(pa_cor_pred_ssp585_30 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_pa_cor_pos_ssp585_50_low <- extract(pa_cor_pred_ssp585_50 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_ssp585_50_mod <- extract(pa_cor_pred_ssp585_50 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_ssp585_50_high <- extract(pa_cor_pred_ssp585_50 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

index_pa_cor_pos_ssp585_70_low <- extract(pa_cor_pred_ssp585_70 > corPred_threshold_1, occThin_cor) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_cor_pos_ssp585_70_mod <- extract(pa_cor_pred_ssp585_70 > corPred_threshold_10, occThin_cor) %>% filter(lyr1 == T)
index_pa_cor_pos_ssp585_70_high <- extract(pa_cor_pred_ssp585_70 > corPred_threshold_50, occThin_cor) %>% filter(lyr1 == T)

# M. fusca
# Historical
index_pa_fus_pos_hist_low <- extract(pa_fus_pred_hist > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_hist_mod <- extract(pa_fus_pred_hist > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_hist_high <- extract(pa_fus_pred_hist > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

# SSP245
index_pa_fus_pos_ssp245_30_low <- extract(pa_fus_pred_ssp245_30 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_ssp245_30_mod <- extract(pa_fus_pred_ssp245_30 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_ssp245_30_high <- extract(pa_fus_pred_ssp245_30 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_pa_fus_pos_ssp245_50_low <- extract(pa_fus_pred_ssp245_50 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_ssp245_50_mod <- extract(pa_fus_pred_ssp245_50 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_ssp245_50_high <- extract(pa_fus_pred_ssp245_50 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_pa_fus_pos_ssp245_70_low <- extract(pa_fus_pred_ssp245_70 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_ssp245_70_mod <- extract(pa_fus_pred_ssp245_70 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_ssp245_70_high <- extract(pa_fus_pred_ssp245_70 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

# SSP585
index_pa_fus_pos_ssp585_30_low <- extract(pa_fus_pred_ssp585_30 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_ssp585_30_mod <- extract(pa_fus_pred_ssp585_30 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_ssp585_30_high <- extract(pa_fus_pred_ssp585_30 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_pa_fus_pos_ssp585_50_low <- extract(pa_fus_pred_ssp585_50 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_ssp585_50_mod <- extract(pa_fus_pred_ssp585_50 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_ssp585_50_high <- extract(pa_fus_pred_ssp585_50 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

index_pa_fus_pos_ssp585_70_low <- extract(pa_fus_pred_ssp585_70 > fusPred_threshold_1, occThin_fus) %>% filter(lyr1 == T) # filter only the true cases # filter only the true cases
index_pa_fus_pos_ssp585_70_mod <- extract(pa_fus_pred_ssp585_70 > fusPred_threshold_10, occThin_fus) %>% filter(lyr1 == T)
index_pa_fus_pos_ssp585_70_high <- extract(pa_fus_pred_ssp585_70 > fusPred_threshold_50, occThin_fus) %>% filter(lyr1 == T)

# return count of the number of occurrences
# positive occurrences (overall)
# list of obj names from above
occ_names <- list(
  index_cor_pos_hist_low,
  index_cor_pos_hist_mod,
  index_cor_pos_hist_high,
  index_cor_pos_ssp245_30_low,
  index_cor_pos_ssp245_30_mod,
  index_cor_pos_ssp245_30_high,
  index_cor_pos_ssp245_50_low,
  index_cor_pos_ssp245_50_mod,
  index_cor_pos_ssp245_50_high,
  index_cor_pos_ssp245_70_low,
  index_cor_pos_ssp245_70_mod,
  index_cor_pos_ssp245_70_high,
  index_cor_pos_ssp585_30_low,
  index_cor_pos_ssp585_30_mod,
  index_cor_pos_ssp585_30_high,
  index_cor_pos_ssp585_50_low,
  index_cor_pos_ssp585_50_mod,
  index_cor_pos_ssp585_50_high,
  index_cor_pos_ssp585_70_low,
  index_cor_pos_ssp585_70_mod,
  index_cor_pos_ssp585_70_high,
  index_fus_pos_hist_low,
  index_fus_pos_hist_mod,
  index_fus_pos_hist_high,
  index_fus_pos_ssp245_30_low,
  index_fus_pos_ssp245_30_mod,
  index_fus_pos_ssp245_30_high,
  index_fus_pos_ssp245_50_low,
  index_fus_pos_ssp245_50_mod,
  index_fus_pos_ssp245_50_high,
  index_fus_pos_ssp245_70_low,
  index_fus_pos_ssp245_70_mod,
  index_fus_pos_ssp245_70_high,
  index_fus_pos_ssp585_30_low,
  index_fus_pos_ssp585_30_mod,
  index_fus_pos_ssp585_30_high,
  index_fus_pos_ssp585_50_low,
  index_fus_pos_ssp585_50_mod,
  index_fus_pos_ssp585_50_high,
  index_fus_pos_ssp585_70_low,
  index_fus_pos_ssp585_70_mod,
  index_fus_pos_ssp585_70_high
)

# return counts in a vector
SRSin_occ <- occ_names %>% 
  lapply(., nrow) %>% 
  unlist()

# positive occurences in PAs
# list of obj names from above
pa_names <- list(
  index_pa_cor_pos_hist_low,
  index_pa_cor_pos_hist_mod,
  index_pa_cor_pos_hist_high,
  index_pa_cor_pos_ssp245_30_low,
  index_pa_cor_pos_ssp245_30_mod,
  index_pa_cor_pos_ssp245_30_high,
  index_pa_cor_pos_ssp245_50_low,
  index_pa_cor_pos_ssp245_50_mod,
  index_pa_cor_pos_ssp245_50_high,
  index_pa_cor_pos_ssp245_70_low,
  index_pa_cor_pos_ssp245_70_mod,
  index_pa_cor_pos_ssp245_70_high,
  index_pa_cor_pos_ssp585_30_low,
  index_pa_cor_pos_ssp585_30_mod,
  index_pa_cor_pos_ssp585_30_high,
  index_pa_cor_pos_ssp585_50_low,
  index_pa_cor_pos_ssp585_50_mod,
  index_pa_cor_pos_ssp585_50_high,
  index_pa_cor_pos_ssp585_70_low,
  index_pa_cor_pos_ssp585_70_mod,
  index_pa_cor_pos_ssp585_70_high,
  index_pa_fus_pos_hist_low,
  index_pa_fus_pos_hist_mod,
  index_pa_fus_pos_hist_high,
  index_pa_fus_pos_ssp245_30_low,
  index_pa_fus_pos_ssp245_30_mod,
  index_pa_fus_pos_ssp245_30_high,
  index_pa_fus_pos_ssp245_50_low,
  index_pa_fus_pos_ssp245_50_mod,
  index_pa_fus_pos_ssp245_50_high,
  index_pa_fus_pos_ssp245_70_low,
  index_pa_fus_pos_ssp245_70_mod,
  index_pa_fus_pos_ssp245_70_high,
  index_pa_fus_pos_ssp585_30_low,
  index_pa_fus_pos_ssp585_30_mod,
  index_pa_fus_pos_ssp585_30_high,
  index_pa_fus_pos_ssp585_50_low,
  index_pa_fus_pos_ssp585_50_mod,
  index_pa_fus_pos_ssp585_50_high,
  index_pa_fus_pos_ssp585_70_low,
  index_pa_fus_pos_ssp585_70_mod,
  index_pa_fus_pos_ssp585_70_high
)
# return counts in a vector
SRSin_pa <- pa_names %>% 
            lapply(., nrow) %>% 
            unlist() 

# Calculate SRSin score

SRSin_score <- ((SRSin_pa/SRSin_occ) * 100) %>% replace_na(., 0) # some values are NAs as there are 0 in the denominator


# create dataframe to store results and work with in ggplot downstream
# can mutate new columns as needed

in_situ <- data.frame(
  species = c(rep('Malus coronaria', times = 21), rep('Malus fusca', times = 21)),
  ssp = c(rep('historical', times = 3), rep(245, times = 9), rep(585, times = 9), rep('historical', times = 3), rep(245, times = 9), rep(585, times = 9)),
  suitability = rep(c('low', 'moderate', 'high'), times = 14),
  period = c(rep(2000, times = 3), rep(2030, times = 3), rep(2050, times = 3), rep(2070, times = 3), 
             rep(2030, times = 3), rep(2050, times = 3), rep(2070, times = 3), 
             rep(2000, times = 3), rep(2030, times = 3), rep(2050, times = 3), rep(2070, times = 3), 
             rep(2030, times = 3), rep(2050, times = 3), rep(2070, times = 3)),
  SRSin = SRSin_score
)


# GRSin score -------------------------------------------------------------
# Need to begin by calculating the area of the total suitable habitat and that of the protected area containing suitable habitat 

# Overal suitable area
# create a list of rasters to calculate total suitable area for each threshold, time period and SSP

suitable_area_rasts <- list(
can_cor_pred_hist > corPred_threshold_1,
can_cor_pred_hist > corPred_threshold_10,
can_cor_pred_hist > corPred_threshold_50,
can_cor_pred_ssp245_30 > corPred_threshold_1,
can_cor_pred_ssp245_30 > corPred_threshold_10,
can_cor_pred_ssp245_30 > corPred_threshold_50,
can_cor_pred_ssp245_50 > corPred_threshold_1,
can_cor_pred_ssp245_50 > corPred_threshold_10,
can_cor_pred_ssp245_50 > corPred_threshold_50,
can_cor_pred_ssp245_70 > corPred_threshold_1,
can_cor_pred_ssp245_70 > corPred_threshold_10,
can_cor_pred_ssp245_70 > corPred_threshold_50,
can_cor_pred_ssp585_30 > corPred_threshold_1,
can_cor_pred_ssp585_30 > corPred_threshold_10,
can_cor_pred_ssp585_30 > corPred_threshold_50,
can_cor_pred_ssp585_50 > corPred_threshold_1,
can_cor_pred_ssp585_50 > corPred_threshold_10,
can_cor_pred_ssp585_50 > corPred_threshold_50,
can_cor_pred_ssp585_70 > corPred_threshold_1,
can_cor_pred_ssp585_70 > corPred_threshold_10,
can_cor_pred_ssp585_70 > corPred_threshold_50,
can_fus_pred_hist > fusPred_threshold_1,
can_fus_pred_hist > fusPred_threshold_10, 
can_fus_pred_hist > fusPred_threshold_50,
can_fus_pred_ssp245_30 > fusPred_threshold_1,
can_fus_pred_ssp245_30 > fusPred_threshold_10,
can_fus_pred_ssp245_30 > fusPred_threshold_50,
can_fus_pred_ssp245_50 > fusPred_threshold_1,
can_fus_pred_ssp245_50 > fusPred_threshold_10,
can_fus_pred_ssp245_50 > fusPred_threshold_50,
can_fus_pred_ssp245_70 > fusPred_threshold_1,
can_fus_pred_ssp245_70 > fusPred_threshold_10,
can_fus_pred_ssp245_70 > fusPred_threshold_50,
can_fus_pred_ssp585_30 > fusPred_threshold_1,
can_fus_pred_ssp585_30 > fusPred_threshold_10,
can_fus_pred_ssp585_30 > fusPred_threshold_50,
can_fus_pred_ssp585_50 > fusPred_threshold_1,
can_fus_pred_ssp585_50 > fusPred_threshold_10,
can_fus_pred_ssp585_50 > fusPred_threshold_50,
can_fus_pred_ssp585_70 > fusPred_threshold_1,
can_fus_pred_ssp585_70 > fusPred_threshold_10,
can_fus_pred_ssp585_70 > fusPred_threshold_50
)


# the following will return a vector of the total suitable habitat in km^2
# lapply allows you to interate functions over a list
# terra::extract calculates the total area of a raster
# I then initiate a anaymous function to create a dataframe of the listed values so that
# they can be filtered such that only TRUE suitable areas (value = 1) are selected
# and then the area is returned using pull (extracts single column) for all cases that match value = 1
# finally I unlist the temporary df so that the area is returned in a vector

total_area <- lapply(suitable_area_rasts, terra::expanse, byValue = T, unit = 'km') %>% 
  lapply(., function(df) {
    df %>% filter(value == 1) %>% pull(area)
  }) %>% 
  unlist()

# Area of suitable habitat in protected areas

pa_area_rasts <- list(
  pa_cor_pred_hist > corPred_threshold_1,
  pa_cor_pred_hist > corPred_threshold_10,
  pa_cor_pred_hist > corPred_threshold_50,
  pa_cor_pred_ssp245_30 > corPred_threshold_1,
  pa_cor_pred_ssp245_30 > corPred_threshold_10,
  pa_cor_pred_ssp245_30 > corPred_threshold_50,
  pa_cor_pred_ssp245_50 > corPred_threshold_1,
  pa_cor_pred_ssp245_50 > corPred_threshold_10,
  pa_cor_pred_ssp245_50 > corPred_threshold_50,
  pa_cor_pred_ssp245_70 > corPred_threshold_1,
  pa_cor_pred_ssp245_70 > corPred_threshold_10,
  pa_cor_pred_ssp245_70 > corPred_threshold_50,
  pa_cor_pred_ssp585_30 > corPred_threshold_1,
  pa_cor_pred_ssp585_30 > corPred_threshold_10,
  pa_cor_pred_ssp585_30 > corPred_threshold_50,
  pa_cor_pred_ssp585_50 > corPred_threshold_1,
  pa_cor_pred_ssp585_50 > corPred_threshold_10,
  pa_cor_pred_ssp585_50 > corPred_threshold_50,
  pa_cor_pred_ssp585_70 > corPred_threshold_1,
  pa_cor_pred_ssp585_70 > corPred_threshold_10,
  pa_cor_pred_ssp585_70 > corPred_threshold_50,
  pa_fus_pred_hist > fusPred_threshold_1,
  pa_fus_pred_hist > fusPred_threshold_10, 
  pa_fus_pred_hist > fusPred_threshold_50,
  pa_fus_pred_ssp245_30 > fusPred_threshold_1,
  pa_fus_pred_ssp245_30 > fusPred_threshold_10,
  pa_fus_pred_ssp245_30 > fusPred_threshold_50,
  pa_fus_pred_ssp245_50 > fusPred_threshold_1,
  pa_fus_pred_ssp245_50 > fusPred_threshold_10,
  pa_fus_pred_ssp245_50 > fusPred_threshold_50,
  pa_fus_pred_ssp245_70 > fusPred_threshold_1,
  pa_fus_pred_ssp245_70 > fusPred_threshold_10,
  pa_fus_pred_ssp245_70 > fusPred_threshold_50,
  pa_fus_pred_ssp585_30 > fusPred_threshold_1,
  pa_fus_pred_ssp585_30 > fusPred_threshold_10,
  pa_fus_pred_ssp585_30 > fusPred_threshold_50,
  pa_fus_pred_ssp585_50 > fusPred_threshold_1,
  pa_fus_pred_ssp585_50 > fusPred_threshold_10,
  pa_fus_pred_ssp585_50 > fusPred_threshold_50,
  pa_fus_pred_ssp585_70 > fusPred_threshold_1,
  pa_fus_pred_ssp585_70 > fusPred_threshold_10,
  pa_fus_pred_ssp585_70 > fusPred_threshold_50
)

pa_area <- lapply(pa_area_rasts, terra::expanse, byValue = T, unit = 'km') %>% 
  lapply(., function(df) {
    df %>% filter(value == 1) %>% pull(area)
  }) %>% 
  unlist()

# calculate the GRSin score
GRSin_score <- (pa_area/total_area) * 100 

# add it to the data frame from abouve

in_situ <- in_situ %>% mutate(GRSin = GRSin_score)


# ERSin -------------------------------------------------------------------
# We will use the same ecoregions we used for defining the background

# Download NA Ecoregion shapefile from: https://www.epa.gov/eco-research/ecoregions-north-america
# Load shapefile from local files
ecoNA <- vect(x = "C:/Users/terre/Documents/Acadia/Malus Project/maps/eco regions/na_cec_eco_l2/", layer = 'NA_CEC_Eco_Level2')
ecoNA <- project(ecoNA, 'WGS84') # project ecoregion vector to same coords ref as basemap



# suitable habitat eco regions
suitable_area_points <- lapply(suitable_area_rasts, terra::as.points) %>% 
  lapply(., function (x) {
    x[x$lyr1 == 1, ]
  }) 


start <- Sys.time()
suitable_area_eco <- lapply(suitable_area_points, function(point) {
  terra::extract(x = ecoNA, y = point)
})
print( Sys.time() - start )

# took 12.93733 hours to complete

# protected area ecoregions
pa_area_points <- lapply(pa_area_rasts, terra::as.points) %>% 
        lapply(., function (x) {
          x[x$lyr1 == 1, ]
        }) 


start <- Sys.time()
pa_area_eco <- lapply(pa_area_points, function(point) {
  terra::extract(x = ecoNA, y = point)
})
print( Sys.time() - start )

# took 3.679289 hours to complete


# SAVE/LOAD eco regiion results 
getwd()
setwd('../gap_analysis')
saveRDS(suitable_area_eco, file = 'suitable_area_eco.Rdata')
saveRDS(pa_area_eco, file = 'pa_area_eco.Rdata')

suitable_area_eco <- readRDS(file = 'suitable_area_eco.Rdata')
pa_area_eco <- readRDS('pa_area_eco.Rdata')

# return unqiue names and number of eco regions

suit_area_eco_num <- suitable_area_eco %>% 
                      bind_rows(.id = 'source_list') %>% # bind the lists of lists together 
                      unnest(cols = c(NA_L2NAME)) %>% # unnest the lists into columns of a dataframe
                      distinct(source_list, NA_L2CODE) %>% # return only the unique eco regions grouped by source list (aka the raster layers)
                      filter(!NA_L2CODE %in% c(NA, '0.0')) %>% # remove NA and 0.0 (water) valyes
                      add_count(source_list, name = 'num_eco') %>%  # count the number of nrow for each source list (or number of eco regions)
                      distinct(source_list, num_eco) # return the number of eco regions per source list



pa_area_eco_num <- pa_area_eco %>% 
                    bind_rows(.id = 'source_list') %>% # bind the lists of lists together 
                    unnest(cols = c(NA_L2NAME)) %>% # unnest the lists into columns of a dataframe
                    distinct(source_list, NA_L2CODE) %>% # return only the unique eco regions grouped by source list (aka the raster layers)
                    filter(!NA_L2CODE %in% c(NA, '0.0')) %>% # remove NA and 0.0 (water) valyes
                    add_count(source_list, name = 'num_eco') %>% # count the number of nrow for each source list (or number of eco regions)
                    distinct(source_list, num_eco) # return the number of eco regions per source list



# Now calculate the ERSin score  
ERSin_score <- (pa_area_eco_num$num_eco / suit_area_eco_num$num_eco) * 100

# Mutate the ERSin score to the master df
in_situ <- in_situ %>% mutate(ERSin = ERSin_score)


# FCSin -------------------------------------------------------------------
# Calculate the final conservation score for in situ conservation 

in_situ <- in_situ %>% mutate(FCSin = ((SRSin + GRSin + ERSin)/3))


# Plot data ---------------------------------------------------------------

# Let's create some new dfs to make plotting a bit easier

cor_ssp245_in_situ <- in_situ %>% dplyr::filter(species == 'Malus coronaria' & ssp %in% c('historical', '245'))
cor_ssp585_in_situ <- in_situ %>% dplyr::filter(species == 'Malus coronaria' & ssp %in% c('historical', '585'))

fus_ssp245_in_situ <- in_situ %>% dplyr::filter(species == 'Malus fusca' & ssp %in% c('historical', '245'))
fus_ssp585_in_situ <- in_situ %>% dplyr::filter(species == 'Malus fusac' & ssp %in% c('historical', '585'))


# Plotting

# M. coronaria ssp245
# GRSin
cor_ssp245_grsin <- ggplot(data = cor_ssp245_in_situ, aes(x = period, y = GRSin, colour = suitability)) + 
  geom_point(size = 5) +
  geom_line() +
  ylim(0, 40) +
  scale_x_continuous(limits = c(2000, 2070),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070))

# SRSin  
cor_ssp245_srsin <- ggplot(data = cor_ssp245_in_situ, aes(x = period, y = SRSin, colour = suitability)) + 
  geom_point(size = 5) +
  geom_line() +
  ylim(0, 50) +
  scale_x_continuous(limits = c(2000, 2070),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070))

# ERSin 
cor_ssp245_ersin <- ggplot(data = cor_ssp245_in_situ, aes(x = period, y = ERSin, colour = suitability)) + 
  geom_point(size = 5) +
  geom_line() +
  ylim(0, 100) +
  scale_x_continuous(limits = c(2000, 2070),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070))

# FCSin
cor_ssp245_fcsin <- ggplot(data = cor_ssp245_in_situ, aes(x = period, y = FCSin, colour = suitability)) + 
  geom_point(size = 5) +
  geom_line() +
  ylim(0, 70) +
  scale_x_continuous(limits = c(2000, 2070),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070))

ggarrange(cor_ssp245_grsin, cor_ssp245_srsin, cor_ssp245_ersin, cor_ssp245_fcsin,
          nrow = 4, ncol = 1,
          legend = "bottom",
          common.legend = T)

# bar plots 
fill_cols <- c("#EDF8B1", "#7FCDBB", "#2C7FB8")

# M. coronaria
cor_srs_combined <- in_situ %>% filter(species == 'Malus coronaria' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = SRSin, fill = ssp)) + 
  geom_col(position = position_dodge()) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

cor_grs_combined <- in_situ %>% filter(species == 'Malus coronaria' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = GRSin, fill = ssp)) + 
  geom_col(position = position_dodge(15), width = 15) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),                      
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

cor_ers_combined <- in_situ %>% filter(species == 'Malus coronaria' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = ERSin, fill = ssp)) + 
  geom_col(position = position_dodge()) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),                     
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

cor_fcs_combined <- in_situ %>% filter(species == 'Malus coronaria' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = FCSin, fill = ssp)) + 
  geom_col(position = position_dodge()) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),                      
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

ggarrange(cor_srs_combined, cor_grs_combined, cor_ers_combined, cor_fcs_combined,
          nrow = 1, ncol = 4,
          legend = "top",
          common.legend = T)

# M. fusca
fus_srs_combined <- in_situ %>% filter(species == 'Malus fusca' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = SRSin, fill = ssp)) + 
  geom_col(position = position_dodge()) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

fus_grs_combined <- in_situ %>% filter(species == 'Malus fusca' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = GRSin, fill = ssp)) + 
  geom_col(position = position_dodge(15), width = 15) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),                      
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

fus_ers_combined <- in_situ %>% filter(species == 'Malus fusca' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = ERSin, fill = ssp)) + 
  geom_col(position = position_dodge()) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),                     
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

fus_fcs_combined <- in_situ %>% filter(species == 'Malus fusca' & suitability == 'high') %>% 
  ggplot(aes(x = period, y = FCSin, fill = ssp)) + 
  geom_col(position = position_dodge()) +
  scale_x_continuous(limits = c(1990, 2080),
                     breaks = c(2000, 2030, 2050, 2070),
                     labels = c("Historical", 2030, 2050, 2070)) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 20, 40, 60, 80, 100),                      
                     expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = fill_cols, 
                    breaks = c('historical', '245', '585'),
                    labels = c('Historical', 'SSP245', 'SSP585')) +
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = 'black'),
        axis.title.x=element_blank(),
        legend.title = element_blank())

ggarrange(fus_srs_combined, fus_grs_combined, fus_ers_combined, fus_fcs_combined,
          nrow = 1, ncol = 4,
          legend = "top",
          common.legend = T)

