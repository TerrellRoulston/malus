# Top ---------------------------------------------------------------------
# Cropping Future World Clim data to ecoregions for SDMs
# Terrell Roulston
# Started Mar 22, 2024

library(tidyverse)
library(terra)
library(geodata)
library(here)


# Load ecoregion rasters --------------------------------------------------
getwd()
setwd('../occ_data/eco_regions/')
ecoNA_cor <- readRDS(file = 'ecoNA_cor.Rdata')
ecoNA_fus <- readRDS(file = 'ecoNA_fus.Rdata')


# Download/load WorldClim data under future climate scenarios -------------
# WARNING DO NOT PUSH WORLDCLIM DATA
setwd('../../wclim_data/')
# Historical climate 1970-2000
wclim <- geodata::worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "../wclim_data/")

# SSP (Shared social-economic pathway) 2.45 
# middle of the road projection, business as usual, high climate adaptation, low climate mitigation
ssp245_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/")

ssp245_2050 <- cmip6_world(model = "CanESM5",
                             ssp = "245",
                             time = "2041-2060",
                             var = "bioc",
                             res = 2.5,
                             path = "../wclim_data/")

ssp245_2070 <- cmip6_world(model = "CanESM5",
                             ssp = "245",
                             time = "2061-2080",
                             var = "bioc",
                             res = 2.5,
                             path = "../wclim_data/")

# SPP 5.85 
# low regard for enviromental sustainability, increased fossil fuel reliance, this is the current tracking projection
ssp585_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/")

ssp585_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/")

ssp585_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/")

# Crop all climate SpatRasters to ecoregions ------------------------------
# SSP 245 
cor_ssp245_2030 <- terra::crop(ssp245_2030, ecoNA_cor, mask = T) # M coronaria
fus_ssp245_2030 <- terra::crop(ssp245_2030, ecoNA_fus, mask = T) # M fusca

cor_ssp245_2050 <- terra::crop(ssp245_2050, ecoNA_cor, mask = T) # M coronaria
fus_ssp245_2050 <- terra::crop(ssp245_2050, ecoNA_fus, mask = T) # M fusca

cor_ssp245_2070 <- terra::crop(ssp245_2070, ecoNA_cor, mask = T) # M coronaria
fus_ssp245_2070 <- terra::crop(ssp245_2070, ecoNA_fus, mask = T) # M fusca

# SSP 585
cor_ssp585_2030 <- terra::crop(ssp585_2030, ecoNA_cor, mask = T) # M coronaria
fus_ssp585_2030 <- terra::crop(ssp585_2030, ecoNA_fus, mask = T) # M fusca

cor_ssp585_2050 <- terra::crop(ssp585_2050, ecoNA_cor, mask = T) # M coronaria
fus_ssp585_2050 <- terra::crop(ssp585_2050, ecoNA_fus, mask = T) # M fusca

cor_ssp585_2070 <- terra::crop(ssp585_2070, ecoNA_cor, mask = T) # M coronaria
fus_ssp585_2070 <- terra::crop(ssp585_2070, ecoNA_fus, mask = T) # M fusca


# Save cropped climate rasters for SDM predictions ------------------------
setwd('../wclim_data/')

# SSP 245
saveRDS(cor_ssp245_2030, file = 'cor_ssp245_2030.Rdata')
saveRDS(fus_ssp245_2030, file = 'fus_ssp245_2030.Rdata')

saveRDS(cor_ssp245_2050, file = 'cor_ssp245_2050.Rdata')
saveRDS(fus_ssp245_2050, file = 'fus_ssp245_2050.Rdata')

saveRDS(cor_ssp245_2070, file = 'cor_ssp245_2070.Rdata')
saveRDS(fus_ssp245_2070, file = 'fus_ssp245_2070.Rdata')

# SSP 585
saveRDS(cor_ssp585_2030, file = 'cor_ssp585_2030.Rdata')
saveRDS(fus_ssp585_2030, file = 'fus_ssp585_2030.Rdata')

saveRDS(cor_ssp585_2050, file = 'cor_ssp585_2050.Rdata')
saveRDS(fus_ssp585_2050, file = 'fus_ssp585_2050.Rdata')

saveRDS(cor_ssp585_2070, file = 'cor_ssp585_2070.Rdata')
saveRDS(fus_ssp585_2070, file = 'fus_ssp585_2070.Rdata')




