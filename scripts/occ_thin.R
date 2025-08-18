# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 20, 2024

## make sure the data is in sync:
source("scripts/occ_clean.R") ## this loads the maps as well

message("** Thinning Records: ", date())

#########################################################################
## The code below will load the saved data. This may be problematic if ##
## the saved data is no longer in sync with the code, although it will ##
## save a 20 seconds over re-cleaning the data via the `source` line   ##
## above.                                                              ##
#########################################################################

## library(tidyverse) # grammar and data management 
## library(terra) # working with spatial data
## library(geodata) # basemaps and climate data


## # Load cleaned occurrence data ---------------------------------------------

## ## occ_cor <- readRDS(file = "./occ_data/cor/occ_cor.Rdata") # GBIF + Husband
## ## occ_fus <- readRDS(file = './occ_data/fus/occ_fus.Rdata') # GBIF + Armstrong + Wickham + Obr. + Fit
## ## occ_ion <- readRDS(file = './occ_data/ion/occ_ion.Rdata') # GBIF
## ## occ_ang <- readRDS(file = './occ_data/ang/occ_ang.Rdata') # GBIF

## occ_cor_orig <- readRDS(file = "./occ_data/cor/occ_cor.Rdata") # GBIF + Husband
## occ_fus_orig <- readRDS(file = './occ_data/fus/occ_fus.Rdata') # GBIF + Armstrong + Wickham + Obr. + Fit
## occ_ion_orig <- readRDS(file = './occ_data/ion/occ_ion.Rdata') # GBIF
## occ_ang_orig <- readRDS(file = './occ_data/ang/occ_ang.Rdata') # GBIF

## occ_cor <- read.table(file = "./occ_data/cor/occ_cor.csv") # GBIF + Husband
## occ_fus <- read.table(file = './occ_data/fus/occ_fus.csv') # GBIF + Armstrong + Wickham + Obr. + Fit
## occ_ion <- read.table(file = './occ_data/ion/occ_ion.csv') # GBIF
## occ_ang <- read.table(file = './occ_data/ang/occ_ang.csv') # GBIF

## ## M. coronaria: note one coastal record from New York is excluded in the
## ## new version:

## dim(occ_cor)
## dim(occ_cor_orig) 

## ## M. fusca: 226 fewer records in the latest version. Looks like duplicate
## ## records from the Wickham, Obrits and Fitsp datasets:

## dim(occ_fus)
## dim(occ_fus_orig)

## ## M. ioensis and M. angustifolia unchanged between previous and current
## ## version: 
## dim(occ_ion)
## dim(occ_ion_orig)

## dim(occ_ang)
## dim(occ_ang_orig)

# Combine occ data for the 3 Chloromeles species
occ_chl <- occ_cor %>% 
  full_join(occ_ion, by = c('species', 'source', 'decimalLongitude', 'decimalLatitude')) %>% 
  full_join(occ_ang, by = c('species', 'source', 'decimalLongitude', 'decimalLatitude')) %>% 
  dplyr::select('species', 'source', 'decimalLongitude', 'decimalLatitude')

# vectorize occurrence df to coordinates for below
occ_cor_vect <- vect(occ_cor, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_fus_vect <- vect(occ_fus, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_ion_vect <- vect(occ_ion, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_ang_vect <- vect(occ_ang, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_chl_vect <- vect(occ_chl, geom = c('decimalLongitude', 'decimalLatitude'),
                           crs = "+proj=longlat +datum=WGS84")

# Thin data using sampler -------------------------------------------------

set.seed(1337) # set random generator seed to get reproducible results
# M. coronaria thinning
occThin_cor <- spatSample(occ_cor_vect, size = 1, 
                      strata = wclim_subs,  #sample one occurrence from each climatic cell
                      method = "random") 

# M. fusca thinning
occThin_fus <- spatSample(occ_fus_vect, size = 1,
                      strata = wclim_subs,
                      method = "random")

#M. ionesis thinning
occThin_ion <- spatSample(occ_ion_vect, size = 1,
                          strata = wclim_subs,
                          method = "random")

#M. angustifolia thinning
occThin_ang <- spatSample(occ_ang_vect, size = 1,
                          strata = wclim_subs,
                          method = "random")
#Sect. Chloromeles thinning
occThin_chl <- spatSample(occ_chl_vect, size = 1,
                          strata = wclim_subs,
                          method = "random")
# Save thinned occurrence points for further analysis ----------------------

##saveRDS(occThin_cor, file = './occ_data/cor/occThin_cor.Rdata')
##saveRDS(occThin_fus, file = './occ_data/fus/occThin_fus.Rdata')
##saveRDS(occThin_ion, file = './occ_data/ion/occThin_ion.Rdata')
##saveRDS(occThin_ang, file = './occ_data/ang/occThin_ang.Rdata')
##saveRDS(occThin_chl, file = './occ_data/chl/occThin_chl.Rdata')

## write.table(occThin_cor, file = './occ_data/cor/occThin_cor.csv')
## write.table(occThin_fus, file = './occ_data/fus/occThin_fus.csv')
## write.table(occThin_ion, file = './occ_data/ion/occThin_ion.csv')
## write.table(occThin_ang, file = './occ_data/ang/occThin_ang.csv')
## write.table(occThin_chl, file = './occ_data/chl/occThin_chl.csv')

# Extract environmental data from wclim ------------------------------------
# Extract and add wclim raster data to spatial occurrence data

#occ_cor <- cbind(occ_cor, extract(wclim, occ_cor)) # M. coronaria
#occ_fus <- cbind(occ_fus, extract(wclim, occ_fus)) # M. fusca

## moved actual plots to occ_plots.R
message("** Records Thinned: ", date())
