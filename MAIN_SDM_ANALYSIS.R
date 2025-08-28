# Top ---------------------------------------------------------------------
# # Code Authors: Terrell T. Roulston, Tyler W. Smith
# Last Updated: August 1st, 2025

# This script contains the complete analysis for 
#### Citation ####
# Roulston, T. T.,  Armstrong, C., Husband, B., Moreau, T., Ulrich, J., Wickham, S., Migicovsky, Z., Smith, T. W. (in prep). Modeling Reveals Conservation Shortfalls for Native North American Apple (Malus) Species.
####

### Native Malus species
# Malus fusca = "fus"
# Malus coronaria = "cor"
# Malus ioensis = "ion"
# Malus angustifolia = "ang"
# Sect. Chloromeles = "chl" (combination of cor, ion and ang)
###

# The goal of this analysis is to map the recent historical and future distribution of suitable habitats for native apple (Malus) species in North America to understand how climate change is predicted to impact these species. We build species distribution models (SDMs) using Maxent.jar algorithum via the ENMeval and predicts packages.
# We further this goal by completing an in situ gap analysis (following Carver et al. (2021) Ecography), although coded herein as the package 'GapAnalysis' is currently being updated.
# Additionally, as a case study into how taxanomic uncertainty impacts the conclusions of SDMs, we model the morphological species (M. coronaria, M. ioensis, and M. angustifolia) as a combined taxon, i.e. Sect. Chloromeles. 

# 000 Source Required Libraries -------------------------------------------
# Load required libraries
source("sdm_analysis/libraries.R")


# 001 Load and Clean Species Occurrence Data ------------------------------
message("** Loading & Cleaning Occurrence Data ", date())

# Datasets from GBIF should always be manually verified and cleaned to removed cultivated and/or spurious records that are mistakenly included.
# Note special attention was payed to outlying and edge occurrence records as these are most likely to miss-influence SDM predictions.
# Additional supplemental occurrence data was added for M. coronaria (Ontario Only) and M. fusca (entire range, but especaially central coast BC)

#### Load GBIF Occ ####
# See scripts/gbif_occ.R
gbif_cor <- read.csv(file = "./occ_data/cor/occ_coronaria.csv") # load csv data
gbif_fus <- read.csv(file = "./occ_data/fus/occ_fusca.csv") # load csv data
gbif_ion <- read.csv(file = "./occ_data/ion/occ_ioensis.csv") # load csv data
gbif_ang <- read.csv(file = "./occ_data/ang/occ_angustifolia.csv") # load csv data

#### Clean GBIF and add Suppl. Occ ####
# See scripts/occ_clean.R

##### M. fusca #####
# GBIF
occ_fus <- gbif_fus %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot cardens 
  cc_sea(ref = seaRef) %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>% 
  filter(!(decimalLongitude > -113)) %>%  # remove record from Idaho herbarium 
  filter(!(decimalLatitude > 64)) %>% # remove record from Interior Alaska
  filter(!(stateProvince == 	'Idaho')) %>% # remove phishy record from Idaho
  filter(!(gbifID == '4908101518')) %>% # questionble outliying Inaturalist observation.
  filter(!(decimalLatitude == 37.27013)) %>% # remove bad record in mountains inland California
  filter(!(decimalLatitude < 37)) %>% # remove unlikely records south of San Jose
  dplyr::select(species, countryCode, decimalLatitude, 
                decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )

# Load cleaned M. fusca data from CG Amrstrong
occ_armstrong <- read.csv(file = "./occ_data/fus/malus_fusca_armstrong.csv") %>% 
  mutate(species = 'Malus fusca') %>% # add species
  mutate(source = 'Armstrong') %>% # add source
  dplyr::select(Latitude, Longitude, source, species) %>% 
  rename(decimalLatitude = Latitude, decimalLongitude = Longitude)

# Load data from Wickham and Obrits and Fitsp
occ_wick_orb_fit <- read.csv(file = "./occ_data/fus/malus_fusca_wickham_orbits_titzpatrick.csv") %>% 
  dplyr::select(latitude, longitude, source, species) %>% 
  rename(decimalLatitude = latitude, decimalLongitude = longitude)

# Combine Datasets
occ_fus <- occ_fus %>% mutate(source = 'GBIF')
occ_fus <- occ_fus %>% full_join(occ_armstrong, by = c("decimalLatitude", "decimalLongitude", "source", "species"))
occ_fus <- occ_fus %>% full_join(occ_wick_orb_fit, by = c("decimalLatitude", "decimalLongitude", "source", "species"))

##### M. coronaria #####
# GBIF DATA
occ_cor <- gbif_cor %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, botanical gardens, and other institutions 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>%
  filter(decimalLongitude >= -100) %>% # remove some records on west coast, two from bot gardens
  filter(!(decimalLatitude < 35 & decimalLongitude < -86)) %>% # remove records from Texas, Oklahoma, Louisiana
  filter(!(decimalLatitude > 45)) %>% # remove record from New Brunswick
  filter(!(decimalLongitude < -98)) %>% # remove iNat record from Kansas, northern Kansas record verified by taxonomist
  filter(!(catalogNumber == 174554608)) %>% # remove inaccurate record 	https://www.inaturalist.org/observations/174554608
  filter(!(catalogNumber == 181065977)) %>% # remove inaccurate record 	https://www.inaturalist.org/observations/181065977
  dplyr::select(species, countryCode, decimalLatitude, 
                decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )

# Brain Husband's Supplementary M. coronaria Data
husband <- read.table("occ_data/cor/malus_coronaria_husband.csv", header = TRUE, sep = "\t")
husband <- husband[substr(husband$tree.location, 0, 1) == "4", ] ## remove records that don't have valid coordinates
husband$tmp <- strsplit(husband$tree.location, "°|'|\"| ")
husband$lat <- sapply(husband$tmp, FUN = function(x) (as.numeric(x[1]) + as.numeric(x[2])/60 + as.numeric(x[3])/3600))
husband$lon <- sapply(husband$tmp, FUN = function(x) (-1 * as.numeric(x[5]) - as.numeric(x[6])/60 - as.numeric(x[7])/3600))
husband$tmp <- NULL

# Combine Datasets
husband_coords <- husband %>% dplyr::select(Accession, lat, lon) %>%
  mutate(source = 'Husband') %>%  
  mutate(species = 'Malus coronaria') %>% 
  rename(decimalLatitude = lat, decimalLongitude = lon) #rename lat lon to match for inner join

occ_cor <- occ_cor %>%  mutate(source = 'GBIF') # add source for tracking and mapping purposes
occ_cor <- occ_cor %>% full_join(husband_coords, by = c("decimalLatitude", "decimalLongitude", "source", "species"))

##### M. ioensis #####
# GBIF DATA
occ_ion <- gbif_ion %>% 
  filter(countryCode %in% c('US')) %>% #limit to US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot cardens 
  cc_sea(ref = seaRef) %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>% 
  filter(!(decimalLongitude < -99)) %>%  # remove records from west of Kansas/Texas occ
  filter(!(gbifID %in% c('1899381516', '5067955163'))) %>% #remove questionable GBIF records
  filter(!(decimalLongitude > -84)) %>%  # remove sporadic east coast records. likely not ioensis based on classical distribtion
  dplyr::select(species, countryCode, decimalLatitude, 
                decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )
occ_ion <- occ_ion %>% mutate(source = 'GBIF')

##### M. angustifolia #####
# GBIF DATA
occ_ang <- gbif_ang %>% 
  filter(countryCode %in% c('US')) %>% #limit to US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, bot cardens 
  cc_sea(ref = seaRef) %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>% 
  filter(decimalLatitude <= 41.99575 | decimalLongitude <= -72.60754) %>% #filter occ northeast of Connecticut
  filter(!(gbifID %in% c('3865090173', '4919535276', '1302641832', '3091187987', '4855490969', '4854957347'))) %>%  #remove a record from greenland (???) and bad inat obs
  dplyr::select(species, countryCode, decimalLatitude, 
                decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )
occ_ang <- occ_ang %>% mutate(source = 'GBIF')

message("** Data Cleaned ", date())


# 002 Thin Occurrence Data ------------------------------------------------
# See scripts/occ_thin.R

# Need to load some map data to use as strata for thinning occurrence records
source("scripts/load_maps.R")
# This script also sources layers for plotting including country boundaries
# And ecoregion data used below in background data

message("** Thinning Records: ", date())

# Combine the cleanned occurrence records for all three species in Chloromeles
# M. coronaria, M. ioensis, M. angustifolia
# This combined set will be thinned seperately to avoid overlaping species occurrences

occ_chl <- occ_cor %>% 
  full_join(occ_ion, by = c('species', 'source', 'decimalLongitude', 'decimalLatitude')) %>% 
  full_join(occ_ang, by = c('species', 'source', 'decimalLongitude', 'decimalLatitude')) %>% 
  dplyr::select('species', 'source', 'decimalLongitude', 'decimalLatitude')

# Vectorize occurrence df to coordinates for below
# This stores the other coloumns as objects in layers
# Make sure projection matches the source (GBIF = WGS84)

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

# Set a random generator seed to get reproducible results
set.seed(1337) 


#### M. fusca
occThin_fus <- spatSample(occ_fus_vect, size = 1,
                          strata = wclim_subs, # Sample one occurrence from each climatic cell
                          method = "random")

#### M. coronaria
occThin_cor <- spatSample(occ_cor_vect, size = 1, 
                          strata = wclim_subs,  
                          method = "random") 

#### M. ioensis
occThin_ion <- spatSample(occ_ion_vect, size = 1,
                          strata = wclim_subs,
                          method = "random")

#### M. angustifolia
occThin_ang <- spatSample(occ_ang_vect, size = 1,
                          strata = wclim_subs,
                          method = "random")

#### Sect. Chloromeles 
occThin_chl <- spatSample(occ_chl_vect, size = 1,
                          strata = wclim_subs,
                          method = "random")

# Finished thinning
message("** Records Thinned: ", date())


# 003 Background Extent and Data ------------------------------------------
# See malus_bg.R script
# The background extent is the area which is a defined study area, where the SDM will randomly sample background points
# We define background using Ecoregions which contain occurrence data for each species/taxon

message("** Calculating Background Data: ", date())

#### M. fusca ####
eco_fus <- intersect(occThin_fus, ecoNA) %>% # extract ecoregion names from points
  as.data.frame() %>% # convert to df
  group_by(NA_L2CODE) %>% # sample 1 point of each eco region
  sample_n(1)

# return vector of eco region codes of the polygons that contain occurrences
eco_fus_code <- eco_fus$NA_L2CODE %>% unique() 
eco_fus_code <- eco_fus_code[eco_fus_code != '0.0'] # remove NA value

#CODES "7.1""6.2"  "10.1" "11.1" "10.2"

ecoNA_fus <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_fus_code)

#### M. coronaria ####
eco_cor <- intersect(occThin_cor, ecoNA) %>% 
  as.data.frame() %>% 
  group_by(NA_L2CODE) %>% 
  sample_n(1)

eco_cor_code <- eco_cor$NA_L2CODE %>% unique() 
eco_cor_code <- eco_cor_code[eco_cor_code != '0.0'] 

#CODES: "8.1" "8.2" "5.3" "8.4" "8.3" "8.5" "9.2" "9.4" "5.2"

ecoNA_cor <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_cor_code) 

#### M. ioensis ####
eco_ion <- intersect(occThin_ion, ecoNA) %>% 
as.data.frame() %>% 
  group_by(NA_L2CODE) %>% 
  sample_n(1)

eco_ion_code <- eco_ion$NA_L2CODE %>% unique() 
eco_ion_code <- eco_ion_code[eco_ion_code != '0.0'] 

#CODES "5.2" "8.1" "8.2" "8.3" "8.4" "8.5" "9.2" "9.4"

ecoNA_ion <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_ion_code)

#### M. angustifolia ####
eco_ang <- intersect(occThin_ang, ecoNA) %>% 
  as.data.frame() %>% 
  group_by(NA_L2CODE) %>% 
  sample_n(1)

eco_ang_code <- eco_ang$NA_L2CODE %>% unique() 
eco_ang_code <- eco_ang_code[eco_ang_code != '0.0']

#CODES "5.3" "8.1" "8.2" "8.3" "8.4" "8.5" "9.5"

ecoNA_ang <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_ang_code)

#### Sect. Chloromeles ####
eco_chl <- intersect(occThin_chl, ecoNA) %>% 
  as.data.frame() %>% 
  group_by(NA_L2CODE) %>% 
  sample_n(1)

eco_chl_code <- eco_chl$NA_L2CODE %>% unique() 
eco_chl_code <- eco_chl_code[eco_chl_code != '0.0'] 

#CODES "5.3" "8.1" "8.2" "8.3" "8.4" "8.5" "9.5"

ecoNA_chl <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_chl_code)


#### Crop Rasters to Ecoregions ####
# Crop extent of WorldClim data to the Malus ecoregions

wclim_cor_subs <- terra::crop(wclim_subs, ecoNA_cor, mask = T)
wclim_fus_subs <- terra::crop(wclim_subs, ecoNA_fus, mask = T)
wclim_ion_subs <- terra::crop(wclim_subs, ecoNA_ion, mask = T)
wclim_ang_subs <- terra::crop(wclim_subs, ecoNA_ang, mask = T)
wclim_chl_subs <- terra::crop(wclim_subs, ecoNA_chl, mask = T)



#### Sample for Boyce Index ####
# Need to sample background points for the Boyce Index calculation

# M. fusca
fus_bg_vec <- spatSample(wclim_fus_subs, 20000, 'random', na.rm = T, as.points = T)

# M. coronaria
cor_bg_vec <- spatSample(wclim_cor_subs, 20000, 'random', na.rm = T, as.points = T)

# M. ioensis
ion_bg_vec <- spatSample(wclim_ion_subs, 20000, 'random', na.rm = T, as.points = T)

# M. angustifolia
ang_bg_vec <- spatSample(wclim_ang_subs, 20000, 'random', na.rm = T, as.points = T)

# Sect. Chloromeles
chl_bg_vec <- spatSample(wclim_chl_subs, 20000, 'random', na.rm = T, as.points = T)

message("** Background Data Calculated: ", date())

# Spatial partitioning preparation. Needed for SDMs and Boyce Index
# M. fusca
occ_fus_coords <- as.data.frame(geom(occThin_fus)[,3:4]) # extract longitude, lattitude from occurence points
bg_fus_coords <- as.data.frame(geom(fus_bg_vec)[,3:4]) # extract longitude, lattitude from background points

# M. coronaria
occ_cor_coords <- as.data.frame(geom(occThin_cor)[,3:4]) 
bg_cor_coords <- as.data.frame(geom(cor_bg_vec)[,3:4]) 

# M. ioensis
occ_ion_coords <- as.data.frame(geom(occThin_ion)[,3:4]) 
bg_ion_coords <- as.data.frame(geom(ion_bg_vec)[,3:4]) 

# M. angustifolia
occ_ang_coords <- as.data.frame(geom(occThin_ang)[,3:4]) 
bg_ang_coords <- as.data.frame(geom(ang_bg_vec)[,3:4]) 

# Sect. Chloromeles
occ_chl_coords <- as.data.frame(geom(occThin_chl)[,3:4]) 
bg_chl_coords <- as.data.frame(geom(chl_bg_vec)[,3:4]) 

# 004 Species Distribution Modeling ---------------------------------------
# See malus_sdm.R for more details. Herein is streamlined and skips data exploration and saving steps included in the base script.
# Relies on sourcing load_maps.R ABOVE to load the future climate data that has been cropped and subsetted for the selected environmental predictors (6 in total)

# Note that this is a computationally heavy step and could take even a few hours to run

message("** Starting SDMs: ", date())

# Run prediction in a parallel using 'socket' clusters to help speed up computation
# <ENMeval> implements parallel functions natively
# But load <parallel> library for additional functions like <decectCores()>
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337) # Double check seed is set to get reproducible results

#### M. fusca ####
# Takes about 10-15 mins to run each SDM with 32 gb ram
message("** Running M. fusca SDM: ", date())
fus_maxent <- ENMevaluate(occ_fus_coords, # occurrence records
                          envs = wclim_fus_subs, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5), # Regularization 0.5-4
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")), # Feature Class: L, H, Q and P
                          partition.settings =    
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                          partitions = 'checkerboard',
                          parallel = TRUE, # run in parellel
                          numCores = cn - 1, # leave one core available for other apps
                          algorithm = 'maxent.jar') # maxent jar application requires Java Development Kit

message("** Finished M. fusca SDM: ", date())

# Select the best performing model based on delta AICc - returns data frame object
best_fus_maxent <- subset(fus_maxent@results, delta.AICc == 0) 
# Extracts the best model - returns MaxEnt object
mod.best_fus_maxent <- eval.models(fus_maxent)[[best_fus_maxent$tune.args]]
# Best = rm.1_fc.LQHPT

eval.variable.importance(fus_maxent)[best_fus_maxent["tune.args"][1, 1]]
head(fus_maxent@results[order(fus_maxent@results$delta.AICc),
                        c("rm", "fc", "delta.AICc")]) 

message("** Running M. fusca Habitat Predictions: ", date())
# Historical
fus_pred_hist <- terra::predict(wclim_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
# SSP245
fus_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
# SSP585
fus_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
message("** Finished M. fusca Habitat Predictions: ", date())

#### M. coronaria ####
message("** Running M. coronaria SDM: ", date())

cor_maxent <- ENMevaluate(occ_cor_coords, 
                          envs = wclim_cor_subs, 
                          n.bg = 20000, 
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5), 
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), 
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, 
                          algorithm = 'maxent.jar')

message("** Finished M. coronaria SDM: ", date())

best_cor_maxent <- subset(cor_maxent@results, delta.AICc == 0) 
mod.best_cor_maxent <- eval.models(cor_maxent)[[best_cor_maxent$tune.args]]
eval.variable.importance(cor_maxent)[best_cor_maxent["tune.args"][1, 1]]
head(cor_maxent@results[order(cor_maxent@results$delta.AICc),
                   c("rm", "fc", "delta.AICc")]) 

message("** Running M. coronaria Habitat Predictions: ", date())
# Historical
cor_pred_hist <- terra::predict(wclim_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
# SSP245
cor_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
# SSP585
cor_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)

message("** Finished M. coronaria Habitat Predictions: ", date())

#### M. ioensis ####
message("** Running M. ioensis SDM: ", date())

ion_maxent <- ENMevaluate(occ_ion_coords, 
                          envs = wclim_ion_subs, 
                          n.bg = 20000, 
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5), 
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), 
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, 
                          algorithm = 'maxent.jar')

message("** Finished M. ioensis SDM: ", date())

best_ion_maxent <- subset(ion_maxent@results, delta.AICc == 0) 
mod.best_ion_maxent <- eval.models(ion_maxent)[[best_ion_maxent$tune.args]]
eval.variable.importance(ion_maxent)[best_ion_maxent["tune.args"][1, 1]]
head(ion_maxent@results[order(ion_maxent@results$delta.AICc),
                        c("rm", "fc", "delta.AICc")]) 

message("** Running M. ioensis Habitat Predictions: ", date())
# Historical
ion_pred_hist <- terra::predict(wclim_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
# SSP245
ion_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
# SSP585
ion_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)

message("** Finished M. ioensis Habitat Predictions: ", date())

#### M. angustifolia ####
message("** Running M. angustifolia SDM: ", date())

ang_maxent <- ENMevaluate(occ_ang_coords, 
                          envs = wclim_ang_subs, 
                          n.bg = 20000, 
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5), 
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), 
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, 
                          algorithm = 'maxent.jar')

message("** Finished M. angustifolia SDM: ", date())

best_ang_maxent <- subset(ang_maxent@results, delta.AICc == 0) 
mod.best_ang_maxent <- eval.models(ang_maxent)[[best_ang_maxent$tune.args]]
eval.variable.importance(ang_maxent)[best_ang_maxent["tune.args"][1, 1]]
head(ang_maxent@results[order(ang_maxent@results$delta.AICc),
                        c("rm", "fc", "delta.AICc")]) 

message("** Running M. angustifolia Habitat Predictions: ", date())
# Historical
ang_pred_hist <- terra::predict(wclim_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
# SSP245
ang_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
# SSP585
ang_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)

message("** Finished M. angustifolia Habitat Predictions: ", date())

#### Sect. Chloromeles ####
message("** Running Sect. Chloromeles SDM: ", date())

chl_maxent <- ENMevaluate(occ_chl_coords, 
                          envs = wclim_chl_subs, 
                          n.bg = 20000, 
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5), 
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), 
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, 
                          algorithm = 'maxent.jar')

message("** Finished Sect. Chloromeles SDM: ", date())

best_chl_maxent <- subset(chl_maxent@results, delta.AICc == 0) 
mod.best_chl_maxent <- eval.models(chl_maxent)[[best_chl_maxent$tune.args]]
eval.variable.importance(chl_maxent)[best_chl_maxent["tune.args"][1, 1]]
head(chl_maxent@results[order(chl_maxent@results$delta.AICc),
                        c("rm", "fc", "delta.AICc")]) 

message("** Running Sect. Chloromeles Habitat Predictions: ", date())
# Historical
chl_pred_hist <- terra::predict(wclim_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
# SSP245
chl_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
# SSP585
chl_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)

message("** Finished Sect. Chloromeles Habitat Predictions: ", date())

# See sdm_plot for code for producing the maps that are seen in the paper
# Note that these maps were processed and arranged using the GIMP software outside of R

message("** Finished Running SDMs: ", date())

# 005 Realized Niche PCA --------------------------------------------------
message("** Realized Niche PCA: ", date())
# See malus_pca.R for plotting
# This analysis builds off function created by authors of the ecospat and ade4 packages
# I source these custom functions to save space herein
source("./scripts/pca_functions.R")

# Convert BG Rasters to Matrices for PCA
bg_mat_full <- values(wclim_chl_subs) %>% na.omit() # remove NA Values
bg_mat_cor <- values(wclim_cor_subs) %>% na.omit()
bg_mat_ion <- values(wclim_ion_subs) %>% na.omit()
bg_mat_ang <- values(wclim_ang_subs) %>% na.omit()

# Extract the wclim raster values from occurrence points then bind them with the SpatVector points
wclim_cor_occ <- cbind(crds(occThin_cor), extract(wclim_cor_subs, occThin_cor, ID = FALSE))
wclim_ion_occ <- cbind(crds(occThin_ion), extract(wclim_ion_subs, occThin_ion, ID = FALSE))
wclim_ang_occ <- cbind(crds(occThin_ang), extract(wclim_ang_subs, occThin_ang, ID = FALSE))

wclim_cor_occ <- wclim_cor_occ[complete.cases(wclim_cor_occ), ]
wclim_ion_occ <- wclim_ion_occ[complete.cases(wclim_ion_occ), ]
wclim_ang_occ <- wclim_ang_occ[complete.cases(wclim_ang_occ), ]

#### PCA ####
# Now create a PCA using the FULL BG EXTENT matrix
# We want to generate components of the entire enviromental variability
# Make sure to center (subtract from mean) and scale (divide by SD to make it 0-1)
# Set scannf to false to skip selecting the number of axes and set it equal to 2 with nf

pca_full <- dudi.pca(bg_mat_full, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
pca_score <- pca_full$li # extract PCA scores

# Extract Eigenvalues and compute statistics
eig <- pca_full$eig
prop_var <- eig / sum(eig) * 100    # Proportion of variance in percentage
cum_var <- cumsum(prop_var)         # Cumulative variance

# Build a summary table for principal components
pca_summary <- data.frame(
  PC = paste0("PC", 1:length(eig)),
  Eigenvalue = round(eig, 2),
  Proportion = round(prop_var, 2),
  Cumulative = round(cum_var, 2)
)

# Display the PCA summary table
print(pca_summary)

# Extract the variable loadings for the retained PCs 
loadings <- pca_full$c1
print(loadings)

#### Niche Equivalency ####
# Niche Equivalency Test that permuates scores to generate a pseudo p-value
####
# NOTE: this is computationally exspensive and will take ~1.1hr to run each test!!!
####
# I ran it remotely on a HPC (see "malus_gap_v2.R" and "scripts/remote_jobs/run_malus_pca_equiv.R"). I comment them out here and load them to save time. Note this is only needed to generate psuedo p-values, so just use the ecospat ecospat.niche.overlap function to do the same test without p values.

# M. coronaria vs. M. ioensis
# cor_ion_test <- niche_equivalency_test(
#                         sp1_scores = cor_occ_score,
#                         sp2_scores = ion_occ_score,
#                         sp1_bg_scores = cor_bg_score,
#                         sp2_bg_scores = ion_bg_score,
#                         bg_scores = pca_score,
#                         reps = 999,
#                         R = 100,
#                         parallel = TRUE,
#                         ncores = 15,
#                         verbose = TRUE
#                       )

# M. coronaria vs. M. angustifolia
# cor_ang_test <- niche_equivalency_test(
#                         sp1_scores = cor_occ_score,
#                         sp2_scores = ang_occ_score,
#                         sp1_bg_scores = cor_bg_score,
#                         sp2_bg_scores = ang_bg_score,
#                         bg_scores = pca_score,
#                         reps = 999,
#                         R = 100,
#                         parallel = TRUE,
#                         ncores = 15,
#                         verbose = TRUE
#                       )

# M. ioensis vs. angustifolia
# ion_ang_test <- niche_equivalency_test(
#                         sp1_scores = ion_occ_score,
#                         sp2_scores = ang_occ_score,
#                         sp1_bg_scores = ion_bg_score,
#                         sp2_bg_scores = ang_bg_score,
#                         bg_scores = pca_score,
#                         reps = 999,
#                         R = 100,
#                         parallel = TRUE,
#                         ncores = 15,
#                         verbose = TRUE
#                       )

# Load Niche Equivalency results from HPC
files <- list.files("./", pattern = "niche_equiv_.*\\.rds", full.names = TRUE)

results <- map_dfr(files, readRDS, .id = "pair")
results$pair <- gsub("niche_equiv_|\\.rds", "", basename(results$pair))

print(results)

#### Niche Overlap Plot ####
# NOTE: I am only showing one pair of species as an example
par(mar = c(7, 8, 4, 2), mgp = c(5, 2, 0))

ecospat.plot.niche.pair(
  z1 = grids[["ion"]], 
  z2 = grids[["ang"]],
  use.zcor = F,        # Toggle between using z.cor (TRUE) and z.uncor (FALSE). Z.cor = z/Z = niche preference, z.uncor = z = realized niche
  quant = 0.1,         # occupancy threshold
  quant_occ = 0.5,     # niche (z.cor) outline threshold
  quant_bg = 0.01,     # background (Z) extent threshold
  drawKD = F,          # Draw kernel density shading
  drawOccOutline = F,
  drawBGextent = T,
  col.unique_i = "#E88E00",  # Color for niche unique to species i (z1)
  col.unique_j = "#007CBE",# Color for niche unique to species j (z2)
  col.bg_i = "#E88E00",      # Color for background outline from z1's Z component
  col.bg_j = "#007CBE",    # Color for background outline from z2's Z component
  name.axis1 = "PC 1",    # X-axis label
  name.axis2 = "PC 2",    # Y-axis label
  title = NULL,  # Plot title
  show.legendOcc = F,
  max.alpha = 0.6,
  cex.lab = 3.5,
  cex.axis = 3.5,
)


####
# Compute overlap metrics for displaying
# Note excluded D metric, only reporting Warren's I in the MS
### 
metrics <- ecospat.niche.overlap(grids[["ion"]], grids[["ang"]], cor = F)
legend("topright", legend = paste0(
  #"D = ", round(metrics$D, 2), "***",
  "I = ", round(metrics$I, 2), "***"),
  bty = "n", cex = 3.5, inset = -0.01)

#### PCA Bi-plot ####
# Compute global axis limits using "Axis1" and "Axis2" across all species
all_scores <- do.call(rbind, occ_scores_list)
x_min <- min(all_scores$Axis1, na.rm = TRUE)
x_max <- max(all_scores$Axis1, na.rm = TRUE)
y_min <- min(all_scores$Axis2, na.rm = TRUE)
y_max <- max(all_scores$Axis2, na.rm = TRUE)

# Axis labels
xlab_text <- paste0("PC 1 (", signif(pca_summary$Proportion[1], 3), "%)")
ylab_text <- paste0("PC 2 (", signif(pca_summary$Proportion[2], 3), "%)")

# Scale factor for loadings (optional tweak to improve visibility)
scale_factor <- 1.75
par(mar = c(7, 8, 6, 2), mgp = c(5, 2, 0))

# Plot occurrence points using Axis1 and Axis2
# Plot M. coronaria
plot(cor_occ_score$Axis1, cor_occ_score$Axis2,
     xlim = c(x_min, x_max), ylim = c(y_min, y_max),
     xlab = xlab_text, ylab = ylab_text,
     main = NULL,
     pch = 19, col = adjustcolor("magenta", alpha.f = 0.6),
     cex.axis = 3.5, cex.lab = 3.5, cex.main = 5, cex = 2.5)


# Plot M. ioensis
points(ion_occ_score$Axis1, ion_occ_score$Axis2,
       pch = 19, col = adjustcolor("#E88E00", alpha.f = 0.6), cex = 2.5)

# Plot M. angustifolia
points(ang_occ_score$Axis1, ang_occ_score$Axis2,
       pch = 19, col = adjustcolor("#007CBE", alpha.f = 0.6), cex = 2.5)

# Draw reference lines at 0
abline(h = 0, v = 0, col = "gray60", lty = 2, lwd = 2)

# Overlay PCA loadings (from your global PCA stored in pca_full$c1)
# Multiply by scale_factor to enhance visibility if needed
arrows(0, 0,
       pca_full$c1[,1] * scale_factor,
       pca_full$c1[,2] * scale_factor,
       length = 0.1, col = "red", lwd = 5)

# Create simplified labels from the original rownames
simplified_labels <- sub("wc2.1_2.5m_bio_", "Bio ", rownames(pca_full$c1))

# Add text labels for each variable loading
text(pca_full$c1[,1] * scale_factor,
     pca_full$c1[,2] * scale_factor,
     labels = simplified_labels,
     col = "black", pos = 3, cex = 3.5, font = 2)

message("** Finished PCA Analysis: ", date())

# 006 Gap Analysis --------------------------------------------------------
# The following gap analysis is referencing Carver et al. (2021)
# https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.05430
# We complete only the in situ analysis

# NOTE: The thresholded objects are from the SDM script (see malus_sdm.R)
# Also see malus_gap_v2. R and/or pca_functions for the source of the metric calculations

# Species list
species_list <- c("cor", "fus", "ion", "ang", "chl")
species_names <- c("Malus coronaria", "Malus fusca", "Malus ioensis", "Malus angustifolia", "Sect. Chloromeles")


#### Run Gap Analysis ####

full_result <- list() #initialize list

for (index in seq_along(species_list)) {
  sp_code <- species_list[index]
  sp_name <- species_names[index]
  
  cat("\n==== Running gap analysis for:", sp_code, "====\n")
  Start <- Sys.time()
  
  base_path <- file.path("./sdm_output", sp_code, "subs")
  thresh_path <- file.path(base_path, "threshold")
  occ_path <- file.path("./occ_data", sp_code)
  
  pa_raster <- terra::rast("./gap_analysis/pa_raster_us_can.tif")
  eco_vec <- readRDS(file.path("./maps/eco_regions", paste0("ecoNA_", sp_code, ".Rdata")))
  
  occ <- readRDS(file.path(occ_path, paste0("occThin_", sp_code, ".Rdata")))
  if (!inherits(occ, "SpatVector")) occ <- terra::vect(occ)
  
  eco_mask_pts <- terra::intersect(eco_vec, occ)
  eco_occ <- unique(eco_mask_pts$NA_L2CODE)
  eco_vec_crop <- eco_vec[eco_vec$NA_L2CODE %in% eco_occ, ]
  
  preds <- list(
    hist = readRDS(file.path(base_path, paste0(sp_code, "_pred_hist_subs.Rdata"))),
    ssp245_30 = readRDS(file.path(base_path, paste0(sp_code, "_pred_ssp245_30_subs.Rdata"))),
    ssp245_50 = readRDS(file.path(base_path, paste0(sp_code, "_pred_ssp245_50_subs.Rdata"))),
    ssp245_70 = readRDS(file.path(base_path, paste0(sp_code, "_pred_ssp245_70_subs.Rdata"))),
    ssp585_30 = readRDS(file.path(base_path, paste0(sp_code, "_pred_ssp585_30_subs.Rdata"))),
    ssp585_50 = readRDS(file.path(base_path, paste0(sp_code, "_pred_ssp585_50_subs.Rdata"))),
    ssp585_70 = readRDS(file.path(base_path, paste0(sp_code, "_pred_ssp585_70_subs.Rdata")))
  )
  
  thresholds <- list(
    low = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_1_subs.Rdata"))),
    mod = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_10_subs.Rdata"))),
    high = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_50_subs.Rdata")))
  )
  
  preds <- lapply(preds, function(x) terra::crop(x, eco_vec_crop, mask = TRUE))
  template <- preds[[1]]
  pa_mask_resampled <- terra::resample(pa_raster == 1, template, method = "near")
  hist_masked <- preds[["hist"]] > thresholds[["high"]]
  hist_masked <- terra::resample(hist_masked, template, method = "near")
  
  eco_rast <- terra::rasterize(eco_vec_crop, template, field = "NA_L2CODE")
  
  out <- list()
  for (pname in names(preds)) {
    for (thresh_name in "high") { # Only analysis high suitability areas
      th <- thresholds[[thresh_name]]
      pr <- preds[[pname]]
      
      cat("[INFO] Processing:", sp_code, pname, thresh_name, "\n")
      
      SRSin <- calculate_srsin(pr, pa_mask_resampled, occ, th)
      GRSin <- calculate_grsin(pr, th, hist_masked, pa_mask_resampled)
      ERSin <- calculate_ersin(pr, th, eco_rast, hist_masked, pa_mask_resampled)
      FCSin <- mean(c(SRSin, GRSin, ERSin))
      
      out[[paste(pname, thresh_name, sep = "_")]] <- tibble(
        species = sp_name,
        sp_code = sp_code,
        ssp = ifelse(pname == 'hist', 'historical', ifelse(grepl('245', pname), '245', '585')),
        period = case_when(
          pname == 'hist' ~ 2000,
          pname == 'ssp245_30' ~ 2030,
          pname == 'ssp245_50' ~ 2050,
          pname == 'ssp245_70' ~ 2070,
          pname == 'ssp585_30' ~ 2030,
          pname == 'ssp585_50' ~ 2050,
          pname == 'ssp585_70' ~ 2070
        ),
        suitability = thresh_name,
        SRSin = SRSin,
        GRSin = GRSin,
        ERSin = ERSin,
        FCSin = FCSin
      )
    }
  }
  
  full_result[[sp_code]] <- bind_rows(out)
}

all_species_results <- bind_rows(full_result)
print(all_species_results)

message("** Finished running gap analysis: ", date())

# END ---------------------------------------------------------------------
# See other exploratory and data cleaning/prep scripts in the "scripts" folder