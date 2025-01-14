## Malus coronaria in a single file


library(tidyverse) #grammar, data management
library(CoordinateCleaner) #helpful functions to clean data
library(terra) #working with vector/raster data
library(geodata) #download basemaps
library(scales) #alpha adjust colours
library(predicts)
library(ENMTools)
library(plotly) # 3D surface Kernel bivariate plots
library(MASS)
library(ENMeval) # Another modeling package, useful for data partitioning 
library(parallel) # speed up computation by running in parallel
library(doParallel) # added functionality to parallel


gadmPath <- "~/data"

#################
## occ_clean.R ##
#################

gbif_cor <- read.csv(file = "./occ_data/occ_coronaria.csv") # load coronaria data

# filter M. corornia data
occ_cor <- gbif_cor %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, botanical gardens, and other institutions 
  cc_sea() %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>%
  filter(decimalLongitude >= -100) %>% # remove some records on west coast, two from bot gardens
  filter(!(decimalLatitude < 35 & decimalLongitude < -86)) %>% # remove records from Texas, Oklahoma, Louisiana
  filter(!(decimalLatitude > 45)) %>% # remove record from New Brunswick
  filter(!(decimalLongitude < -98)) %>% # remove iNat record from Kansas, northern Kansas record verified by taxonomist
  dplyr::select(species, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
         )


## Brian Husband's verified Ontario records
husband <- read.table("occ_data/malus_coronaria_husband.csv",
                      header = TRUE, sep = "\t")

## remove records that don't have valid coordinates
husband <- husband[substr(husband$tree.location, 0, 1) == "4", ]

husband$tmp <- strsplit(husband$tree.location, "°|'|\"| ")
husband$lat <- sapply(husband$tmp, FUN = function(x) (as.numeric(x[1])
  + as.numeric(x[2])/60 + as.numeric(x[3])/3600))

husband$lon <- sapply(husband$tmp, FUN = function(x) (-1 * as.numeric(x[5])
  - as.numeric(x[6])/60 - as.numeric(x[7])/3600))

husband$tmp <- NULL

us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = gadmPath) #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = gadmPath) #Canada basemap w. 

canUS_map <- rbind(us_map, ca_map) #combine US and Canada vector map

## Swap in Brian Husband's records for Ontario:

occ_hus <- vect(occ_cor, geom = c("decimalLongitude",
                                  "decimalLatitude"))
ONT <- ca_map[ca_map$NAME_1 == "Ontario"]

occ_hus <- occ_hus[!relate(occ_hus, ONT, "intersects")]

husband <- vect(husband, geom = c("lon", "lat"))

occ_hus <- rbind(occ_hus, husband)

################
## occ_thin.R ##
################

occ_cor <- vect(occ_cor, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

wclim <- worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = gadmPath)

# M. coronaria thinning
occ_cor <- spatSample(occ_cor, size = 1, 
                      strata = wclim) #sample one occurrence from each climatic cell

occ_hus <- spatSample(occ_hus, size = 1, 
                      strata = wclim) #sample one occurrence from each 

################
## malus_bg.R ##
################

ecoNA <- vect(x = "~/data/maps/NA_CEC_Eco_Level2.shp",
              layer = 'NA_CEC_Eco_Level2') 
ecoNA <- project(ecoNA, 'WGS84') # project ecoregion vector to same coords ref as basemap

mapPath <- "~/data"

mex_map <-gadm(country = 'MX', level = 1, resolution = 2,
               path = mapPath) # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map

eco_cor <- extract(ecoNA, occ_hus) # extract what polygons contained

eco_cor_code <- eco_cor$NA_L2CODE %>% unique() 
eco_cor_code <- eco_cor_code[eco_cor_code != '0.0']  #remove the 'water'
                                        #'0.0' ecoregion 

ecoNA_cor <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_cor_code) # subset eco region spat vector by the codes

wclim_cor <- terra::crop(wclim, ecoNA_cor, mask = T)

cor_bg_vec <- spatSample(wclim_cor, 20000, 'random', na.rm = T, as.points = T) #ignore NA values

cor_predvals <- extract(wclim_cor, occ_cor) # M. coronaria
cor_predvals <- cor_predvals[-1] # drop ID column

hus_predvals <- extract(wclim_cor, occ_hus) # M. coronaria
hus_predvals <- hus_predvals[-1] # drop ID column

#################
## malus_sdm.R ##
#################

great_lakes <- vect('~/data/maps/greatlakes.shp')
great_lakes <- project(great_lakes, crs(wclim))

NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in
                                 # Western Hemisphere

bg_cor_coords <- as.data.frame(geom(cor_bg_vec)[,3:4]) # extract longitude, lattitude from background points

occ_cor_coords <- as.data.frame(geom(occ_cor)[,3:4]) # extract longitude, lattitude from occurence points

occ_hus_coords <- as.data.frame(geom(occ_hus)[,3:4]) # extract longitude, lattitude from occurence points

cn <- detectCores(logical = F) # logical = F, is number of physical RAM

date()
cor_maxent <- ENMevaluate(occ_cor_coords, # occurrence records
                          envs = wclim_cor,
                          bg = bg_cor_coords,
                          tune.args =
                            list(
                                        #rm = seq(0.5, 8, 0.5),
                              rm = seq(0.5, 4, 0.5),
                              fc = c("L", "LQ", "H",
                                     "LQH", "LQHP")),
                          ## dropping Threshold, Phillips et al. 2017
                          partition.settings =
                            list(aggregation.factor = c(9, 9),
                                 gridSampleN = 20000), # 9,9 agg 
                            partitions = 'checkerboard2',
                            parallel = TRUE,
                            numCores = 6,
                          ##cn - 1, # leave one core available for other apps
                            parallelType = "doParallel", # use doParrallel on Windows - socket cluster  
                            algorithm = 'maxent.jar')
date()


date()
hus_maxent <- ENMevaluate(occ_hus_coords, # occurrence records
                          envs = wclim_cor,
                          bg = bg_cor_coords,
                          tune.args =
                            list(
                                        #rm = seq(0.5, 8, 0.5),
                              rm = seq(0.5, 4, 0.5),
                              fc = c("L", "LQ", "H",
                                     "LQH", "LQHP")),
                          ## dropping Threshold, Phillips et al. 2017
                          partition.settings =
                            list(aggregation.factor = c(9, 9),
                                 gridSampleN = 20000), # 9,9 agg 
                            partitions = 'checkerboard2',
                            parallel = TRUE,
                            numCores = 6,
                          ##cn - 1, # leave one core available for other apps
                            parallelType = "doParallel", # use doParrallel on Windows - socket cluster  
                            algorithm = 'maxent.jar')
date()

best_cor_maxent <- subset(cor_maxent@results, delta.AICc == 0)
best_hus_maxent <- subset(hus_maxent@results, delta.AICc == 0)

mod.best_cor_maxent <- eval.models(cor_maxent)[[best_cor_maxent$tune.args]]
                                        # extracts the best
mod.best_hus_maxent <- eval.models(hus_maxent)[[best_hus_maxent$tune.args]]
                                        # extracts the best

eval.variable.importance(cor_maxent)$rm.0.5_fc.LQ
eval.variable.importance(hus_maxent)$rm.0.5_fc.LQ

wclim_na <- crop(wclim, NA_ext)

date()
cor_pred_hist <- terra::predict(wclim_na, mod.best_cor_maxent, cores = 6,
                                na.rm = T) 
date()

date()
hus_pred_hist <- terra::predict(wclim_na, mod.best_hus_maxent, cores = 4,
                                na.rm = T) 
date()
