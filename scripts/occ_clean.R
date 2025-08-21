# Top ---------------------------------------------------------------------
# cleaning occurrence data of M. coronaria, M. fusca, M. ioensis, and M. angustifolia
# Terrell Roulston
# Started Feb 16, 2024

source("./scripts/libraries.R")
source("./scripts/load_maps.R")
message("** Loading & Cleaning Occurrence Data ", date())

# load occurrence csv files ------------------------------------------------
gbif_cor <- read.csv(file = "./occ_data/cor/occ_coronaria.csv") # load coronaria data
gbif_fus <- read.csv(file = "./occ_data/fus/occ_fusca.csv") # load fusca data
gbif_ion <- read.csv(file = "./occ_data/ion/occ_ioensis.csv") # load ioensis data
gbif_ang <- read.csv(file = "./occ_data/ang/occ_angustifolia.csv") # load angustifolia data

# Clean Coronaria ---------------------------------------------------------
# filter M. corornia data
occ_cor <- gbif_cor %>% 
  filter(countryCode %in% c('US', 'CA')) %>% #limit to CA and US
  filter(!is.na(decimalLongitude)) %>% # remove records w/o coords
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>% 
  cc_cen(buffer = 2000) %>% # remove records within 2km of country centroids
  cc_inst(buffer = 2000) %>% # remove records within 2km of herbariums, botanical gardens, and other institutions 
  cc_sea(ref = seaRef) %>% 
  distinct(decimalLatitude, decimalLongitude, speciesKey, datasetKey, .keep_all = T) %>%
  filter(decimalLongitude >= -100) %>% # remove some records on west coast, two from bot gardens
  filter(!(decimalLatitude < 35 & decimalLongitude < -86)) %>% # remove records from Texas, Oklahoma, Louisiana
  filter(!(decimalLatitude > 45)) %>% # remove record from New Brunswick
  filter(!(decimalLongitude < -98)) %>% # remove iNat record from Kansas, northern Kansas record verified by taxonomist
  filter(!(catalogNumber == 174554608)) %>% # remove inaccurate record 	https://www.inaturalist.org/observations/174554608
  filter(!(catalogNumber == 181065977)) %>% # remove inaccurate record 	https://www.inaturalist.org/observations/181065977
  dplyr::select(species, gbifID, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
         )

#write.table(occ_cor, file = "./occ_data/cor/occ_cor_gbif.csv")
##saveRDS(occ_cor, file = "./occ_data/cor/occ_cor_gbif.Rdata") # Note that this copy of occurrence data is used for the sythesis paper in figure 3.

# Note it is helpful to plot the occurrences bellow, and then add more conditions to clean inaccurate points
# Pay special attention to points at the edge of the range of occurrences, as these are most likely to be suspicious
# as well as influence the model in strange ways, in comparison to inaccurate points well within the other occurrences


# Import Brian Husband's data ---------------------------------------------
## Brian's verified Ontario records
husband <- read.table("occ_data/cor/malus_coronaria_husband.csv",
                      header = TRUE, sep = "\t")

## remove records that don't have valid coordinates
husband <- husband[substr(husband$tree.location, 0, 1) == "4", ]

husband$tmp <- strsplit(husband$tree.location, "°|'|\"| ")
husband$lat <- sapply(husband$tmp, FUN = function(x) (as.numeric(x[1])
  + as.numeric(x[2])/60 + as.numeric(x[3])/3600))

husband$lon <- sapply(husband$tmp, FUN = function(x) (-1 * as.numeric(x[5])
  - as.numeric(x[6])/60 - as.numeric(x[7])/3600))

husband$tmp <- NULL

## Swap in Brian Husband's records for Ontario:

#occ_husband <- vect(occ_cor, geom = c("decimalLongitude", "decimalLatitude"))
#ont <- ca_map[ca_map$NAME_1 == "Ontario"]

#occ_husband <- occ_husband[!relate(occ_husband, ont, "intersects")]

#husband <- vect(husband, geom = c("lon", "lat"))

#occ_husband <- rbind(occ_husband, husband)

#saveRDS(occ_husband, file = "occ_data/cor/occ_husband.Rdata") # this has US reocrds from GBIF and CA recrods from Brian

# Keep all GBIF records but add Brian's data

husband_coords <- husband %>% dplyr::select(Accession, lat, lon) %>%
                              mutate(source = 'Husband') %>%  
                              mutate(species = 'Malus coronaria') %>% 
                              rename(decimalLatitude = lat, decimalLongitude = lon) #rename lat lon to match for inner join
                
occ_cor <- occ_cor %>%  mutate(source = 'GBIF') # add source for tracking and mapping purposes

occ_cor <- occ_cor %>% full_join(husband_coords, by = c("decimalLatitude", "decimalLongitude", "source", "species"))

# Save M. coronaria occurrence dataframe
write.table(occ_cor, file = "./occ_data/cor/occ_cor.csv") # Note this copy of occurrence data 
saveRDS(occ_cor, file = "./occ_data/cor/occ_cor.Rdata") # Note this copy of occurrence data will be used in downstream SDM work


# Clean Fusca -------------------------------------------------------------
# filter M. fusca data in df
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
  dplyr::select(species, gbifID, countryCode, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )

#write.table(occ_fus, file = "./occ_data/fus/occ_fus_gbif.csv") # Note that this copy of 
##saveRDS(occ_fus, file = "./occ_data/fus/occ_fus_gbif.Rdata") # Note that this copy of occurrence data is used for the sythesis paper in figure 3.

# load cleaned fusca data from CG Amrstrong
occ_armstrong <- read.csv(file = "./occ_data/fus/malus_fusca_armstrong.csv") %>% #Note that this is only Armstrong data (not a combination of GBIF data)
                  mutate(species = 'Malus fusca') %>% # add species
                  mutate(source = 'Armstrong') %>% # add source
                  dplyr::select(Latitude, Longitude, source, species) %>% 
                  rename(decimalLatitude = Latitude, decimalLongitude = Longitude)
                  
occ_fus <- occ_fus %>% mutate(source = 'GBIF')

occ_fus <- occ_fus %>% full_join(occ_armstrong, by = c("decimalLatitude", "decimalLongitude", "source", "species"))

#load data from Wickham and Obrits and Fitsp
occ_wick_orb_fit <- read.csv(file = "./occ_data/fus/malus_fusca_wickham_orbits_titzpatrick.csv") %>% 
  dplyr::select(latitude, longitude, source, species) %>% 
  rename(decimalLatitude = latitude, decimalLongitude = longitude)

occ_fus <- occ_fus %>% full_join(occ_wick_orb_fit, by = c("decimalLatitude", "decimalLongitude", "source", "species"))

write.table(occ_fus, file = "./occ_data/fus/occ_fus.csv") 
saveRDS(occ_fus, file = "./occ_data/fus/occ_fus.Rdata") # Note this copy of occurrence data will be used in downstream SDM work

# Clean Ioensis -----------------------------------------------------------
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
  dplyr::select(species, gbifID, countryCode, decimalLatitude, 
                decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )

occ_ion <- occ_ion %>% mutate(source = 'GBIF')

write.table(occ_ion, file = "./occ_data/ion/occ_ion.csv")
saveRDS(occ_ion, file = "./occ_data/ion/occ_ion.Rdata")


# Clean Angustifolia ------------------------------------------------------
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
  dplyr::select(species, gbifID, countryCode, decimalLatitude, 
                decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord
  )
occ_ang <- occ_ang %>% mutate(source = 'GBIF')

write.table(occ_ang,  file = "./occ_data/ang/occ_ang.csv")
saveRDS(occ_ang,  file = "./occ_data/ang/occ_ang.Rdata")

message("** Data Cleaned ", date())
