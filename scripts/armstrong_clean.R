# Top ---------------------------------------------------------------------
# cleaning Malus fusca records from Chelsey Geralda Armstrong
# Terrell Roulston
# Started Feb 26, 2025

#Issues in the coordinates for records converted from UTM measurments incorreclty need to be adjusted
#this data is corrected in the occurrence data folder (malus_fusca_armstrong.csv)

library(tidyverse)
library(geodata)
library(terra)
library(oce)


# Load basemaps -----------------------------------------------------------

us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "./maps/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = './maps/base_maps') #Canada basemap w. Provinces

mex_map <-gadm(country = 'MX', level = 1, resolution = 2,
               path = './maps/base_maps') # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map



# Load Malus fusca data from CG Armstrong ---------------------------------
armstrong <- readxl::read_xlsx("C:/Users/terre/Documents/Acadia/Malus Project/ex situ data/Core MF geodata.xlsx") 
armstrong <- armstrong  %>% filter(Easting != 'NA') # drop one record that does have any coordinates

# There apprears to be an issues with some of the coordiates, I suspect that the UTMs were converted incorrectly.
plot(canUSMex_map, xlim = c(-180, -100), ylim = c(30, 70))
points(x=armstrong$Longitude, y=armstrong$Latitude, col = 'red')

# Seperate UTM into easting and northing values
armstrong <- armstrong %>% separate(UTM, into = c("Easting", "Northing"), sep = ", ", convert = TRUE, remove = FALSE)
armstrong$Easting <- as.numeric(armstrong$Easting)
armstrong$Northing <- as.numeric(armstrong$Northing)
# filter and convert occ form UTM.ZONE 9

armstrong_utm_9 <- armstrong %>% filter(`UTM ZONE` == 9)
armstrong_latLon_9 <- as.data.frame(utm2lonlat(easting = armstrong_utm_9$Easting, northing = armstrong_utm_9$Northing, zone = 9, hemisphere = 'N', km = F)) %>% rename(newLongitude = longitude, newLatitude = latitude)
#Join the new lat/lon values with the utm df to rejoin with main df
armstrong_utm_9 <- cbind(armstrong_utm_9, armstrong_latLon_9)

# filter and convert occ form UTM.ZONE 10
armstrong_utm_10 <- armstrong %>% filter(`UTM ZONE` == 10)
armstrong_latLon_10 <- as.data.frame(utm2lonlat(easting = armstrong_utm_10$Easting, northing = armstrong_utm_10$Northing, zone = 10, hemisphere = 'N')) %>% rename(newLongitude = longitude, newLatitude = latitude)
#Join the new lat/lon values with the utm df to rejoin with main df
armstrong_utm_10 <- cbind(armstrong_utm_10, armstrong_latLon_10)

#Now join the utm df with the main df

armstrong <- left_join(armstrong, armstrong_utm_9)
armstrong <- left_join(armstrong, armstrong_utm_10)

#Now replot to check if fixed

# There apprears to be an issues with some of the coordiates, I suspect that the UTMs were converted incorrectly.
plot(canUSMex_map, xlim = c(-160, -110), ylim = c(40, 70), main = 'Armstong M. fusca data')
points(x=armstrong$Longitude, y=armstrong$Latitude, col = 'red', pch = 16) #old
points(x=armstrong_utm_9$newLongitude, y=armstrong_utm_9$newLatitude, col = 'blue',pch = 16) #new
points(x=armstrong_utm_10$newLongitude, y=armstrong_utm_10$newLatitude, col = 'blue', pch = 16) #new

legend(x= -150,
       y = 50,
       title = 'Recalculated UTM Occurrences',
       legend = c('Originial Points', 'Adjusted Points'),
       fill = c('red', 'blue'))



# Overlay GBIF and Armstrong data -----------------------------------------
occ_fus <- readRDS(file = './occ_data/fus/occ_fus.Rdata')
occ_armstrong <- read.csv(file = "./occ_data/fus/malus_fusca_armstrong.csv") #Note that this is only Armstrong data (not a combination of GBIF data)


occ_fus <- vect(occ_fus, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")
occ_armstrong <- vect(occ_armstrong, geom = c('Longitude', 'Latitude'),
                crs = "+proj=longlat +datum=WGS84")


plot(canUSMex_map, xlim = c(-160, -110), ylim = c(40, 70), main = 'Armstong M. fusca data')
points(x=occ_fus, col = 'red', pch = 16) #old
points(occ_armstrong, col = 'blue',pch = 16) #new



