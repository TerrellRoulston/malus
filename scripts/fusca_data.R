library(tidyverse)
library(geodata)
library(terra)

us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "./maps/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = './maps/base_maps') #Canada basemap w. Provinces

mex_map <-gadm(country = 'MX', level = 1, resolution = 2,
               path = './maps/base_maps') # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map


occ_fus <- readRDS(file = './occ_data/fus/occ_fus_gbif.Rdata')
occ_armstrong <- read.csv(file = "./occ_data/fus/malus_fusca_armstrong.csv") #Note that this is only Armstrong data (not a combination of GBIF data)
occ_orb_fit <- readxl::read_xlsx("C:/Users/terre/Documents/Acadia/Malus Project/Obrist_2017_Fitzpatrick_2020_M_fusca.xlsx") 
occ_wickham <- readxl::read_xlsx("C:/Users/terre/Documents/Acadia/Malus Project/Wickham_Malusfusca.xlsx") 

occ_fus <- vect(occ_fus, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")
occ_armstrong <- vect(occ_armstrong, geom = c('Longitude', 'Latitude'),
                      crs = "+proj=longlat +datum=WGS84")
occ_orb_fit <- vect(occ_orb_fit, geom = c('longitude', 'latitude'),
                    crs = "+proj=longlat +datum=WGS84")
occ_wickham <- vect(occ_wickham, geom = c('Longitude', 'Latitude'),
                    crs = "+proj=longlat +datum=WGS84")

plot(canUSMex_map, xlim = c(-152, -118), ylim = c(40, 63))
points(x=occ_fus, col = 'red', pch = 1, cex = 1) #old
points(occ_armstrong, col = 'blue',pch = 1, cex = 1) #new
points(occ_orb_fit, col = 'orange',pch = 1, cex = 1) #new new
points(occ_wickham, col = 'purple', pch = 1, cex = 1) #new new new
#rect(xleft = -129, ybottom = 51.25, xright = -126.5, ytop = 52.1, border = "green", lwd = 2.25)

legend(-150, 50, legend = c("GBIF", "Armstrong", "Wickham", "Obr. and Fit."), 
       col = c("red", "blue", "purple", "orange"), pch = 16,  bty = "y")


arcmin_res <- 2.5 / 60  # Convert 2.5 arc minutes to degrees
extent <- ext((-128.4684) - 0.1, (-126.8941) + 0.1, (51.45909) - 0.1, (52.06380) + 0.1)

# Calculate number of rows and columns based on the extent and resolution
ncols <- ceiling((extent[2] - extent[1]) / arcmin_res)  # Number of columns
nrows <- ceiling((extent[4] - extent[3]) / arcmin_res)  # Number of row

# Create a raster with a 2.5 arc minute resolution
raster_grid <- rast(extent, ncol = ncols, nrow = nrows, crs = "EPSG:4326")
values(raster_grid) <- 1

plot(ca_map, xlim = c(-128.4684, -126.8941), ylim = c(51.45909, 52.06380))

wclim_polygons <- as.polygons(wclim)
#plot(wclim, add = TRUE, col = NA, border = "black")
plot(wclim_polygons, add = TRUE, border = "black")

grid_x <- seq(extent[1], extent[2], by = arcmin_res)
grid_y <- seq(extent[3], extent[4], by = arcmin_res)

# Add grid lines to the plot
abline(h = grid_y, v = grid_x, col = "black", lwd = 1)

terra::points(occ_orb_fit, col = 'orange', pch=1, cex = 1)
terra::points(occ_wickham, col = 'purple', pch = 1, cex = 1)
points(x=occ_fus, col = 'red', pch = 1, cex = 1) #old

legend(-127.3, 52, legend = c("Wickham", "Orb. and Fit."), 
       col = c("purple", "orange"), pch = 16,  bty = "y")

