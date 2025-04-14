# Top ---------------------------------------------------------------------
# Plotting of M. coronaria and M. fusca occurrences with Protected Areas
# Terrell Roulston
# Started Jan 21 2025

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package
library(geodata)

# Load occurrence data ----------------------------------------------------
#NOTE: This is not the thinned occurrence data

occ_cor <- readRDS(file = "./occ_data/cor/occ_cor_gbif.Rdata") # note this is cleaned occurrence data from GBIF ONLY
occ_fus <- readRDS(file = './occ_data/fus/occ_fus_gbif.Rdata') # note this is cleaned occurrence data from GBIF ONLY

occ_cor <- vect(occ_cor, geom = c('decimalLongitude', 'decimalLatitude'), crs = "+proj=longlat +datum=WGS84")
occ_fus <- vect(occ_fus, geom = c('decimalLongitude', 'decimalLatitude'), crs = "+proj=longlat +datum=WGS84")

# Load raster and geoms for plotting --------------------------------------
#International and national  basemaps
us_map <- gadm(country = 'USA', level = 1, resolution = 2, path = "./maps/base_maps") #USA w. States
us_map_0 <- gadm(country = 'USA', level = 0, resolution = 2, path = "./maps/base_maps") #USA w.o States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2, path = './maps/base_maps') #Canada w. Provinces
ca_map_0 <- gadm(country = 'CA', level = 0, resolution = 2, path = './maps/base_maps') #Canada w.o Provinces

mex_map <-gadm(country = 'MX', level = 1, resolution = 2, path = './maps/base_maps') # Mexico w. States
mex_map_0 <-gadm(country = 'MX', level = 0, resolution = 2, path = './maps/base_maps') # Mexico w.o States


gl_map_0 <- gadm(country = 'GL', level = 0, resolution = 2, path = './maps/base_maps') # Greenland w.o States
# Shape file for Canada/US/Mexico borders

can_us_mex_border <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/can_us_mex_border')

# Shape files downloaded from the USGS (https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')


# Protected Areas of Canada
pro_area <- terra::vect("C:/Users/terre/Documents/Acadia/Malus Project/maps/canada_pa/ProtectedConservedArea_2023/ProtectedConservedArea_2023/ProtectedConservedArea_2023.gdb")
pro_area <- pro_area[pro_area$BIOME == 'T'] # filter only the terrestrial protected areas and OECMs
pro_area <- pro_area[pro_area$PA_OECM_DF == '1'] # filter only protected areas (remove OECMs)

#extract basemaps of provinces for crop PA raster to speed up plotting
on_basemap <- ca_map[ca_map$NAME_1 == 'Ontario']
on_basemap_proj <- project(on_basemap, crs(pro_area))
on_pa <- intersect(pro_area, on_basemap_proj)

qc_basemap <- ca_map[ca_map$NAME_1 == 'Québec']
qc_basemap_proj <- project(qc_basemap, crs(pro_area))
qc_pa <- intersect(pro_area, on_basemap_proj)

bc_basemap <- ca_map[ca_map$NAME_1 == 'British Columbia']
bc_basemap_proj <- project(bc_basemap, crs(pro_area))
bc_pa <- intersect(pro_area, bc_basemap_proj)

wa_basemap <- us_map[us_map$NAME_1 == 'Washington']


# Map projection ----------------------------------------------------------
# Lambert projection
projLam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

us_map_0.lcc <- project(us_map_0, projLam)
ca_map_0.lcc <- project(ca_map_0, projLam)
mex_map_0.lcc <- project(mex_map_0, projLam)
gl_map_0.lcc <- project(gl_map_0, projLam)
on_basemap.lcc <- project(on_basemap, projLam)
qc_basemap.lcc <- project(qc_basemap, projLam)
bc_basemap.lcc <- project(bc_basemap, projLam)
wa_basemap.lcc <- project(wa_basemap, projLam)
occ_cor.lcc <- project(occ_cor, projLam)
occ_fus.lcc <- project(occ_fus, projLam)
can_us_mex_border.lcc <- project(can_us_mex_border, projLam)
great_lakes.lcc <- project(great_lakes, projLam)
on_pa.lcc <- project(on_pa, projLam)
qc_pa.lcc <- project(qc_pa, projLam)
bc_pa.lcc <-project(bc_pa, projLam)

# Some prep for adding names to maps
na_coords <- data.frame(
  name = c("Canada", "U.S.A.", "Mexico", 'Greenland', 'Hudson Bay', 'Atlantic Ocean', 'Pacific Ocean'),
  lon = c(-111.93, -101.53, -107.95, -39, -85.23, -66.58, -131.19),
  lat = c(60.10, 39.64, 29.61, 72, 59.64, 38.54, 48.48)
)
na_points <- vect(na_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
na_points.lcc <- project(na_points, projLam) # Reproject to Lambert Conformal Conic
na_points_coords <- data.frame(crds(na_points.lcc)) #extract coordinates
na_coords <- cbind(na_coords, na_points_coords) #mutate together

on_coords <- data.frame(
  name = c("Toronto", "Lake Erie", "Lake Ontario", "Lake Huron"),
  lon = c(-79.3832, -81.2, -77.9, -82.5),
  lat = c(43.6532, 42.2, 43.7, 44.5)
)
on_points <- vect(on_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
on_points.lcc <- project(on_points, projLam) # Reproject to Lambert Conformal Conic
on_points_coords <- data.frame(crds(on_points.lcc)) #extract coordinates
on_coords <- cbind(on_coords, on_points_coords) #mutate together

bc_coords <- data.frame(
  name = c("Vancouver", "Pacific Ocean", "Salish Sea", "Straigt of Georgia"),
  lon = c(-123.1215, -125.25, -123.12, -123.46),
  lat = c(49.2820, 48.44, 48.27, 49.21)
)
bc_points <- vect(bc_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
bc_points.lcc <- project(bc_points, projLam) # Reproject to Lambert Conformal Conic
bc_points_coords <- data.frame(crds(bc_points.lcc)) #extract coordinates
bc_coords <- cbind(bc_coords, bc_points_coords) #mutate together

# Plot occurrences and protected areas ------------------------------------
# North America plot with occurrences for both species
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/maps/malus_occ_and_pa/occ_north_america_v2.jpeg", width = 10000, height = 6666, res = 300)

na.xlim <- c(-3.5*10^6, 3.1*10^6)
na.ylim <- c(-1.8*10^6, 2.6*10^6)
terra::plot(ca_map_0.lcc , col = c('white'),
            background = 'lightblue',
            border = 'transparent',
            legend = F, 
            xlim = na.xlim, ylim = na.ylim, 
            main = "",
            cex.main = 3,
            axes = F,
            box = F,
            mar = c(1,1,1,1))
terra::plot(us_map_0.lcc, col = '#E8E8E8', add = T, border = 'transparent')
terra::plot(mex_map_0.lcc, col = 'white', add = T, border = 'transparent')
terra::plot(gl_map_0.lcc, col = '#E8E8E8', add = T, border = 'transparent')
terra::plot(great_lakes.lcc, box = F, add = T, col = 'lightblue', border = 'grey')
terra::plot(can_us_mex_border.lcc, box = F,  add = T, col = 'grey')
terra::points(occ_cor.lcc, pch = 16, col = alpha("magenta", 1))
terra::points(occ_fus.lcc, pch = 16, col = alpha("#39FF14", 1))
terra::add_box(col = 'grey')
legend(
  x = 1.8e6,
  y = 2.1e6,
  title = c(expression(underline('GBIF Occurrence Data'))),
  legend = c(expression(italic("Malus coronaria")), expression(italic("Malus fusca"))),
  fill = c("magenta", "#39FF14"),
  col = "black",
  box.col = "transparent",  # No border around legend
  bg = "transparent",
  text.col = 'black',
  cex = 3.5,
  xjust = 0,              
  yjust = 1,
  title.adj = 0.25
)
terra::text(533792.2, 1206373.8, labels = "Hudson\nBay", cex = 2.2, col = "steelblue")  
terra::text(2525623.3, -617100.0, labels = "Atlantic\nOcean", cex = 2.5, col = "steelblue")  
terra::text(-2838692.4, 684470.5, labels = "Pacific\nOcean", cex = 2.5, col = "steelblue")  
terra::text(-585874.4, 1236224, labels = 'Canada', cex = 3, col = "black")
terra::text(-585874.4, -1032531, labels = 'U.S.A.', cex = 3, col = "black")
terra::text(-2186283.9, 2489474.75, labels = 'U.S.A.', cex = 3, col = "black") # Alaska
terra::text(2035568.47, 2509411.76, labels = 'Greenland', cex = 2.5, col = 'black')

dev.off()

# Plot Malus coronaria in Southern Ontario
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/maps/malus_occ_and_pa/malus_coronaria_ontario_v3.jpeg",  width = 5000 , height = 6667, res = 300)

cor.xlim <- c(0.9e6, 1.5e6)
cor.ylim <- c(-0.872e6, -0.178e6)

terra::plot(on_basemap.lcc , col = c('white', 'white'),
            background = '#E8E8E8',
            border = 'grey',
            legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = "",
            cex.main = 3,
            axes = F,
            box = F,
            mar = c(1,1,1,1))
terra::plot(qc_basemap.lcc, col = 'white', xlim = cor.xlim, ylim = cor.ylim, axes = F, box = F, add = T, border = 'grey')
terra::plot(great_lakes.lcc, box = F, add = T, col = 'lightblue', border = 'grey')
terra::points(occ_cor.lcc, pch = 16, col = "magenta", cex = 1.5)
terra::plot(on_pa.lcc, box = F, add = T, border = 'black', lwd = 1, col = '#00000040')
terra::add_box(col = 'grey')
#terra::points(1277388 , -443572.4, pch = 9, cex = 2)
#terra::text(1297388 , -443572.4, labels = "Toronto", cex = 1.5, col = "black") 
terra::text(1167004, -640048.9, labels = "Lake Erie", cex = 2.5, col = "steelblue")  
terra::text(1394488 , -413268.6, labels = "Lake Ontario", cex = 2.5, col = "steelblue")  
terra::text(1007380 , -405710.4, labels = "Lake Huron", cex = 2.5, col = "steelblue")  
terra::text(1404488, -570048.9, labels = 'U.S.A.', cex = 3, col = "black")
terra::text(947380, -590048.9, labels = 'U.S.A.', cex = 3, col = "black")


dev.off()

# Plot Malus fusca in BC
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/maps/malus_occ_and_pa/malus_fucas_bc_v3.jpeg", width = 5000 , height = 6667, res = 300)

fus.xlim <- c(-2.25e6, -1.8e6)
fus.ylim <- c(0.3e6, 0.82e6)

terra::plot(bc_basemap.lcc , col = c('white', 'white'),
            background = 'lightblue',
            border = 'grey',
            legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = "",
            cex.main = 3,
            axes = F,
            box = F,
            mar = c(1,1,1,1))
terra::plot(wa_basemap.lcc, col = '#E8E8E8', add = T, border = 'grey')
terra::points(occ_fus.lcc, pch = 16, col = "#39FF14", cex = 1.5)
terra::plot(bc_pa.lcc, box = F, add = T, border = 'black', lwd = 1, col = '#00000040')
terra::add_box(col = 'grey')
#terra::points(-1977905 , 475844.1, pch = 9, cex = 3, col = 'black')
#terra::text(-1958000 , 489844.1, labels = "Vancouver", cex = 1.5, col = "black")  # Vacouver
terra::text(-2199490, 461068.4, labels = "Pacific Ocean", cex = 2.5, col = "steelblue")  # Lake Erie
terra::text(-2025998 , 378987.6, labels = "Salish Sea", cex = 2.2, col = "steelblue")  # Lake Ontario
terra::text(-2003575 , 479216.9, labels = "Georgia\nStrait", cex = 2, col = "steelblue")  # Lake Huron
terra::text(-1885998 , 335987.6, labels = 'U.S.A.', cex = 3, col = "black")

dev.off()

# Number of occ in ea Prov ------------------------------------------------
occ_cor_on <- crop(occ_cor, on_basemap)
length(occ_cor_on) #n=135
length(occ_cor) #n=1030

occ_fus_bc <- crop(occ_fus, bc_basemap)
length(occ_fus_bc) #n=696
length(occ_fus)#n=1461

occ_cor_on_pa <- crop(occ_cor.lcc, on_pa.lcc) # make sure the projections are the same! the pa is not WGS84
length(occ_cor_on_pa) #16

occ_fus_bc_pa <- crop(occ_fus.lcc, bc_pa.lcc)
length(occ_fus_bc_pa) #238
