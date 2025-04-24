# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 20, 2024

library(tidyverse) # grammar and data management 
library(terra) # working with spatial data
library(geodata) # basemaps and climate data


# Load cleaned occurrence data ---------------------------------------------

## occ_cor <- readRDS(file = "./occ_data/cor/occ_cor.Rdata") # GBIF + Husband
## occ_fus <- readRDS(file = './occ_data/fus/occ_fus.Rdata') # GBIF + Armstrong + Wickham + Obr. + Fit
## occ_ion <- readRDS(file = './occ_data/ion/occ_ion.Rdata') # GBIF
## occ_ang <- readRDS(file = './occ_data/ang/occ_ang.Rdata') # GBIF

occ_cor_orig <- readRDS(file = "./occ_data/cor/occ_cor.Rdata") # GBIF + Husband
occ_fus_orig <- readRDS(file = './occ_data/fus/occ_fus.Rdata') # GBIF + Armstrong + Wickham + Obr. + Fit
occ_ion_orig <- readRDS(file = './occ_data/ion/occ_ion.Rdata') # GBIF
occ_ang_orig <- readRDS(file = './occ_data/ang/occ_ang.Rdata') # GBIF

occ_cor <- read.table(file = "./occ_data/cor/occ_cor.csv") # GBIF + Husband
occ_fus <- read.table(file = './occ_data/fus/occ_fus.csv') # GBIF + Armstrong + Wickham + Obr. + Fit
occ_ion <- read.table(file = './occ_data/ion/occ_ion.csv') # GBIF
occ_ang <- read.table(file = './occ_data/ang/occ_ang.csv') # GBIF

## M. coronaria: note one coastal record from New York is excluded in the
## new version:

dim(occ_cor)
dim(occ_cor_orig) 

## M. fusca: 226 fewer records in the latest version. Looks like duplicate
## records from the Wickham, Obrits and Fitsp datasets:

dim(occ_fus)
dim(occ_fus_orig)

## M. ioensis and M. angustifolia unchanged between previous and current
## version: 
dim(occ_ion)
dim(occ_ion_orig)

dim(occ_ang)
dim(occ_ang_orig)

# Combine occ data for the 3 Chloromeles species
occ_chl <- occ_cor %>% 
  full_join(occ_ion, by = c('species', 'source', 'decimalLongitude', 'decimalLatitude')) %>% 
  full_join(occ_ang, by = c('species', 'source', 'decimalLongitude', 'decimalLatitude')) %>% 
  dplyr::select('species', 'source', 'decimalLongitude', 'decimalLatitude')

# vectorize occurrence df to coordinates for below
occ_cor <- vect(occ_cor, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_fus <- vect(occ_fus, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_ion <- vect(occ_ion, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_ang <- vect(occ_ang, geom = c('decimalLongitude', 'decimalLatitude'),
                crs = "+proj=longlat +datum=WGS84")

occ_chl <- vect(occ_chl, geom = c('decimalLongitude', 'decimalLatitude'),
                           crs = "+proj=longlat +datum=WGS84")
# Download WorldClim Bioclimatic raster -----------------------------------

# Note DO NOT PUSH wclim data**
wclim <- worldclim_global(var = 'bio', res = 2.5, version = '2.1', path = "./wclim_data/")

# plot a raster to check it downloaded properly
#plot(wclim$wc2.1_2.5m_bio_1, main = 'Annual Mean Temperature')

# Thin data using sampler -------------------------------------------------

set.seed(1337) # set random generator seed to get reproducible results
# M. coronaria thinning
occThin_cor <- spatSample(occ_cor, size = 1, 
                      strata = wclim, #sample one occurrence from each climatic cell
                      method = "random") 

# M. fusca thinning
occThin_fus <- spatSample(occ_fus, size = 1,
                      strata = wclim,
                      method = "random")

#M. ionesis thinning
occThin_ion <- spatSample(occ_ion, size = 1,
                          strata = wclim,
                          method = "random")

#M. angustifolia thinning
occThin_ang <- spatSample(occ_ang, size = 1,
                          strata = wclim,
                          method = "random")
#Sect. Chloromeles thinning
occThin_chl <- spatSample(occ_chl, size = 1,
                          strata = wclim,
                          method = "random")
# Save thinned occurrence points for further analysis ----------------------

##saveRDS(occThin_cor, file = './occ_data/cor/occThin_cor.Rdata')
##saveRDS(occThin_fus, file = './occ_data/fus/occThin_fus.Rdata')
##saveRDS(occThin_ion, file = './occ_data/ion/occThin_ion.Rdata')
##saveRDS(occThin_ang, file = './occ_data/ang/occThin_ang.Rdata')
##saveRDS(occThin_chl, file = './occ_data/chl/occThin_chl.Rdata')

write.table(occThin_cor, file = './occ_data/cor/occThin_cor.csv')
write.table(occThin_fus, file = './occ_data/fus/occThin_fus.csv')
write.table(occThin_ion, file = './occ_data/ion/occThin_ion.csv')
write.table(occThin_ang, file = './occ_data/ang/occThin_ang.csv')
write.table(occThin_chl, file = './occ_data/chl/occThin_chl.csv')

# Extract environmental data from wclim ------------------------------------
# Extract and add wclim raster data to spatial occurrence data

#occ_cor <- cbind(occ_cor, extract(wclim, occ_cor)) # M. coronaria
#occ_fus <- cbind(occ_fus, extract(wclim, occ_fus)) # M. fusca

# Plot thinned occurrence data --------------------------------------------
us_map_0 <- gadm(country = 'USA', level = 0, resolution = 2, path = "./maps/base_maps") #USA w.o States
ca_map_0 <- gadm(country = 'CA', level = 0, resolution = 2, path = './maps/base_maps') #Canada w.o Provinces
mex_map_0 <-gadm(country = 'MX', level = 0, resolution = 2, path = './maps/base_maps') # Mexico w.o States
gl_map_0  <-gadm(country = 'GL', level = 0, resolution = 2, path = './maps/base_maps') # Mexico w.o States

caribbean_codes <- c("BS", "CU", "JM", "HT", "DO", "PR", "BM", "TC", "KY") # Caribbean codes
gadm_list <- lapply(caribbean_codes, function(code) {
  tryCatch(
    gadm(country = code, level = 0, path = './maps/base_maps'),  # Save to specified path
    error = function(e) NULL
  )
})
gadm_list <- gadm_list[!sapply(gadm_list, is.null)]
# Combine all downloaded boundaries into a single spatial object
car_map_0 <- do.call(rbind, gadm_list) # Spatvertor of Caribbean Islands

# Shape file for Canada/US/Mexico borders

###########################################
## TERRELL: this path is not accessible! ##
###########################################

## can_us_mex_border <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/can_us_mex_border')

can_us_mex_border <- rbind(us_map_0, ca_map_0, mex_map_0)

# Shape files downloaded from the USGS (https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
erie <- vect("~/data/maps/hydro_p_LakeErie.shp")
ontario <- vect("~/data/maps/hydro_p_LakeOntario.shp")
michigan <- vect("~/data/maps/hydro_p_LakeMichigan.shp")
huron <- vect("~/data/maps/hydro_p_LakeHuron.shp")
superior <- vect("~/data/maps/hydro_p_LakeSuperior.shp")
stclair <- vect("~/data/maps/hydro_p_LakeStClair.shp")
great_lakes <- rbind(erie, ontario, michigan, huron, superior, stclair)

## great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')

# separate GBIF from other data sources
occThin_cor_gbif <- subset(occThin_cor,occThin_cor$source == 'GBIF')
occThin_cor_hus <- subset(occThin_cor,occThin_cor$source == 'Husband')

occThin_fus_gbif <- subset(occThin_fus, occThin_fus$source == 'GBIF')
occThin_fus_arm <- subset(occThin_fus, occThin_fus$source == 'Armstrong')
occThin_fus_obr_fit <- subset(occThin_fus, occThin_fus$source %in% c('Obrist', 'Fitzpatrick'))
occThin_fus_wick <- subset(occThin_fus, occThin_fus$source == 'Wickham')

#Project data
projLam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

us_map_0.lcc <- project(us_map_0, projLam)
ca_map_0.lcc <- project(ca_map_0, projLam)
mex_map_0.lcc <- project(mex_map_0, projLam)
gl_map_0.lcc <- project(gl_map_0, projLam)
car_map_0.lcc <- project(car_map_0, projLam)
occThin_cor_gbif.lcc <- project(occThin_cor_gbif, projLam)
occThin_cor_hus.lcc <- project(occThin_cor_hus, projLam)
occThin_fus_gbif.lcc <- project(occThin_fus_gbif, projLam)
occThin_fus_arm.lcc <- project(occThin_fus_arm, projLam)
occThin_fus_obr_fit.lcc <- project(occThin_fus_obr_fit, projLam)
occThin_fus_wick.lcc <- project(occThin_fus_wick, projLam)
occThin_ion.lcc <- project(occThin_ion, projLam)
occThin_ang.lcc <- project(occThin_ang, projLam)
occThin_chl.lcc <- project(occThin_chl, projLam)
can_us_mex_border.lcc <- project(can_us_mex_border, projLam)
great_lakes.lcc <- project(great_lakes, projLam)


# Plot thinned occurrence points ------------------------------------------

#Start plotting
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plots/occ_plot/malus_occurrence_data_sources_v3.jpeg", width = 10000, height = 6666, res = 300)

# Plot M. fusca and M. coronaria, M. ionesis and M. angustifolia 
na.xlim <- c(-3.5*10^6, 3.1*10^6) #plot lims dependent on CSR
na.ylim <- c(-2.3*10^6, 2.6*10^6)
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
terra::plot(car_map_0.lcc, col = 'white', add = T, border = 'transparent')
terra::plot(great_lakes.lcc, box = F, add = T, col = 'lightblue', border = 'grey')
terra::plot(can_us_mex_border.lcc, box = F,  add = T, border = 'grey')
terra::points(occThin_cor_gbif.lcc, pch = 16, col = alpha("magenta", 1), cex = 1.3) # COR GBIF
terra::points(occThin_fus_gbif.lcc, pch = 16, col = alpha("#228B22", 1), cex = 1.3) # FUS GBIF
terra::points(occThin_cor_hus.lcc, pch = 16, col = alpha("#C90076", 1), cex = 1.3) # COR HUSBAND
terra::points(occThin_fus_arm.lcc, pch = 16, col = alpha("#333f07", 1), cex = 1.3) # FUS ARMSTRONG
terra::points(occThin_fus_wick.lcc, pch = 16, col = alpha("#333f07", 1), cex = 1.3) # FUS WICKHAM
terra::points(occThin_fus_obr_fit.lcc, pch = 16, col = alpha("#333f07", 1), cex = 1.3) # FUS OBRITS ET AL. AND FITZPATRICK ET AL.
terra::points(occThin_ang.lcc, pch = 16, col = alpha("#007CBE", 1), cex = 1.3) # ANG GBIF
terra::points(occThin_ion.lcc, pch = 16, col = alpha("#E88E00", 1), cex = 1.3) # ION GBIF
terra::add_box(col = 'grey')
legend( # legend for data
  x = 1.7e6,
  y = 2.3e6,
  title = c(expression(underline('Thinned Occurrence Data'))),
  legend = c(expression(italic("Malus coronaria")*"—GBIF"), 
             expression(italic("Malus coronaria")* "—Other Data"), 
             expression(italic("Malus fusca")*"—GBIF"), 
             expression(italic("Malus fusca")* "—Other Data"), 
             expression(italic("Malus angustifolium")*"—GBIF"), 
             expression(italic("Malus ioensis")*"—GBIF")),
  fill = c("magenta", "#C90076","#228B22", "#333f07", "#007CBE", "#E88E00"),
  col = "black",
  box.col = "black",  # No border around legend
  bg = "white",
  text.col = 'black',
  cex = 2,
  xjust = 0,              
  yjust = 1,
  title.adj = 0.25
)
terra::text(533792.2, 1206373.8, labels = "Hudson\nBay", cex = 1.2, col = "steelblue")  
terra::text(2525623.3, -617100.0, labels = "North\nAtlantic\nOcean", cex = 1.5, col = "steelblue")  
terra::text(-2838692.4, 684470.5, labels = "North\nPacific\nOcean", cex = 1.5, col = "steelblue")  
terra::text(-585874.4, 1236224, labels = 'Canada', cex = 2, col = "black")
terra::text(-585874.4, -1032531, labels = 'U.S.A.', cex = 2, col = "black")

dev.off()



# Plot ONLY sect. Chloromeles

#Start plotting
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/occ_plot/chloromeles_occurrence_data_v1.jpeg", width = 10000, height = 6666, res = 300)

na.xlim <- c(-3.5*10^6, 3.1*10^6) #plot lims dependent on CSR
na.ylim <- c(-2.3*10^6, 2.6*10^6)
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
terra::plot(car_map_0.lcc, col = 'white', add = T, border = 'transparent')
terra::plot(great_lakes.lcc, box = F, add = T, col = 'lightblue', border = 'grey')
terra::plot(can_us_mex_border.lcc, box = F,  add = T, border = 'grey')
terra::points(occThin_chl.lcc, pch = 16, col = alpha("magenta", 1), cex = 1.3) # Sect. Chloromeles
terra::add_box(col = 'grey')
legend( # legend for data
  x = 1.7e6,
  y = 2.3e6,
  title = c(expression(underline('Thinned Occurrence Data'))),
  legend = c(expression("Sect. Chloromeles")),
  fill = c("magenta"),
  col = "black",
  box.col = "black",  # No border around legend
  bg = "white",
  text.col = 'black',
  cex = 2,
  xjust = 0,              
  yjust = 1,
  title.adj = 0.25
)
terra::text(533792.2, 1206373.8, labels = "Hudson\nBay", cex = 1.2, col = "steelblue")  
terra::text(2525623.3, -617100.0, labels = "North\nAtlantic\nOcean", cex = 1.5, col = "steelblue")  
terra::text(-2838692.4, 684470.5, labels = "North\nPacific\nOcean", cex = 1.5, col = "steelblue")  
terra::text(-585874.4, 1236224, labels = 'Canada', cex = 2, col = "black")
terra::text(-585874.4, -1032531, labels = 'U.S.A.', cex = 2, col = "black")

dev.off()
