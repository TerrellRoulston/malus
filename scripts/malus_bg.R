# Top ---------------------------------------------------------------------
# Setting Background area and predictor vars for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024

source("scripts/occ_thin.R") ## also loads data, maps; takes about 20 seconds

message("** Calculating Background Data: ", date())
## # plot basemap
## plot(canUSMex_map, xlim = c(-180, -50))
## # plot ecoregions 
## lines(ecoNA, col = 'red')

## plot(ecoNA)
## points(occThin_cor, col = 'red')

## dev.off()

# M. coronaria ecoregions -------------------------------------------------
# extract ecoregion polygon that contain M. coronaria occurrence points
eco_cor <- intersect(occThin_cor, ecoNA) %>% # extract ecoregion names from points
  as.data.frame() %>% # convert to df
  group_by(NA_L2CODE) %>% # sample 1 point of each eco region
  sample_n(1)

# return vector of eco region codes of the polygons that contain occurrences
eco_cor_code <- eco_cor$NA_L2CODE %>% unique() 
eco_cor_code <- eco_cor_code[eco_cor_code != '0.0']  #remove the 'water' '0.0' ecoregion

#CODES: "8.1" "8.2" "5.3" "8.4" "8.3" "8.5" "9.2" "9.4" "5.2"

ecoNA_cor <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_cor_code) # subset eco region spat vector by the codes

## plot(ecoNA_cor) # plot the subseted M. coronaria eco regions
## points(occThin_cor, pch = 3, col = 'red') # plot M. coronaria points

# M. fusca eco regions ----------------------------------------------------
# # extract eco region polygon that contain M. fusca occurrence points
eco_fus <- intersect(occThin_fus, ecoNA) %>% # extract ecoregion names from points
  as.data.frame() %>% # convert to df
  group_by(NA_L2CODE) %>% # sample 1 point of each eco region
  sample_n(1)

# return vector of eco region codes of the polygons that contain occurrences
eco_fus_code <- eco_fus$NA_L2CODE %>% unique() 
eco_fus_code <- eco_fus_code[eco_fus_code != '0.0'] # remove NA value

#CODES "7.1""6.2"  "10.1" "11.1" "10.2"

ecoNA_fus <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_fus_code)

## plot(ecoNA_fus)
## points(occThin_fus, pch = 3, col = 'red') # plot M. coronaria points


# M. ionesis eco regions --------------------------------------------------
# # extract eco region polygon that contain M. ionesis occurrence points
eco_ion <- intersect(occThin_ion, ecoNA) %>% # extract ecoregion names from points
  as.data.frame() %>% # convert to df
  group_by(NA_L2CODE) %>% # sample 1 point of each eco region
  sample_n(1)

# return vector of eco region codes of the polygons that contain occurrences
eco_ion_code <- eco_ion$NA_L2CODE %>% unique() 
eco_ion_code <- eco_ion_code[eco_ion_code != '0.0'] # remove NA value

#CODES "5.2" "8.1" "8.2" "8.3" "8.4" "8.5" "9.2" "9.4"

ecoNA_ion <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_ion_code)

## plot(ecoNA_ion)
## points(occThin_ion, pch = 3, col = 'red') 


# M. angustifolia eco regions ---------------------------------------------
# extract eco region polygon that contain M. angustifolia occurrence points
eco_ang <- intersect(occThin_ang, ecoNA) %>% # extract ecoregion names from points
  as.data.frame() %>% # convert to df
  group_by(NA_L2CODE) %>% # sample 1 point of each eco region
  sample_n(1)

# return vector of eco region codes of the polygons that contain occurrences
eco_ang_code <- eco_ang$NA_L2CODE %>% unique() 
eco_ang_code <- eco_ang_code[eco_ang_code != '0.0'] # remove NA value

#CODES "5.3" "8.1" "8.2" "8.3" "8.4" "8.5" "9.5"

ecoNA_ang <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_ang_code)

## plot(ecoNA_ang)
## points(occThin_ang, pch = 3, col = 'red') 


# Sect. Chloromeles eco region --------------------------------------------
# extract eco region polygon that contain sect. Chloromeles occurrence points
eco_chl <- intersect(occThin_chl, ecoNA) %>% # extract ecoregion names from points
  as.data.frame() %>% # convert to df
  group_by(NA_L2CODE) %>% # sample 1 point of each eco region
  sample_n(1)

# return vector of eco region codes of the polygons that contain occurrences
eco_chl_code <- eco_chl$NA_L2CODE %>% unique() 
eco_chl_code <- eco_chl_code[eco_chl_code != '0.0'] # remove NA value

#CODES "5.3" "8.1" "8.2" "8.3" "8.4" "8.5" "9.5"

ecoNA_chl <- subset(ecoNA, ecoNA$NA_L2CODE %in% eco_chl_code)

## plot(ecoNA_chl)
## points(occThin_chl, pch = 3, col = 'red') 


## NB: could save this as shapefiles instead!
# SAVE
## saveRDS(ecoNA_cor, file = './maps/eco_regions/ecoNA_cor.Rdata')
## saveRDS(ecoNA_fus, file = './maps/eco_regions/ecoNA_fus.Rdata')
## saveRDS(ecoNA_ion, file = './maps/eco_regions/ecoNA_ion.Rdata')
## saveRDS(ecoNA_ang, file = './maps/eco_regions/ecoNA_ang.Rdata')
## saveRDS(ecoNA_chl, file = './maps/eco_regions/ecoNA_chl.Rdata')

# Crop WorldClim to Ecoregions and Create Background ----------------------

# crop+mask extent of WorldClim data to the Malus ecoregions

wclim_cor_subs <- terra::crop(wclim_subs, ecoNA_cor, mask = T) ##
wclim_fus_subs <- terra::crop(wclim_subs, ecoNA_fus, mask = T) ##
wclim_ion_subs <- terra::crop(wclim_subs, ecoNA_ion, mask = T) ##
wclim_ang_subs <- terra::crop(wclim_subs, ecoNA_ang, mask = T) ##
wclim_chl_subs <- terra::crop(wclim_subs, ecoNA_chl, mask = T) ##
# Save cropped wclim data for downsteam SDM workflow

## writeRaster(wclim_cor, file = './wclim_data/wclim_cor.tif')
## writeRaster(wclim_fus, file = './wclim_data/wclim_fus.tif')
## writeRaster(wclim_ion, file = './wclim_data/wclim_ion.tif')
## writeRaster(wclim_ang, file = './wclim_data/wclim_ang.tif')
## writeRaster(wclim_chl, file = './wclim_data/wclim_chl.tif')

# Sample background points ------------------------------------------------


set.seed(1337) # set a seed to ensure reproducible results

# NOTE: Set arguments as.raster = T to return raster
# OR as.points to return spatvector = T to return spatvector

# M. coronaria background
# SpatVector

# Note upped bg points from 5000 to 20000

########################################################################
## TERRELL: You are not using any of the objects created past here, I ##
## think this code is redundant? Maybe useful for inspecting data,    ##
## but used in the actual analysis anywhere                           ##
########################################################################

## cor_bg_vec <- spatSample(wclim_cor, 20000, 'random', na.rm = T,
##                          as.points = T) #ignore NA values      

## plot(wclim_cor[[1]])

## points(cor_bg_vec, cex = 0.01)

## expanse(wclim_cor[[1]], unit = 'km') # total area of raster in km^2
# 5683684 km^2
##20000/5683684
# 0.00035 samples/km

# M. fusca background
# SpatVector
##fus_bg_vec <- spatSample(wclim_fus, 20000, 'random', na.rm = T, as.points = T) #ignore NA values
## plot(wclim_fus[[1]])
## points(fus_bg_vec, cex = 0.01)

## expanse(wclim_fus[[1]], unit = 'km') # total area of raster in km^2
# 4659175 km^2
# 5000/4659175 = 0.00107 samples/km

# M. ionesis background
# SpatVector
##ion_bg_vec <- spatSample(wclim_ion, 20000, 'random', na.rm = T, as.points = T) #ignore NA values
## plot(wclim_ion[[1]])
## points(ion_bg_vec, cex = 0.01)

## expanse(wclim_ion[[1]], unit = 'km') # total area of raster in km^2

# M. angustifolia background
# SpatVector
##ang_bg_vec <- spatSample(wclim_ang, 20000, 'random', na.rm = T, as.points = T) #ignore NA values
## plot(wclim_ang[[1]])
## points(ang_bg_vec, cex = 0.01)

## expanse(wclim_ang[[1]], unit = 'km') # total area of raster in km^2


# Sect. Chloromeles background
# SpatVector
##chl_bg_vec <- spatSample(wclim_chl, 20000, 'random', na.rm = T, as.points = T) #ignore NA values
## plot(wclim_chl[[1]])
## points(chl_bg_vec, cex = 0.01)

## expanse(wclim_chl[[1]], unit = 'km') # total area of raster in km^2


## # Save background SpatVectors
## saveRDS(cor_bg_vec, file = './occ_data/cor/cor_bg_vec.Rdata')
## saveRDS(fus_bg_vec, file = './occ_data/fus/fus_bg_vec.Rdata')
## saveRDS(ion_bg_vec, file = './occ_data/ion/ion_bg_vec.Rdata')
## saveRDS(ang_bg_vec, file = './occ_data/ang/ang_bg_vec.Rdata')
## saveRDS(chl_bg_vec, file = './occ_data/chl/chl_bg_vec.Rdata')

## # Load background SpatVectors

## cor_bg_vec <- readRDS(file = './occ_data/cor/cor_bg_vec.Rdata')
## fus_bg_vec <- readRDS(file = './occ_data/fus/fus_bg_vec.Rdata')
## ion_bg_vec <- readRDS(file = './occ_data/ion/ion_bg_vec.Rdata')
## ang_bg_vec <- readRDS(file = './occ_data/ang/ang_bg_vec.Rdata')
## chl_bg_vec <- readRDS(file = './occ_data/chl/chl_bg_vec.Rdata')

# Extracting presence-background raster values ----------------------------

## cor_predvals <- extract(wclim_cor, occThin_cor, ID = FALSE) # M. coronaria

## fus_predvals <- extract(wclim_fus, occThin_fus, ID = FALSE) # M. fusca

## ion_predvals <- extract(wclim_ion, occThin_ion, ID = FALSE) # M. fusca

## ang_predvals <- extract(wclim_ang, occThin_ang, ID = FALSE) # M. fusca

## chl_predvals <- extract(wclim_chl, occThin_chl, ID = FALSE) # M. fusca

## cor_bgvals <- values(cor_bg_vec) # Extract raster values for bg points
## fus_bgvals <- values(fus_bg_vec) # Extract raster values for bg points
## ion_bgvals <- values(ion_bg_vec) # Extract raster values for bg points
## ang_bgvals <- values(ang_bg_vec) # Extract raster values for bg points
## chl_bgvals <- values(chl_bg_vec) # Extract raster values for bg points


# Create a df for presence-background raster values for SDM ---------------

## cor_pb <- c(rep(1, nrow(cor_predvals)), rep(0, nrow(cor_bgvals))) #T/F presence or background string
## fus_pb <- c(rep(1, nrow(fus_predvals)), rep(0, nrow(fus_bgvals))) #T/F presence or background string
## ion_pb <- c(rep(1, nrow(ion_predvals)), rep(0, nrow(ion_bgvals))) #T/F presence or background string
## ang_pb <- c(rep(1, nrow(ang_predvals)), rep(0, nrow(ang_bgvals))) #T/F presence or background string
## chl_pb <- c(rep(1, nrow(chl_predvals)), rep(0, nrow(chl_bgvals))) #T/F presence or background string

# combine presence and background dataframes for SDM

## cor_sdmData <- data.frame(cbind(cor_pb, rbind(cor_predvals, cor_bgvals)))
## fus_sdmData <- data.frame(cbind(fus_pb, rbind(fus_predvals, cor_bgvals)))
## ion_sdmData <- data.frame(cbind(ion_pb, rbind(ion_predvals, cor_bgvals)))
## ang_sdmData <- data.frame(cbind(ang_pb, rbind(ang_predvals, cor_bgvals)))
## chl_sdmData <- data.frame(cbind(chl_pb, rbind(chl_predvals, cor_bgvals)))

#Save df for downstream SDM work

## saveRDS(cor_sdmData, file = './occ_data/cor/cor_sdmData.Rdata')
## saveRDS(fus_sdmData, file = './occ_data/fus/fus_sdmData.Rdata')
## saveRDS(ion_sdmData, file = './occ_data/ion/ion_sdmData.Rdata')
## saveRDS(ang_sdmData, file = './occ_data/ang/ang_sdmData.Rdata')
## saveRDS(chl_sdmData, file = './occ_data/chl/chl_sdmData.Rdata')

## cor_sdmData <- readRDS(file = 'cor_sdmData.Rdata')

# Check for colinearity of predictor variables for presence-bg -----------
# Although Maxent is equipped to handle variable selection through its built in variable regularization,
# it is still important to under the relationships between predictors.
# Dendograms useful for identifying groupings or clusters of colinear variables

## # M. coronaria
## pairs(cor_sdmData[,-1]) # drop the first column of 0/1

## # Dendogram cluster of predictor colinearity
## threshold <- 0.7 # set the threshold for colinearity 
## cor_cors <- raster.cor.matrix(wclim_cor) # pearson correlation
## cor_dist <- as.dist(1 - abs(cor_cors)) # calculate distance of predictors

## cor_clust <- hclust(cor_dist, method = 'single') # calculate cluster dendogram
## cor_groups <- cutree(cor_clust, h = 1 - threshold) #calculate groupings of variables

## plot(cor_clust, hang = -1, main = "M. coronaria Predictors")
## rect.hclust(cor_clust, h = 1 - threshold)


## #M. fusca
## pairs(fus_sdmData[,2:5])
## pairs(fus_sdmData[,6:9])
## pairs(fus_sdmData[,10:13])
## pairs(fus_sdmData[,14:17])
## pairs(fus_sdmData[,17:20])

## pairs(fus_sdmData[,-1])

## # Dendogram cluster of predictor colinearity
## threshold <- 0.7 # set the threshold for colinearity 
## fus_cors <- raster.cor.matrix(wclim_fus)
## fus_dist <- as.dist(1 - abs(fus_cors)) # calculate distance of predictors

## fus_clust <- hclust(fus_dist, method = 'single') # calculate cluster dendogram
## fus_groups <- cutree(fus_clust, h = 1 - threshold) #calculate groupings of variables

## plot(fus_clust, hang = -1, main = "M. fusca Predictors")
## rect.hclust(fus_clust, h = 1 - threshold)



# Predictor Kernel Density Plots -------------------------------------------------
# It is helpful to visualize two predictors pairs of
# Presence points and background points


## cor_occ.temp <- cor_sdmData %>% filter(cor_pb == 1) %>% # Presence points
##   dplyr::select(wc2.1_2.5m_bio_1) %>% # Mean annual temp
##   drop_na() %>% 
##   unlist()
## cor_bg.temp <- cor_sdmData %>% filter(cor_pb == 0) %>% # Background points
##   dplyr::select(wc2.1_2.5m_bio_1) %>% # Mean annual temp
##   drop_na() %>% 
##   unlist()

## cor_occ.perc <- cor_sdmData %>% filter(cor_pb == 1) %>% # Presence points
##   dplyr::select(wc2.1_2.5m_bio_12) %>% # Annual precipitation
##   drop_na() %>% 
##   unlist()

## cor_bg.perc <- cor_sdmData %>% filter(cor_pb == 0) %>% # Background points
##   dplyr::select(wc2.1_2.5m_bio_12) %>% # Annual precipitation
##   drop_na() %>% 
##   unlist()


## cor_occ.3d <- kde2d(cor_occ.temp, cor_occ.perc)
## cor_bg.3d <- kde2d(cor_bg.temp, cor_bg.perc)

## #Plot 3D surface Kernel density estimation

## plot_cor.occ_3d <- plot_ly(x=cor_occ.3d$x, y=cor_occ.3d$y, z=cor_occ.3d$z) %>% 
##   add_surface() %>% 
##   layout(scene = list(xaxis = list(title = 'Mean Annual Temp (C)', autotick = F, nticks = 5, tickvals = list(0,5,10,15,20)), 
##                       yaxis = list(title = 'Mean Annual Percip. (mm)', tick0=0, tick1=2000, dtick=200), 
##                       zaxis = list(title = 'Kernel Density', tick0=0, tick1=0.001, dtick=0.0002)),
##          title = list(text = "<i>M. coronaria<i> Occurrence Points", 
##                       y = 0.95, x = 0.5, 
##                       xanchor = 'center', 
##                       yanchor = 'top'))

## plot_cor.occ_3d # run to view

## plot_cor.bg_3d <- plot_ly(x=cor_bg.3d$x, y=cor_bg.3d$y, z=cor_bg.3d$z) %>% 
##   add_surface() %>% 
##   layout(scene = list(xaxis = list(title = 'Mean Annual Temp (C)', tick0=0, tick1=20, dtick=5), 
##                       yaxis = list(title = 'Mean Annual Percip. (mm)', tick0=0, tick1=2000, dtick=200), 
##                       zaxis = list(title = 'Kernel Density')),
##          title = list(text = "<i>M. coronaria<i> Background Points", 
##                       y = 0.95, x = 0.5, 
##                       xanchor = 'center', 
##                       yanchor = 'top'))

## plot_cor.bg_3d

message("** Background Data Calculated: ", date())

