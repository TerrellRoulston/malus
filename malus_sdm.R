# Top ---------------------------------------------------------------------
# MaxEnt Species Distribution Modeling (SDM) for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package
library(predicts) # SDM package
library(geodata) # basemaps
library(rJava) # MaxEnt models are dependant on JDK
library(ENMeval) # Another modeling package, useful for data partitioning (Checkerboarding)
library(raster) # RasterStack dependancy (a now deprecated function)
library(ecospat) # Useful spatial ecology tools
library(parallel) # speed up computation by running in parallel
library(doParallel) # added functionality to parallel


# Load occurrence data and basemaps -------------------------------------------------------
getwd() # check you directory location

# Background points in SpatVectors
setwd("./occ_data/")
cor_bg_vec <- readRDS(file = 'cor_bg_vec.Rdata')
fus_bg_vec <- readRDS(file = 'fus_bg_vec.Rdata')

# Occurrence Points in SpatVectors
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = 'occThin_fus.Rdata') # M. fusca

# Download/load basemaps
us_map <- gadm(country = 'USA', level = 0, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 0, resolution = 2,
               path = '../occ_data/base_maps') # Canada basemap w. Provinces

mex_map <-gadm(country = 'MX', level = 0, resolution = 2,
               path = '../occ_data/base_maps') # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map


NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in Western Hemisphere

canUSMex_map <- crop(canUSMex_map, NA_ext) # crop to Western Hemisphere

plot(canUSMex_map) # plot basemap


# Great Lakes shapefiles for making pretty maps and cropping
great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')
great_lakes <- crop(great_lakes, NA_ext)


# Download/load WorldClim data under future climate scenarios -------------
# WARNING DO NOT PUSH WORLDCLIM DATA
setwd('../wclim_data/')
# Historical climate 1970-2000
wclim <- geodata::worldclim_global(var = 'bio',
                                   res = 2.5, 
                                   version = '2.1', 
                                   path = "../wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

# SSP (Shared social-economic pathway) 2.45 
# middle of the road projection, high climate adaptation, low climate mitigation
ssp245_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp245_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp245_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

# SPP 5.85 
# low regard for enviromental sustainability, increased fossil fuel reliance, this is the current tracking projection
ssp585_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp585_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp585_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "../wclim_data/")%>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

# Load cropped climate Rasters --------------------------------------------
# These Rasters are useful for sampling spatial checkerboards 
# and making habitat suitability predictions (Historical and under future SSPs climate scenarios)

setwd('../wclim_data')
# Historical (1970-2000)
wclim_cor <- readRDS(file = 'wclim_cor.Rdata') 
#wclim_cor_stack <- raster::stack(wclim_cor) # covert SpatRaster to RasterStack for dependency in ENMeval checkboarding

wclim_fus <- readRDS(file = 'wclim_fus.Rdata')
#wclim_fus_stack <- raster::stack(wclim_fus) # covert SpatRaster to RasterStack for dependency in ENMeval checkboarding

climate_predictors <- names(wclim_cor) # extract climate predictor names, to rename layers in the rasters below
# This is important to do for making predictions once the SDMs have been made on future climate data
# Note that the names of the layers still correspond to the same environmental variables

# Future SSPs
# Do not need to create RasterStacks
# SSP 245
names(ssp245_2030) <- climate_predictors #rename raster layers for downsteam analysis
names(ssp245_2050) <- climate_predictors 
names(ssp245_2070) <- climate_predictors 

# SSP 585
names(ssp585_2030) <- climate_predictors #rename raster layers for downsteam analysis
names(ssp585_2050) <- climate_predictors 
names(ssp585_2070) <- climate_predictors 



# Coronaria - MaxEnt Model  ------------------------------------------------

# Spatial partitioning preparation
occ_cor_coords <- as.data.frame(geom(occThin_cor)[,3:4]) # extract longitude, lattitude from occurence points
bg_cor_coords <- as.data.frame(geom(cor_bg_vec)[,3:4]) # extract longitude, lattitude from background points

occ_fus_coords <- as.data.frame(geom(occThin_fus)[,3:4]) # extract longitude, lattitude from occurence points
bg_fus_coords <- as.data.frame(geom(fus_bg_vec)[,3:4]) # extract longitude, lattitude from background points

# Build Species Distribution Model using MaxEnt from the <ENMeval> package

# Run prediction in a parallel using 'socket' clusters to help speed up computation
# <ENMeval> implements parallel functions natively
# But load <parallel> library for additional functions like <decectCores()>
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

# current version of maxent.jar =  v3.4.4

cor_maxent <- ENMevaluate(occ_cor_coords, # occurrence records
                            envs = wclim_cor, # environment from background training area
                            n.bg = 20000, # 20000 bg points
                            tune.args =
                              list(rm = seq(0.5, 8, 0.5),
                                   fc = c("L", "LQ", "H",
                                          "LQH", "LQHP", "LQHPT")),
                            partition.settings =
                              list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                            partitions = 'checkerboard2',
                            parallel = TRUE,
                            numCores = cn - 1, # leave one core available for other apps
                            parallelType = "doParallel", # use doParrallel on Windows - socket cluster  
                            algorithm = 'maxent.jar')


# Save the MaxEnt model so you do not have to waste time re-running the model
setwd('../sdm_output')
saveRDS(cor_maxent, file = 'cor_maxent.Rdata') # save
cor_maxent <- readRDS(file = 'cor_maxent.Rdata') # load 


# M. coronaria Model Selection --------------------------------------------
# Note that maxent results provide Continuous Boyce Index (cbi)
best_cor_maxent <- subset(cor_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_cor_maxent <- eval.models(cor_maxent)[[best_cor_maxent$tune.args]] # extracts the best model - returns MaxEnt object


# M. coronaria predictions ------------------------------------------------
# Now use the <terra> package to plot the SDM prediction.
# Wclim is the historical climatic conditions (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
cor_pred_hist <- terra::predict(wclim, mod.best_cor_maxent, cores = cn - 1, na.rm = T)

plot(cor_pred_hist)
points(occThin_cor, cex = 0.05)

# Evaluate predictions using Boyce Index
# the number of true presences should decline with suitability groups 100-91, 90-81, etc. 
# First extract suitability values for the background and presence points, make sure to omit NA values
corPred_bg_val <- terra::extract(cor_pred_hist, bg_cor_coords)$lyr1 %>% 
  na.omit()

corPred_val_na <- terra::extract(cor_pred_hist, occ_cor_coords)$lyr1 %>% 
  na.omit()

# Evaluate predictions using Boyce Index
ecospat.boyce(fit = corPred_bg_val, # vector of predicted habitat suitability of bg points
                           obs = corPred_val_na, # vector of 
                           nclass = 0, 
                           PEplot = TRUE,
                           method = 'spearman')

# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
corPred_val <- terra::extract(cor_pred_hist, occ_cor_coords)$lyr1
corPred_threshold_1 <- quantile(corPred_val, 0.01, na.rm = T) # Low suitability
corPred_threshold_10 <- quantile(corPred_val, 0.1, na.rm = T) # Moderate suitability
corPred_threshold_50 <- quantile(corPred_val, 0.5, na.rm = T) # High suitability

# Plotting the prediction
legend_labs <- c('Low Suitability', 'Moderate Suitability', 'High Suitability')

# brewer.pal(3, "YlOrBr") # Brewer palettes are helpful for mapping colour gradients
fill_cols <- c("#FFF7BC", "#FEC44F", "#D95F0E")
#fill_cols <- c('#D81B60', '#1E88E5', '#FFC107') # old colours

dev.off()
par(mar = c(5, 5, 5, 5))
terra::plot(cor_pred_hist > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-100, -50), ylim = c(30, 60), main = expression(atop(italic('Malus coronaria'), " Historical Suitability (1970-2000)")), background = 'lightskyblue1')
terra::plot(cor_pred_hist > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_hist > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T, cex = .1)
#points(occThin_cor, col = 'black', cex = 0.75, pch = 4)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)


# Future Climate predictions
# SSP 245
cor_pred_ssp245_30 <- terra::predict(ssp245_2030, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp245_50 <- terra::predict(ssp245_2050, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp245_70 <- terra::predict(ssp245_2070, mod.best_cor_maxent, cores = cn - 1, na.rm = T)

# SSP 585
cor_pred_ssp585_30 <- terra::predict(ssp585_2030, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_50 <- terra::predict(ssp585_2050, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_70 <- terra::predict(ssp585_2070, mod.best_cor_maxent, cores = cn - 1, na.rm = T)

# Plot SSP 585 2030
par(mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp585_30 > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-100, -50), ylim = c(30, 60), main = expression(atop(italic('Malus coronaria'), " SSP5-8.5 Prediction: Early Century (2020-2040)")), background = 'lightskyblue1')
terra::plot(cor_pred_ssp585_30 > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_30 > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

# Save/Load M cor. SDM predictions ----------------------------------------
setwd('../sdm_output')

# Save
saveRDS(cor_pred_hist, file = 'cor_pred_hist.Rdata')

saveRDS(cor_pred_ssp245_30, file = 'cor_pred_ssp245_30.Rdata')
saveRDS(cor_pred_ssp245_50, file = 'cor_pred_ssp245_50.Rdata')
saveRDS(cor_pred_ssp245_70, file = 'cor_pred_ssp245_70.Rdata')

saveRDS(cor_pred_ssp585_30, file = 'cor_pred_ssp585_30.Rdata')
saveRDS(cor_pred_ssp585_50, file = 'cor_pred_ssp585_50.Rdata')
saveRDS(cor_pred_ssp585_70, file = 'cor_pred_ssp585_70.Rdata')

# Load
cor_pred_hist <- readRDS(file = 'cor_pred_hist.Rdata')

cor_pred_ssp245_30 <- readRDS(file = 'cor_pred_ssp245_30.Rdata')
cor_pred_ssp245_50 <- readRDS(file = 'cor_pred_ssp245_50.Rdata')
cor_pred_ssp245_70 <- readRDS(file = 'cor_pred_ssp245_70.Rdata')

cor_pred_ssp585_30 <- readRDS(file = 'cor_pred_ssp585_30.Rdata')
cor_pred_ssp585_50 <- readRDS(file = 'cor_pred_ssp585_50.Rdata')
cor_pred_ssp585_70 <- readRDS(file = 'cor_pred_ssp585_70.Rdata')



# Save coronaria thresholds -----------------------------------------------
setwd('../sdm_output/thresholds')
saveRDS(corPred_threshold_1, file = 'corPred_threshold_1.Rdata')
saveRDS(corPred_threshold_10, file = 'corPred_threshold_10.Rdata')
saveRDS(corPred_threshold_50, file = 'corPred_threshold_50.Rdata')

corPred_threshold_1 <- readRDS(file = 'corPred_threshold_1.Rdata')
corPred_threshold_10 <- readRDS(file = 'corPred_threshold_10.Rdata')
corPred_threshold_50 <- readRDS(file = 'corPred_threshold_50.Rdata')


# M coronaria Habitat predictions -----------------------------------------
# Categorical habitat suitability
# Historical
cor_pred_high_hist <- cor_pred_hist > corPred_threshold_50
cor_pred_mod_hist <- cor_pred_hist > corPred_threshold_10
cor_pred_low_hist <- cor_pred_hist > corPred_threshold_1

#SSP245 
cor_pred_high_ssp245_30 <- cor_pred_ssp245_30 > corPred_threshold_50
cor_pred_mod_ssp245_30 <- cor_pred_ssp245_30 > corPred_threshold_10
cor_pred_low_ssp245_30 <- cor_pred_ssp245_30 > corPred_threshold_1

cor_pred_high_ssp245_50 <- cor_pred_ssp245_50 > corPred_threshold_50
cor_pred_mod_ssp245_50 <- cor_pred_ssp245_50 > corPred_threshold_10
cor_pred_low_ssp245_50 <- cor_pred_ssp245_50 > corPred_threshold_1

cor_pred_high_ssp245_70 <- cor_pred_ssp245_70 > corPred_threshold_50
cor_pred_mod_ssp245_70 <- cor_pred_ssp245_70 > corPred_threshold_10
cor_pred_low_ssp245_70 <- cor_pred_ssp245_70 > corPred_threshold_1

#SSP585
cor_pred_high_ssp585_30 <- cor_pred_ssp585_30 > corPred_threshold_50
cor_pred_mod_ssp585_30 <- cor_pred_ssp585_30 > corPred_threshold_10
cor_pred_low_ssp585_30 <- cor_pred_ssp585_30 > corPred_threshold_1

cor_pred_high_ssp585_50 <- cor_pred_ssp585_50 > corPred_threshold_50
cor_pred_mod_ssp585_50 <- cor_pred_ssp585_50 > corPred_threshold_10
cor_pred_low_ssp585_50 <- cor_pred_ssp585_50 > corPred_threshold_1

cor_pred_high_ssp585_70 <- cor_pred_ssp585_70 > corPred_threshold_50
cor_pred_mod_ssp585_70 <- cor_pred_ssp585_70 > corPred_threshold_10
cor_pred_low_ssp585_70 <- cor_pred_ssp585_70 > corPred_threshold_1

# Save
getwd()
setwd('../sdm_output/habitat_predictions/high_moderate_low_predictions')

# Historical
saveRDS(cor_pred_high_hist, file = 'cor_pred_high_hist.Rdata')
saveRDS(cor_pred_mod_hist, file = 'cor_pred_mod_hist.Rdata')
saveRDS(cor_pred_low_hist, file = 'cor_pred_low_hist.Rdata')

# SSP245
saveRDS(cor_pred_high_ssp245_30, file = 'cor_pred_high_ssp245_30.Rdata')
saveRDS(cor_pred_mod_ssp245_30, file = 'cor_pred_mod_ssp245_30.Rdata')
saveRDS(cor_pred_low_ssp245_30, file = 'cor_pred_low_ssp245_30.Rdata')

saveRDS(cor_pred_high_ssp245_50, file = 'cor_pred_high_ssp245_50.Rdata')
saveRDS(cor_pred_mod_ssp245_50, file = 'cor_pred_mod_ssp245_50.Rdata')
saveRDS(cor_pred_low_ssp245_50, file = 'cor_pred_low_ssp245_50.Rdata')

saveRDS(cor_pred_high_ssp245_70, file = 'cor_pred_high_ssp245_70.Rdata')
saveRDS(cor_pred_mod_ssp245_70, file = 'cor_pred_mod_ssp245_70.Rdata')
saveRDS(cor_pred_low_ssp245_70, file = 'cor_pred_low_ssp245_70.Rdata')

# SSP585
saveRDS(cor_pred_high_ssp585_30, file = 'cor_pred_high_ssp585_30.Rdata')
saveRDS(cor_pred_mod_ssp585_30, file = 'cor_pred_mod_ssp585_30.Rdata')
saveRDS(cor_pred_low_ssp585_30, file = 'cor_pred_low_ssp585_30.Rdata')

saveRDS(cor_pred_high_ssp585_50, file = 'cor_pred_high_ssp585_50.Rdata')
saveRDS(cor_pred_mod_ssp585_50, file = 'cor_pred_mod_ssp585_50.Rdata')
saveRDS(cor_pred_low_ssp585_50, file = 'cor_pred_low_ssp585_50.Rdata')

saveRDS(cor_pred_high_ssp585_70, file = 'cor_pred_high_ssp585_70.Rdata')
saveRDS(cor_pred_mod_ssp585_70, file = 'cor_pred_mod_ssp585_70.Rdata')
saveRDS(cor_pred_low_ssp585_70, file = 'cor_pred_low_ssp585_70.Rdata')

# Load



# Binary pred habitat GAP ANALYSIS ----------------------------------------


cor_pa <- predicts::pa_evaluate(p = occ_cor_coords_mat, a = bg_cor_coords_mat, model = cor_maxent, x = wclim_cor)
cor_threshold <- predicts::threshold(cor_pa)

cor_hist_habitat <- cor_pred_hist > cor_threshold$max_spec_sens #the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest


# Fusca - MaxEnt Model ----------------------------------------------------

# Build Species Distribution Model using MaxEnt from the <ENMeval> package

# Run prediction in a parallel using 'socket' clusters to help speed up computation
# <ENMeval> implements parallel functions natively
# But load <parallel> library for additional functions like <decectCores()>
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

fus_maxent <- ENMevaluate(occ_fus_coords, # occurrence records
                          envs = wclim_fus, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 8, 0.5),
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP", "LQHPT")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                          partitions = 'checkerboard2',
                          parallel = TRUE,
                          numCores = cn - 1, # leave one core available for other apps
                          parallelType = "doParallel", # use doParrallel on Windows - socket cluster  
                          algorithm = 'maxent.jar')

# Save the MaxEnt model so you do not have to waste time re-running the model
setwd('../sdm_output')
saveRDS(fus_maxent, file = 'fus_maxent.Rdata') # save
fus_maxent <- readRDS(file = 'fus_maxent.Rdata') # load 



# M. fusca Model Selection ------------------------------------------------
# Note that maxent results provide Continuous Boyce Index (cbi)
# Two models had a delta AIC < 2, rm.1_fc.LQHPT and rm.1.5_fc.LQHPT
best_fus_maxent <- subset(fus_maxent@results, delta.AICc < 2) # selects the best performing model based on delta AICc - returns data frame object
mod.best_fus_maxent <- eval.models(fus_maxent)[[best_fus_maxent$tune.args]] # extracts the best model - returns MaxEnt object


# M. fusca Predictions ----------------------------------------------------
# Now use the <terra> package to plot the SDM prediction.
# Wclim is the historical climatic conditions (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
fus_pred_hist <- terra::predict(wclim, mod.best_fus_maxent, cores = cn - 1, na.rm = T)

plot(fus_pred_hist)
points(occThin_fus, cex = 0.05)

# Evaluate predictions using Boyce Index
# the number of true presences should decline with suitability groups 100-91, 90-81, etc. 
# First extract suitability values for the background and presence points, make sure to omit NA values
fusPred_bg_val <- terra::extract(fus_pred_hist, bg_fus_coords)$lyr1 %>% 
  na.omit()

fusPred_val_na <- terra::extract(fus_pred_hist, occ_fus_coords)$lyr1 %>% 
  na.omit()

# Evaluate predictions using Boyce Index
ecospat.boyce(fit = fusPred_bg_val, # vector of predicted habitat suitability of bg points
              obs = fusPred_val_na, # vector of 
              nclass = 0, 
              PEplot = TRUE,
              method = 'spearman')

# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
fusPred_val <- terra::extract(fus_pred_hist, occ_fus_coords)$lyr1
fusPred_threshold_1 <- quantile(fusPred_val, 0.01, na.rm = T) # Low suitability
fusPred_threshold_10 <- quantile(fusPred_val, 0.1, na.rm = T) # Moderate suitability
fusPred_threshold_50 <- quantile(fusPred_val, 0.5, na.rm = T) # High suitability

legend_labs <- c('Low Suitability', 'Moderate Suitability', 'High Suitability')
fill_cols <- c("#FFF7BC", "#FEC44F", "#D95F0E")

par(mar = c(5, 5, 5, 5))
terra::plot(fus_pred_hist > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " Historical Suitability (1970-2000)")), background = 'lightskyblue1')
terra::plot(fus_pred_hist > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_hist > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
#points(occThin_fus, col = 'black', cex = 0.75, pch = 4)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)


# Future Climate predictions
# SSP 245
fus_pred_ssp245_30 <- terra::predict(ssp245_2030, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp245_50 <- terra::predict(ssp245_2050, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp245_70 <- terra::predict(ssp245_2070, mod.best_fus_maxent, cores = cn - 1, na.rm = T)

# SSP 585
fus_pred_ssp585_30 <- terra::predict(ssp585_2030, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp585_50 <- terra::predict(ssp585_2050, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp585_70 <- terra::predict(ssp585_2070, mod.best_fus_maxent, cores = cn - 1, na.rm = T)

# Plot SSP 585 2030
par(mar = c(5, 5, 5, 5))
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_1, col = c('lightgrey', '#D81B60'), legend = F,  xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus coronaria'), " SSP5-8.5 Prediction: Early Century (2020-2040)")), background = 'lightskyblue1')
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_10, col = c(NA, '#1E88E5'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_50, col = c(NA, '#FFC107'), add = T, legend = F)
terra::plot(canUSMex_map, add = T)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)


# Save/Load M fus. SDM predictions ----------------------------------------
setwd('../sdm_output/../../')

# Save
saveRDS(fus_pred_hist, file = 'fus_pred_hist.Rdata')

saveRDS(fus_pred_ssp245_30, file = 'fus_pred_ssp245_30.Rdata')
saveRDS(fus_pred_ssp245_50, file = 'fus_pred_ssp245_50.Rdata')
saveRDS(fus_pred_ssp245_70, file = 'fus_pred_ssp245_70.Rdata')

saveRDS(fus_pred_ssp585_30, file = 'fus_pred_ssp585_30.Rdata')
saveRDS(fus_pred_ssp585_50, file = 'fus_pred_ssp585_50.Rdata')
saveRDS(fus_pred_ssp585_70, file = 'fus_pred_ssp585_70.Rdata')

# Load
fus_pred_hist <- readRDS(file = 'fus_pred_hist.Rdata')

fus_pred_ssp245_30 <- readRDS(file = 'fus_pred_ssp245_30.Rdata')
fus_pred_ssp245_50 <- readRDS(file = 'fus_pred_ssp245_50.Rdata')
fus_pred_ssp245_70 <- readRDS(file = 'fus_pred_ssp245_70.Rdata')

fus_pred_ssp585_30 <- readRDS(file = 'fus_pred_ssp585_30.Rdata')
fus_pred_ssp585_50 <- readRDS(file = 'fus_pred_ssp585_50.Rdata')
fus_pred_ssp585_70 <- readRDS(file = 'fus_pred_ssp585_70.Rdata')


# Save M. fusca thresholds ------------------------------------------------
setwd('../sdm_output/thresholds')
saveRDS(fusPred_threshold_1, file = 'fusPred_threshold_1.Rdata')
saveRDS(fusPred_threshold_10, file = 'fusPred_threshold_10.Rdata')
saveRDS(fusPred_threshold_50, file = 'fusPred_threshold_50.Rdata')

# M fusca Habitat predictions ---------------------------------------------
# Categorical habitat suitability
# Historical
fus_pred_high_hist <- fus_pred_hist > fusPred_threshold_50
fus_pred_mod_hist <- fus_pred_hist > fusPred_threshold_10
fus_pred_low_hist <- fus_pred_hist > fusPred_threshold_1

#SSP245 
fus_pred_high_ssp245_30 <- fus_pred_ssp245_30 > fusPred_threshold_50
fus_pred_mod_ssp245_30 <- fus_pred_ssp245_30 > fusPred_threshold_10
fus_pred_low_ssp245_30 <- fus_pred_ssp245_30 > fusPred_threshold_1

fus_pred_high_ssp245_50 <- fus_pred_ssp245_50 > fusPred_threshold_50
fus_pred_mod_ssp245_50 <- fus_pred_ssp245_50 > fusPred_threshold_10
fus_pred_low_ssp245_50 <- fus_pred_ssp245_50 > fusPred_threshold_1

fus_pred_high_ssp245_70 <- fus_pred_ssp245_70 > fusPred_threshold_50
fus_pred_mod_ssp245_70 <- fus_pred_ssp245_70 > fusPred_threshold_10
fus_pred_low_ssp245_70 <- fus_pred_ssp245_70 > fusPred_threshold_1

#SSP585
fus_pred_high_ssp585_30 <- fus_pred_ssp585_30 > fusPred_threshold_50
fus_pred_mod_ssp585_30 <- fus_pred_ssp585_30 > fusPred_threshold_10
fus_pred_low_ssp585_30 <- fus_pred_ssp585_30 > fusPred_threshold_1

fus_pred_high_ssp585_50 <- fus_pred_ssp585_50 > fusPred_threshold_50
fus_pred_mod_ssp585_50 <- fus_pred_ssp585_50 > fusPred_threshold_10
fus_pred_low_ssp585_50 <- fus_pred_ssp585_50 > fusPred_threshold_1

fus_pred_high_ssp585_70 <- fus_pred_ssp585_70 > fusPred_threshold_50
fus_pred_mod_ssp585_70 <- fus_pred_ssp585_70 > fusPred_threshold_10
fus_pred_low_ssp585_70 <- fus_pred_ssp585_70 > fusPred_threshold_1

# Save
getwd()
setwd('../sdm_output/habitat_predictions/high_moderate_low_predictions')

# Historical
saveRDS(fus_pred_high_hist, file = 'fus_pred_high_hist.Rdata')
saveRDS(fus_pred_mod_hist, file = 'fus_pred_mod_hist.Rdata')
saveRDS(fus_pred_low_hist, file = 'fus_pred_low_hist.Rdata')

# SSP245
saveRDS(fus_pred_high_ssp245_30, file = 'fus_pred_high_ssp245_30.Rdata')
saveRDS(fus_pred_mod_ssp245_30, file = 'fus_pred_mod_ssp245_30.Rdata')
saveRDS(fus_pred_low_ssp245_30, file = 'fus_pred_low_ssp245_30.Rdata')

saveRDS(fus_pred_high_ssp245_50, file = 'fus_pred_high_ssp245_50.Rdata')
saveRDS(fus_pred_mod_ssp245_50, file = 'fus_pred_mod_ssp245_50.Rdata')
saveRDS(fus_pred_low_ssp245_50, file = 'fus_pred_low_ssp245_50.Rdata')

saveRDS(fus_pred_high_ssp245_70, file = 'fus_pred_high_ssp245_70.Rdata')
saveRDS(fus_pred_mod_ssp245_70, file = 'fus_pred_mod_ssp245_70.Rdata')
saveRDS(fus_pred_low_ssp245_70, file = 'fus_pred_low_ssp245_70.Rdata')

# SSP585
saveRDS(fus_pred_high_ssp585_30, file = 'fus_pred_high_ssp585_30.Rdata')
saveRDS(fus_pred_mod_ssp585_30, file = 'fus_pred_mod_ssp585_30.Rdata')
saveRDS(fus_pred_low_ssp585_30, file = 'fus_pred_low_ssp585_30.Rdata')

saveRDS(fus_pred_high_ssp585_50, file = 'fus_pred_high_ssp585_50.Rdata')
saveRDS(fus_pred_mod_ssp585_50, file = 'fus_pred_mod_ssp585_50.Rdata')
saveRDS(fus_pred_low_ssp585_50, file = 'fus_pred_low_ssp585_50.Rdata')

saveRDS(fus_pred_high_ssp585_70, file = 'fus_pred_high_ssp585_70.Rdata')
saveRDS(fus_pred_mod_ssp585_70, file = 'fus_pred_mod_ssp585_70.Rdata')
saveRDS(fus_pred_low_ssp585_70, file = 'fus_pred_low_ssp585_70.Rdata')

# Load
