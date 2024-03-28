# Top ---------------------------------------------------------------------
# MaxEnt Species Distribution Modeling (SDM) for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024

library(tidyverse) # Grammar and data management
library(terra)# Spatial Data package
library(predicts) # SDM package
library(geodata) # basemaps
library(rJava) # MaxEnt models are dependant on JDK
library(ENMeval) # Another modeling package, useful for data partitioning (Checkerboarding)
library(raster) # RasterStack dependancy (a now deprecated function)
library(ecospat) # Useful spatial ecology tools


# Load occurrence data and basemaps -------------------------------------------------------
getwd() # check you directory location

# Background points in SpatVectors
setwd("../occ_data/")
cor_bg_vec <- readRDS(file = 'cor_bg_vec.Rdata')
fus_bg_vec <- readRDS(file = 'fus_bg_vec.Rdata')

# Occurrence Points in SpatVectors
setwd("../occ_data/")
occThin_cor <- readRDS(file = 'occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = 'occThin_fus.Rdata') # M. fusca

# Download/load basemaps
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "../occ_data/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = '../occ_data/base_maps') # Canada basemap w. Provinces

mex_map <-gadm(country = 'MX', level = 1, resolution = 2,
               path = '../occ_data/base_maps') # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map


NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in Western Hemisphere

canUSMex_map <- crop(canUSMex_map, NA_ext) # crop to Western Hemisphere

plot(canUSMex_map) # plot basemap

# Load cropped climate Rasters --------------------------------------------
# These Rasters are useful for sampling spatial checkerboards 
# and making habitat suitability predictions (Historical and under future SSPs climate scenarios)

setwd('../wclim_data')
# Historical (1970-2000)
wclim_cor <- readRDS(file = 'wclim_cor.Rdata')
wclim_cor_stack <- raster::stack(wclim_cor) # covert SpatRaster to RasterStack for depandacy in ENMeval checkboarding

wclim_fus <- readRDS(file = 'wclim_fus.Rdata')
wclim_fus_stack <- raster::stack(wclim_fus) # covert SpatRaster to RasterStack for depandacy in ENMeval checkboarding

climate_predictors <- names(wclim_cor) # extract climate predictor names, to rename layers in the rasters below
# This is important to do for making predictions once the SDMs have been made on future climate data
# Note that the names of the layers still correspond to the same enviromental variables

# Future SSPs
# Do not need to create RasterStacks
# SSP 245
cor_ssp245_2030 <- readRDS(file = 'cor_ssp245_2030.Rdata')
names(cor_ssp245_2030) <- climate_predictors #rename raster layers for downsteam anaylsis
fus_ssp245_2030 <- readRDS(file = 'fus_ssp245_2030.Rdata')
names(fus_ssp245_2030) <- climate_predictors 

cor_ssp245_2050 <- readRDS(file = 'cor_ssp245_2050.Rdata')
names(cor_ssp245_2050) <- climate_predictors
fus_ssp245_2050 <- readRDS(file = 'fus_ssp245_2050.Rdata')
names(fus_ssp245_2050) <- climate_predictors

cor_ssp245_2070 <- readRDS(file = 'cor_ssp245_2070.Rdata')
names(cor_ssp245_2070) <- climate_predictors
fus_ssp245_2070 <- readRDS(file = 'fus_ssp245_2070.Rdata')
names(fus_ssp245_2070) <- climate_predictors

# SSP 585
cor_ssp585_2030 <- readRDS(file = 'cor_ssp585_2030.Rdata')
names(cor_ssp585_2030) <- climate_predictors
fus_ssp585_2030 <- readRDS(file = 'fus_ssp585_2030.Rdata')
names(fus_ssp585_2030) <- climate_predictors

cor_ssp585_2050 <- readRDS(file = 'cor_ssp585_2050.Rdata')
names(cor_ssp585_2050) <- climate_predictors
fus_ssp585_2050 <- readRDS(file = 'fus_ssp585_2050.Rdata')
names(fus_ssp585_2050) <- climate_predictors

cor_ssp585_2070 <- readRDS(file = 'cor_ssp585_2070.Rdata')
names(cor_ssp585_2070) <- climate_predictors
fus_ssp585_2070 <- readRDS(file = 'fus_ssp585_2070.Rdata')
names(fus_ssp585_2070) <- climate_predictors



# Sub-sample Test and Training Occurrences --------------------------------
# Simple random paritioning of occurrance datapoints into k-folds ('k' number of groups)
# This does not always capture properly spattially replicated variation in the data 
# See spatial partitioning bellow using spatial checker boarding
# M. coronaria
folds_cor <- folds(occThin_cor, k = 5)
test_cor <- occThin_cor[folds_cor == 1, ]
train_cor <- occThin_cor[folds_cor != 1, ]

#take a look at the spatial distribution
plot(canUSMex_map)
points(train_cor, col = 'blue')
points(test_cor, col = 'red')


# M. fusca
folds_fus <- folds(occThin_fus, k = 5)
test_fus <- occThin_fus[folds_fus == 1, ]
train_fus <- occThin_fus[folds_fus != 1, ]


#take a look at the spatial distribution
plot(canUSMex_map)
points(train_fus, col = 'blue')
points(test_fus, col = 'red')



# Spatial partitioning preparation ----------------------------------------

occ_cor_coords <- as.data.frame(geom(occThin_cor)[,3:4]) # extract longitude, lattitude from occurence points
bg_cor_coords <- as.data.frame(geom(cor_bg_vec)[,3:4]) # extract longitude, lattitude from background points

occ_fus_coords <- as.data.frame(geom(occThin_fus)[,3:4]) # extract longitude, lattitude from occurence points
bg_fus_coords <- as.data.frame(geom(fus_bg_vec)[,3:4]) # extract longitude, lattitude from background points


# Spatial Checkerboard sampling -------------------------------------------
# Spatial Checkerboard partitioning of occurrence and background points cross-fold validation
set.seed(1337)

agg_factor <- c(9,9) # this defines how many adjacent cells are aggregated into the same checkboard 'sqaure'
# I visually tested different aggregations and found 9,9 was the best fit given the spatial extent


cor_cb <- get.checkerboard2(occs = occ_cor_coords, 
                            bg = bg_cor_coords, 
                            envs = wclim_cor_stack, 
                            aggregation.factor = agg_factor
                            )

evalplot.grps(pts = occ_cor_coords, pts.grp = cor_cb$occs.grp, envs = wclim_cor_stack, pts.size = .75) # plot the checkerboard partitions of occurrence points
evalplot.grps(pts = bg_cor_coords, pts.grp = cor_cb$bg.grp, envs = wclim_cor_stack, pts.size = .75) # plot the checkerboard partitions of occurrence bg points
# background points may return an error if the checkboard sampling, this can happen due to the aggregation factor removing some obs
# It is fine as we do not need to visualize groups. If you want to, must remove the difference in observations.

# Separate testing and training data
cor_occ_cbGroup <- cor_cb$occs.grp # extract group patitions into a vector
cor_occ_train <- occThin_cor[cor_occ_cbGroup != 1, ] # 3 groups used for training the model
cor_occ_test <- occThin_cor[cor_occ_cbGroup == 1, ] # 1 group withheld for testing the model


fus_cb <- get.checkerboard2(occs = occ_fus_coords, 
                            bg = bg_fus_coords, 
                            envs = wclim_fus_stack, 
                            aggregation.factor = agg_factor,
                            )

evalplot.grps(pts = occ_fus_coords, pts.grp = fus_cb$occs.grp, envs = wclim_fus_stack, pts.size = .75) # plot the checkerboard partitions of occurrence points
evalplot.grps(pts = bg_fus_coords, pts.grp = fus_cb$bg.grp, envs = wclim_fus_stack, pts.size = .75) # plot the checkerboard partitions of occurrence bg points

# Seperate testing and training data
fus_cbGroup <- fus_cb$occs.grp # extract group patitions into a vector
fus_occ_train <- occThin_fus[fus_cbGroup != 1, ] # 3 groups used for training the model
fus_occ_test <- occThin_fus[fus_cbGroup == 1, ] # 1 group withheld for testing the model


# Model from Predictors in Raster -------------------------------------------------------
# Build Species Distribution Model using MaxEnt from the <predicts> package

# M. coronaria
cor_trainModel <- MaxEnt(x = wclim_cor, p = cor_occ_train, a = cor_bg_vec) # build a training model
cor_testModel <- MaxEnt(x = wclim_cor, p = cor_occ_test, a = cor_bg_vec) #build a testing model
# Compare AUC values between the training and test models to see how well they both perform.  Ideally have similar AUC values




# M. fusca
fus_trainModel <- MaxEnt(x = wclim_fus, p = fus_occ_train, a = fus_bg_vec) # build a training model
fus_testModel <- MaxEnt(x = wclim_fus, p = fus_occ_test, a = fus_bg_vec) #build a testing model



# Habitat Suitability Predictions (Histocal and Future Projections) -------
# Run prediction in a parallel 'socket' cluster to help speed up computation
# Terra implements parallel functions natively
# But load <parallel> library for additional functions like <decectCores()>
library(parallel)
cn <- detectCores(logical = F) # number of physical RAM cores in your computure

# M. coronaria predictions ------------------------------------------------

# 'Historical Habitat Suitability'
cor_pred_hist <- terra::predict(wclim_cor, cor_trainModel, cores = cn - 1) # leave one core for other processes, browser etc.
plot(cor_pred_hist, main = 'Malus coronaria (Historical)') #Plot 
points(occThin_cor, cex = 0.5) # Overlay thinned occurrences

# Boyce Index, useful for assessing model performance
#NOT WORKING?
ecospat.boyce(fit = cor_pred_hist, obs = occ_cor_coords_mat)

# SSP 245
cor_pred_ssp245_30 <- terra::predict(cor_ssp245_2030, cor_trainModel, cores = cn - 1)
cor_pred_ssp245_50 <- terra::predict(cor_ssp245_2050, cor_trainModel, cores = cn - 1)
cor_pred_ssp245_70 <- terra::predict(cor_ssp245_2070, cor_trainModel, cores = cn - 1)

# SSP585
cor_pred_ssp585_30 <- terra::predict(cor_ssp585_2030, cor_trainModel, cores = cn - 1)
cor_pred_ssp585_50 <- terra::predict(cor_ssp585_2050, cor_trainModel, cores = cn - 1)
cor_pred_ssp585_70 <- terra::predict(cor_ssp585_2070, cor_trainModel, cores = cn - 1)

# Save/Load M cor. SDM predictions ----------------------------------------
setwd('../sdm_output')
saveRDS(cor_pred_ssp245_30, file = 'cor_pred_ssp245_30.Rdata')
saveRDS(cor_pred_ssp245_50, file = 'cor_pred_ssp245_50.Rdata')
saveRDS(cor_pred_ssp245_70, file = 'cor_pred_ssp245_70.Rdata')

saveRDS(cor_pred_ssp585_30, file = 'cor_pred_ssp585_30.Rdata')
saveRDS(cor_pred_ssp585_50, file = 'cor_pred_ssp585_50.Rdata')
saveRDS(cor_pred_ssp585_70, file = 'cor_pred_ssp585_70.Rdata')

# Plot M. coronaria historical, early/mid/late century projections of habitat suitability
# Example: SSP245
cor_par <- par(mfrow = c(2,2), mar=c(3,3,1,1), oma=c(0,0,3,1))  # oma creates space
plot(cor_pred_hist, main = 'Historical (1970-2000)')
plot(cor_pred_ssp245_30, main = 'Early-Century (2021-2040)')
plot(cor_pred_ssp245_50, main = 'Mid-Century (2041-2060)')
plot(cor_pred_ssp245_70, main = 'Late-Century (2061-2080)')
mtext("Malus coronaria Probability of Habitat Suitability", outer = T)


# M. coronaria model eval -------------------------------------------------
# for evaluating the p/bg performance of models
# pa_evaluate expects occurrences in matric form
occ_cor_coords_mat <- as.matrix(occ_cor_coords)
bg_cor_coords_mat <- as.matrix(bg_cor_coords)

# Note make sure to change the predicted climate Raster for each model evaluation
cor_pa <- predicts::pa_evaluate(p = occ_cor_coords_mat, a = bg_cor_coords_mat, model = cor_trainModel, x = wclim_cor)
cor_threshold <- predicts::threshold(cor_pa)

# M. coronaria Habititat Suitability Maps ---------------------------------
#Historical and predicted maps
cor_hist_habitat <- cor_pred_hist > cor_threshold$max_spec_sens #the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest
plot(cor_hist_habitat, col = c('white', 'blue'))
plot(canUSMex_map, add = T)
points(occThin_cor, cex = 0.5, pch = 4, col = 'red')

cor_ssp245_2030_habitat <- cor_pred_ssp245_30 > cor_threshold$max_spec_sens
cor_ssp245_2050_habitat <- cor_pred_ssp245_50 > cor_threshold$max_spec_sens
cor_ssp245_2070_habitat <- cor_pred_ssp245_70 > cor_threshold$max_spec_sens

cor_ssp585_2030_habitat <- cor_pred_ssp585_30 > cor_threshold$max_spec_sens
cor_ssp585_2050_habitat <- cor_pred_ssp585_50 > cor_threshold$max_spec_sens
cor_ssp585_2070_habitat <- cor_pred_ssp585_70 > cor_threshold$max_spec_sens


# Save M. cor. Predicted Suitable Habitat  --------------------------------
setwd('../sdm_output')
saveRDS(cor_ssp245_2030_habitat, file = 'cor_ssp245_2030_habitat.Rdata')
saveRDS(cor_ssp245_2050_habitat, file = 'cor_ssp245_2050_habitat.Rdata')
saveRDS(cor_ssp245_2070_habitat, file = 'cor_ssp245_2070_habitat.Rdata')

saveRDS(cor_ssp585_2030_habitat, file = 'cor_ssp585_2030_habitat.Rdata')
saveRDS(cor_ssp585_2030_habitat, file = 'cor_ssp585_2050_habitat.Rdata')
saveRDS(cor_ssp585_2030_habitat, file = 'cor_ssp585_2070_habitat.Rdata')

# Plot predicted Habitat suitability 
# SSP245
cor_par <- par(mfrow = c(2,2), mar=c(3,3,1,1), oma=c(0,0,3,1))  # oma creates space
plot(cor_hist_habitat, main = 'Historical (1970-2000)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(cor_ssp245_2030_habitat, main = 'Early-Century (2021-2040)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(cor_ssp245_2050_habitat, main = 'Mid-Century (2041-2060)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(cor_ssp245_2070_habitat, main = 'Late-Century (2061-2080)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
mtext(text = expression(paste("SSP 2-4.5 ", italic("Malus coronaria "),"Predicted Suitable Habitat")), outer = T, cex = 2)

#SSP585
# Plot predicted Habitat suitability 
cor_par <- par(mfrow = c(2,2), mar=c(3,3,1,1), oma=c(0,0,3,1))  # oma creates space
plot(cor_hist_habitat, main = 'Historical (1970-2000)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(cor_ssp585_2030_habitat, main = 'Early-Century (2021-2040)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(cor_ssp585_2050_habitat, main = 'Mid-Century (2041-2060)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(cor_ssp585_2070_habitat, main = 'Late-Century (2061-2080)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
mtext(text = expression(paste("SSP 5- 8.5 ", italic("Malus coronaria "),"Predicted Suitable Habitat")), outer = T, cex = 2)



# M. fusca predictions ----------------------------------------------------
# M. fusca
fus_pred_hist <- terra::predict(wclim_fus, fus_trainModel, cores = cn - 1) # leave one core for other processes, browser etc.
#plot
plot(fus_pred, main = 'Malus fusca (Historical)')
points(occThin_fus, cex = 0.5) # Overlay thinned occurrences

# SSP 245
fus_pred_ssp245_30 <- terra::predict(fus_ssp245_2030, fus_trainModel, cores = cn - 1)
fus_pred_ssp245_50 <- terra::predict(fus_ssp245_2050, fus_trainModel, cores = cn - 1)
fus_pred_ssp245_70 <- terra::predict(fus_ssp245_2070, fus_trainModel, cores = cn - 1)

# SSP585
fus_pred_ssp585_30 <- terra::predict(fus_ssp585_2030, fus_trainModel, cores = cn - 1)
fus_pred_ssp585_50 <- terra::predict(fus_ssp585_2050, fus_trainModel, cores = cn - 1)
fus_pred_ssp585_70 <- terra::predict(fus_ssp585_2070, fus_trainModel, cores = cn - 1)


# Save/Load M fus. SDM predictions ----------------------------------------
setwd('../sdm_output')
saveRDS(fus_pred_ssp245_30, file = 'fus_pred_ssp245_30.Rdata')
saveRDS(fus_pred_ssp245_50, file = 'fus_pred_ssp245_50.Rdata')
saveRDS(fus_pred_ssp245_70, file = 'fus_pred_ssp245_70.Rdata')

saveRDS(fus_pred_ssp585_30, file = 'fus_pred_ssp585_30.Rdata')
saveRDS(fus_pred_ssp585_50, file = 'fus_pred_ssp585_50.Rdata')
saveRDS(fus_pred_ssp585_70, file = 'fus_pred_ssp585_70.Rdata')

# Plot M. fusca historical, early/mid/late century projections of habitat suitability
# Example: SSP245
fus_par <- par(mfrow = c(2,2), mar=c(3,3,1,1), oma=c(0,0,3,1))  # oma creates space
plot(fus_pred_hist, main = 'Historical (1970-2000)')
plot(fus_pred_ssp245_30, main = 'Early-Century (2021-2040)')
plot(fus_pred_ssp245_50, main = 'Mid-Century (2041-2060)')
plot(fus_pred_ssp245_70, main = 'Late-Century (2061-2080)')
mtext("Malus coronaria Probability of Habitat Suitability", outer = T)

# M. fusca model eval -----------------------------------------------------
# for evaluating the p/bg performance of models
# pa_evaluate expects occurrences in matric form
occ_fus_coords_mat <- as.matrix(occ_fus_coords)
bg_fus_coords_mat <- as.matrix(bg_fus_coords)

# Note make sure to change the predicted climate Raster for each model evaluation
fus_pa <- predicts::pa_evaluate(p = occ_fus_coords_mat, a = bg_fus_coords_mat, model = fus_trainModel, x = wclim_fus)
fus_threshold <- predicts::threshold(fus_pa)


# M. fusca Habititat Suitability Maps -------------------------------------
#Historical and predicted maps
fus_hist_habitat <- fus_pred_hist > fus_threshold$max_spec_sens #the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest
plot(fus_hist_habitat, col = c('white', 'blue'))
plot(canUSMex_map, add = T)
points(occThin_cor, cex = 0.5, pch = 4, col = 'red')

fus_ssp245_2030_habitat <- fus_pred_ssp245_30 > fus_threshold$max_spec_sens
fus_ssp245_2050_habitat <- fus_pred_ssp245_50 > fus_threshold$max_spec_sens
fus_ssp245_2070_habitat <- fus_pred_ssp245_70 > fus_threshold$max_spec_sens

fus_ssp585_2030_habitat <- fus_pred_ssp585_30 > fus_threshold$max_spec_sens
fus_ssp585_2050_habitat <- fus_pred_ssp585_50 > fus_threshold$max_spec_sens
fus_ssp585_2070_habitat <- fus_pred_ssp585_70 > fus_threshold$max_spec_sens


# Save M. fus. Predicted Suitable Habitat  --------------------------------
setwd('../sdm_output')
saveRDS(fus_ssp245_2030_habitat, file = 'fus_ssp245_2030_habitat.Rdata')
saveRDS(fus_ssp245_2050_habitat, file = 'fus_ssp245_2050_habitat.Rdata')
saveRDS(fus_ssp245_2070_habitat, file = 'fus_ssp245_2070_habitat.Rdata')

saveRDS(fus_ssp585_2030_habitat, file = 'fus_ssp585_2030_habitat.Rdata')
saveRDS(fus_ssp585_2030_habitat, file = 'fus_ssp585_2050_habitat.Rdata')
saveRDS(fus_ssp585_2030_habitat, file = 'fus_ssp585_2070_habitat.Rdata')

# Plot predicted Habitat suitability

# SSP 245
fus_par <- par(mfrow = c(2,2), mar=c(3,3,1,1), oma=c(0,0,3,1))  # oma creates space
plot(fus_hist_habitat, main = 'Historical (1970-2000)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(fus_ssp245_2030_habitat, main = 'Early-Century (2021-2040)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(fus_ssp245_2050_habitat, main = 'Mid-Century (2041-2060)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(fus_ssp245_2070_habitat, main = 'Late-Century (2061-2080)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
mtext(text = expression(paste('SSP2-4.5',italic(" Malus fusca "),"Predicted Suitable Habitat")), outer = T, cex = 1.5)

# SSP 585
fus_par <- par(mfrow = c(2,2), mar=c(3,3,1,1), oma=c(0,0,3,1))  # oma creates space
plot(fus_hist_habitat, main = 'Historical (1970-2000)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(fus_ssp585_2030_habitat, main = 'Early-Century (2021-2040)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(fus_ssp585_2050_habitat, main = 'Mid-Century (2041-2060)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
plot(fus_ssp585_2070_habitat, main = 'Late-Century (2061-2080)', col = c('white', 'blue'))
plot(canUSMex_map, add = T)
mtext(text = expression(paste('SSP5-8.5',italic(" Malus fusca "),"Predicted Suitable Habitat")), outer = T, cex = 1.5)


