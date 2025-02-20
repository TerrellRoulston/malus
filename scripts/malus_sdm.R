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

# Background points in SpatVectors
cor_bg_vec <- readRDS(file = './occ_data/cor_bg_vec.Rdata')
fus_bg_vec <- readRDS(file = './occ_data/fus_bg_vec.Rdata')

# Occurrence Points in SpatVectors
occThin_cor <- readRDS(file = './occ_data/occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = './occ_data/occThin_fus.Rdata') # M. fusca

# Great Lakes shapefiles for making pretty maps and cropping
great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')

NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in Western Hemisphere

# Download/load WorldClim data under future climate scenarios -------------
# WARNING DO NOT PUSH WORLDCLIM DATA
# Historical climate 1970-2000
wclim <- geodata::worldclim_global(var = 'bio',
                                   res = 2.5, 
                                   version = '2.1', 
                                  path = "./wclim_data/") %>% 
  terra::crop(NA_ext)  %>% #crop raster to NA 
  terra::mask(great_lakes, inverse = T) # cut out the great lakes

# SSP (Shared social-economic pathway) 2.45 
# middle of the road projection, high climate adaptation, low climate mitigation
ssp245_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp245_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp245_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

# SPP 5.85 
# low regard for enviromental sustainability, increased fossil fuel reliance, this is the current tracking projection
ssp585_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp585_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

ssp585_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/")%>% 
  crop(NA_ext) %>% #crop raster to NA 
  mask(great_lakes, inverse = T) # cut out the great lakes

# Load cropped climate Rasters --------------------------------------------
# These Rasters are useful for sampling spatial checkerboards 
# and making habitat suitability predictions (Historical and under future SSPs climate scenarios)
# Historical (1970-2000)
wclim_cor <- readRDS(file = './wclim_data/wclim_cor.Rdata') 
#wclim_cor_stack <- raster::stack(wclim_cor) # covert SpatRaster to RasterStack for dependency in ENMeval checkboarding

wclim_fus <- readRDS(file = './wclim_data/wclim_fus.Rdata')
#wclim_fus_stack <- raster::stack(wclim_fus) # covert SpatRaster to RasterStack for dependency in ENMeval checkboarding

climate_predictors <- names(wclim_cor) # extract climate predictor names, to rename layers in the rasters below
# This is important to do for making predictions once the SDMs have been made on future climate data
# Note that the names of the layers still correspond to the same environmental variables

# Future SSPs
# Do not need to create RasterStacks
# SSP 245
names(wclim) <- climate_predictors
names(ssp245_2030) <- climate_predictors #rename raster layers for downsteam analysis
names(ssp245_2050) <- climate_predictors 
names(ssp245_2070) <- climate_predictors 

# SSP 585
names(ssp585_2030) <- climate_predictors #rename raster layers for downsteam analysis
names(ssp585_2050) <- climate_predictors 
names(ssp585_2070) <- climate_predictors 


# Subset climate variables for SDM analysis -------------------------------
wclim_subs <- wclim %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
ssp245_2030_subs <- ssp245_2030 %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
ssp245_2050_subs <- ssp245_2050 %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
ssp245_2070_subs <- ssp245_2070 %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))

ssp585_2030_subs <- ssp585_2030 %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
ssp585_2050_subs <- ssp585_2050 %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
ssp585_2070_subs <- ssp585_2070 %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))

wclim_cor_subs <- wclim_cor %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_fus_subs <- wclim_fus %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))

# Coronaria - MaxEnt Model  ------------------------------------------------

# Spatial partitioning preparation
occ_cor_coords <- as.data.frame(geom(occThin_cor)[,3:4]) # extract longitude, lattitude from occurence points
bg_cor_coords <- as.data.frame(geom(cor_bg_vec)[,3:4]) # extract longitude, lattitude from background points


# Build Species Distribution Model using MaxEnt from the <ENMeval> package

# Run prediction in a parallel using 'socket' clusters to help speed up computation
# <ENMeval> implements parallel functions natively
# But load <parallel> library for additional functions like <decectCores()>
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

# current version of maxent.jar =  v3.4.4

cor_maxent <- ENMevaluate(occ_cor_coords, # occurrence records
                            envs = wclim_cor_subs, # NOTE CHANGE THE ENVS inputed, environment from background training area
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

saveRDS(cor_maxent, file = './sdm_output/cor_maxent_subs.Rdata') # save
cor_maxent <- readRDS(file = './sdm_output/cor_maxent.Rdata') # load 

#subsetted model
cor_maxent <- readRDS(file = './sdm_output/cor_maxent.Rdata') # load 


# M. coronaria Model Selection --------------------------------------------
# Note that maxent results provide Continuous Boyce Index (cbi)
best_cor_maxent <- subset(cor_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_cor_maxent <- eval.models(cor_maxent)[[best_cor_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.1_fc.LQHPT

eval.variable.importance(cor_maxent)[best_cor_maxent["tune.args"][1, 1]]

cor_maxent@results[order(cor_maxent@results$delta.AICc), c("rm", "fc", "delta.AICc")]

# BIO5 = Max Temperature of Warmest Month - 26.69% contribution
# BIO6 = Min Temperature of Coldest Month - 15.95% contribution
# BIO1 = Annual Mean Temperature - 15.898% contribution
# BIO15 = Precipitation Seasonality (Coefficient of Variation) - 12.27% contribution



# M. coronaria predictions ------------------------------------------------
# Now use the <terra> package to plot the SDM prediction.
# Wclim is the historical climatic conditions (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
cor_pred_hist <- terra::predict(wclim_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)

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

dev.off()
dev.new()
par(mar = c(4, 4, 4, 4), mfcol = c(1, 2))
terra::plot(cor_pred_hist > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            xlim = c(-100, -50), ylim = c(30, 60), 
            main = expression(atop(italic('Malus coronaria'), " Historical Suitability (1970-2000)")), 
            background = 'lightskyblue1', box = 'black')
terra::plot(cor_pred_hist > corPred_threshold_10, add = T, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), legend = F, xlim = c(-100, -50), ylim = c(30, 60))
terra::plot(cor_pred_hist > corPred_threshold_50, add = T, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), legend = F, xlim = c(-100, -50), ylim = c(30, 60))
#terra::plot(canUSMex_map, add = T, cex = .1)
#points(occThin_cor, col = 'black', cex = 0.75, pch = 4)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)


# Subsetted historical
terra::plot(cor_pred_hist_subs > corPred_threshold_1_subs, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            xlim = c(-100, -50), ylim = c(30, 60), 
            main = expression(atop(italic('Malus coronaria'), " Historical Suitability (1970-2000) (SUBSETTED MODEL)")), 
            background = 'lightskyblue1', box = 'black')
terra::plot(cor_pred_hist_subs > corPred_threshold_10_subs, add = T, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), legend = F, xlim = c(-100, -50), ylim = c(30, 60))
terra::plot(cor_pred_hist_subs > corPred_threshold_50_subs, add = T, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), legend = F, xlim = c(-100, -50), ylim = c(30, 60))
#terra::plot(canUSMex_map, add = T, cex = .1)
#points(occThin_cor, col = 'black', cex = 0.75, pch = 4)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)


# Future Climate predictions
# SSP 245
cor_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)

# SSP 585
cor_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)


dev.off()
dev.new()
par(mar = c(4, 4, 4, 4), mfcol = c(1, 2))
# Plot SSP 585 2030
terra::plot(cor_pred_ssp585_70 > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-100, -50), ylim = c(30, 60), main = expression(atop(italic('Malus coronaria'), " SSP5-8.5 Prediction: Late Century (2060-2080)")), background = 'lightskyblue1')
terra::plot(cor_pred_ssp585_70 > corPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_70 > corPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

# plot 
par(mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp585_30_subs > corPred_threshold_1_subs, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-100, -50), ylim = c(30, 60), main = expression(atop(italic('Malus coronaria'), " SSP5-8.5 Prediction: Late Century (2060-2080) (SUBSETTED MODEL)")), background = 'lightskyblue1')
terra::plot(cor_pred_ssp585_30_subs > corPred_threshold_10_subs, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_30_subs > corPred_threshold_50_subs, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
points(occThin_cor, col = 'black', cex = 0.75, pch = 4)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

# Save/Load M cor. SDM predictions ----------------------------------------
# Save
saveRDS(cor_pred_hist, file = './sdm_output/cor_pred_hist_subs.Rdata')

saveRDS(cor_pred_ssp245_30, file = './sdm_output/cor_pred_ssp245_30.Rdata')
saveRDS(cor_pred_ssp245_50, file = './sdm_output/cor_pred_ssp245_50.Rdata')
saveRDS(cor_pred_ssp245_70, file = './sdm_output/cor_pred_ssp245_70.Rdata')

saveRDS(cor_pred_ssp585_30, file = './sdm_output/cor_pred_ssp585_30.Rdata')
saveRDS(cor_pred_ssp585_50, file = './sdm_output/cor_pred_ssp585_50.Rdata')
saveRDS(cor_pred_ssp585_70, file = './sdm_output/cor_pred_ssp585_70.Rdata')

# Save
saveRDS(cor_pred_hist_subs, file = './sdm_output/cor_pred_hist_subs.Rdata')

saveRDS(cor_pred_ssp245_30_subs, file = './sdm_output/cor_pred_ssp245_30_subs.Rdata')
saveRDS(cor_pred_ssp245_50_subs, file = './sdm_output/cor_pred_ssp245_50_subs.Rdata')
saveRDS(cor_pred_ssp245_70_subs, file = './sdm_output/cor_pred_ssp245_70_subs.Rdata')

saveRDS(cor_pred_ssp585_30_subs, file = './sdm_output/cor_pred_ssp585_30_subs.Rdata')
saveRDS(cor_pred_ssp585_50_subs, file = './sdm_output/cor_pred_ssp585_50_subs.Rdata')
saveRDS(cor_pred_ssp585_70_subs, file = './sdm_output/cor_pred_ssp585_70_subs.Rdata')

# Load
cor_pred_hist <- readRDS(file = './sdm_output/cor_pred_hist.Rdata')

cor_pred_ssp245_30 <- readRDS(file = './sdm_output/cor_pred_ssp245_30.Rdata')
cor_pred_ssp245_50 <- readRDS(file = './sdm_output/cor_pred_ssp245_50.Rdata')
cor_pred_ssp245_70 <- readRDS(file = './sdm_output/cor_pred_ssp245_70.Rdata')

cor_pred_ssp585_30 <- readRDS(file = './sdm_output/cor_pred_ssp585_30.Rdata')
cor_pred_ssp585_50 <- readRDS(file = './sdm_output/cor_pred_ssp585_50.Rdata')
cor_pred_ssp585_70 <- readRDS(file = './sdm_output/cor_pred_ssp585_70.Rdata')

#subsetted climate varaible rasters
cor_pred_hist_subs <- readRDS(file = './sdm_output/cor_pred_hist_subs.Rdata')

cor_pred_ssp245_30_subs <- readRDS(file = './sdm_output/cor_pred_ssp245_30_subs.Rdata')
cor_pred_ssp245_50_subs <- readRDS(file = './sdm_output/cor_pred_ssp245_50_subs.Rdata')
cor_pred_ssp245_70_subs <- readRDS(file = './sdm_output/cor_pred_ssp245_70_subs.Rdata')

cor_pred_ssp585_30_subs <- readRDS(file = './sdm_output/cor_pred_ssp585_30_subs.Rdata')
cor_pred_ssp585_50_subs <- readRDS(file = './sdm_output/cor_pred_ssp585_50_subs.Rdata')
cor_pred_ssp585_70_subs <- readRDS(file = './sdm_output/cor_pred_ssp585_70_subs.Rdata')


# M coronaria thresholds --------------------------------------------------

saveRDS(corPred_threshold_1, file = './sdm_output/thresholds/corPred_threshold_1_subs.Rdata')
saveRDS(corPred_threshold_10, file = './sdm_output/thresholds/corPred_threshold_10_subs.Rdata')
saveRDS(corPred_threshold_50, file = './sdm_output/thresholds/corPred_threshold_50_subs.Rdata')

# Load
corPred_threshold_1 <- readRDS(file = './sdm_output/thresholds/corPred_threshold_1.Rdata')
corPred_threshold_10 <- readRDS(file = './sdm_output/thresholds/corPred_threshold_10.Rdata')
corPred_threshold_50 <- readRDS(file = './sdm_output/thresholds/corPred_threshold_50.Rdata')

corPred_threshold_1_subs <- readRDS(file = './sdm_output/thresholds/corPred_threshold_1_subs.Rdata')
corPred_threshold_10_subs <- readRDS(file = './sdm_output/thresholds/corPred_threshold_10_subs.Rdata')
corPred_threshold_50_subs <- readRDS(file = './sdm_output/thresholds/corPred_threshold_50_subs.Rdata')

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

# Historical
saveRDS(cor_pred_high_hist, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_hist_subs.Rdata')
saveRDS(cor_pred_mod_hist, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_hist_subs.Rdata')
saveRDS(cor_pred_low_hist, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(cor_pred_high_ssp245_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_ssp245_30_subs.Rdata')
saveRDS(cor_pred_mod_ssp245_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_ssp245_30_subs.Rdata')
saveRDS(cor_pred_low_ssp245_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_ssp245_30_subs.Rdata')

saveRDS(cor_pred_high_ssp245_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_ssp245_50_subs.Rdata')
saveRDS(cor_pred_mod_ssp245_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_ssp245_50_subs.Rdata')
saveRDS(cor_pred_low_ssp245_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_ssp245_50_subs.Rdata')

saveRDS(cor_pred_high_ssp245_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_ssp245_70_subs.Rdata')
saveRDS(cor_pred_mod_ssp245_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_ssp245_70_subs.Rdata')
saveRDS(cor_pred_low_ssp245_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(cor_pred_high_ssp585_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_ssp585_30_subs.Rdata')
saveRDS(cor_pred_mod_ssp585_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_ssp585_30_subs.Rdata')
saveRDS(cor_pred_low_ssp585_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_ssp585_30_subs.Rdata')

saveRDS(cor_pred_high_ssp585_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_ssp585_50_subs.Rdata')
saveRDS(cor_pred_mod_ssp585_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_ssp585_50_subs.Rdata')
saveRDS(cor_pred_low_ssp585_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_ssp585_50_subs.Rdata')

saveRDS(cor_pred_high_ssp585_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_high_ssp585_70_subs.Rdata')
saveRDS(cor_pred_mod_ssp585_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_mod_ssp585_70_subs.Rdata')
saveRDS(cor_pred_low_ssp585_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/cor_pred_low_ssp585_70_subs.Rdata')

# Load
setwd('../sdm_output/habitat_predictions/high_moderate_low_predictions')

# Historical
cor_pred_high_hist <- readRDS(file = 'cor_pred_high_hist.Rdata')
cor_pred_mod_hist <- readRDS(file = 'cor_pred_mod_hist.Rdata')
cor_pred_low_hist <- readRDS(file = 'cor_pred_low_hist.Rdata')

# SSP245
cor_pred_high_ssp245_30 <- readRDS(file = 'cor_pred_high_ssp245_30.Rdata')
cor_pred_mod_ssp245_30 <- readRDS(file = 'cor_pred_mod_ssp245_30.Rdata')
cor_pred_low_ssp245_30 <- readRDS(file = 'cor_pred_low_ssp245_30.Rdata')

cor_pred_high_ssp245_50 <- readRDS(file = 'cor_pred_high_ssp245_50.Rdata')
cor_pred_mod_ssp245_50 <- readRDS(file = 'cor_pred_mod_ssp245_50.Rdata')
cor_pred_low_ssp245_50 <- readRDS(file = 'cor_pred_low_ssp245_50.Rdata')

cor_pred_high_ssp245_70 <- readRDS(file = 'cor_pred_high_ssp245_70.Rdata')
cor_pred_mod_ssp245_70 <- readRDS(file = 'cor_pred_mod_ssp245_70.Rdata')
cor_pred_low_ssp245_70 <- readRDS(file = 'cor_pred_low_ssp245_70.Rdata')

# SSP585
cor_pred_high_ssp585_30 <- readRDS(file = 'cor_pred_high_ssp585_30.Rdata')
cor_pred_mod_ssp585_30 <- readRDS(file = 'cor_pred_mod_ssp585_30.Rdata')
cor_pred_low_ssp585_30 <- readRDS(file = 'cor_pred_low_ssp585_30.Rdata')

cor_pred_high_ssp585_50 <- readRDS(file = 'cor_pred_high_ssp585_50.Rdata')
cor_pred_mod_ssp585_50 <- readRDS(file = 'cor_pred_mod_ssp585_50.Rdata')
cor_pred_low_ssp585_50 <- readRDS(file = 'cor_pred_low_ssp585_50.Rdata')

cor_pred_high_ssp585_70 <- readRDS(file = 'cor_pred_high_ssp585_70.Rdata')
cor_pred_mod_ssp585_70 <- readRDS(file = 'cor_pred_mod_ssp585_70.Rdata')
cor_pred_low_ssp585_70 <- readRDS(file = 'cor_pred_low_ssp585_70.Rdata')

# Crop categorical layers to restrict the extent of predicted suitability.
# In this case the model is making predictions outside of what ecologicaly makes sense!
cor_ext <- ext(-100, -52, 30, 70)

# Historical
cor_pred_high_hist_crop <- crop(cor_pred_high_hist, cor_ext)
cor_pred_mod_hist_crop <- crop(cor_pred_mod_hist, cor_ext)
cor_pred_low_hist_crop <- crop(cor_pred_low_hist, cor_ext)

#SSP245 
cor_pred_high_ssp245_30_crop <- crop(cor_pred_high_ssp245_30, cor_ext)
cor_pred_mod_ssp245_30_crop <- crop(cor_pred_mod_ssp245_30, cor_ext)
cor_pred_low_ssp245_30_crop <- crop(cor_pred_low_ssp245_30, cor_ext)

cor_pred_high_ssp245_50_crop <- crop(cor_pred_high_ssp245_50, cor_ext)
cor_pred_mod_ssp245_50_crop <- crop(cor_pred_mod_ssp245_50, cor_ext)
cor_pred_low_ssp245_50_crop <- crop(cor_pred_low_ssp245_50, cor_ext)

cor_pred_high_ssp245_70_crop <- crop(cor_pred_high_ssp245_70, cor_ext)
cor_pred_mod_ssp245_70_crop <- crop(cor_pred_mod_ssp245_70, cor_ext)
cor_pred_low_ssp245_70_crop <- crop(cor_pred_low_ssp245_70, cor_ext)

#SSP585
cor_pred_high_ssp585_30_crop <- crop(cor_pred_high_ssp585_30, cor_ext)
cor_pred_mod_ssp585_30_crop <- crop(cor_pred_mod_ssp585_30, cor_ext)
cor_pred_low_ssp585_30_crop <- crop(cor_pred_low_ssp585_30, cor_ext)

cor_pred_high_ssp585_50_crop <- crop(cor_pred_high_ssp585_50, cor_ext)
cor_pred_mod_ssp585_50_crop <- crop(cor_pred_mod_ssp585_50, cor_ext)
cor_pred_low_ssp585_50_crop <- crop(cor_pred_low_ssp585_50, cor_ext)

cor_pred_high_ssp585_70_crop <- crop(cor_pred_high_ssp585_70, cor_ext)
cor_pred_mod_ssp585_70_crop <- crop(cor_pred_mod_ssp585_70, cor_ext)
cor_pred_low_ssp585_70_crop <- crop(cor_pred_low_ssp585_70, cor_ext)


# Save Cropped rasters as Tiffs -------------------------------------------
setwd('../..')
getwd()
setwd('./sdm_output/habitat_predictions/high_moderate_low_predictions/cropped_predictions')

# Historical 
terra::writeRaster(cor_pred_high_hist_crop, "cor_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_hist_crop, "cor_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_hist_crop, "cor_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(cor_pred_high_ssp245_30_crop, "cor_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp245_30_crop, "cor_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp245_30_crop, "cor_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp245_50_crop, "cor_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp245_50_crop, "cor_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp245_50_crop, "cor_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp245_70_crop, "cor_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp245_70_crop, "cor_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp245_70_crop, "cor_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(cor_pred_high_ssp585_30_crop, "cor_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp585_30_crop, "cor_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp585_30_crop, "cor_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp585_50_crop, "cor_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp585_50_crop, "cor_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp585_50_crop, "cor_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp585_70_crop, "cor_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp585_70_crop, "cor_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp585_70_crop, "cor_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)


# Fusca - MaxEnt Model ----------------------------------------------------

occ_fus_coords <- as.data.frame(geom(occThin_fus)[,3:4]) # extract longitude, lattitude from occurence points
bg_fus_coords <- as.data.frame(geom(fus_bg_vec)[,3:4]) # extract longitude, lattitude from background points

# Build Species Distribution Model using MaxEnt from the <ENMeval> package

# Run prediction in a parallel using 'socket' clusters to help speed up computation
# <ENMeval> implements parallel functions natively
# But load <parallel> library for additional functions like <decectCores()>
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

fus_maxent <- ENMevaluate(occ_fus_coords, # occurrence records
                          envs = wclim_fus_subs, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5),
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
saveRDS(fus_maxent, file = './sdm_output/fus_maxent_subs.Rdata') # save
fus_maxent <- readRDS(file = './sdm_output/fus_maxent.Rdata') # load 
fus_maxent_subs <- readRDS(file = './sdm_output/fus_maxent.Rdata') # subset model 



# M. fusca Model Selection ------------------------------------------------
# Note that maxent results provide Continuous Boyce Index (cbi)
best_fus_maxent <- subset(fus_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_fus_maxent <- eval.models(fus_maxent)[[best_fus_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.2_fc.LQHPT

eval.variable.importance(fus_maxent)[best_fus_maxent["tune.args"][1, 1]]
# BIO19 = Precipitation of Coldest Quarter - 41.73% Contribution
# BIO7 = Temperature Annual Range (BIO5-BIO6) (Max Temp Warmest Month - Min Temp of Coldest Month) - 20.219% Contribution
# BIO10 = Mean Temperature of Warmest Quarter - 10.782% Contribution
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) - 10.54% Contribution

# M. fusca Predictions ----------------------------------------------------
# Now use the <terra> package to plot the SDM prediction.
# Wclim is the historical climatic conditions (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
fus_pred_hist <- terra::predict(wclim_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)

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



dev.off()
dev.new()
par(mar = c(4, 4, 4, 4), mfcol = c(1, 2))
terra::plot(fus_pred_hist > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " Historical Suitability (1970-2000)")), background = 'lightskyblue1')
terra::plot(fus_pred_hist > fusPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_hist > fusPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
#points(occThin_fus, col = 'black', cex = 0.75, pch = 4)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

terra::plot(fus_pred_hist_subs > fusPred_threshold_1_subs, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " Historical Suitability (1970-2000) (SUBSETTED MODEL)")), background = 'lightskyblue1')
terra::plot(fus_pred_hist_subs > fusPred_threshold_10_subs, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_hist_subs > fusPred_threshold_50_subs, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
#points(occThin_fus, col = 'black', cex = 0.75, pch = 4)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)



# Future Climate predictions
# SSP 245
fus_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)

# SSP 585
fus_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)
fus_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_fus_maxent, cores = cn - 1, na.rm = T)

# Plot SSP 585 2030
dev.off()
dev.new()
par(mar = c(4, 4, 4, 4), mfcol = c(1, 2))
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " SSP5-8.5 Prediction: Early Century (2020-2040)")), background = 'lightskyblue1')
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
#points(occThin_fus, col = 'black', cex = 0.75, pch = 4)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

terra::plot(fus_pred_ssp585_30_subs > fusPred_threshold_1_subs, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " SSP5-8.5 Prediction: Early Century (2020-2040) (SUBSETTED MODEL)")), background = 'lightskyblue1')
terra::plot(fus_pred_ssp585_30_subs > fusPred_threshold_10_subs, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30_subs > fusPred_threshold_50_subs, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
#terra::plot(canUSMex_map, add = T)
#points(occThin_fus, col = 'black', cex = 0.75, pch = 4)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)



# Save/Load M fus. SDM predictions ----------------------------------------

# Save
saveRDS(fus_pred_hist, file = './sdm_output/fus_pred_hist_subs.Rdata')

saveRDS(fus_pred_ssp245_30, file = './sdm_output/fus_pred_ssp245_30_subs.Rdata')
saveRDS(fus_pred_ssp245_50, file = './sdm_output/fus_pred_ssp245_50_subs.Rdata')
saveRDS(fus_pred_ssp245_70, file = './sdm_output/fus_pred_ssp245_70_subs.Rdata')

saveRDS(fus_pred_ssp585_30, file = './sdm_output/fus_pred_ssp585_30_subs.Rdata')
saveRDS(fus_pred_ssp585_50, file = './sdm_output/fus_pred_ssp585_50_subs.Rdata')
saveRDS(fus_pred_ssp585_70, file = './sdm_output/fus_pred_ssp585_70_subs.Rdata')

# Load
fus_pred_hist <- readRDS(file = './sdm_output/fus_pred_hist.Rdata')

fus_pred_ssp245_30 <- readRDS(file = './sdm_output/fus_pred_ssp245_30.Rdata')
fus_pred_ssp245_50 <- readRDS(file = './sdm_output/fus_pred_ssp245_50.Rdata')
fus_pred_ssp245_70 <- readRDS(file = './sdm_output/fus_pred_ssp245_70.Rdata')

fus_pred_ssp585_30 <- readRDS(file = './sdm_output/fus_pred_ssp585_30.Rdata')
fus_pred_ssp585_50 <- readRDS(file = './sdm_output/fus_pred_ssp585_50.Rdata')
fus_pred_ssp585_70 <- readRDS(file = './sdm_output/fus_pred_ssp585_70.Rdata')

#Subsetted models
fus_pred_hist_subs <- readRDS(file = './sdm_output/fus_pred_hist_subs.Rdata')

fus_pred_ssp245_30_subs <- readRDS(file = './sdm_output/fus_pred_ssp245_30_subs.Rdata')
fus_pred_ssp245_50_subs <- readRDS(file = './sdm_output/fus_pred_ssp245_50_subs.Rdata')
fus_pred_ssp245_70_subs <- readRDS(file = './sdm_output/fus_pred_ssp245_70_subs.Rdata')

fus_pred_ssp585_30_subs <- readRDS(file = './sdm_output/fus_pred_ssp585_30_subs.Rdata')
fus_pred_ssp585_50_subs <- readRDS(file = './sdm_output/fus_pred_ssp585_50_subs.Rdata')
fus_pred_ssp585_70_subs <- readRDS(file = './sdm_output/fus_pred_ssp585_70_subs.Rdata')

# Save M. fusca thresholds ------------------------------------------------

saveRDS(fusPred_threshold_1, file = './sdm_output/thresholds/fusPred_threshold_1_subs.Rdata')
saveRDS(fusPred_threshold_10, file = './sdm_output/thresholds/fusPred_threshold_10_subs.Rdata')
saveRDS(fusPred_threshold_50, file = './sdm_output/thresholds/fusPred_threshold_50_subs.Rdata')

# Load
fusPred_threshold_1 <- readRDS(file = './sdm_output/thresholds/fusPred_threshold_1.Rdata')
fusPred_threshold_10 <- readRDS(file = './sdm_output/thresholds/fusPred_threshold_10.Rdata')
fusPred_threshold_50 <- readRDS(file = './sdm_output/thresholds/fusPred_threshold_50.Rdata')

# Load subsetted
fusPred_threshold_1_subs <- readRDS(file = './sdm_output/thresholds/fusPred_threshold_1_subs.Rdata')
fusPred_threshold_10_subs <- readRDS(file = './sdm_output/thresholds/fusPred_threshold_10_subs.Rdata')
fusPred_threshold_50_subs <- readRDS(file = './sdm_output/thresholds/fusPred_threshold_50_subs.Rdata')


# Prediction plotting -----------------------------------------------------



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

# Historical
saveRDS(fus_pred_high_hist, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_hist_subs.Rdata')
saveRDS(fus_pred_mod_hist, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_hist_subs.Rdata')
saveRDS(fus_pred_low_hist, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(fus_pred_high_ssp245_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_ssp245_30_subs.Rdata')
saveRDS(fus_pred_mod_ssp245_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_ssp245_30_subs.Rdata')
saveRDS(fus_pred_low_ssp245_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_ssp245_30_subs.Rdata')

saveRDS(fus_pred_high_ssp245_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_ssp245_50_subs.Rdata')
saveRDS(fus_pred_mod_ssp245_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_ssp245_50_subs.Rdata')
saveRDS(fus_pred_low_ssp245_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_ssp245_50_subs.Rdata')

saveRDS(fus_pred_high_ssp245_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_ssp245_70_subs.Rdata')
saveRDS(fus_pred_mod_ssp245_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_ssp245_70_subs.Rdata')
saveRDS(fus_pred_low_ssp245_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(fus_pred_high_ssp585_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_ssp585_30_subs.Rdata')
saveRDS(fus_pred_mod_ssp585_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_ssp585_30_subs.Rdata')
saveRDS(fus_pred_low_ssp585_30, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_ssp585_30_subs.Rdata')

saveRDS(fus_pred_high_ssp585_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_ssp585_50_subs.Rdata')
saveRDS(fus_pred_mod_ssp585_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_ssp585_50_subs.Rdata')
saveRDS(fus_pred_low_ssp585_50, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_ssp585_50_subs.Rdata')

saveRDS(fus_pred_high_ssp585_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_high_ssp585_70_subs.Rdata')
saveRDS(fus_pred_mod_ssp585_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_mod_ssp585_70_subs.Rdata')
saveRDS(fus_pred_low_ssp585_70, file = './sdm_output/habitat_predictions/high_moderate_low_predictions/fus_pred_low_ssp585_70_subs.Rdata')

# Load
# Load
setwd('./sdm_output/habitat_predictions/high_moderate_low_predictions/')

# Historical
fus_pred_high_hist <- readRDS(file = 'fus_pred_high_hist.Rdata')
fus_pred_mod_hist <- readRDS(file = 'fus_pred_mod_hist.Rdata')
fus_pred_low_hist <- readRDS(file = 'fus_pred_low_hist.Rdata')

# SSP245
fus_pred_high_ssp245_30 <- readRDS(file = 'fus_pred_high_ssp245_30.Rdata')
fus_pred_mod_ssp245_30 <- readRDS(file = 'fus_pred_mod_ssp245_30.Rdata')
fus_pred_low_ssp245_30 <- readRDS(file = 'fus_pred_low_ssp245_30.Rdata')

fus_pred_high_ssp245_50 <- readRDS(file = 'fus_pred_high_ssp245_50.Rdata')
fus_pred_mod_ssp245_50 <- readRDS(file = 'fus_pred_mod_ssp245_50.Rdata')
fus_pred_low_ssp245_50 <- readRDS(file = 'fus_pred_low_ssp245_50.Rdata')

fus_pred_high_ssp245_70 <- readRDS(file = 'fus_pred_high_ssp245_70.Rdata')
fus_pred_mod_ssp245_70 <- readRDS(file = 'fus_pred_mod_ssp245_70.Rdata')
fus_pred_low_ssp245_70 <- readRDS(file = 'fus_pred_low_ssp245_70.Rdata')

# SSP585
fus_pred_high_ssp585_30 <- readRDS(file = 'fus_pred_high_ssp585_30.Rdata')
fus_pred_mod_ssp585_30 <- readRDS(file = 'fus_pred_mod_ssp585_30.Rdata')
fus_pred_low_ssp585_30 <- readRDS(file = 'fus_pred_low_ssp585_30.Rdata')

fus_pred_high_ssp585_50 <- readRDS(file = 'fus_pred_high_ssp585_50.Rdata')
fus_pred_mod_ssp585_50 <- readRDS(file = 'fus_pred_mod_ssp585_50.Rdata')
fus_pred_low_ssp585_50 <- readRDS(file = 'fus_pred_low_ssp585_50.Rdata')

fus_pred_high_ssp585_70 <- readRDS(file = 'fus_pred_high_ssp585_70.Rdata')
fus_pred_mod_ssp585_70 <- readRDS(file = 'fus_pred_mod_ssp585_70.Rdata')
fus_pred_low_ssp585_70 <- readRDS(file = 'fus_pred_low_ssp585_70.Rdata')


# Crop categorical layers to restrict the extent of predicted suitability.
# In this case the model is making predictions outside of what ecologicaly makes sense!
fus_ext <- ext(-180, -100, 30, 70)

# Historical
fus_pred_high_hist_crop <- crop(fus_pred_high_hist, fus_ext)
fus_pred_mod_hist_crop <- crop(fus_pred_mod_hist, fus_ext)
fus_pred_low_hist_crop <- crop(fus_pred_low_hist, fus_ext)

#SSP245 
fus_pred_high_ssp245_30_crop <- crop(fus_pred_high_ssp245_30, fus_ext)
fus_pred_mod_ssp245_30_crop <- crop(fus_pred_mod_ssp245_30, fus_ext)
fus_pred_low_ssp245_30_crop <- crop(fus_pred_low_ssp245_30, fus_ext)

fus_pred_high_ssp245_50_crop <- crop(fus_pred_high_ssp245_50, fus_ext)
fus_pred_mod_ssp245_50_crop <- crop(fus_pred_mod_ssp245_50, fus_ext)
fus_pred_low_ssp245_50_crop <- crop(fus_pred_low_ssp245_50, fus_ext)

fus_pred_high_ssp245_70_crop <- crop(fus_pred_high_ssp245_70, fus_ext)
fus_pred_mod_ssp245_70_crop <- crop(fus_pred_mod_ssp245_70, fus_ext)
fus_pred_low_ssp245_70_crop <- crop(fus_pred_low_ssp245_70, fus_ext)

#SSP585
fus_pred_high_ssp585_30_crop <- crop(fus_pred_high_ssp585_30, fus_ext)
fus_pred_mod_ssp585_30_crop <- crop(fus_pred_mod_ssp585_30, fus_ext)
fus_pred_low_ssp585_30_crop <- crop(fus_pred_low_ssp585_30, fus_ext)

fus_pred_high_ssp585_50_crop <- crop(fus_pred_high_ssp585_50, fus_ext)
fus_pred_mod_ssp585_50_crop <- crop(fus_pred_mod_ssp585_50, fus_ext)
fus_pred_low_ssp585_50_crop <- crop(fus_pred_low_ssp585_50, fus_ext)

fus_pred_high_ssp585_70_crop <- crop(fus_pred_high_ssp585_70, fus_ext)
fus_pred_mod_ssp585_70_crop <- crop(fus_pred_mod_ssp585_70, fus_ext)
fus_pred_low_ssp585_70_crop <- crop(fus_pred_low_ssp585_70, fus_ext)

# Save Cropped rasters as Tiffs -------------------------------------------
getwd()
setwd('./cropped_predictions')

# Historical 
terra::writeRaster(fus_pred_high_hist_crop, "fus_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_hist_crop, "fus_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_hist_crop, "fus_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(fus_pred_high_ssp245_30_crop, "fus_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp245_30_crop, "fus_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp245_30_crop, "fus_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp245_50_crop, "fus_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp245_50_crop, "fus_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp245_50_crop, "fus_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp245_70_crop, "fus_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp245_70_crop, "fus_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp245_70_crop, "fus_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(fus_pred_high_ssp585_30_crop, "fus_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp585_30_crop, "fus_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp585_30_crop, "fus_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp585_50_crop, "fus_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp585_50_crop, "fus_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp585_50_crop, "fus_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp585_70_crop, "fus_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp585_70_crop, "fus_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp585_70_crop, "fus_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)


# Binary Thresholds for Gap Analysis --------------------------------------
# For the purposed of the Gap Analysis it is best to use a binary threshold of habitat suitability.
# You could use the moderate, or high suitability thresholds from above, but that is not very well accepted in the literature.
# In this case there are several binary thresholds that exist, including several in the <predicts> package
# See ?predicts::threshold


cor_pa <- predicts::pa_evaluate(p = occ_cor_coords_mat, a = bg_cor_coords_mat, model = cor_maxent, x = wclim_cor)
cor_binary_threshold <- predicts::threshold(cor_pa)

cor_hist_habitat <- cor_pred_hist > cor_threshold$max_spec_sens #the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest


# Area calculations -------------------------------------------------------
area_cor_high_hist_subs <- expanse(cor_pred_high_hist_subs, byValue = T, unit = 'km') %>% filter(value == 1) %>% pull(area)
area_cor_high_ssp585_30_subs <- expanse(cor_pred_high_585_30_subs, byValue = T, unit = 'km') %>% filter(value == 1) %>% pull(area)

((area_cor_high_ssp585_30_subs - area_cor_high_hist_subs)/(area_cor_high_hist_subs))*100


# Occurrences in suitability ----------------------------------------------
cor_suit_df <- extract(cor_pred_low_ssp585_50, occThin_cor)

cor_suit_df_f <- cor_suit_df %>% filter(lyr1 == 'FALSE')
cor_suit_df_t <- cor_suit_df %>% filter(lyr1 == 'TRUE')

nrow(cor_suit_df_f)/nrow(cor_suit_df)*100
