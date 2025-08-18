# Top ---------------------------------------------------------------------
# MaxEnt Species Distribution Modeling (SDM) for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024
START <- date()
#source("scripts/malus_bg.R") 
END <- date()

source("scripts/functions.R") ## for twsBoyce

message("Loaded & prepped all data, starting at:\n     ", START,
        "\n and ending at: \n     ", END) 

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


# Install Tyler's hotfix for Ecospat Boyce Index --------------------------

library(devtools)
install_github("plantarum/ecospat", ref = "boyce", subdir = "ecospat")

library(ecospat)

# Load occurrence data and basemaps -------------------------------------------------------

## # Background points in SpatVectors
cor_bg_vec <- readRDS(file = './occ_data/cor/cor_bg_vec.Rdata')
fus_bg_vec <- readRDS(file = './occ_data/fus/fus_bg_vec.Rdata')
ion_bg_vec <- readRDS(file = './occ_data/ion/ion_bg_vec.Rdata')
ang_bg_vec <- readRDS(file = './occ_data/ang/ang_bg_vec.Rdata')
chl_bg_vec <- readRDS(file = './occ_data/chl/chl_bg_vec.Rdata')

## # Occurrence Points in SpatVectors
occThin_cor <- readRDS(file = './occ_data/cor/occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = './occ_data/fus/occThin_fus.Rdata') # M. fusca
occThin_ion <- readRDS(file = './occ_data/ion/occThin_ion.Rdata') # M. fusca
occThin_ang <- readRDS(file = './occ_data/ang/occThin_ang.Rdata') # M. fusca
occThin_chl <- readRDS(file = './occ_data/chl/occThin_chl.Rdata') # M. fusca

## # Great Lakes shapefiles for making pretty maps and cropping
## great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')


# Climate Data
great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')

NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in Western Hemisphere

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
#wclim_cor_stack <- raster::stack(wclim_cor) # covert SpatRaster to RasterStack for dependency in ENMeval checkboarding


wclim_cor <- readRDS(file = './wclim_data/wclim_cor.Rdata')
wclim_fus <- readRDS(file = './wclim_data/wclim_fus.Rdata')
wclim_ion <- readRDS(file = './wclim_data/wclim_ion.Rdata')
wclim_ang <- readRDS(file = './wclim_data/wclim_ang.Rdata')
wclim_chl <- readRDS(file = './wclim_data/wclim_chl.Rdata')

wclim <- geodata::worldclim_global(var = 'bio',
                                   res = 2.5, 
                                   version = '2.1', 
                                   path = "./wclim_data/") %>% 
  terra::crop(NA_ext)  %>% #crop raster to NA 
  terra::mask(great_lakes, inverse = T) # cut out the great lakes

# SSP (Shared social-economic pathway) 2.45 
# middle of the road projection, high climate adaptation, low climate mitigation
climate_predictors <- names(wclim_cor) # extract climate predictor names, to ren
# Future SSPs
# Do not need to create RasterStacks
# SSP 245
names(wclim) <- climate_predictors
names(ssp245_2030) <- climate_predictors #rename raster layers for downsteam analysis
names(ssp245_2050) <- climate_predictors 
names(ssp245_2070) <- climate_predictors 
names(ssp585_2030) <- climate_predictors
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
wclim_ion_subs <- wclim_ion %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_ang_subs <- wclim_ang %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_chl_subs <- wclim_chl %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))

# M. coronaria - MaxEnt Model ---------------------------------------------

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

date()
cor_maxent <- ENMevaluate(occ_cor_coords, # occurrence records
                          envs = wclim_cor_subs, # NOTE CHANGE THE ENVS
                                        # inputed, environment from background training area
                          n.bg = 20000, # 
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5), # Regularization 0.5-4
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                            partitions = 'checkerboard',
                            parallel = TRUE,
                            numCores = cn - 1, # leave one core available for other apps
                          algorithm = 'maxent.jar')
date()

# Save the MaxEnt model so you do not have to waste time re-running the model
##saveRDS(cor_maxent, file = './sdm_output/cor/subs/cor_maxent_subs.Rdata') # save
# Load Maxent model
saveRDS(cor_maxent, file = './sdm_output/cor/subs/cor_maxent_subs.Rdata')
cor_maxent <- readRDS(file = './sdm_output/cor/subs/cor_maxent_subs.Rdata')

##cor_maxent <- readRDS(file = './sdm_output/cor/subs/cor_maxent_subs.Rdata') # load #subsetted model


# M. coronaria Model Selection --------------------------------------------
# Note that maxent results provide Continuous Boyce Index (cbi)

best_cor_maxent <- subset(cor_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_cor_maxent <- eval.models(cor_maxent)[[best_cor_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.1_fc.LQHPT

eval.variable.importance(cor_maxent)[best_cor_maxent["tune.args"][1, 1]]

cor_maxent@results[order(cor_maxent@results$delta.AICc),
                   c("rm", "fc", "delta.AICc")] 

plot(mod.best_cor_maxent)

partialResponse(model = mod.best_cor_maxent, var = "wc2.1_2.5m_bio_10",
                main = "Bio 10")

# M. coronaria historical prediction --------------------------------------
# Now use the <terra> package to plot the SDM prediction.
# Wclim is the historical climatic conditions (1970-2000)
set.seed(1337)
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


# M. coronaria Boyce Index ----------------------------------------------------
# Evaluate predictions using Boyce Index
png('C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/boyce_index/corboyce_plot_malus_coronaria.png', width = 1600, height = 1200, res = 300)
corBoyce <- ecospat.boyce(fit = corPred_bg_val, # vector of predicted habitat suitability of bg points
              obs = corPred_val_na, # vector of 
              nclass = 0, 
              PEplot = TRUE,
              method = 'spearman')

title(main = bquote0(italic("Malus coronaria") ~ ", Boyce Cor." ==
                      .(corBoyce$cor))) 
dev.off()

## I thought limiting the prediction to eastern NA might speed it up. If it
## does, it's not by a lot?

# corPredExt <- ext(c(-100, -50, 30, 60))
# corPredWclim <- crop(wclim_subs, corPredExt)
# 
# cor_pred_hist <- terra::predict(corPredWclim, mod.best_cor_maxent,
#                                 cores = cn - 1, na.rm = T)
# 
# saveRDS(cor_pred_hist, file = './sdm_output/cor/subs/cor_pred_hist_2025-04-30.Rdata')
# cor_pred_hist <- readRDS(file = './sdm_output/cor/subs/cor_pred_hist_2025-04-30.Rdata')
# 
# plot(cor_pred_hist)
# points(occThin_cor, cex = 0.5, col = 1, pch = 21, bg = "white" )



# M. coronaria Boyce Index ------------------------------------------------
# Evaluate predictions using Boyce Index
# the number of true presences should decline with suitability groups 100-91, 90-81, etc. 

## Note there is a bug in ecospat.boyce,
## see https://github.com/ecospat/ecospat/issues/99
## Until it is fixed, use my own corrected version, in functions.R:

# Evaluate predictions using Boyce Index
## png('C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/boyce_index/corboyce_plot_malus_coronaria.png', width = 1600, height = 1200, res = 300)
# 
# corBoyce <- twsBoyce(fit = cor_pred_hist, # suitability raster
#                      obs = occ_cor_coords, 
#                      PEplot = TRUE, method = 'spearman')
# 
# title(main = bquote(italic("Malus coronaria") ~ ",Boyce Cor." ==
#                       .(corBoyce$cor))) 
# 
# 
# 
# ## dev.off()



# M. coronaria - Thresholds -----------------------------------------------
# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
corPred_val <- terra::extract(cor_pred_hist, occ_cor_coords)$lyr1
corPred_threshold_1 <- quantile(corPred_val, 0.01, na.rm = T) # Low suitability
corPred_threshold_10 <- quantile(corPred_val, 0.1, na.rm = T) # Moderate suitability
corPred_threshold_50 <- quantile(corPred_val, 0.5, na.rm = T) # High suitability

# Plotting the prediction
legend_labs <- c('Low Suitability', 'Moderate Suitability', 'High Suitability')

# brewer.pal(3, "YlOrBr") # Brewer palettes are helpful for mapping colour gradients
fill_cols <- c("#FFF7BC", "#FEC44F", "#D95F0E")

# Subsetted historical
terra::plot(cor_pred_hist > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            xlim = c(-100, -50), ylim = c(30, 60), 
            main = expression(atop(italic('Malus coronaria'),
                                   " Historical Suitability (1970-2000) (SUBSETTED MODEL)")), 
            background = 'lightskyblue1', box = 'black')
terra::plot(cor_pred_hist> corPred_threshold_10, add = T,
            col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), legend = F,
            xlim = c(-100, -50), ylim = c(30, 60))
terra::plot(cor_pred_hist > corPred_threshold_50, add = T,
            col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'),
            legend = F, xlim = c(-100, -50), ylim = c(30, 60)) 
#points(occThin_cor, col = 'black', cex = 0.75, pch = 4)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)


######################################################################
## Tyler got this far. Moving on to other parts, as the projections ##
## should not present any new issues?                               ##
######################################################################

# M. coronaria climate predictions ----------------------------------------
# SSP 245
cor_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_cor_maxent,
                                     cores = cn - 1, na.rm = T)
cor_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_cor_maxent,
                                     cores = cn - 1, na.rm = T)
cor_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_cor_maxent,
                                     cores = cn - 1, na.rm = T)

# SSP 585
cor_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)
cor_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_cor_maxent, cores = cn - 1, na.rm = T)


# Quick comparison of historical and 2030s SSP585 (today's climate)
# Plot historical 

# Plotting the prediction
legend_labs <- c('Low Suitability', 'Moderate Suitability', 'High Suitability')

# brewer.pal(3, "YlOrBr") # Brewer palettes are helpful for mapping colour gradients
fill_cols <- c("#FFF7BC", "#FEC44F", "#D95F0E")

par(mar = c(4, 4, 4, 4), mfcol = c(1, 2))
terra::plot(cor_pred_hist > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-100, -50), ylim = c(30, 60), main = expression(atop(italic('Malus coronaria'), " SSP5-8.5 Prediction: Late Century (2060-2080)")), background = 'lightskyblue1')
terra::plot(cor_pred_hist > corPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_hist > corPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

# Plot SSP585 2030
par(mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp585_30 > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-100, -50), ylim = c(30, 60), main = expression(atop(italic('Malus coronaria'), " SSP5-8.5 Prediction: Late Century (2060-2080) (SUBSETTED MODEL)")), background = 'lightskyblue1')
terra::plot(cor_pred_ssp585_30 > corPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_30 > corPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
points(occThin_cor, col = 'black', cex = 0.75, pch = 4)
legend(x = -72, y = 40, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

# Save/Load M cor. SDM predictions ----------------------------------------
# Save
saveRDS(cor_pred_hist, file = './sdm_output/cor/subs/cor_pred_hist_subs.Rdata')

saveRDS(cor_pred_ssp245_30, file = './sdm_output/cor/subs/cor_pred_ssp245_30_subs.Rdata')
saveRDS(cor_pred_ssp245_50, file = './sdm_output/cor/subs/cor_pred_ssp245_50_subs.Rdata')
saveRDS(cor_pred_ssp245_70, file = './sdm_output/cor/subs/cor_pred_ssp245_70_subs.Rdata')

saveRDS(cor_pred_ssp585_30, file = './sdm_output/cor/subs/cor_pred_ssp585_30_subs.Rdata')
saveRDS(cor_pred_ssp585_50, file = './sdm_output/cor/subs/cor_pred_ssp585_50_subs.Rdata')
saveRDS(cor_pred_ssp585_70, file = './sdm_output/cor/subs/cor_pred_ssp585_70_subs.Rdata')

# Load
cor_pred_hist <- readRDS(file = './sdm_output/cor/subs/cor_pred_hist_subs.Rdata')

cor_pred_ssp245_30 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp245_30_subs.Rdata')
cor_pred_ssp245_50 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp245_50_subs.Rdata')
cor_pred_ssp245_70 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp245_70_subs.Rdata')

cor_pred_ssp585_30 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_30_subs.Rdata')
cor_pred_ssp585_50 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_50_subs.Rdata')
cor_pred_ssp585_70 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_70_subs.Rdata')

# M coronaria thresholds --------------------------------------------------

saveRDS(corPred_threshold_1, file = './sdm_output/cor/subs/threshold/corPred_threshold_1_subs.Rdata')
saveRDS(corPred_threshold_10, file = './sdm_output/cor/subs/threshold/corPred_threshold_10_subs.Rdata')
saveRDS(corPred_threshold_50, file = './sdm_output/cor/subs/threshold/corPred_threshold_50_subs.Rdata')

# Load
corPred_threshold_1 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_1_subs.Rdata')
corPred_threshold_10 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_10_subs.Rdata')
corPred_threshold_50 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_50_subs.Rdata')

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
saveRDS(cor_pred_high_hist, file = './sdm_output/cor/subs/habitat_pred/hist/cor_pred_high_hist_subs.Rdata')
saveRDS(cor_pred_mod_hist, file = './sdm_output/cor/subs/habitat_pred/hist/cor_pred_mod_hist_subs.Rdata')
saveRDS(cor_pred_low_hist, file = './sdm_output/cor/subs/habitat_pred/hist/cor_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(cor_pred_high_ssp245_30, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_high_ssp245_30_subs.Rdata')
saveRDS(cor_pred_mod_ssp245_30, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_mod_ssp245_30_subs.Rdata')
saveRDS(cor_pred_low_ssp245_30, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_low_ssp245_30_subs.Rdata')

saveRDS(cor_pred_high_ssp245_50, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_high_ssp245_50_subs.Rdata')
saveRDS(cor_pred_mod_ssp245_50, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_mod_ssp245_50_subs.Rdata')
saveRDS(cor_pred_low_ssp245_50, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_low_ssp245_50_subs.Rdata')

saveRDS(cor_pred_high_ssp245_70, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_high_ssp245_70_subs.Rdata')
saveRDS(cor_pred_mod_ssp245_70, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_mod_ssp245_70_subs.Rdata')
saveRDS(cor_pred_low_ssp245_70, file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(cor_pred_high_ssp585_30, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_high_ssp585_30_subs.Rdata')
saveRDS(cor_pred_mod_ssp585_30, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_mod_ssp585_30_subs.Rdata')
saveRDS(cor_pred_low_ssp585_30, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_low_ssp585_30_subs.Rdata')

saveRDS(cor_pred_high_ssp585_50, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_high_ssp585_50_subs.Rdata')
saveRDS(cor_pred_mod_ssp585_50, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_mod_ssp585_50_subs.Rdata')
saveRDS(cor_pred_low_ssp585_50, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_low_ssp585_50_subs.Rdata')

saveRDS(cor_pred_high_ssp585_70, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_high_ssp585_70_subs.Rdata')
saveRDS(cor_pred_mod_ssp585_70, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_mod_ssp585_70_subs.Rdata')
saveRDS(cor_pred_low_ssp585_70, file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_low_ssp585_70_subs.Rdata')

# Load

# Historical
cor_pred_high_hist <- readRDS(file = './sdm_output/cor/subs/habitat_pred/hist/cor_pred_high_hist_subs.Rdata')
cor_pred_mod_hist <- readRDS(file = './sdm_output/cor/subs/habitat_pred/hist/cor_pred_mod_hist_subs.Rdata')
cor_pred_low_hist <- readRDS(file = './sdm_output/cor/subs/habitat_pred/hist/cor_pred_low_hist_subs.Rdata')

# SSP245
cor_pred_high_ssp245_30 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_high_ssp245_30_subs.Rdata')
cor_pred_mod_ssp245_30 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_mod_ssp245_30_subsRdata')
cor_pred_low_ssp245_30 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_low_ssp245_30_subs.Rdata')

cor_pred_high_ssp245_50 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_high_ssp245_50_subs.Rdata')
cor_pred_mod_ssp245_50 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_mod_ssp245_50_subs.Rdata')
cor_pred_low_ssp245_50 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_low_ssp245_50_subs.Rdata')

cor_pred_high_ssp245_70 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_high_ssp245_70_subs.Rdata')
cor_pred_mod_ssp245_70 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_mod_ssp245_70_subs.Rdata')
cor_pred_low_ssp245_70 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp245/cor_pred_low_ssp245_70_subs.Rdata')

# SSP585
cor_pred_high_ssp585_30 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_high_ssp585_30_subs.Rdata')
cor_pred_mod_ssp585_30 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_mod_ssp585_30_subs.Rdata')
cor_pred_low_ssp585_30 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_low_ssp585_30_subs.Rdata')

cor_pred_high_ssp585_50 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_high_ssp585_50_subs.Rdata')
cor_pred_mod_ssp585_50 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_mod_ssp585_50_subs.Rdata')
cor_pred_low_ssp585_50 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_low_ssp585_50_subs.Rdata')

cor_pred_high_ssp585_70 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_high_ssp585_70_subs.Rdata')
cor_pred_mod_ssp585_70 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_mod_ssp585_70_subs.Rdata')
cor_pred_low_ssp585_70 <- readRDS(file = './sdm_output/cor/subs/habitat_pred/ssp585/cor_pred_low_ssp585_70_subs.Rdata')


# M. coronaria - Crop rasters to study extent - TIFF ----------------------
# Crop categorical layers to restrict the extent of predicted suitability.
# In this case the model is making predictions outside of what ecologicaly makes sense!
cor_ext <- ext(-110, -52, 20, 70)

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

# Historical 
terra::writeRaster(cor_pred_high_hist_crop, "./sdm_output/cor/subs/cropped/hist/cor_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_hist_crop, "./sdm_output/cor/subs/cropped/hist/cor_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_hist_crop, "./sdm_output/cor/subs/cropped/hist/cor_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(cor_pred_high_ssp245_30_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp245_30_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp245_30_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp245_50_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp245_50_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp245_50_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp245_70_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp245_70_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp245_70_crop, "./sdm_output/cor/subs/cropped/ssp245/cor_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(cor_pred_high_ssp585_30_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp585_30_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp585_30_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp585_50_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp585_50_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp585_50_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(cor_pred_high_ssp585_70_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_mod_ssp585_70_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(cor_pred_low_ssp585_70_crop, "./sdm_output/cor/subs/cropped/ssp585/cor_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)



# M. fusca - MaxEnt Model -------------------------------------------------

occ_fus_coords <- as.data.frame(geom(occThin_fus)[,3:4]) # extract longitude, lattitude from occurence points
bg_fus_coords <- as.data.frame(geom(fus_bg_vec)[,3:4]) # extract longitude, lattitude from background points

cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

fus_maxent <- ENMevaluate(occ_fus_coords, # occurrence records
                          envs = wclim_fus_subs, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5),
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =    
                          list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, # leave one core available for other apps
                          algorithm = 'maxent.jar')

# Save the MaxEnt model so you do not have to waste time re-running the model
saveRDS(fus_maxent, file = './sdm_output/fus/subs/fus_maxent_subs.Rdata') # save

# Load
fus_maxent <- readRDS(file = './sdm_output/fus/subs/fus_maxent_subs.Rdata') # load 


# M. fusca Model Selection ------------------------------------------------
# Note that maxent results provide Continuous Boyce Index (cbi)
best_fus_maxent <- subset(fus_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_fus_maxent <- eval.models(fus_maxent)[[best_fus_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.2_fc.LQHPT

eval.variable.importance(fus_maxent)[best_fus_maxent["tune.args"][1, 1]]


# M. fusca historical prediction ------------------------------------------
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


# M. fusca Boyce Index ----------------------------------------------------
# Evaluate predictions using Boyce Index
png('C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/boyce_index/corboyce_plot_malus_fusca.png', width = 1600, height = 1200, res = 300)
fusBoyce <- ecospat.boyce(fit = fusPred_bg_val, # vector of predicted habitat suitability of bg points
              obs = fusPred_val_na, # vector of 
              nclass = 0, 
              PEplot = TRUE,
              method = 'spearman')

title(main = bquote(italic("Malus fusca") ~ ", Boyce Cor." ==
                       .(fusBoyce$cor))) 
dev.off()

# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
fusPred_val <- terra::extract(fus_pred_hist, occ_fus_coords)$lyr1
fusPred_threshold_1 <- quantile(fusPred_val, 0.01, na.rm = T) # Low suitability
fusPred_threshold_10 <- quantile(fusPred_val, 0.1, na.rm = T) # Moderate suitability
fusPred_threshold_50 <- quantile(fusPred_val, 0.5, na.rm = T) # High suitability


# M. fusca - Climate Predictions ------------------------------------------
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
terra::plot(fus_pred_hist > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " SSP5-8.5 Prediction: Early Century (2020-2040)")), background = 'lightskyblue1')
terra::plot(fus_pred_hist > fusPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_hist > fusPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
legend(x = -165, y = 45, xpd = NA, inset = c(1, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

terra::plot(fus_pred_ssp585_30 > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus fusca'), " SSP5-8.5 Prediction: Early Century (2020-2040) (SUBSETTED MODEL)")), background = 'lightskyblue1')
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30 > fusPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)



# Save/Load M fus. SDM predictions ----------------------------------------

# Save
saveRDS(fus_pred_hist, file = './sdm_output/fus/subs/fus_pred_hist_subs.Rdata')

saveRDS(fus_pred_ssp245_30, file = './sdm_output/fus/subs/fus_pred_ssp245_30_subs.Rdata')
saveRDS(fus_pred_ssp245_50, file = './sdm_output/fus/subs/fus_pred_ssp245_50_subs.Rdata')
saveRDS(fus_pred_ssp245_70, file = './sdm_output/fus/subs/fus_pred_ssp245_70_subs.Rdata')

saveRDS(fus_pred_ssp585_30, file = './sdm_output/fus/subs/fus_pred_ssp585_30_subs.Rdata')
saveRDS(fus_pred_ssp585_50, file = './sdm_output/fus/subs/fus_pred_ssp585_50_subs.Rdata')
saveRDS(fus_pred_ssp585_70, file = './sdm_output/fus/subs/fus_pred_ssp585_70_subs.Rdata')

# Load
#Subsetted models
fus_pred_hist <- readRDS(file = './sdm_output/fus/subs/fus_pred_hist_subs.Rdata')

fus_pred_ssp245_30 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp245_30_subs.Rdata')
fus_pred_ssp245_50 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp245_50_subs.Rdata')
fus_pred_ssp245_70 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp245_70_subs.Rdata')

fus_pred_ssp585_30 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_30_subs.Rdata')
fus_pred_ssp585_50 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_50_subs.Rdata')
fus_pred_ssp585_70 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_70_subs.Rdata')

# Save M. fusca thresholds ------------------------------------------------

saveRDS(fusPred_threshold_1, file = './sdm_output/fus/subs/threshold/fusPred_threshold_1_subs.Rdata')
saveRDS(fusPred_threshold_10, file = './sdm_output/fus/subs/threshold/fusPred_threshold_10_subs.Rdata')
saveRDS(fusPred_threshold_50, file = './sdm_output/fus/subs/threshold/fusPred_threshold_50_subs.Rdata')


# Load subsetted
fusPred_threshold_1 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_1_subs.Rdata')
fusPred_threshold_10 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_10_subs.Rdata')
fusPred_threshold_50 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_50_subs.Rdata')

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
saveRDS(fus_pred_high_hist, file = './sdm_output/fus/subs/habitat_pred/hist/fus_pred_high_hist_subs.Rdata')
saveRDS(fus_pred_mod_hist, file = './sdm_output/fus/subs/habitat_pred/hist/fus_pred_mod_hist_subs.Rdata')
saveRDS(fus_pred_low_hist, file = './sdm_output/fus/subs/habitat_pred/hist/fus_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(fus_pred_high_ssp245_30, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_high_ssp245_30_subs.Rdata')
saveRDS(fus_pred_mod_ssp245_30, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_mod_ssp245_30_subs.Rdata')
saveRDS(fus_pred_low_ssp245_30, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_low_ssp245_30_subs.Rdata')

saveRDS(fus_pred_high_ssp245_50, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_high_ssp245_50_subs.Rdata')
saveRDS(fus_pred_mod_ssp245_50, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_mod_ssp245_50_subs.Rdata')
saveRDS(fus_pred_low_ssp245_50, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_low_ssp245_50_subs.Rdata')

saveRDS(fus_pred_high_ssp245_70, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_high_ssp245_70_subs.Rdata')
saveRDS(fus_pred_mod_ssp245_70, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_mod_ssp245_70_subs.Rdata')
saveRDS(fus_pred_low_ssp245_70, file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(fus_pred_high_ssp585_30, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_high_ssp585_30_subs.Rdata')
saveRDS(fus_pred_mod_ssp585_30, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_mod_ssp585_30_subs.Rdata')
saveRDS(fus_pred_low_ssp585_30, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_low_ssp585_30_subs.Rdata')

saveRDS(fus_pred_high_ssp585_50, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_high_ssp585_50_subs.Rdata')
saveRDS(fus_pred_mod_ssp585_50, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_mod_ssp585_50_subs.Rdata')
saveRDS(fus_pred_low_ssp585_50, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_low_ssp585_50_subs.Rdata')

saveRDS(fus_pred_high_ssp585_70, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_high_ssp585_70_subs.Rdata')
saveRDS(fus_pred_mod_ssp585_70, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_mod_ssp585_70_subs.Rdata')
saveRDS(fus_pred_low_ssp585_70, file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_low_ssp585_70_subs.Rdata')

# Load


# Historical
fus_pred_high_hist <- readRDS(file = './sdm_output/fus/subs/habitat_pred/hist/fus_pred_high_hist_subs.Rdata')
fus_pred_mod_hist <- readRDS(file = './sdm_output/fus/subs/habitat_pred/hist/fus_pred_mod_hist_subs.Rdata')
fus_pred_low_hist <- readRDS(file = './sdm_output/fus/subs/habitat_pred/hist/fus_pred_low_hist_subs.Rdata')

# SSP245
fus_pred_high_ssp245_30 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_high_ssp245_30_subs.Rdata')
fus_pred_mod_ssp245_30 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_mod_ssp245_30_subs.Rdata')
fus_pred_low_ssp245_30 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_low_ssp245_30_subs.Rdata')

fus_pred_high_ssp245_50 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_high_ssp245_50_subs.Rdata')
fus_pred_mod_ssp245_50 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_mod_ssp245_50_subs.Rdata')
fus_pred_low_ssp245_50 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_low_ssp245_50_subs.Rdata')

fus_pred_high_ssp245_70 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_high_ssp245_70_subs.Rdata')
fus_pred_mod_ssp245_70 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_mod_ssp245_70_subs.Rdata')
fus_pred_low_ssp245_70 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp245/fus_pred_low_ssp245_70_subs.Rdata')

# SSP585
fus_pred_high_ssp585_30 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_high_ssp585_30_subs.Rdata')
fus_pred_mod_ssp585_30 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_mod_ssp585_30_subs.Rdata')
fus_pred_low_ssp585_30 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_low_ssp585_30_subs.Rdata')

fus_pred_high_ssp585_50 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_high_ssp585_50_subs.Rdata')
fus_pred_mod_ssp585_50 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_mod_ssp585_50_subs.Rdata')
fus_pred_low_ssp585_50 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_low_ssp585_50_subs.Rdata')

fus_pred_high_ssp585_70 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_high_ssp585_70_subs.Rdata')
fus_pred_mod_ssp585_70 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_mod_ssp585_70_subs.Rdata')
fus_pred_low_ssp585_70 <- readRDS(file = './sdm_output/fus/subs/habitat_pred/ssp585/fus_pred_low_ssp585_70_subs.Rdata')


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


# M. fusca - Crop rasters - TIFF ------------------------------------------


# Historical 
terra::writeRaster(fus_pred_high_hist_crop, "./sdm_output/fus/subs/cropped/hist/fus_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_hist_crop, "./sdm_output/fus/subs/cropped/hist/fus_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_hist_crop, "./sdm_output/fus/subs/cropped/hist/fus_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(fus_pred_high_ssp245_30_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp245_30_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp245_30_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp245_50_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp245_50_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp245_50_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp245_70_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp245_70_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp245_70_crop, "./sdm_output/fus/subs/cropped/ssp245/fus_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(fus_pred_high_ssp585_30_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp585_30_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp585_30_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp585_50_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp585_50_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp585_50_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(fus_pred_high_ssp585_70_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_mod_ssp585_70_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(fus_pred_low_ssp585_70_crop, "./sdm_output/fus/subs/cropped/ssp585/fus_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

# M. ioensis - Maxent Model -----------------------------------------------

occ_ion_coords <- as.data.frame(geom(occThin_ion)[,3:4]) # extract longitude, lattitude from occurence points
bg_ion_coords <- as.data.frame(geom(ion_bg_vec)[,3:4]) # extract longitude, lattitude from background points

set.seed(1337)

ion_maxent <- ENMevaluate(occ_ion_coords, # occurrence records
                          envs = wclim_ion_subs, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5),
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, # leave one core available for other apps
                          algorithm = 'maxent.jar')

# Save the MaxEnt model so you do not have to waste time re-running the model
saveRDS(ion_maxent, file = './sdm_output/ion/subs/ion_maxent_subs.Rdata') # save

#Load
ion_maxent <- readRDS(file = './sdm_output/ion/subs/ion_maxent_subs.Rdata') # subset model


# M. ionesis Model Selection ----------------------------------------------

# Note that maxent results provide Continuous Boyce Index (cbi)
best_ion_maxent <- subset(ion_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_ion_maxent <- eval.models(ion_maxent)[[best_ion_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.2_fc.LQHPT
  
eval.variable.importance(ion_maxent)[best_ion_maxent["tune.args"][1, 1]]


# M. ioensis historical prediction ----------------------------------------
# Now use the <terra> package to plot the SDM prediction.
# Wclim is the historical climatic conditions (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
ion_pred_hist <- terra::predict(wclim_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)

plot(ion_pred_hist)
points(occThin_ion, cex = 0.05)

# Evaluate predictions using Boyce Index
# the number of true presences should decline with suitability groups 100-91, 90-81, etc. 
# First extract suitability values for the background and presence points, make sure to omit NA values
ionPred_bg_val <- terra::extract(ion_pred_hist, bg_ion_coords)$lyr1 %>% 
  na.omit()

ionPred_val_na <- terra::extract(ion_pred_hist, occ_ion_coords)$lyr1 %>% 
  na.omit()



# M. ioensis Boyce Index --------------------------------------------------
# Evaluate predictions using Boyce Index
png('C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/boyce_index/corboyce_plot_malus_ioensis.png', width = 1600, height = 1200, res = 300)
ionBoyce <- ecospat.boyce(fit = ionPred_bg_val, # vector of predicted habitat suitability of bg points
              obs = ionPred_val_na, # vector of 
              nclass = 0, 
              PEplot = TRUE,
              method = 'spearman')

title(main = bquote(italic("Malus ioensis") ~ ", Boyce Cor." ==
                      .(ionBoyce$cor))) 

dev.off()

# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
ionPred_val <- terra::extract(ion_pred_hist, occ_ion_coords)$lyr1
ionPred_threshold_1 <- quantile(ionPred_val, 0.01, na.rm = T) # Low suitability
ionPred_threshold_10 <- quantile(ionPred_val, 0.1, na.rm = T) # Moderate suitability
ionPred_threshold_50 <- quantile(ionPred_val, 0.5, na.rm = T) # High suitability



# M. ioensis- Climate Predictions -----------------------------------------


# Future Climate predictions
# SSP 245
ion_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)

# SSP 585
ion_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)
ion_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_ion_maxent, cores = cn - 1, na.rm = T)

# Plot SSP 585 2030
dev.off()
dev.new()
par(mar = c(4, 4, 4, 4), mfcol = c(1, 2))
terra::plot(ion_pred_hist > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus ionca'), " SSP5-8.5 Prediction: Early Century (2020-2040)")), background = 'lightskyblue1')
terra::plot(ion_pred_hist > ionPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_hist > ionPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
legend(x = -165, y = 45, xpd = NA, inset = c(1, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)

terra::plot(ion_pred_ssp585_30 > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, xlim = c(-170, -110), ylim = c(30, 65), main = expression(atop(italic('Malus ionca'), " SSP5-8.5 Prediction: Early Century (2020-2040) (SUBSETTED MODEL)")), background = 'lightskyblue1')
terra::plot(ion_pred_ssp585_30 > ionPred_threshold_10, col = c(rgb(1, 1, 1, alpha=0), '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp585_30 > ionPred_threshold_50, col = c(rgb(1, 1, 1, alpha=0), '#D95F0E'), add = T, legend = F)
legend(x = -165, y = 45, xpd = NA, inset = c(5, 0), 
       title = 'Habitat Suitability', 
       legend = legend_labs,
       fill = fill_cols)



# Save/Load M ion. SDM predictions ----------------------------------------

# Save
saveRDS(ion_pred_hist, file = './sdm_output/ion/subs/ion_pred_hist_subs.Rdata')

saveRDS(ion_pred_ssp245_30, file = './sdm_output/ion/subs/ion_pred_ssp245_30_subs.Rdata')
saveRDS(ion_pred_ssp245_50, file = './sdm_output/ion/subs/ion_pred_ssp245_50_subs.Rdata')
saveRDS(ion_pred_ssp245_70, file = './sdm_output/ion/subs/ion_pred_ssp245_70_subs.Rdata')

saveRDS(ion_pred_ssp585_30, file = './sdm_output/ion/subs/ion_pred_ssp585_30_subs.Rdata')
saveRDS(ion_pred_ssp585_50, file = './sdm_output/ion/subs/ion_pred_ssp585_50_subs.Rdata')
saveRDS(ion_pred_ssp585_70, file = './sdm_output/ion/subs/ion_pred_ssp585_70_subs.Rdata')

# Load
#Subsetted models
ion_pred_hist <- readRDS(file = './sdm_output/ion/subs/ion_pred_hist_subs.Rdata')

ion_pred_ssp245_30 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp245_30_subs.Rdata')
ion_pred_ssp245_50 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp245_50_subs.Rdata')
ion_pred_ssp245_70 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp245_70_subs.Rdata')

ion_pred_ssp585_30 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_30_subs.Rdata')
ion_pred_ssp585_50 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_50_subs.Rdata')
ion_pred_ssp585_70 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_70_subs.Rdata')


# Save M. ioensis thresholds ----------------------------------------------

saveRDS(ionPred_threshold_1, file = './sdm_output/ion/subs/threshold/ionPred_threshold_1_subs.Rdata')
saveRDS(ionPred_threshold_10, file = './sdm_output/ion/subs/threshold/ionPred_threshold_10_subs.Rdata')
saveRDS(ionPred_threshold_50, file = './sdm_output/ion/subs/threshold/ionPred_threshold_50_subs.Rdata')


# Load subsetted
ionPred_threshold_1 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_1_subs.Rdata')
ionPred_threshold_10 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_10_subs.Rdata')
ionPred_threshold_50 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_50_subs.Rdata')


# M ioensis Habitat predictions -------------------------------------------
# Categorical habitat suitability
# Historical
ion_pred_high_hist <- ion_pred_hist > ionPred_threshold_50
ion_pred_mod_hist <- ion_pred_hist > ionPred_threshold_10
ion_pred_low_hist <- ion_pred_hist > ionPred_threshold_1

#SSP245 
ion_pred_high_ssp245_30 <- ion_pred_ssp245_30 > ionPred_threshold_50
ion_pred_mod_ssp245_30 <- ion_pred_ssp245_30 > ionPred_threshold_10
ion_pred_low_ssp245_30 <- ion_pred_ssp245_30 > ionPred_threshold_1

ion_pred_high_ssp245_50 <- ion_pred_ssp245_50 > ionPred_threshold_50
ion_pred_mod_ssp245_50 <- ion_pred_ssp245_50 > ionPred_threshold_10
ion_pred_low_ssp245_50 <- ion_pred_ssp245_50 > ionPred_threshold_1

ion_pred_high_ssp245_70 <- ion_pred_ssp245_70 > ionPred_threshold_50
ion_pred_mod_ssp245_70 <- ion_pred_ssp245_70 > ionPred_threshold_10
ion_pred_low_ssp245_70 <- ion_pred_ssp245_70 > ionPred_threshold_1

#SSP585
ion_pred_high_ssp585_30 <- ion_pred_ssp585_30 > ionPred_threshold_50
ion_pred_mod_ssp585_30 <- ion_pred_ssp585_30 > ionPred_threshold_10
ion_pred_low_ssp585_30 <- ion_pred_ssp585_30 > ionPred_threshold_1

ion_pred_high_ssp585_50 <- ion_pred_ssp585_50 > ionPred_threshold_50
ion_pred_mod_ssp585_50 <- ion_pred_ssp585_50 > ionPred_threshold_10
ion_pred_low_ssp585_50 <- ion_pred_ssp585_50 > ionPred_threshold_1

ion_pred_high_ssp585_70 <- ion_pred_ssp585_70 > ionPred_threshold_50
ion_pred_mod_ssp585_70 <- ion_pred_ssp585_70 > ionPred_threshold_10
ion_pred_low_ssp585_70 <- ion_pred_ssp585_70 > ionPred_threshold_1

# Save

# Historical
saveRDS(ion_pred_high_hist, file = './sdm_output/ion/subs/habitat_pred/hist/ion_pred_high_hist_subs.Rdata')
saveRDS(ion_pred_mod_hist, file = './sdm_output/ion/subs/habitat_pred/hist/ion_pred_mod_hist_subs.Rdata')
saveRDS(ion_pred_low_hist, file = './sdm_output/ion/subs/habitat_pred/hist/ion_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(ion_pred_high_ssp245_30, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_high_ssp245_30_subs.Rdata')
saveRDS(ion_pred_mod_ssp245_30, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_mod_ssp245_30_subs.Rdata')
saveRDS(ion_pred_low_ssp245_30, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_low_ssp245_30_subs.Rdata')

saveRDS(ion_pred_high_ssp245_50, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_high_ssp245_50_subs.Rdata')
saveRDS(ion_pred_mod_ssp245_50, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_mod_ssp245_50_subs.Rdata')
saveRDS(ion_pred_low_ssp245_50, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_low_ssp245_50_subs.Rdata')

saveRDS(ion_pred_high_ssp245_70, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_high_ssp245_70_subs.Rdata')
saveRDS(ion_pred_mod_ssp245_70, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_mod_ssp245_70_subs.Rdata')
saveRDS(ion_pred_low_ssp245_70, file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(ion_pred_high_ssp585_30, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_high_ssp585_30_subs.Rdata')
saveRDS(ion_pred_mod_ssp585_30, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_mod_ssp585_30_subs.Rdata')
saveRDS(ion_pred_low_ssp585_30, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_low_ssp585_30_subs.Rdata')

saveRDS(ion_pred_high_ssp585_50, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_high_ssp585_50_subs.Rdata')
saveRDS(ion_pred_mod_ssp585_50, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_mod_ssp585_50_subs.Rdata')
saveRDS(ion_pred_low_ssp585_50, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_low_ssp585_50_subs.Rdata')

saveRDS(ion_pred_high_ssp585_70, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_high_ssp585_70_subs.Rdata')
saveRDS(ion_pred_mod_ssp585_70, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_mod_ssp585_70_subs.Rdata')
saveRDS(ion_pred_low_ssp585_70, file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_low_ssp585_70_subs.Rdata')

# Load

# Historical
ion_pred_high_hist <- readRDS(file = './sdm_output/ion/subs/habitat_pred/hist/ion_pred_high_hist_subs.Rdata')
ion_pred_mod_hist <- readRDS(file = './sdm_output/ion/subs/habitat_pred/hist/ion_pred_mod_hist_subs.Rdata')
ion_pred_low_hist <- readRDS(file = './sdm_output/ion/subs/habitat_pred/hist/ion_pred_low_hist_subs.Rdata')

# SSP245
ion_pred_high_ssp245_30 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_high_ssp245_30_subs.Rdata')
ion_pred_mod_ssp245_30 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_mod_ssp245_30_subs.Rdata')
ion_pred_low_ssp245_30 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_low_ssp245_30_subs.Rdata')

ion_pred_high_ssp245_50 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_high_ssp245_50_subs.Rdata')
ion_pred_mod_ssp245_50 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_mod_ssp245_50_subs.Rdata')
ion_pred_low_ssp245_50 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_low_ssp245_50_subs.Rdata')

ion_pred_high_ssp245_70 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_high_ssp245_70_subs.Rdata')
ion_pred_mod_ssp245_70 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_mod_ssp245_70_subs.Rdata')
ion_pred_low_ssp245_70 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp245/ion_pred_low_ssp245_70_subs.Rdata')

# SSP585
ion_pred_high_ssp585_30 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_high_ssp585_30_subs.Rdata')
ion_pred_mod_ssp585_30 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_mod_ssp585_30_subs.Rdata')
ion_pred_low_ssp585_30 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_low_ssp585_30_subs.Rdata')

ion_pred_high_ssp585_50 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_high_ssp585_50_subs.Rdata')
ion_pred_mod_ssp585_50 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_mod_ssp585_50_subs.Rdata')
ion_pred_low_ssp585_50 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_low_ssp585_50_subs.Rdata')

ion_pred_high_ssp585_70 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_high_ssp585_70_subs.Rdata')
ion_pred_mod_ssp585_70 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_mod_ssp585_70_subs.Rdata')
ion_pred_low_ssp585_70 <- readRDS(file = './sdm_output/ion/subs/habitat_pred/ssp585/ion_pred_low_ssp585_70_subs.Rdata')


# Crop categorical layers to restrict the extent of predicted suitability.
# In this case the model is making predictions outside of what ecologicaly makes sense!
ion_ext <- ext(-110, -52, 20, 70)

# Historical
ion_pred_high_hist_crop <- crop(ion_pred_high_hist, ion_ext)
ion_pred_mod_hist_crop <- crop(ion_pred_mod_hist, ion_ext)
ion_pred_low_hist_crop <- crop(ion_pred_low_hist, ion_ext)

#SSP245 
ion_pred_high_ssp245_30_crop <- crop(ion_pred_high_ssp245_30, ion_ext)
ion_pred_mod_ssp245_30_crop <- crop(ion_pred_mod_ssp245_30, ion_ext)
ion_pred_low_ssp245_30_crop <- crop(ion_pred_low_ssp245_30, ion_ext)

ion_pred_high_ssp245_50_crop <- crop(ion_pred_high_ssp245_50, ion_ext)
ion_pred_mod_ssp245_50_crop <- crop(ion_pred_mod_ssp245_50, ion_ext)
ion_pred_low_ssp245_50_crop <- crop(ion_pred_low_ssp245_50, ion_ext)

ion_pred_high_ssp245_70_crop <- crop(ion_pred_high_ssp245_70, ion_ext)
ion_pred_mod_ssp245_70_crop <- crop(ion_pred_mod_ssp245_70, ion_ext)
ion_pred_low_ssp245_70_crop <- crop(ion_pred_low_ssp245_70, ion_ext)

#SSP585
ion_pred_high_ssp585_30_crop <- crop(ion_pred_high_ssp585_30, ion_ext)
ion_pred_mod_ssp585_30_crop <- crop(ion_pred_mod_ssp585_30, ion_ext)
ion_pred_low_ssp585_30_crop <- crop(ion_pred_low_ssp585_30, ion_ext)

ion_pred_high_ssp585_50_crop <- crop(ion_pred_high_ssp585_50, ion_ext)
ion_pred_mod_ssp585_50_crop <- crop(ion_pred_mod_ssp585_50, ion_ext)
ion_pred_low_ssp585_50_crop <- crop(ion_pred_low_ssp585_50, ion_ext)

ion_pred_high_ssp585_70_crop <- crop(ion_pred_high_ssp585_70, ion_ext)
ion_pred_mod_ssp585_70_crop <- crop(ion_pred_mod_ssp585_70, ion_ext)
ion_pred_low_ssp585_70_crop <- crop(ion_pred_low_ssp585_70, ion_ext)


# M. ionesis- Crop rasters - TIFF -----------------------------------------

# Historical 
terra::writeRaster(ion_pred_high_hist_crop, "./sdm_output/ion/subs/cropped/hist/ion_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_hist_crop, "./sdm_output/ion/subs/cropped/hist/ion_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_hist_crop, "./sdm_output/ion/subs/cropped/hist/ion_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(ion_pred_high_ssp245_30_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_ssp245_30_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_ssp245_30_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ion_pred_high_ssp245_50_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_ssp245_50_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_ssp245_50_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ion_pred_high_ssp245_70_crop$lyr1, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_ssp245_70_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_ssp245_70_crop, "./sdm_output/ion/subs/cropped/ssp245/ion_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(ion_pred_high_ssp585_30_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_ssp585_30_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_ssp585_30_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ion_pred_high_ssp585_50_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_ssp585_50_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_ssp585_50_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ion_pred_high_ssp585_70_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_mod_ssp585_70_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ion_pred_low_ssp585_70_crop, "./sdm_output/ion/subs/cropped/ssp585/ion_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

# M. angustifolia - Maxent Model ------------------------------------------

occ_ang_coords <- as.data.frame(geom(occThin_ang)[,3:4]) # extract longitude, lattitude from occurence points
bg_ang_coords <- as.data.frame(geom(ang_bg_vec)[,3:4]) # extract longitude, lattitude from background points

cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

ang_maxent <- ENMevaluate(occ_ang_coords, # occurrence records
                          envs = wclim_ang_subs, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5),
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, # leave one core available for other apps
                          algorithm = 'maxent.jar')

# Save the MaxEnt model so you do not have to waste time re-running the model
saveRDS(ang_maxent, file = './sdm_output/ang/subs/ang_maxent_subs.Rdata') # save


ang_maxent <- readRDS(file = './sdm_output/ang/subs/ang_maxent_subs.Rdata') # subset model 


# M. angustifolia Model Selecting -----------------------------------------

# Note that maxent results provide Continuous Boyce Index (cbi)
best_ang_maxent <- subset(ang_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_ang_maxent <- eval.models(ang_maxent)[[best_ang_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.2_fc.LQHPT

eval.variable.importance(ang_maxent)[best_ang_maxent["tune.args"][1, 1]]



# M. angustifolia historical prediction -----------------------------------

# Now use the <terra> package to plot the SDM predictang.
# Wclim is the historical climatic conditangs (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
ang_pred_hist <- terra::predict(wclim_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)

plot(ang_pred_hist)
points(occThin_ang, cex = 0.05)

# Evaluate predictangs using Boyce Index
# the number of true presences should decline with suitability groups 100-91, 90-81, etc. 
# First extract suitability values for the background and presence points, make sure to omit NA values
angPred_bg_val <- terra::extract(ang_pred_hist, bg_ang_coords)$lyr1 %>% 
  na.omit()

angPred_val_na <- terra::extract(ang_pred_hist, occ_ang_coords)$lyr1 %>% 
  na.omit()

# M. angustifolia Boyce Index ---------------------------------------------

# Evaluate predictangs using Boyce Index
png('C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/boyce_index/corboyce_plot_malus_angustifolia.png', width = 1600, height = 1200, res = 300)
angBoyce <- ecospat.boyce(fit = angPred_bg_val, # vector of predicted habitat suitability of bg points
              obs = angPred_val_na, # vector of 
              nclass = 0, 
              PEplot = TRUE,
              method = 'spearman')

title(main = bquote(italic("Malus angustifolia") ~ ", Boyce Cor." ==
                      .(angBoyce$cor))) 
dev.off()

# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
angPred_val <- terra::extract(ang_pred_hist, occ_ang_coords)$lyr1
angPred_threshold_1 <- quantile(angPred_val, 0.01, na.rm = T) # Low suitability
angPred_threshold_10 <- quantile(angPred_val, 0.1, na.rm = T) # Moderate suitability
angPred_threshold_50 <- quantile(angPred_val, 0.5, na.rm = T) # High suitability

# M. angustifolia - Climate Predictions -----------------------------------


# Future Climate predictangs
# SSP 245
ang_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)

# SSP 585
ang_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)
ang_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_ang_maxent, cores = cn - 1, na.rm = T)


# Save/Load M ang. SDM predictangs ----------------------------------------

# Save
saveRDS(ang_pred_hist, file = './sdm_output/ang/subs/ang_pred_hist_subs.Rdata')

saveRDS(ang_pred_ssp245_30, file = './sdm_output/ang/subs/ang_pred_ssp245_30_subs.Rdata')
saveRDS(ang_pred_ssp245_50, file = './sdm_output/ang/subs/ang_pred_ssp245_50_subs.Rdata')
saveRDS(ang_pred_ssp245_70, file = './sdm_output/ang/subs/ang_pred_ssp245_70_subs.Rdata')

saveRDS(ang_pred_ssp585_30, file = './sdm_output/ang/subs/ang_pred_ssp585_30_subs.Rdata')
saveRDS(ang_pred_ssp585_50, file = './sdm_output/ang/subs/ang_pred_ssp585_50_subs.Rdata')
saveRDS(ang_pred_ssp585_70, file = './sdm_output/ang/subs/ang_pred_ssp585_70_subs.Rdata')

# Load
#Subsetted models
ang_pred_hist <- readRDS(file = './sdm_output/ang/subs/ang_pred_hist_subs.Rdata')

ang_pred_ssp245_30 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp245_30_subs.Rdata')
ang_pred_ssp245_50 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp245_50_subs.Rdata')
ang_pred_ssp245_70 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp245_70_subs.Rdata')

ang_pred_ssp585_30 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_30_subs.Rdata')
ang_pred_ssp585_50 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_50_subs.Rdata')
ang_pred_ssp585_70 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_70_subs.Rdata')



# Save M. angustifolia thresholds -----------------------------------------


saveRDS(angPred_threshold_1, file = './sdm_output/ang/subs/threshold/angPred_threshold_1_subs.Rdata')
saveRDS(angPred_threshold_10, file = './sdm_output/ang/subs/threshold/angPred_threshold_10_subs.Rdata')
saveRDS(angPred_threshold_50, file = './sdm_output/ang/subs/threshold/angPred_threshold_50_subs.Rdata')


# Load subsetted
angPred_threshold_1 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_1_subs.Rdata')
angPred_threshold_10 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_10_subs.Rdata')
angPred_threshold_50 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_50_subs.Rdata')


# M angustifolia Habitat predictions --------------------------------------
# Categorical habitat suitability
# Historical
ang_pred_high_hist <- ang_pred_hist > angPred_threshold_50
ang_pred_mod_hist <- ang_pred_hist > angPred_threshold_10
ang_pred_low_hist <- ang_pred_hist > angPred_threshold_1

#SSP245 
ang_pred_high_ssp245_30 <- ang_pred_ssp245_30 > angPred_threshold_50
ang_pred_mod_ssp245_30 <- ang_pred_ssp245_30 > angPred_threshold_10
ang_pred_low_ssp245_30 <- ang_pred_ssp245_30 > angPred_threshold_1

ang_pred_high_ssp245_50 <- ang_pred_ssp245_50 > angPred_threshold_50
ang_pred_mod_ssp245_50 <- ang_pred_ssp245_50 > angPred_threshold_10
ang_pred_low_ssp245_50 <- ang_pred_ssp245_50 > angPred_threshold_1

ang_pred_high_ssp245_70 <- ang_pred_ssp245_70 > angPred_threshold_50
ang_pred_mod_ssp245_70 <- ang_pred_ssp245_70 > angPred_threshold_10
ang_pred_low_ssp245_70 <- ang_pred_ssp245_70 > angPred_threshold_1

#SSP585
ang_pred_high_ssp585_30 <- ang_pred_ssp585_30 > angPred_threshold_50
ang_pred_mod_ssp585_30 <- ang_pred_ssp585_30 > angPred_threshold_10
ang_pred_low_ssp585_30 <- ang_pred_ssp585_30 > angPred_threshold_1

ang_pred_high_ssp585_50 <- ang_pred_ssp585_50 > angPred_threshold_50
ang_pred_mod_ssp585_50 <- ang_pred_ssp585_50 > angPred_threshold_10
ang_pred_low_ssp585_50 <- ang_pred_ssp585_50 > angPred_threshold_1

ang_pred_high_ssp585_70 <- ang_pred_ssp585_70 > angPred_threshold_50
ang_pred_mod_ssp585_70 <- ang_pred_ssp585_70 > angPred_threshold_10
ang_pred_low_ssp585_70 <- ang_pred_ssp585_70 > angPred_threshold_1

# Save

# Historical
saveRDS(ang_pred_high_hist, file = './sdm_output/ang/subs/habitat_pred/hist/ang_pred_high_hist_subs.Rdata')
saveRDS(ang_pred_mod_hist, file = './sdm_output/ang/subs/habitat_pred/hist/ang_pred_mod_hist_subs.Rdata')
saveRDS(ang_pred_low_hist, file = './sdm_output/ang/subs/habitat_pred/hist/ang_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(ang_pred_high_ssp245_30, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_high_ssp245_30_subs.Rdata')
saveRDS(ang_pred_mod_ssp245_30, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_mod_ssp245_30_subs.Rdata')
saveRDS(ang_pred_low_ssp245_30, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_low_ssp245_30_subs.Rdata')

saveRDS(ang_pred_high_ssp245_50, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_high_ssp245_50_subs.Rdata')
saveRDS(ang_pred_mod_ssp245_50, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_mod_ssp245_50_subs.Rdata')
saveRDS(ang_pred_low_ssp245_50, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_low_ssp245_50_subs.Rdata')

saveRDS(ang_pred_high_ssp245_70, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_high_ssp245_70_subs.Rdata')
saveRDS(ang_pred_mod_ssp245_70, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_mod_ssp245_70_subs.Rdata')
saveRDS(ang_pred_low_ssp245_70, file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(ang_pred_high_ssp585_30, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_high_ssp585_30_subs.Rdata')
saveRDS(ang_pred_mod_ssp585_30, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_mod_ssp585_30_subs.Rdata')
saveRDS(ang_pred_low_ssp585_30, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_low_ssp585_30_subs.Rdata')

saveRDS(ang_pred_high_ssp585_50, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_high_ssp585_50_subs.Rdata')
saveRDS(ang_pred_mod_ssp585_50, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_mod_ssp585_50_subs.Rdata')
saveRDS(ang_pred_low_ssp585_50, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_low_ssp585_50_subs.Rdata')

saveRDS(ang_pred_high_ssp585_70, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_high_ssp585_70_subs.Rdata')
saveRDS(ang_pred_mod_ssp585_70, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_mod_ssp585_70_subs.Rdata')
saveRDS(ang_pred_low_ssp585_70, file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_low_ssp585_70_subs.Rdata')

# Load

# Historical
ang_pred_high_hist <- readRDS(file = './sdm_output/ang/subs/habitat_pred/hist/ang_pred_high_hist_subs.Rdata')
ang_pred_mod_hist <- readRDS(file = './sdm_output/ang/subs/habitat_pred/hist/ang_pred_mod_hist_subs.Rdata')
ang_pred_low_hist <- readRDS(file = './sdm_output/ang/subs/habitat_pred/hist/ang_pred_low_hist_subs.Rdata')

# SSP245
ang_pred_high_ssp245_30 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_high_ssp245_30_subs.Rdata')
ang_pred_mod_ssp245_30 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_mod_ssp245_30_subs.Rdata')
ang_pred_low_ssp245_30 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_low_ssp245_30_subs.Rdata')

ang_pred_high_ssp245_50 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_high_ssp245_50_subs.Rdata')
ang_pred_mod_ssp245_50 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_mod_ssp245_50_subs.Rdata')
ang_pred_low_ssp245_50 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_low_ssp245_50_subs.Rdata')

ang_pred_high_ssp245_70 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_high_ssp245_70_subs.Rdata')
ang_pred_mod_ssp245_70 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_mod_ssp245_70_subs.Rdata')
ang_pred_low_ssp245_70 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp245/ang_pred_low_ssp245_70_subs.Rdata')

# SSP585
ang_pred_high_ssp585_30 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_high_ssp585_30_subs.Rdata')
ang_pred_mod_ssp585_30 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_mod_ssp585_30_subs.Rdata')
ang_pred_low_ssp585_30 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_low_ssp585_30_subs.Rdata')

ang_pred_high_ssp585_50 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_high_ssp585_50_subs.Rdata')
ang_pred_mod_ssp585_50 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_mod_ssp585_50_subs.Rdata')
ang_pred_low_ssp585_50 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_low_ssp585_50_subs.Rdata')

ang_pred_high_ssp585_70 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_high_ssp585_70_subs.Rdata')
ang_pred_mod_ssp585_70 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_mod_ssp585_70_subs.Rdata')
ang_pred_low_ssp585_70 <- readRDS(file = './sdm_output/ang/subs/habitat_pred/ssp585/ang_pred_low_ssp585_70_subs.Rdata')


# Crop categorical layers to restrict the extent of predicted suitability.
# In this case the model is making predictangs outside of what ecologicaly makes sense!
ang_ext <- ext(-110, -52, 20, 70)

# Historical
ang_pred_high_hist_crop <- crop(ang_pred_high_hist, ang_ext)
ang_pred_mod_hist_crop <- crop(ang_pred_mod_hist, ang_ext)
ang_pred_low_hist_crop <- crop(ang_pred_low_hist, ang_ext)

#SSP245 
ang_pred_high_ssp245_30_crop <- crop(ang_pred_high_ssp245_30, ang_ext)
ang_pred_mod_ssp245_30_crop <- crop(ang_pred_mod_ssp245_30, ang_ext)
ang_pred_low_ssp245_30_crop <- crop(ang_pred_low_ssp245_30, ang_ext)

ang_pred_high_ssp245_50_crop <- crop(ang_pred_high_ssp245_50, ang_ext)
ang_pred_mod_ssp245_50_crop <- crop(ang_pred_mod_ssp245_50, ang_ext)
ang_pred_low_ssp245_50_crop <- crop(ang_pred_low_ssp245_50, ang_ext)

ang_pred_high_ssp245_70_crop <- crop(ang_pred_high_ssp245_70, ang_ext)
ang_pred_mod_ssp245_70_crop <- crop(ang_pred_mod_ssp245_70, ang_ext)
ang_pred_low_ssp245_70_crop <- crop(ang_pred_low_ssp245_70, ang_ext)

#SSP585
ang_pred_high_ssp585_30_crop <- crop(ang_pred_high_ssp585_30, ang_ext)
ang_pred_mod_ssp585_30_crop <- crop(ang_pred_mod_ssp585_30, ang_ext)
ang_pred_low_ssp585_30_crop <- crop(ang_pred_low_ssp585_30, ang_ext)

ang_pred_high_ssp585_50_crop <- crop(ang_pred_high_ssp585_50, ang_ext)
ang_pred_mod_ssp585_50_crop <- crop(ang_pred_mod_ssp585_50, ang_ext)
ang_pred_low_ssp585_50_crop <- crop(ang_pred_low_ssp585_50, ang_ext)

ang_pred_high_ssp585_70_crop <- crop(ang_pred_high_ssp585_70, ang_ext)
ang_pred_mod_ssp585_70_crop <- crop(ang_pred_mod_ssp585_70, ang_ext)
ang_pred_low_ssp585_70_crop <- crop(ang_pred_low_ssp585_70, ang_ext)


#  M. angustifolia - Crop rasters - TIFF  ---------------------------------

# Historical 
terra::writeRaster(ang_pred_high_hist_crop, "./sdm_output/ang/subs/cropped/hist/ang_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_hist_crop, "./sdm_output/ang/subs/cropped/hist/ang_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_hist_crop, "./sdm_output/ang/subs/cropped/hist/ang_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(ang_pred_high_ssp245_30_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_ssp245_30_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_ssp245_30_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ang_pred_high_ssp245_50_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_ssp245_50_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_ssp245_50_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ang_pred_high_ssp245_70_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_ssp245_70_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_ssp245_70_crop, "./sdm_output/ang/subs/cropped/ssp245/ang_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(ang_pred_high_ssp585_30_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_ssp585_30_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_ssp585_30_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ang_pred_high_ssp585_50_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_ssp585_50_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_ssp585_50_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(ang_pred_high_ssp585_70_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_mod_ssp585_70_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(ang_pred_low_ssp585_70_crop, "./sdm_output/ang/subs/cropped/ssp585/ang_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)



# Chloromeles - Maxent Model ----------------------------------------------
occ_chl_coords <- as.data.frame(geom(occThin_chl)[,3:4]) # extract longitude, lattitude from occurence points
bg_chl_coords <- as.data.frame(geom(chl_bg_vec)[,3:4]) # extract longitude, lattitude from background points

cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
set.seed(1337)

chl_maxent <- ENMevaluate(occ_chl_coords, # occurrence records
                          envs = wclim_chl_subs, # environment from background training area
                          n.bg = 20000, # 20000 bg points
                          tune.args =
                            list(rm = seq(0.5, 4, 0.5),
                                 fc = c("L", "LQ", "H",
                                        "LQH", "LQHP")),
                          partition.settings =
                            list(aggregation.factor = c(9, 9), gridSampleN = 20000), # 9,9 agg
                          partitions = 'checkerboard',
                          parallel = TRUE,
                          numCores = cn - 1, # leave one core available for other apps
                          algorithm = 'maxent.jar')

# Save the MaxEnt model so you do not have to waste time re-running the model
saveRDS(chl_maxent, file = './sdm_output/chl/subs/chl_maxent_subs.Rdata') # save


chl_maxent <- readRDS(file = './sdm_output/chl/subs/chl_maxent_subs.Rdata') # load subset model 
# Sect. Chloromeles Model Selecting -----------------------------------------

# Note that maxent results provide Continuous Boyce Index (cbi)
best_chl_maxent <- subset(chl_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_chl_maxent <- eval.models(chl_maxent)[[best_chl_maxent$tune.args]] # extracts the best model - returns MaxEnt object
# Best = rm.2_fc.LQHPT

eval.variable.importance(chl_maxent)[best_chl_maxent["tune.args"][1, 1]]



# Sect. Chloromeles historical prediction -----------------------------------

# Now use the <terra> package to plot the SDM predictchl.
# Wclim is the historical climatic conditchls (1970-2000)
cn <- detectCores(logical = F) # logical = F, is number of physical RAM cores in your computer
chl_pred_hist <- terra::predict(wclim_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)

plot(chl_pred_hist)
points(occThin_chl, cex = 0.05)

# Evaluate predictchls using Boyce Index
# the number of true presences should decline with suitability groups 100-91, 90-81, etc. 
# First extract suitability values for the background and presence points, make sure to omit NA values
chlPred_bg_val <- terra::extract(chl_pred_hist, bg_chl_coords)$lyr1 %>% 
  na.omit()

chlPred_val_na <- terra::extract(chl_pred_hist, occ_chl_coords)$lyr1 %>% 
  na.omit()

# Sect. Chloromeles Boyce Index ---------------------------------------------

# Evaluate predictchls using Boyce Indextea
png('C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/boyce_index/corboyce_plot_sect_chloromeles.png', width = 1600, height = 1200, res = 300)
chlBoyce <- ecospat.boyce(fit = chlPred_bg_val, # vector of predicted habitat suitability of bg points
              obs = chlPred_val_na, # vector of 
              nclass = 0, 
              PEplot = TRUE,
              method = 'spearman')

title(main = bquote("Sect." ~ italic("Chloromeles") ~ ", Boyce Cor." ==
                      .(chlBoyce$cor)))
dev.off()

# Gradients can be hard to understand at a glance, so lets create categorical bins of high suitability, moderate suitability, low suitability using thresholds
chlPred_val <- terra::extract(chl_pred_hist, occ_chl_coords)$lyr1
chlPred_threshold_1 <- quantile(chlPred_val, 0.01, na.rm = T) # Low suitability
chlPred_threshold_10 <- quantile(chlPred_val, 0.1, na.rm = T) # Moderate suitability
chlPred_threshold_50 <- quantile(chlPred_val, 0.5, na.rm = T) # High suitability

# Sect. Chloromeles - Climate Predictions -----------------------------------


# Future Climate predictchls
# SSP 245
chl_pred_ssp245_30 <- terra::predict(ssp245_2030_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp245_50 <- terra::predict(ssp245_2050_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp245_70 <- terra::predict(ssp245_2070_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)

# SSP 585
chl_pred_ssp585_30 <- terra::predict(ssp585_2030_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp585_50 <- terra::predict(ssp585_2050_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)
chl_pred_ssp585_70 <- terra::predict(ssp585_2070_subs, mod.best_chl_maxent, cores = cn - 1, na.rm = T)


# Save/Load Sect. Chloromeles SDM predictions ----------------------------------------

# Save
saveRDS(chl_pred_hist, file = './sdm_output/chl/subs/chl_pred_hist_subs.Rdata')

saveRDS(chl_pred_ssp245_30, file = './sdm_output/chl/subs/chl_pred_ssp245_30_subs.Rdata')
saveRDS(chl_pred_ssp245_50, file = './sdm_output/chl/subs/chl_pred_ssp245_50_subs.Rdata')
saveRDS(chl_pred_ssp245_70, file = './sdm_output/chl/subs/chl_pred_ssp245_70_subs.Rdata')

saveRDS(chl_pred_ssp585_30, file = './sdm_output/chl/subs/chl_pred_ssp585_30_subs.Rdata')
saveRDS(chl_pred_ssp585_50, file = './sdm_output/chl/subs/chl_pred_ssp585_50_subs.Rdata')
saveRDS(chl_pred_ssp585_70, file = './sdm_output/chl/subs/chl_pred_ssp585_70_subs.Rdata')

# Load
#Subsetted models
chl_pred_hist <- readRDS(file = './sdm_output/chl/subs/chl_pred_hist_subs.Rdata')

chl_pred_ssp245_30 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp245_30_subs.Rdata')
chl_pred_ssp245_50 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp245_50_subs.Rdata')
chl_pred_ssp245_70 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp245_70_subs.Rdata')

chl_pred_ssp585_30 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_30_subs.Rdata')
chl_pred_ssp585_50 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_50_subs.Rdata')
chl_pred_ssp585_70 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_70_subs.Rdata')



# Save Sect. Chloromeles thresholds -----------------------------------------


saveRDS(chlPred_threshold_1, file = './sdm_output/chl/subs/threshold/chlPred_threshold_1_subs.Rdata')
saveRDS(chlPred_threshold_10, file = './sdm_output/chl/subs/threshold/chlPred_threshold_10_subs.Rdata')
saveRDS(chlPred_threshold_50, file = './sdm_output/chl/subs/threshold/chlPred_threshold_50_subs.Rdata')


# Load subsetted
chlPred_threshold_1 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_1_subs.Rdata')
chlPred_threshold_10 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_10_subs.Rdata')
chlPred_threshold_50 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_50_subs.Rdata')


# Sect. Chloromeles Habitat predictions --------------------------------------
# Categorical habitat suitability
# Historical
chl_pred_high_hist <- chl_pred_hist > chlPred_threshold_50
chl_pred_mod_hist <- chl_pred_hist > chlPred_threshold_10
chl_pred_low_hist <- chl_pred_hist > chlPred_threshold_1

#SSP245 
chl_pred_high_ssp245_30 <- chl_pred_ssp245_30 > chlPred_threshold_50
chl_pred_mod_ssp245_30 <- chl_pred_ssp245_30 > chlPred_threshold_10
chl_pred_low_ssp245_30 <- chl_pred_ssp245_30 > chlPred_threshold_1

chl_pred_high_ssp245_50 <- chl_pred_ssp245_50 > chlPred_threshold_50
chl_pred_mod_ssp245_50 <- chl_pred_ssp245_50 > chlPred_threshold_10
chl_pred_low_ssp245_50 <- chl_pred_ssp245_50 > chlPred_threshold_1

chl_pred_high_ssp245_70 <- chl_pred_ssp245_70 > chlPred_threshold_50
chl_pred_mod_ssp245_70 <- chl_pred_ssp245_70 > chlPred_threshold_10
chl_pred_low_ssp245_70 <- chl_pred_ssp245_70 > chlPred_threshold_1

#SSP585
chl_pred_high_ssp585_30 <- chl_pred_ssp585_30 > chlPred_threshold_50
chl_pred_mod_ssp585_30 <- chl_pred_ssp585_30 > chlPred_threshold_10
chl_pred_low_ssp585_30 <- chl_pred_ssp585_30 > chlPred_threshold_1

chl_pred_high_ssp585_50 <- chl_pred_ssp585_50 > chlPred_threshold_50
chl_pred_mod_ssp585_50 <- chl_pred_ssp585_50 > chlPred_threshold_10
chl_pred_low_ssp585_50 <- chl_pred_ssp585_50 > chlPred_threshold_1

chl_pred_high_ssp585_70 <- chl_pred_ssp585_70 > chlPred_threshold_50
chl_pred_mod_ssp585_70 <- chl_pred_ssp585_70 > chlPred_threshold_10
chl_pred_low_ssp585_70 <- chl_pred_ssp585_70 > chlPred_threshold_1

# Save

# Historical
saveRDS(chl_pred_high_hist, file = './sdm_output/chl/subs/habitat_pred/hist/chl_pred_high_hist_subs.Rdata')
saveRDS(chl_pred_mod_hist, file = './sdm_output/chl/subs/habitat_pred/hist/chl_pred_mod_hist_subs.Rdata')
saveRDS(chl_pred_low_hist, file = './sdm_output/chl/subs/habitat_pred/hist/chl_pred_low_hist_subs.Rdata')

# SSP245
saveRDS(chl_pred_high_ssp245_30, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_high_ssp245_30_subs.Rdata')
saveRDS(chl_pred_mod_ssp245_30, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_mod_ssp245_30_subs.Rdata')
saveRDS(chl_pred_low_ssp245_30, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_low_ssp245_30_subs.Rdata')

saveRDS(chl_pred_high_ssp245_50, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_high_ssp245_50_subs.Rdata')
saveRDS(chl_pred_mod_ssp245_50, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_mod_ssp245_50_subs.Rdata')
saveRDS(chl_pred_low_ssp245_50, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_low_ssp245_50_subs.Rdata')

saveRDS(chl_pred_high_ssp245_70, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_high_ssp245_70_subs.Rdata')
saveRDS(chl_pred_mod_ssp245_70, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_mod_ssp245_70_subs.Rdata')
saveRDS(chl_pred_low_ssp245_70, file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_low_ssp245_70_subs.Rdata')

# SSP585
saveRDS(chl_pred_high_ssp585_30, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_high_ssp585_30_subs.Rdata')
saveRDS(chl_pred_mod_ssp585_30, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_mod_ssp585_30_subs.Rdata')
saveRDS(chl_pred_low_ssp585_30, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_low_ssp585_30_subs.Rdata')

saveRDS(chl_pred_high_ssp585_50, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_high_ssp585_50_subs.Rdata')
saveRDS(chl_pred_mod_ssp585_50, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_mod_ssp585_50_subs.Rdata')
saveRDS(chl_pred_low_ssp585_50, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_low_ssp585_50_subs.Rdata')

saveRDS(chl_pred_high_ssp585_70, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_high_ssp585_70_subs.Rdata')
saveRDS(chl_pred_mod_ssp585_70, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_mod_ssp585_70_subs.Rdata')
saveRDS(chl_pred_low_ssp585_70, file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_low_ssp585_70_subs.Rdata')

# Load

# Historical
chl_pred_high_hist <- readRDS(file = './sdm_output/chl/subs/habitat_pred/hist/chl_pred_high_hist_subs.Rdata')
chl_pred_mod_hist <- readRDS(file = './sdm_output/chl/subs/habitat_pred/hist/chl_pred_mod_hist_subs.Rdata')
chl_pred_low_hist <- readRDS(file = './sdm_output/chl/subs/habitat_pred/hist/chl_pred_low_hist_subs.Rdata')

# SSP245
chl_pred_high_ssp245_30 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_high_ssp245_30_subs.Rdata')
chl_pred_mod_ssp245_30 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_mod_ssp245_30_subs.Rdata')
chl_pred_low_ssp245_30 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_low_ssp245_30_subs.Rdata')

chl_pred_high_ssp245_50 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_high_ssp245_50_subs.Rdata')
chl_pred_mod_ssp245_50 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_mod_ssp245_50_subs.Rdata')
chl_pred_low_ssp245_50 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_low_ssp245_50_subs.Rdata')

chl_pred_high_ssp245_70 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_high_ssp245_70_subs.Rdata')
chl_pred_mod_ssp245_70 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_mod_ssp245_70_subs.Rdata')
chl_pred_low_ssp245_70 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp245/chl_pred_low_ssp245_70_subs.Rdata')

# SSP585
chl_pred_high_ssp585_30 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_high_ssp585_30_subs.Rdata')
chl_pred_mod_ssp585_30 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_mod_ssp585_30_subs.Rdata')
chl_pred_low_ssp585_30 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_low_ssp585_30_subs.Rdata')

chl_pred_high_ssp585_50 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_high_ssp585_50_subs.Rdata')
chl_pred_mod_ssp585_50 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_mod_ssp585_50_subs.Rdata')
chl_pred_low_ssp585_50 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_low_ssp585_50_subs.Rdata')

chl_pred_high_ssp585_70 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_high_ssp585_70_subs.Rdata')
chl_pred_mod_ssp585_70 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_mod_ssp585_70_subs.Rdata')
chl_pred_low_ssp585_70 <- readRDS(file = './sdm_output/chl/subs/habitat_pred/ssp585/chl_pred_low_ssp585_70_subs.Rdata')


# Crop categorical layers to restrict the extent of predicted suitability.
# In this case the model is making predictchls outside of what ecologicaly makes sense!
chl_ext <- ext(-170, -50, 5, 85)

# Historical
chl_pred_high_hist_crop <- crop(chl_pred_high_hist, chl_ext)
chl_pred_mod_hist_crop <- crop(chl_pred_mod_hist, chl_ext)
chl_pred_low_hist_crop <- crop(chl_pred_low_hist, chl_ext)

#SSP245 
chl_pred_high_ssp245_30_crop <- crop(chl_pred_high_ssp245_30, chl_ext)
chl_pred_mod_ssp245_30_crop <- crop(chl_pred_mod_ssp245_30, chl_ext)
chl_pred_low_ssp245_30_crop <- crop(chl_pred_low_ssp245_30, chl_ext)

chl_pred_high_ssp245_50_crop <- crop(chl_pred_high_ssp245_50, chl_ext)
chl_pred_mod_ssp245_50_crop <- crop(chl_pred_mod_ssp245_50, chl_ext)
chl_pred_low_ssp245_50_crop <- crop(chl_pred_low_ssp245_50, chl_ext)

chl_pred_high_ssp245_70_crop <- crop(chl_pred_high_ssp245_70, chl_ext)
chl_pred_mod_ssp245_70_crop <- crop(chl_pred_mod_ssp245_70, chl_ext)
chl_pred_low_ssp245_70_crop <- crop(chl_pred_low_ssp245_70, chl_ext)

#SSP585
chl_pred_high_ssp585_30_crop <- crop(chl_pred_high_ssp585_30, chl_ext)
chl_pred_mod_ssp585_30_crop <- crop(chl_pred_mod_ssp585_30, chl_ext)
chl_pred_low_ssp585_30_crop <- crop(chl_pred_low_ssp585_30, chl_ext)

chl_pred_high_ssp585_50_crop <- crop(chl_pred_high_ssp585_50, chl_ext)
chl_pred_mod_ssp585_50_crop <- crop(chl_pred_mod_ssp585_50, chl_ext)
chl_pred_low_ssp585_50_crop <- crop(chl_pred_low_ssp585_50, chl_ext)

chl_pred_high_ssp585_70_crop <- crop(chl_pred_high_ssp585_70, chl_ext)
chl_pred_mod_ssp585_70_crop <- crop(chl_pred_mod_ssp585_70, chl_ext)
chl_pred_low_ssp585_70_crop <- crop(chl_pred_low_ssp585_70, chl_ext)


# Sect. Chloromeles. chlustifolia - Crop rasters - TIFF  ---------------------------------

# Historical 
terra::writeRaster(chl_pred_high_hist_crop, "./sdm_output/chl/subs/cropped/hist/chl_pred_high_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_hist_crop, "./sdm_output/chl/subs/cropped/hist/chl_pred_mod_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_hist_crop, "./sdm_output/chl/subs/cropped/hist/chl_pred_low_hist_crop.tif", filetype = "GTiff", overwrite = TRUE)

# SSP245 
terra::writeRaster(chl_pred_high_ssp245_30_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_high_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_ssp245_30_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_mod_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_ssp245_30_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_low_ssp245_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(chl_pred_high_ssp245_50_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_high_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_ssp245_50_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_mod_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_ssp245_50_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_low_ssp245_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(chl_pred_high_ssp245_70_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_high_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_ssp245_70_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_mod_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_ssp245_70_crop, "./sdm_output/chl/subs/cropped/ssp245/chl_pred_low_ssp245_70_crop.tif", filetype = "GTiff", overwrite = TRUE)

#SSP585
terra::writeRaster(chl_pred_high_ssp585_30_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_high_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_ssp585_30_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_mod_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_ssp585_30_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_low_ssp585_30_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(chl_pred_high_ssp585_50_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_high_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_ssp585_50_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_mod_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_ssp585_50_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_low_ssp585_50_crop.tif", filetype = "GTiff", overwrite = TRUE)

terra::writeRaster(chl_pred_high_ssp585_70_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_high_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_mod_ssp585_70_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_mod_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
terra::writeRaster(chl_pred_low_ssp585_70_crop, "./sdm_output/chl/subs/cropped/ssp585/chl_pred_low_ssp585_70_crop.tif", filetype = "GTiff", overwrite = TRUE)
