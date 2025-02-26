# Top ---------------------------------------------------------------------
# SDM plotting and map making
# Terrell Roulston 
# Started May 7th, 2024

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package


# Load occurrences and  raster/vectors  -----------------------------------
# Occurrence Points in SpatVectors

occThin_cor <- readRDS(file = './occ_data/occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = './occ_data/occThin_fus.Rdata') # M. fusca

# Great Lakes shapefiles for making pretty maps
# Shape files downloaded from the USGS (https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')

# Canada/US Border for showing the international line. Gadm admin boundaries trace the entire country, vs this is just the border
# Much easier to see the SDM results along coastlines where tracing obscures the data
# Downloaded from  https://koordinates.com/layer/111012-canada-and-us-border/
can_us_border <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/can_us border')

# Two line segments are in the water and are not needed in this case, lets remove them to make the maps look prettier
segments_to_remove <- c("Gulf of Maine", "Straits of Georgia and Juan de Fuca")
can_us_border <- can_us_border[!can_us_border$SectionEng %in% segments_to_remove, ]

# Predicted habitat suitability rasters
# M. coronaria
getwd()
setwd('../sdm_output/')

cor_pred_hist <- readRDS(file = 'cor_pred_hist.Rdata')

cor_pred_ssp245_30 <- readRDS(file = 'cor_pred_ssp245_30.Rdata')
cor_pred_ssp245_50 <- readRDS(file = 'cor_pred_ssp245_50.Rdata')
cor_pred_ssp245_70 <- readRDS(file = 'cor_pred_ssp245_70.Rdata')

cor_pred_ssp585_30 <- readRDS(file = 'cor_pred_ssp585_30.Rdata')
cor_pred_ssp585_50 <- readRDS(file = 'cor_pred_ssp585_50.Rdata')
cor_pred_ssp585_70 <- readRDS(file = 'cor_pred_ssp585_70.Rdata')

#M. fusca
fus_pred_hist <- readRDS(file = 'fus_pred_hist.Rdata')

fus_pred_ssp245_30 <- readRDS(file = 'fus_pred_ssp245_30.Rdata')
fus_pred_ssp245_50 <- readRDS(file = 'fus_pred_ssp245_50.Rdata')
fus_pred_ssp245_70 <- readRDS(file = 'fus_pred_ssp245_70.Rdata')

fus_pred_ssp585_30 <- readRDS(file = 'fus_pred_ssp585_30.Rdata')
fus_pred_ssp585_50 <- readRDS(file = 'fus_pred_ssp585_50.Rdata')
fus_pred_ssp585_70 <- readRDS(file = 'fus_pred_ssp585_70.Rdata')

# Thresholds
# M. coronaria
setwd('../sdm_output/thresholds')
corPred_threshold_1 <- readRDS(file = 'corPred_threshold_1.Rdata')
corPred_threshold_10 <- readRDS(file = 'corPred_threshold_10.Rdata')
corPred_threshold_50 <- readRDS(file = 'corPred_threshold_50.Rdata')

#M. fusca
fusPred_threshold_1 <- readRDS(file = 'fusPred_threshold_1.Rdata')
fusPred_threshold_10 <- readRDS(file = 'fusPred_threshold_10.Rdata')
fusPred_threshold_50 <- readRDS(file = 'fusPred_threshold_50.Rdata')



# Project to Lambert Conformal Conic --------------------------------------

projLam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Gratitcules
g.cor <- terra::graticule(
  lon = seq(-110, -55, by = 10),  # Longitude from 130°W to 50°W
  lat = seq(30, 55, by = 10),     # Latitude from 25°N to 75°N
  crs = projLam
)
# Great lakes
great_lakes.lcc <- project(great_lakes, projLam)

# Can/US border
can_us_border.lcc <- project(can_us_border, projLam)

# Occurrences
occThin_cor.lcc <- project(occThin_cor, projLam)
occThin_fus.lcc <- project(occThin_fus, projLam)

# Habitat suitability
cor_pred_hist.lcc <- project(cor_pred_hist, projLam)

cor_pred_ssp245_30.lcc <- project(cor_pred_ssp245_30, projLam)
cor_pred_ssp245_50.lcc <- project(cor_pred_ssp245_50, projLam)
cor_pred_ssp245_70.lcc <- project(cor_pred_ssp245_70, projLam)

cor_pred_ssp585_30.lcc <- project(cor_pred_ssp585_30, projLam)
cor_pred_ssp585_50.lcc <- project(cor_pred_ssp585_50, projLam)
cor_pred_ssp585_70.lcc <- project(cor_pred_ssp585_70, projLam)

#M. fusca
fus_pred_hist.lcc <- project(fus_pred_hist, projLam)

fus_pred_ssp245_30.lcc <- project(fus_pred_ssp245_30, projLam)
fus_pred_ssp245_50.lcc <- project(fus_pred_ssp245_50, projLam)
fus_pred_ssp245_70.lcc <- project(fus_pred_ssp245_70, projLam)

fus_pred_ssp585_30.lcc <- project(fus_pred_ssp585_30, projLam)
fus_pred_ssp585_50.lcc <- project(fus_pred_ssp585_50, projLam)
fus_pred_ssp585_70.lcc <- project(fus_pred_ssp585_70, projLam)


# M. coronaria future habitat plot ----------------------------------------
# Predicted historical distribtuion
legend_labs <- c('Low Suitability', 'Moderate Suitability', 'High Suitability')
fill_cols <- c("#FFF7BC", "#FEC44F", "#D95F0E")

cor.xlim <- c(-5*10^5, 3.1*10^6)
cor.ylim <- c(-2*10^6, 2*10^6)

# Plot a legend that can be saved on its own
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', xpd = NA, title = c(as.expression(bquote(bold('Habitat Suitability')))), legend = legend_labs, fill = fill_cols, cex = 3)

# Graticule lines for the Lambert projection for cornaria
g.cor <- terra::graticule(
  lon = seq(-105, -45, by = 10),  # Longitude from 105°W to 45°W
  lat = seq(30, 75, by = 10),     # Latitude from 30°N to 75°N
  crs = projLam
)

# Plot species occurrences

terra::plot(cor_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = "M. coronaria Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(can_us_border.lcc, add = T)
terra::points(occThin_cor.lcc)



# Plot historical distribtion 
terra::plot(cor_pred_hist.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = " Historical Suitability (1970-2000)",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_hist.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_hist.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)
#terra::plot(g.cor, add = T, retro = T, col = '#999999',lab.lat = c(2,3,4), lab.lon = c(2,3,4,5), lab.cex = 1.5, mar = c(5,5,5,5), lab.loc = 1, off.lat = -0.01)

# SSP245
# 2030

terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2030',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2050

terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2050',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2070

terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2070',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# SSP585
par(mfrow = c(1, 3))
# 2030

terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2030',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2050

terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2050',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2070

terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '2070',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)


# M. fusca future habitat plot --------------------------------------------

# Predicted historical distribtuion

fus.xlim <- c(-4*10^6, -1*10^6)
fus.ylim <- c(-2*10^6, 3.2*10^6)


terra::plot(fus_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = "M. fusca Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5, 5, 5, 5))
terra::plot(can_us_border.lcc, add = T)
terra::points(occThin_fus.lcc)


terra::plot(fus_pred_hist.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = fus.xlim, ylim = fus.ylim, 
            main = NULL,
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_hist.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_hist.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# SSP245
# 2030

terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '2030',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2050

terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '2050',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2070

terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '2070',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# SSP585
# 2030

terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '2030',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2050

terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '2050',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)

# 2070

terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '2070',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(5,5,5,5))
terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_10, col = c(NA, '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_50, col = c(NA, '#D95F0E'), add = T, legend = F)
terra::plot(can_us_border.lcc, add = T)


