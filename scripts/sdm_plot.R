# Top ---------------------------------------------------------------------
# SDM plotting and map making
# Terrell Roulston 
# Started May 7th, 2024

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package


# Load occurrences and  raster/vectors  -----------------------------------
# Occurrence Points in SpatVectors

occThin_cor <- readRDS(file = './occ_data/cor/occThin_cor.Rdata') # M. coronaria
occThin_fus <- readRDS(file = './occ_data/fus/occThin_fus.Rdata') # M. fusca
occThin_ion <- readRDS(file = './occ_data/ion/occThin_ion.Rdata') # M. ioensis
occThin_ang <- readRDS(file = './occ_data/ang/occThin_ang.Rdata') # M. angustifolia
occThin_chl <- readRDS(file = './occ_data/chl/occThin_chl.Rdata') # Chloromeles

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

can_us_mex_border <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/can_us_mex_border')

# Predicted habitat suitability rasters
# M. coronaria
cor_pred_hist <- readRDS(file = './sdm_output/cor/subs/cor_pred_hist_subs.Rdata')

cor_pred_ssp245_30 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp245_30_subs.Rdata')
cor_pred_ssp245_50 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp245_50_subs.Rdata')
cor_pred_ssp245_70 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp245_70_subs.Rdata')

cor_pred_ssp585_30 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_30_subs.Rdata')
cor_pred_ssp585_50 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_50_subs.Rdata')
cor_pred_ssp585_70 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_70_subs.Rdata')

# M. fusca
fus_pred_hist <- readRDS(file = './sdm_output/fus/subs/fus_pred_hist_subs.Rdata')

fus_pred_ssp245_30 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp245_30_subs.Rdata')
fus_pred_ssp245_50 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp245_50_subs.Rdata')
fus_pred_ssp245_70 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp245_70_subs.Rdata')

fus_pred_ssp585_30 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_30_subs.Rdata')
fus_pred_ssp585_50 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_50_subs.Rdata')
fus_pred_ssp585_70 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_70_subs.Rdata')

# M. ioensis
ion_pred_hist <- readRDS(file = './sdm_output/ion/subs/ion_pred_hist_subs.Rdata')

ion_pred_ssp245_30 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp245_30_subs.Rdata')
ion_pred_ssp245_50 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp245_50_subs.Rdata')
ion_pred_ssp245_70 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp245_70_subs.Rdata')

ion_pred_ssp585_30 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_30_subs.Rdata')
ion_pred_ssp585_50 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_50_subs.Rdata')
ion_pred_ssp585_70 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_70_subs.Rdata')

# M. angustifolia
ang_pred_hist <- readRDS(file = './sdm_output/ang/subs/ang_pred_hist_subs.Rdata')

ang_pred_ssp245_30 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp245_30_subs.Rdata')
ang_pred_ssp245_50 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp245_50_subs.Rdata')
ang_pred_ssp245_70 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp245_70_subs.Rdata')

ang_pred_ssp585_30 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_30_subs.Rdata')
ang_pred_ssp585_50 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_50_subs.Rdata')
ang_pred_ssp585_70 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_70_subs.Rdata')

# Chloromeles
chl_pred_hist <- readRDS(file = './sdm_output/chl/subs/chl_pred_hist_subs.Rdata')

chl_pred_ssp245_30 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp245_30_subs.Rdata')
chl_pred_ssp245_50 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp245_50_subs.Rdata')
chl_pred_ssp245_70 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp245_70_subs.Rdata')

chl_pred_ssp585_30 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_30_subs.Rdata')
chl_pred_ssp585_50 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_50_subs.Rdata')
chl_pred_ssp585_70 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_70_subs.Rdata')

# Thresholds
# M. coronaria
corPred_threshold_1 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_1_subs.Rdata')
corPred_threshold_10 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_10_subs.Rdata')
corPred_threshold_50 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_50_subs.Rdata')

#M. fusca
fusPred_threshold_1 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_1_subs.Rdata')
fusPred_threshold_10 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_10_subs.Rdata')
fusPred_threshold_50 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_50_subs.Rdata')

#M. ioensis
ionPred_threshold_1 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_1_subs.Rdata')
ionPred_threshold_10 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_10_subs.Rdata')
ionPred_threshold_50 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_50_subs.Rdata')

#M. angustifolia
angPred_threshold_1 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_1_subs.Rdata')
angPred_threshold_10 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_10_subs.Rdata')
angPred_threshold_50 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_50_subs.Rdata')

#Chloromeles
chlPred_threshold_1 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_1_subs.Rdata')
chlPred_threshold_10 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_10_subs.Rdata')
chlPred_threshold_50 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_50_subs.Rdata')


# Project to Lambert Conformal Conic --------------------------------------

projLam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Gratitcules
#g.cor <- terra::graticule(
#  lon = seq(-110, -55, by = 10),  # Longitude from 130°W to 50°W
#  lat = seq(30, 55, by = 10),     # Latitude from 25°N to 75°N
#  crs = projLam
#)
# Great lakes
great_lakes.lcc <- project(great_lakes, projLam)

# Can/US border
can_us_border.lcc <- project(can_us_border, projLam)
can_us_mex_border.lcc <- project(can_us_mex_border, projLam)

# Occurrences
occThin_cor.lcc <- project(occThin_cor, projLam)
occThin_fus.lcc <- project(occThin_fus, projLam)
occThin_ion.lcc <- project(occThin_ion, projLam)
occThin_ang.lcc <- project(occThin_ang, projLam)
occThin_chl.lcc <- project(occThin_chl, projLam)

# Habitat suitability
# M. coronaria
cor_pred_hist.lcc <- project(cor_pred_hist, projLam)

cor_pred_ssp245_30.lcc <- project(cor_pred_ssp245_30, projLam)
cor_pred_ssp245_50.lcc <- project(cor_pred_ssp245_50, projLam)
cor_pred_ssp245_70.lcc <- project(cor_pred_ssp245_70, projLam)

cor_pred_ssp585_30.lcc <- project(cor_pred_ssp585_30, projLam)
cor_pred_ssp585_50.lcc <- project(cor_pred_ssp585_50, projLam)
cor_pred_ssp585_70.lcc <- project(cor_pred_ssp585_70, projLam)

# M. fusca
fus_pred_hist.lcc <- project(fus_pred_hist, projLam)

fus_pred_ssp245_30.lcc <- project(fus_pred_ssp245_30, projLam)
fus_pred_ssp245_50.lcc <- project(fus_pred_ssp245_50, projLam)
fus_pred_ssp245_70.lcc <- project(fus_pred_ssp245_70, projLam)

fus_pred_ssp585_30.lcc <- project(fus_pred_ssp585_30, projLam)
fus_pred_ssp585_50.lcc <- project(fus_pred_ssp585_50, projLam)
fus_pred_ssp585_70.lcc <- project(fus_pred_ssp585_70, projLam)

# M. ionesis
ion_pred_hist.lcc <- project(ion_pred_hist, projLam)

ion_pred_ssp245_30.lcc <- project(ion_pred_ssp245_30, projLam)
ion_pred_ssp245_50.lcc <- project(ion_pred_ssp245_50, projLam)
ion_pred_ssp245_70.lcc <- project(ion_pred_ssp245_70, projLam)

ion_pred_ssp585_30.lcc <- project(ion_pred_ssp585_30, projLam)
ion_pred_ssp585_50.lcc <- project(ion_pred_ssp585_50, projLam)
ion_pred_ssp585_70.lcc <- project(ion_pred_ssp585_70, projLam)

# M. angustifolia
ang_pred_hist.lcc <- project(ang_pred_hist, projLam)

ang_pred_ssp245_30.lcc <- project(ang_pred_ssp245_30, projLam)
ang_pred_ssp245_50.lcc <- project(ang_pred_ssp245_50, projLam)
ang_pred_ssp245_70.lcc <- project(ang_pred_ssp245_70, projLam)

ang_pred_ssp585_30.lcc <- project(ang_pred_ssp585_30, projLam)
ang_pred_ssp585_50.lcc <- project(ang_pred_ssp585_50, projLam)
ang_pred_ssp585_70.lcc <- project(ang_pred_ssp585_70, projLam)

# Chloromeles
chl_pred_hist.lcc <- project(chl_pred_hist, projLam)

chl_pred_ssp245_30.lcc <- project(chl_pred_ssp245_30, projLam)
chl_pred_ssp245_50.lcc <- project(chl_pred_ssp245_50, projLam)
chl_pred_ssp245_70.lcc <- project(chl_pred_ssp245_70, projLam)

chl_pred_ssp585_30.lcc <- project(chl_pred_ssp585_30, projLam)
chl_pred_ssp585_50.lcc <- project(chl_pred_ssp585_50, projLam)
chl_pred_ssp585_70.lcc <- project(chl_pred_ssp585_70, projLam)


# M. coronaria future habitat plot ----------------------------------------
# Predicted historical distribtuion
legend_labs <- rev(c('Low', 'Moderate', 'High'))
fill_cols <- rev(c("#FFF7BC", "#FEC44F", "#D95F0E"))

jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/vertBIG_suitability_legend.jpeg", width = 9999, height = 6666, res = 300)
# Plot a legend that can be saved on its own
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#legend('center', xpd = NA, title = c(as.expression(bquote(bold('Habitat Suitability')))), legend = legend_labs, fill = fill_cols, cex = 3)
legend('center', xpd = NA, box.lwd = 2, legend = legend_labs, fill = fill_cols, cex = 9, horiz = F, bty = "o", title = 'Habitat Suitability')
dev.off()

cor.xlim <- c(-0.9*10^6, 3.1*10^6)  # Expand westward and eastward
cor.ylim <- c(-2.7*10^6, 2.5*10^6) # Expand southward and northward

# Graticule lines for the Lambert projection for cornaria
#g.cor <- terra::graticule(
#  lon = seq(-105, -45, by = 10),  # Longitude from 105°W to 45°W
#  lat = seq(30, 75, by = 10),     # Latitude from 30°N to 75°N
#  crs = projLam
#)

# Plot species occurrences
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/occ/cor_occ_map.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = "M. coronaria Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(can_us_mex_border.lcc, add = T)
terra::points(occThin_cor.lcc)

dev.off()

# Plot historical distribtion
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_historical.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_hist.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_hist.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_hist.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
#terra::plot(g.cor, add = T, retro = T, col = '#999999',lab.lat = c(2,3,4), lab.lon = c(2,3,4,5), lab.cex = 1.5, mar = c(5,5,5,5), lab.loc = 1, off.lat = -0.01)

# SSP245
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp245_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp245_30.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp245_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp245_50.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp245_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp245_70.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP585
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_30.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_50.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = cor.xlim, ylim = cor.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(cor_pred_ssp585_70.lcc > corPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# M. fusca future habitat plot --------------------------------------------

# Predicted historical distribtuion
fus.xlim <- c(-4*10^6, -1*10^6)
fus.ylim <- c(-2*10^6, 3.2*10^6)

# Occurrences
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/occ/fusca_occ_map.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = fus.xlim, ylim = fus.ylim, 
            main = "M. fusca Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(can_us_mex_border.lcc, add = T)
terra::points(occThin_fus.lcc)

dev.off()
# HISTORICAL 
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_historical.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_hist.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_hist.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_hist.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP245
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp245_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp245_30.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp245_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp245_50.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp245_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp245_70.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP585
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_30.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_50.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = fus.xlim, ylim = fus.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(fus_pred_ssp585_70.lcc > fusPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# M. ionesis future habitat plot ------------------------------------------

# Predicted historical distribtuion
ion.xlim <- c(-0.9*10^6, 3.1*10^6)  # Expand westward and eastward
ion.ylim <- c(-2.7*10^6, 2.5*10^6) # Expand southward and northward

# OCCURRENCE
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/occ/ionesis_occ_map.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = ion.xlim, ylim = ion.ylim, 
            main = "M. ioensis Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(can_us_mex_border.lcc, add = T)
terra::points(occThin_ion.lcc)

dev.off()
# HISTORICAL
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_historical.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_hist.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_hist.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_hist.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP245
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp245_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_ssp245_30.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_ssp245_30.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp245_30.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp245_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_ssp245_50.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_ssp245_50.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp245_50.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp245_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_ssp245_70.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_ssp245_70.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp245_70.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP585
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_ssp585_30.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_ssp585_30.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp585_30.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_ssp585_50.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_ssp585_50.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp585_50.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ion_pred_ssp585_70.lcc > ionPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ion.xlim, ylim = ion.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ion_pred_ssp585_70.lcc > ionPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ion_pred_ssp585_70.lcc > ionPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# M. angustifolia future habitat plot -------------------------------------

# Predicted historical distribtuion

ang.xlim <- c(-0.9*10^6, 3.1*10^6)  # Expand westward and eastward
ang.ylim <- c(-2.7*10^6, 2.5*10^6) # Expand southward and northward

# OCCURRENCE
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/occ/angustifolia_occ_map.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = ang.xlim, ylim = ang.ylim, 
            main = "M. angustifolia Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(can_us_mex_border.lcc, add = T)
terra::points(occThin_ang.lcc)

dev.off()
# HISTORICAL
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_historical.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_hist.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_hist.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_hist.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP245
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp245_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_ssp245_30.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_ssp245_30.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_ssp245_30.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp245_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_ssp245_50.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_ssp245_50.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_ssp245_50.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp245_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_ssp245_70.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_ssp245_70.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_ssp245_70.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP585
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_ssp585_30.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_ssp585_30.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_ssp585_30.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_ssp585_50.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_ssp585_50.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_ssp585_50.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(ang_pred_ssp585_70.lcc > angPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = ang.xlim, ylim = ang.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(ang_pred_ssp585_70.lcc > angPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(ang_pred_ssp585_70.lcc > angPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()

# Chloromeles future habitat plot -----------------------------------------

# Predicted historical distribtuion

chl.xlim <- c(-0.9*10^6, 3.1*10^6)  # Expand westward and eastward
chl.ylim <- c(-2.7*10^6, 2.5*10^6) # Expand southward and northward

# OCCURRENCE
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/occ/chloromeles_occ_map.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_hist.lcc , col = c('#E8E8E8', '#E8E8E8'),
            background = 'lightskyblue1',
            legend = F, 
            xlim = chl.xlim, ylim = chl.ylim, 
            main = "Sect. Chloromeles Occurrences",
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(can_us_mex_border.lcc, add = T)
terra::points(occThin_chl.lcc)

dev.off()
# HISTORICAL
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_historical.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_hist.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             
            xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_hist.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_hist.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP245
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp245_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_ssp245_30.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_ssp245_30.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_ssp245_30.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp245_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_ssp245_50.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_ssp245_50.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_ssp245_50.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp245_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_ssp245_70.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_ssp245_70.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_ssp245_70.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# SSP585
# 2030
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2030.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_ssp585_30.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_ssp585_30.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_ssp585_30.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2050
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2050.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_ssp585_50.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes = F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_ssp585_50.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_ssp585_50.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()
# 2070
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2070.jpeg", width = 3333, height = 6666, res = 300)

terra::plot(chl_pred_ssp585_70.lcc > chlPred_threshold_1, col = c('#E8E8E8', '#FFF7BC'), legend = F, 
            background = 'lightskyblue1',             xlim = chl.xlim, ylim = chl.ylim, 
            main = '',
            cex.main = 3,
            axes= F,
            box = T,
            mar = c(1, 1, 1, 1))
terra::plot(chl_pred_ssp585_70.lcc > chlPred_threshold_10, col = c("#FFFFFF00", '#FEC44F'), add = T, legend = F)
terra::plot(chl_pred_ssp585_70.lcc > chlPred_threshold_50, col = c("#FFFFFF00", '#D95F0E'), add = T, legend = F)
terra::plot(can_us_mex_border.lcc, add = T)

dev.off()


# SDM GRID PLOT -----------------------------------------------------------
# Lets load the jpeg images and arrange them using the <magick> package
# Create a 3x5 grid with each species, and recent (historical), 2030 and 2070 SSP585
library(magick)

# Load jpegs
cor_hist <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_historical.jpeg")
cor_2030 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2030.jpeg")
cor_2070 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2070.jpeg")

fus_hist <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_historical.jpeg")
fus_2030 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2030.jpeg")
fus_2070 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2070.jpeg")

ion_hist <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_historical.jpeg")
ion_2030 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2030.jpeg")
ion_2070 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2070.jpeg")

ang_hist <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_historical.jpeg")
ang_2030 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2030.jpeg")
ang_2070 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2070.jpeg")

chl_hist <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_historical.jpeg")
chl_2030 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2030.jpeg")
chl_2070 <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2070.jpeg")

# Stack each species horizontally (time series)
cor_row <- image_append(c(cor_hist, cor_2030, cor_2070))
fus_row <- image_append(c(fus_hist, fus_2030, fus_2070))
ion_row <- image_append(c(ion_hist, ion_2030, ion_2070))
ang_row <- image_append(c(ang_hist, ang_2030, ang_2070))
chl_row <- image_append(c(chl_hist, chl_2030, chl_2070))

# Load and resize the legend to match a row
legend_path <- "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/vertBIG_suitability_legend.jpeg"
legend_img <- image_read(legend_path)

# Create blank images matching other plot size (3333 x 6666)
#blank_img <- image_blank(width = 3333, height = 6666, color = "white")

# Make a row for the legend followed by two blanks
#leg_row <- image_append(c(legend_img, blank_img))

# Combine into final layout: 3 rows × 2 columns
row1 <- image_append(c(fus_row, legend_img))
row2 <- image_append(c(cor_row, ion_row))
row3 <- image_append(c(ang_row, chl_row))

mega_grid <- image_append(c(row1, row2, row3), stack = TRUE)

# Save the final mega figure
image_write(mega_grid, path = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/SDM_mega_figure_V2.jpeg", format = "jpeg")


#row1 <- image_append(c(fus_row, chl_row, legend_img))
#row2 <- image_append(c(cor_row, ion_row, ang_row))

#final_grid <- image_append(c(row1, row2), stack = TRUE)

# Save final figure
#image_write(final_grid, "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/SDM_custom_figure.jpeg")

library(magick)

# -- SETUP --
load_panel <- function(path) {
  if (!file.exists(path)) {
    stop("The file does not exist: ", path)
  }
  image_read(path)
}

img_width <- 3333
img_height <- 6666
font_family <- "Arial"

# -- LOAD PANELS IN YOUR DESIRED ORDER --
fus <- list(
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_historical.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2030.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/fus/fusca_ssp585_2070.jpeg")
)

cor <- list(
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_historical.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2030.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/cor/coronaria_ssp585_2070.jpeg")
)

ion <- list(
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_historical.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2030.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ion/ioensis_ssp585_2070.jpeg")
)

ang <- list(
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_historical.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2030.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/ang/angustifolia_ssp585_2070.jpeg")
)

chl <- list(
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_historical.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2030.jpeg"),
  load_panel("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/chl/chloromeles_ssp585_2070.jpeg")
)

# Legend and blank
legend_img <- image_read("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/vert_suitability_legend.jpeg")
blank_img <- image_blank(width = img_width, height = img_height, color = "white")

# -- LABEL PANELS A–R --
all_panels <- list(fus[[1]], fus[[2]], fus[[3]],
                   cor[[1]], cor[[2]], cor[[3]],
                   ion[[1]], ion[[2]], ion[[3]],
                   ang[[1]], ang[[2]], ang[[3]],
                   chl[[1]], chl[[2]], chl[[3]])

panel_labels <- LETTERS[1:length(all_panels)]

annotated_panels <- mapply(function(img, label) {
  image_annotate(img, label, gravity = "northwest", location = "+100+100",
                 size = 100, font = font_family, weight = 700, color = "black")
}, all_panels, panel_labels, SIMPLIFY = FALSE)

# -- LABEL TIME PERIODS: Recent, 2030, 2070 --
column_titles <- c("Recent (1970–2000)", "2030", "2070")
title_strip <- function(title) {
  blank_img <- image_blank(width = img_width, height = 600, color = "white")
  image_annotate(blank_img, title, size = 120, gravity = "center", 
                 font = font_family, weight = 700, color = "black")
}
header_row <- image_append(do.call(c, header_list))

# -- BUILD ROWS WITH SPECIES TITLES --
build_species_row <- function(images, label) {
  # Create the species title bar using image_blank and image_annotate
  title_bar <- image_blank(width = img_width * 3, height = 500, color = "white")
  title_bar <- image_annotate(title_bar, label, gravity = "center", 
                              size = 120, font = font_family, weight = 700, color = "black")
  
  # Convert the list of images into a proper magick image vector
  row_images_vector <- do.call(c, images)
  # Append the images horizontally to form the row of maps
  row <- image_append(row_images_vector)
  
  # Combine the title bar and the row vertically. 
  # First, combine the images into a magick vector, then append vertically (stack = TRUE)
  combined <- image_append(do.call(c, list(title_bar, row)), stack = TRUE)
  return(combined)
}
# Rows (species name + 3 maps)
row_fus <- build_species_row(annotated_panels[1:3], "Malus fusca")
row_cor <- build_species_row(annotated_panels[4:6], "Malus coronaria")
row_ion <- build_species_row(annotated_panels[7:9], "Malus ioensis")
row_ang <- build_species_row(annotated_panels[10:12], "Malus angustifolia")
row_chl <- build_species_row(annotated_panels[13:15], "Sect. Chloromeles")

# -- ADD BLANKS + LEGEND AS SIDE PANELS TO FORM 2-COLUMN GRID --
# Row 1: fusca + legend
row1 <- image_append(c(row_fus, image_append(c(legend_img, blank_img), stack = TRUE)))

# Row 2: coronaria + ioensis
row2 <- image_append(c(row_cor, row_ion))

# Row 3: angustifolia + chloromeles
row3 <- image_append(c(row_ang, row_chl))

# -- BUILD FINAL GRID --
final_grid <- image_append(c(header_row, row1, row2, row3), stack = TRUE)

# -- SAVE OUTPUT --
out_path <- "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/habitat/SDM_final_labeled_correct_order.jpeg"

image_write(final_grid, path = out_path, format = "jpeg")



# Supplemental SDM Figure.  -----------------------------------------------
# Function for automating the plotting 
plot_sdm_greyscale <- function(r, threshold1, threshold10, threshold50, xlim, ylim, out_path) {
  jpeg(filename = out_path, width = 1600, height = 1200, res = 300)
  terra::plot(r > threshold1, col = c("#F2F2F2", "#CCCCCC"), 
              legend = FALSE, background = "white",
              xlim = xlim, ylim = ylim, main = "", axes = FALSE,
              box = FALSE, mar = c(1, 1, 1, 1))
  terra::plot(r > threshold10, col = c("#FFFFFF00", "#999999"), add = TRUE, legend = FALSE)
  terra::plot(r > threshold50, col = c("#FFFFFF00", "#333333"), add = TRUE, legend = FALSE)
  terra::plot(can_us_mex_border.lcc, add = TRUE, col = "black", lwd = 0.6)
  dev.off()
}

# Plot limits
cor.xlim <- c(-0.9*10^6, 3.1*10^6)  
cor.ylim <- c(-2.7*10^6, 2.5*10^6) 

fus.xlim <- c(-4*10^6, -1*10^6)
fus.ylim <- c(-2*10^6, 3.2*10^6)

ion.xlim <- c(-0.9*10^6, 3.1*10^6)  
ion.ylim <- c(-2.7*10^6, 2.5*10^6) 

ang.xlim <- c(-0.9*10^6, 3.1*10^6) 
ang.ylim <- c(-2.7*10^6, 2.5*10^6) 

chl.xlim <- c(-0.9*10^6, 3.1*10^6)  
chl.ylim <- c(-2.7*10^6, 2.5*10^6) 




# Setup a meta list for all the SDM layers
species_info <- list(
  cor = list(name = "coronaria", rasters = list(
    hist = cor_pred_hist.lcc,
    ssp245_2030 = cor_pred_ssp245_30.lcc,
    ssp245_2050 = cor_pred_ssp245_50.lcc,
    ssp245_2070 = cor_pred_ssp245_70.lcc,
    ssp585_2030 = cor_pred_ssp585_30.lcc,
    ssp585_2050 = cor_pred_ssp585_50.lcc,
    ssp585_2070 = cor_pred_ssp585_70.lcc
  ), thresholds = list(
    t1 = corPred_threshold_1,
    t10 = corPred_threshold_10,
    t50 = corPred_threshold_50
  ), xlim = cor.xlim, ylim = cor.ylim),
  
  fus = list(name = "fusca", rasters = list(
    hist = fus_pred_hist.lcc,
    ssp245_2030 = fus_pred_ssp245_30.lcc,
    ssp245_2050 = fus_pred_ssp245_50.lcc,
    ssp245_2070 = fus_pred_ssp245_70.lcc,
    ssp585_2030 = fus_pred_ssp585_30.lcc,
    ssp585_2050 = fus_pred_ssp585_50.lcc,
    ssp585_2070 = fus_pred_ssp585_70.lcc
  ), thresholds = list(
    t1 = fusPred_threshold_1,
    t10 = fusPred_threshold_10,
    t50 = fusPred_threshold_50
  ), xlim = fus.xlim, ylim = fus.ylim),
  
  ion = list(name = "ioensis", rasters = list(
    hist = ion_pred_hist.lcc,
    ssp245_2030 = ion_pred_ssp245_30.lcc,
    ssp245_2050 = ion_pred_ssp245_50.lcc,
    ssp245_2070 = ion_pred_ssp245_70.lcc,
    ssp585_2030 = ion_pred_ssp585_30.lcc,
    ssp585_2050 = ion_pred_ssp585_50.lcc,
    ssp585_2070 = ion_pred_ssp585_70.lcc
  ), thresholds = list(
    t1 = ionPred_threshold_1,
    t10 = ionPred_threshold_10,
    t50 = ionPred_threshold_50
  ), xlim = ion.xlim, ylim = ion.ylim),
  
  ang = list(name = "angustifolia", rasters = list(
    hist = ang_pred_hist.lcc,
    ssp245_2030 = ang_pred_ssp245_30.lcc,
    ssp245_2050 = ang_pred_ssp245_50.lcc,
    ssp245_2070 = ang_pred_ssp245_70.lcc,
    ssp585_2030 = ang_pred_ssp585_30.lcc,
    ssp585_2050 = ang_pred_ssp585_50.lcc,
    ssp585_2070 = ang_pred_ssp585_70.lcc
  ), thresholds = list(
    t1 = angPred_threshold_1,
    t10 = angPred_threshold_10,
    t50 = angPred_threshold_50
  ), xlim = ang.xlim, ylim = ang.ylim),
  
  chl = list(name = "chloromeles", rasters = list(
    hist = chl_pred_hist.lcc,
    ssp245_2030 = chl_pred_ssp245_30.lcc,
    ssp245_2050 = chl_pred_ssp245_50.lcc,
    ssp245_2070 = chl_pred_ssp245_70.lcc,
    ssp585_2030 = chl_pred_ssp585_30.lcc,
    ssp585_2050 = chl_pred_ssp585_50.lcc,
    ssp585_2070 = chl_pred_ssp585_70.lcc
  ), thresholds = list(
    t1 = chlPred_threshold_1,
    t10 = chlPred_threshold_10,
    t50 = chlPred_threshold_50
  ), xlim = chl.xlim, ylim = chl.ylim)
)

# Loop through plots
scenarios <- c("hist", "ssp245_2030", "ssp245_2050", "ssp245_2070", 
               "ssp585_2030", "ssp585_2050", "ssp585_2070")

for (sp in names(species_info)) {
  s <- species_info[[sp]]
  for (scen in scenarios) {
    out_file <- paste0("C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/supplement/", sp, "_", scen, "_grey.jpeg")
    plot_sdm_greyscale(
      r = s$rasters[[scen]],
      threshold1 = s$thresholds$t1,
      threshold10 = s$thresholds$t10,
      threshold50 = s$thresholds$t50,
      xlim = s$xlim,
      ylim = s$ylim,
      out_path = out_file
    )
  }
}
