# Top ---------------------------------------------------------------------
# Plotting thinned Occurence Data
# Figure 1 in SDM Paper
# Terrell Roulston
# Started May 8 2025

library(tidyverse) # grammar and data management 
library(terra) # working with spatial data
library(geodata) # basemaps and climate data

occThin_cor <- readRDS(file = "./occ_data/cor/occThin_cor.Rdata") # GBIF + Husband
occThin_fus <- readRDS(file = './occ_data/fus/occThin_fus.Rdata') # GBIF + Armstrong + Wickham + Obr. + Fit
occThin_ion <- readRDS(file = './occ_data/ion/occThin_ion.Rdata') # GBIF
occThin_ang <- readRDS(file = './occ_data/ang/occThin_ang.Rdata') # GBIF
occThin_chl <- readRDS(file = './occ_data/chl/occThin_chl.Rdata') 

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
can_us_mex_border <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/can_us_mex_border')

# Shape files downloaded from the USGS (https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
great_lakes <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/great lakes/combined great lakes/')

# seperate GBIF from other data sources
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
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/occ_plot/malus_occurrence_data_sources_NEW.jpeg", width = 10000, height = 6666, res = 300)

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
terra::plot(can_us_mex_border.lcc, box = F,  add = T, col = 'grey')
terra::points(occThin_cor_gbif.lcc, pch = 16, col = alpha("#882255", 1), cex = 1.3) # COR GBIF
terra::points(occThin_fus_gbif.lcc, pch = 16, col = alpha("#228B22", 1), cex = 1.3) # FUS GBIF
terra::points(occThin_cor_hus.lcc, pch = 16, col = alpha("magenta", 1), cex = 1.3) # COR HUSBAND
terra::points(occThin_fus_arm.lcc, pch = 16, col = alpha("#333f07", 1), cex = 1.3) # FUS ARMSTRONG
terra::points(occThin_fus_wick.lcc, pch = 16, col = alpha("#333f07", 1), cex = 1.3) # FUS WICKHAM
terra::points(occThin_fus_obr_fit.lcc, pch = 16, col = alpha("#333f07", 1), cex = 1.3) # FUS OBRITS ET AL. AND FITZPATRICK ET AL.
terra::points(occThin_ion.lcc, pch = 16, col = alpha("#E88E00", 1), cex = 1.3) # ION GBIF
terra::points(occThin_ang.lcc, pch = 16, col = alpha("#007CBE", 1), cex = 1.3) # ANG GBIF
terra::add_box(col = 'grey')
legend( # legend for data
  x = 1e6,
  y = 2.5e6,
  title = c(expression(underline('Thinned Occurrence Data'))),
  legend = c(expression(italic("Malus fusca")*"—GBIF"),
             expression(italic("Malus fusca")* "—Suppl. Data"),
             expression(italic("Malus coronaria")*"—GBIF"), 
             expression(italic("Malus coronaria")* "—Suppl. Data"),
             expression(italic("Malus ioensis")*"—GBIF"),
             expression(italic("Malus angustifolium")*"—GBIF")),
  fill = c("#228B22", "#333f07", "#882255", "magenta", "#E88E00", "#007CBE"),
  col = "black",
  box.col = "black",  # No border around legend
  bg = "white",
  text.col = 'black',
  cex = 3.25,
  xjust = 0,              
  yjust = 1,
  title.adj = 0.25
)
terra::text(533792.2, 1206373.8, labels = "Hudson\nBay", cex = 2.5, col = "steelblue")  
terra::text(2525623.3, -617100.0, labels = "Atlantic\nOcean", cex = 2.75, col = "steelblue")  
terra::text(-2838692.4, 684470.5, labels = "Pacific\nOcean", cex = 2.75, col = "steelblue")  
terra::text(-585874.4, 1236224, labels = 'Canada', cex = 3.5, col = "black")
terra::text(-585874.4, -1032531, labels = 'U.S.A.', cex = 3.5, col = "black")

dev.off()

"#AA4499"

"#CC79A7"

"#117733"

"#882255"
# Plot ONLY sect. Chloromeles

#Start plotting
jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/occ_plot/chloromeles_occurrence_data.jpeg", width = 10000, height = 6666, res = 300)

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
terra::plot(can_us_mex_border.lcc, box = F,  add = T, col = 'grey')
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