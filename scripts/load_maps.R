library(sf)
library(rnaturalearth)
library(terra)

message("**  Loading Maps")
# download/load maps
## with provinces and states:
us_map <- gadm(country = 'USA', level = 1, resolution = 2,
               path = "./maps/base_maps") #USA basemap w. States

ca_map <- gadm(country = 'CA', level = 1, resolution = 2,
               path = './maps/base_maps') #Canada basemap w. Provinces

mex_map <-gadm(country = 'MX', level = 1, resolution = 2,
               path = './maps/base_maps') # Mexico basemap w. States

canUSMex_map <- rbind(us_map, ca_map, mex_map) # Combine Mexico, US and Canada vector map

## Country boundaries:
us_map_0 <- gadm(country = 'USA', level = 0, resolution = 2, path = "./maps/base_maps") 
ca_map_0 <- gadm(country = 'CA', level = 0, resolution = 2, path = './maps/base_maps') 
mex_map_0 <-gadm(country = 'MX', level = 0, resolution = 2, path = './maps/base_maps') 
gl_map_0  <-gadm(country = 'GL', level = 0, resolution = 2, path = './maps/base_maps') 

can_us_mex_border <- rbind(us_map_0, ca_map_0, mex_map_0)

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

# Shape files downloaded from the USGS (https://www.sciencebase.gov/catalog/item/530f8a0ee4b0e7e46bd300dd)
erie <- vect("./maps/base_maps/hydro_p_LakeErie.shp")
ontario <- vect("./maps/base_maps/hydro_p_LakeOntario.shp")
michigan <- vect("./maps/base_maps/hydro_p_LakeMichigan.shp")
huron <- vect("./maps/base_maps/hydro_p_LakeHuron.shp")
superior <- vect("./maps/base_maps/hydro_p_LakeSuperior.shp")
stclair <- vect("./maps/base_maps/hydro_p_LakeStClair.shp")
great_lakes <- rbind(erie, ontario, michigan, huron, superior, stclair)

great_lakes <- project(great_lakes, "WGS84")

## Project to Lambert Conformal Conic for pretty maps:

projLam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

us_map_0.lcc <- project(us_map_0, projLam)
ca_map_0.lcc <- project(ca_map_0, projLam)
mex_map_0.lcc <- project(mex_map_0, projLam)
gl_map_0.lcc <- project(gl_map_0, projLam)
car_map_0.lcc <- project(car_map_0, projLam)
can_us_mex_border.lcc <- project(can_us_mex_border, projLam)
great_lakes.lcc <- project(great_lakes, projLam)

# Ecoregion prep ----------------------------------------------------------
# Download NA Ecoregion shapefile from: https://www.epa.gov/eco-research/ecoregions-north-america
# Load shapefile from local files in ./maps/eco_regions:
ecoNA <- vect(x = "maps/eco_regions/NA_CEC_Eco_Level2.shp")
ecoNA <- project(ecoNA, 'WGS84') # project ecoregion vector to same coords ref as basemap

## cc_sea reference map:
if(file.exists("./maps/seaRef/ne_110m_land.shp")){
  seaRef <- st_read("./maps/seaRef/ne_110m_land.shp")
} else {
  seaRef <- ne_download(scale = 110, type = 'land', category = 'physical',
                        returnclass = 'sf', destdir = "./maps/seaRef/")
}

# Note DO NOT PUSH wclim data**
NA_ext <- ext(-180, -30, 18, 85) # Set spatial extent of analyis to NA in Western Hemisphere

# Download/load WorldClim data under future climate scenarios -------------
# WARNING DO NOT PUSH WORLDCLIM DATA
# Historical climate 1970-2000
wclim <- geodata::worldclim_global(var = 'bio',
                                   res = 2.5, 
                                   version = '2.1', 
                                  path = "./wclim_data/") %>% terra::crop(NA_ext)
wclim <- mask(wclim, great_lakes, inverse = T,
                    wopt = list(names = names(wclim)))

# SSP (Shared social-economic pathway) 2.45 
# middle of the road projection, high climate adaptation, low climate mitigation
ssp245_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% crop(NA_ext)

ssp245_2030 <- mask(ssp245_2030, great_lakes, inverse = T,
                    wopt = list(names =
                                  gsub("bioc_.*(_.*)", "bio\\1",
                                       names(ssp245_2030)))) 

ssp245_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% crop(NA_ext)

ssp245_2050 <- mask(ssp245_2050, great_lakes, inverse = T,
                    wopt = list(names =
                                  gsub("bioc_.*(_.*)", "bio\\1",
                                       names(ssp245_2050)))) 

ssp245_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "245",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% crop(NA_ext)

ssp245_2070 <- mask(ssp245_2070, great_lakes, inverse = T,
                    wopt = list(names =
                                  gsub("bioc_.*(_.*)", "bio\\1",
                                       names(ssp245_2070)))) 

# SPP 5.85 
# low regard for enviromental sustainability, increased fossil fuel reliance, this is the current tracking projection
ssp585_2030 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2021-2040",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% crop(NA_ext) 

ssp585_2030 <- mask(ssp585_2030, great_lakes, inverse = T,
                    wopt = list(names =
                                  gsub("bioc_.*(_.*)", "bio\\1",
                                       names(ssp585_2030)))) 

ssp585_2050 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2041-2060",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% crop(NA_ext)

ssp585_2050 <- mask(ssp585_2050, great_lakes, inverse = T,
                    wopt = list(names =
                                  gsub("bioc_.*(_.*)", "bio\\1",
                                       names(ssp585_2050)))) 


ssp585_2070 <- cmip6_world(model = "CanESM5",
                           ssp = "585",
                           time = "2061-2080",
                           var = "bioc",
                           res = 2.5,
                           path = "./wclim_data/") %>% crop(NA_ext)

ssp585_2070 <- mask(ssp585_2070, great_lakes, inverse = T,
                    wopt = list(names =
                                  gsub("bioc_.*(_.*)", "bio\\1",
                                       names(ssp585_2070)))) 

message("**  Maps Loaded")
