# Script for getting protected area data from the WDPA
# This is more straight forward then getting the data from the Canadian and US databases
# And allows for easier filtering of Marine PAs


library(terra)
library(wdpar)

# Fetch WDPA Data
# Canada
#wdpa_fetch("Canada", download_dir = "./gap_analysis/wdpar_us_can/", wait = TRUE, page_wait = 10, datatype = 'shp')
# US
#wdpa_fetch("United States of America", download_dir = "./gap_analysis/wdpar_us_can/", wait = TRUE, page_wait = 10, datatype = 'shp')

# Load and combine Canadian WDPA shapefiles
ca1 <- vect("./gap_analysis/wdpar_us_can/WDPA_May2025_CAN-shapefile/WDPA_WDOECM_May2025_Public_CAN_shp_0/WDPA_WDOECM_May2025_Public_CAN_shp-polygons.shp")
ca2 <- vect("./gap_analysis/wdpar_us_can/WDPA_May2025_CAN-shapefile/WDPA_WDOECM_May2025_Public_CAN_shp_1/WDPA_WDOECM_May2025_Public_CAN_shp-polygons.shp")
ca3 <- vect("./gap_analysis/wdpar_us_can/WDPA_May2025_CAN-shapefile/WDPA_WDOECM_May2025_Public_CAN_shp_2/WDPA_WDOECM_May2025_Public_CAN_shp-polygons.shp")
ca_all <- rbind(ca1, ca2, ca3)
ca_all_filt <- ca_all[ca_all$MARINE %in% c(0, 2), ] # filter only terrestrial parks (and marine + terrestrial parks)

# Load and filter US shapefile
us1 <- vect("./gap_analysis/wdpar_us_can/WDPA_May2025_USA-shapefile/WDPA_WDOECM_May2025_Public_USA_shp_0/WDPA_WDOECM_May2025_Public_USA_shp-polygons.shp")
us2 <- vect("./gap_analysis/wdpar_us_can/WDPA_May2025_USA-shapefile/WDPA_WDOECM_May2025_Public_USA_shp_1/WDPA_WDOECM_May2025_Public_USA_shp-polygons.shp")
us3 <- vect("./gap_analysis/wdpar_us_can/WDPA_May2025_USA-shapefile/WDPA_WDOECM_May2025_Public_USA_shp_2/WDPA_WDOECM_May2025_Public_USA_shp-polygons.shp")
us_all <- rbind(us1, us2, us3)
us_all_filt <- us_all[us_all$MARINE %in% c(0, 2), ]

# Combine US + Canada and reproject
pa_vect <- rbind(ca_all_filt, us_all_filt)
pa_vect <- project(pa_vect, "EPSG:4326")

# Define 2.5 arc-min raster template
template <- rast(
  xmin = -180, xmax = -50,
  ymin = 25, ymax = 85,
  resolution = c(2.5 / 60, 2.5 / 60),
  crs = "EPSG:4326"
)

# Rasterize vector into template
pa_raster <- terra::rasterize(pa_vect, template, field = 1, background = NA)

# Save output
writeRaster(pa_raster, "./gap_analysis/pa_raster_us_can.tif", overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=DEFLATE"))
