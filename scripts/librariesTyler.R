library(MASS) ## loading first so the geospatial packages will take
## precedence, hopefully?
library(sf)
library(rnaturalearth)
library(tidyverse) #grammar, data management
library(CoordinateCleaner) #helpful functions to clean data
library(terra) #working with vector/raster data
library(geodata) #download basemaps
library(scales) #alpha adjust colours
library(oce)
library(predicts)
library(ENMTools)
library(plotly) # 3D surface Kernel bivariate plots
library(ENMeval) # Another modeling package, useful for data partitioning (Checkerboarding)
library(ecospat) # Useful spatial ecology tools
library(parallel) # speed up computation by running in parallel
library(doParallel) # added functionality to parallel
library(ade4)
library(grid)



