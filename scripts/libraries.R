# List of libraries required for Main Analysis 
# Source this script in the analysis script to load all libs

### Data Management and Grammer
library(tidyverse) # Grammer and fundamental packages
library(readr) # Loading text files

### GIS Data
library(terra) # GIS package, geospatial data managment, spatial analysis, plotting
library(geodata) # Basemaps, and WorldClim Data (Historical and climate change projections) (NOTE: loading WorldClim data from memory via project)
library(wdpar) # World Database for Protected Areas, PA data for gap analysis
library(CoordinateCleaner) # Automated cleaning of occurrence data


### SDM Dependacies
library(rJava) # Dependacy for Maxent.jar via Java
library(predicts) # Base SDM package, Runs via ENMeval
library(ENMeval) # Load for additioanl functions for tuning SDMs
library(parallel) # Speed up SDM computation by running in parallel
library(doParallel) # Added functionality to parallel
library(dismo) # For plotting response curves (accepts Maxent objects) -- another SDM package


### Ecological Niche Modeling
library(ade4) # Multivariate PCA niche analysis
# Note, one function from the Ecospat package has a flaw.
# Tyler Smith has updated this function to work correctly, install this fix via his github
library(devtools) # Github remote access
devtools::install_github("plantarum/ecospat", ref = "boyce", subdir = "ecospat")
library(ecospat) # Ecolgical niche analysis
