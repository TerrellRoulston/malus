# Top ---------------------------------------------------------------------
# thinning occurrence data of M. coronaria and M. fusca
# Terrell Roulston
# Started Feb 20, 2024

library(tidyverse)
library(dismo)
library(terra) 
library(raster)
library(sp)


# Load cleaned occurrence data ---------------------------------------------
occ_cor <- readRDS(file = "occ_cor.Rdata")
occ_fus <- readRDS(file = 'occ_fus.Rdata')


