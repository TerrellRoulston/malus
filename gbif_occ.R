# Top ---------------------------------------------------------------------
# Data collection for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 16, 2024

library(tidyverse) #data management, grammar
library(rgbif) #access GBIF data

# Download occurrence data from GBIF

# GBIF use info
user='terrell_roulston'
pwd='MalusApples#123'
email='terrellroulston@gmail.com'

# Taxon IDs
# M coronaria = 3001166
# M fusca = 3001080

# Coronaria download ------------------------------------------------------

taxonKey <- 3001166
basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE') # excluded living specimens and material samples (germplasm)
hasCoordinates <- TRUE # limit to records with coordinates
years <- seq(1970, 2024, 1) # 1970 to modern day - this is baseline

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred("hasCoordinate", hasCoordinates),
  pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)

setwd("./") 
download_coronaria <- occ_download_get(down_code[1], overwrite = TRUE)
#extract and save as csv

gbif_cor <- read.csv(file = "occ_coronaria.csv") # load csv data



# Fusca download ----------------------------------------------------------

taxonKey <- 3001080
basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE') # excluded living specimens and material samples (germplasm)
hasCoordinates <- TRUE # limit to records with coordinates
years <- seq(1970, 2024, 1) # 1970 to modern day - this is baseline

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred("hasCoordinate", hasCoordinates),
  pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)


download_fusca <- occ_download_get(down_code[1], overwrite = TRUE)
#extract and save as csv

gbif_fusca <- read.csv(file = "occ_fusca.csv") # load csv data

