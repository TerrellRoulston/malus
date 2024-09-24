# Top ---------------------------------------------------------------------
# Data collection for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 16, 2024


library(tidyverse) # data management, grammar
library(rgbif) # access GBIF data

  getwd()
setwd('../malus/')

# Download occurrence data from GBIF

# GBIF user info
user='REDACTED'
pwd='REDACTED'
email='REDACTED'

# Taxon IDs
# M coronaria = 3001166
# M fusca = 3001080

# Coronaria download ------------------------------------------------------

taxonKey <- 3001166
basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE', 'MATERIAL_SAMPLE', 'LIVING_SPECIMEN') 
hasCoordinates <- TRUE # limit to records with coordinates
#years <- seq(1970, 2024, 1) # 1970 to modern day - this is baseline
# Tyler suggest not limiting records to 1970
# Main reason for Malus not occurring in places it did pre-1970 is due to habitat destruction
# Not due to unsuitable climatic niches

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred("hasCoordinate", hasCoordinates),
  #pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)

getwd() # check your working directory (wd)
setwd("./occ_data/") # set wd to a location where you want to save the csv file.
download_coronaria <- occ_download_get(down_code[1], overwrite = TRUE)
# extract csv from zipper folder and save as clearly named csv in excel or equivalent.

gbif_cor <- read.csv(file = "occ_coronaria.csv") # load named csv data



# fusca download ----------------------------------------------------------

taxonKey <- 3001080
basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE', 'MATERIAL_SAMPLE', 'LIVING_SPECIMEN')
hasCoordinates <- TRUE # limit to records with coordinates
# years <- seq(1970, 2024, 1) # 1970 to modern day - this is baseline
# see comments above about years of occurrence

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
down_code = occ_download(
  rgbif::pred("taxonKey", taxonKey),
  pred_in("basisOfRecord", basisOfRecord),
  pred("hasCoordinate", hasCoordinates),
  #pred_in("year", years),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email)


download_fusca <- occ_download_get(down_code[1], overwrite = TRUE)
# extract and save as csv

gbif_fusca <- read.csv(file = "occ_fusca.csv") # load csv data

