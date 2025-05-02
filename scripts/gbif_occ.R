# Top ---------------------------------------------------------------------
# Data collection for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 16, 2024

## Reload GBIF data (after it was downloaded using the code below)

gbif_cor <- read.csv(file = "./occ_data/cor/occ_coronaria.csv") # load named csv data
gbif_fusca <- read.csv(file = "./occ_data/fus/occ_fusca.csv") # load csv data
gbif_ioensis <- read.csv(file = "./occ_data/ion/occ_ioensis.csv") # load csv data
gbif_angustifolia <- read.csv(file = "./occ_data/ang/occ_angustifolia.csv") # load csv data

## The following is needed only the first time the data is downloaded from
## GBIF. If you get the data from this repository, you don't need to run
## this at all.

## library(tidyverse) # data management, grammar
## library(rgbif) # access GBIF data

# GBIF user info
## user='REDACTED'
## pwd='REDACTED'
## email='REDACTED'

# Taxon IDs
# M coronaria = 3001166
# M fusca = 3001080
# M. ioensis = 3001596
# M. angustifolia = 3001548

# Coronaria download ------------------------------------------------------

## taxonKey <- 3001166
## basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE', 'MATERIAL_SAMPLE', 'LIVING_SPECIMEN') 
## hasCoordinates <- TRUE # limit to records with coordinates

#years <- seq(1970, 2024, 1) # 1970 to modern day - this is baseline
# Tyler suggest not limiting records to 1970
# Main reason for Malus not occurring in places it did pre-1970 is due to habitat destruction
# Not due to unsuitable climatic niches

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
## down_code = occ_download(
##   pred("taxonKey", taxonKey),
##   pred_in("basisOfRecord", basisOfRecord),
##   pred("hasCoordinate", hasCoordinates),
##   #pred_in("year", years),
##   format = "SIMPLE_CSV",
##   user=user, pwd=pwd, email=email)


## download_coronaria <- occ_download_get(down_code[1], overwrite = TRUE, path = "./occ_data/")
# extract csv from zipper folder and save as clearly named csv in excel or equivalent.

## gbif_cor <- read.csv(file = "./occ_data/occ_coronaria.csv") # load named csv data



# fusca download ----------------------------------------------------------

## taxonKey <- 3001080
## basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE', 'MATERIAL_SAMPLE', 'LIVING_SPECIMEN')
## hasCoordinates <- TRUE # limit to records with coordinates

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple

## down_code = occ_download(
##   rgbif::pred("taxonKey", taxonKey),
##   pred_in("basisOfRecord", basisOfRecord),
##   pred("hasCoordinate", hasCoordinates),
##   #pred_in("year", years),
##   format = "SIMPLE_CSV",
##   user=user, pwd=pwd, email=email)


## download_fusca <- occ_download_get(down_code[1], overwrite = TRUE, path = "./occ_data/")
# extract and save as csv

## gbif_fusca <- read.csv(file = "./occ_data/occ_fusca.csv") # load csv data


# Ioensis download --------------------------------------------------------

## taxonKey <- 3001596
## basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE', 'MATERIAL_SAMPLE', 'LIVING_SPECIMEN')
## hasCoordinates <- TRUE # limit to records with coordinates

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple

## down_code = occ_download(
##   rgbif::pred("taxonKey", taxonKey),
##   pred_in("basisOfRecord", basisOfRecord),
##   pred("hasCoordinate", hasCoordinates),
##   #pred_in("year", years),
##   format = "SIMPLE_CSV",
##   user=user, pwd=pwd, email=email)


## download_ioensis <- occ_download_get(down_code[1], overwrite = TRUE, path = "./occ_data/")
# extract and save as csv

## gbif_ioensis <- read.csv(file = "./occ_data/occ_ioensis.csv") # load csv data


# Angustifolia download ---------------------------------------------------

## taxonKey <- 3001548
## basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION', 'OCCURRENCE', 'MATERIAL_SAMPLE', 'LIVING_SPECIMEN')
## hasCoordinates <- TRUE # limit to records with coordinates

# Download data
# Use 'pred()' if there is a single argument, or 'pred_in()' if there are multiple
## down_code = occ_download(
##   rgbif::pred("taxonKey", taxonKey),
##   pred_in("basisOfRecord", basisOfRecord),
##   pred("hasCoordinate", hasCoordinates),
##   #pred_in("year", years),
##   format = "SIMPLE_CSV",
##   user=user, pwd=pwd, email=email)


## download_angustifolia <- occ_download_get(down_code[1], overwrite = TRUE, path = "./occ_data/")
# extract and save as csv

## gbif_angustifolia <- read.csv(file = "./occ_data/occ_angustifolia.csv") # load csv data
