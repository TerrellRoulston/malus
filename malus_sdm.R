# Top ---------------------------------------------------------------------
# MaxEnt Species Distribution Modeling (SDM) for Malus CWR (M. cornistica, M. fusca)
# Terrell Roulston
# Started Feb 29, 2024

library(tidyverse)
library(terra) 
library(predicts)
library(geodata)



# Load data for SDM -------------------------------------------------------

setwd("../occ_data/")

# Predictors values 
cor_sdmData <- readRDS(file = 'cor_sdmData.Rdata') # M. coronaria
cor_pb <- cor_sdmData[1] # save object with 0/1 background/presence vector

fus_sdmData <- readRDS(file = 'fus_sdmData.Rdata') # M. fusca
fus_pb <- fus_sdmData[1] # save object with 0/1 background/presence vector



# Select predictor variables for SDM --------------------------------------
# See dendograms (in malus_bg) for reference
# Select variables that are not colinear together, or one var from each colinear group
# See https://www.worldclim.org/data/bioclim.html for definition of variables from codes

# M. coronaria
names(cor_sdmData) # peak at predictor variables for assistance
cor_sdmData.x <- cor_sdmData %>% 
  dplyr::select('wc2.1_2.5m_bio_8', 'wc2.1_2.5m_bio_2', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
#BIO8 = Mean Temperature of Wettest Quarter; BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO12 =  Annual Mean Precipitation; BIO1 = Annual Mean Temperature

# M. fusca
names(fus_sdmData) # peak at predictor variables for assistance
fus_sdmData.x <- fus_sdmData %>% 
  dplyr::select('wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_8', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
#BIO15 = Precipitation Seasonality (Coefficient of Variation); BIO8 = Mean Temperature of Wettest Quarter
#BIO12 =  Annual Mean Precipitation; BIO1 = Annual Mean Temperature



predicts::MaxEnt()