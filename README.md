# Native Canadian *Malus* CWR Species Distribution Modeling for *M. coronaria* and *M. fusca*
Terrell Roulston
Started: Feb 2024
Last Updated: Nov 2024

This repo contains code and data to produce species distribution models (SDM) for *M. coronaria* and *M. fusca*; two native *Malus* crop wild relatives (CWR) of domesticated apples (*Malus domestica*). 

## Workflow
Visit <https://terrellroulston.github.io/malus/> for a detailed reproducible workflow to reproduce the development of species distribution models (SDMs) using the Maxent machine learning algorithm within the `predicts` package, using presence/background occurrence data.

Also see a Shiny webtool to access the results of the SDM developed by Jens Ulrich at <https://julrich.shinyapps.io/CWR_SDM_Shiny/>.

## Citation (In Preparation)
This work is associated with a *in-prep* manuscript on *Malus* crop wild relatives (CWR) that are native to Canada, including *M. coronaria* (Sweet Crabapple) and *M. fusca* (Pacific Crabapple).


Roulston, T. T., Armstrong, C., Batstone, M., Bunsha, D., Ciotir, C., Borda, S. G., Husband, B., Manning, P., Moreau, T., Singh, A., Smith, T. W., Ulrich, J., & Migicovsky, Z. (*in prep*). *Conservation Challenges and Opportunities for Native Malus Species in Canada*.

## Background and funding
This work was funded in part by the UBC Botanical Garden, who is advancing action for adaptation through the development of community-based adaptation plans with the Sustainable Communities Field School. 

Ensuring the long-term preservation and sustainable utilization of apple genetic diversity across Canada is of paramount importance. This collaborative initiative aims to inventory and understand existing apple genetic diversity within the country. Additionally, it seeks to investigate how the suitability of habitats for apples may be affected by changing climatic conditions and emerging challenges. In response to these potential impacts, the project aims to identify urgent and longer-term adaptation strategies necessary for the continued resilience of apple populations in Canada.

## Scripts
All scripts for this analysis can be found in the `scripts` folder. 

Here is a brief overview of the steps and associated scripts for each part of the analysis. 

Step| Script| Purpose
-| ---| --------
1| `gbif_occ.R` | Prepare data request and download from GBIF
2| `occ_clean.R` | Clean occurrence (presence) data points
3| `occ_thin.R` | Thin occurrence points to a single observation per predictor (raster) cell
4| `malus_bg.R` | Download ecoregions, sample background points, predictor cluster, 3D Kernel density plots
5| `malus_sdm.R` | Climate predictor data, Maxent modeling and tuning, habitat suitability predictions
6| `sdm_plot.R` | Produce publication quality plots of SDM predictions
7| `malus_gap` | Conservation gap analysis to assess the level of conservation of species

## Data
###occ_data
The `occ_data` folder contains subfolders for the 'raw' CSV files downloaded from GBIF, as well as formated csv files saved as `occ_<species>.csv`, the cleanded occurrences are saved as Rdata files as `occ_xxx.Rdata`. Thinned occurrence points used for the sdm are saved as `occThin_xxx`. Additional files used only for the some background preliminary analysis is `xxx_sdmData` (combining bg points and predictor variables at occurrence points), which includes `xxx_bg` and `xxx_bg_vect`.

###maps
The `maps` folder contains a subfolder for GADM raster of admin base maps of Canada, Mexico and U.S., as well as extracted ecoregion rasters for each species.

###sdm_output
The `sdm_output` folder contains Rdata files for SpatRasters for the predictions from the SDMs. There are separate files for 2x species (*M. coronaria*, *M. fusca*) X 2x SSPs (SSP245, SSP585) X 4x time periods (historical, 2030, 2050, 2070). In the subfolder `thresholds` there are Rdata files for 10th, 5th and 1st percentile suitability thresholds. In the subfolder `habitat_predictions/high_moderate_low_predictions` there are thresholded rasters saved as Rdata, and the `cropped_predictions` as .tif for the Webtool.

###mess_data
The `mess_data` folder contains objects to compute similarity in the MESS analysis.

