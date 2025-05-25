# Malus Gap Analysis — Restricted by Ecoregions Containing Occurrences Masked by Historical Suitability 
# Terrell Roulston
# Started May 20th 2025

# Load Libraries ----------------------------------------------------------
library(terra)        # For raster operations
library(tidyverse)    # For dplyr, tibble, etc.
library(geodata)      # For loading GADM and other data


# Species Info ------------------------------------------------------------
species_list  <- c("cor", "fus", "ion", "ang", "chl")
species_names <- c("Malus coronaria", "Malus fusca", "Malus ioensis", "Malus angustifolia", "Sect. Chloromeles")


# Start Timer -------------------------------------------------------------
Start_all <- Sys.time()


# Store All Results -------------------------------------------------------
full_result <- list()

# Loop Through Each Species -----------------------------------------------
for (index in seq_along(species_list)) {
  sp_code <- species_list[index]
  sp_name <- species_names[index]

  cat("\n==== Running gap analysis for:", sp_code, "====\n")
  Start <- Sys.time()
  
  # ---- Set Paths ----
  terraOptions(threads = 15)
  base_path   <- file.path("./sdm_output", sp_code)
  pred_path   <- file.path(base_path, "subs")
  thresh_path <- file.path(pred_path, "threshold")
  occ_path    <- file.path("./occ_data", sp_code)
  
  # ---- Load and Prep Maps ----
  ca_bound <- geodata::gadm("CAN", level = 0, resolution = 1, path = "./maps/base_maps")
  us_bound <- geodata::gadm("USA", level = 0, resolution = 1, path = "./maps/base_maps")
  na_bound <- rbind(ca_bound, us_bound)
  pa_raster <- terra::rast("./gap_analysis/pa_raster_us_can.tif")
  
  eco_vec <- readRDS(file.path("./maps/eco_regions", paste0("ecoNA_", sp_code, ".Rdata")))
  eco_vec <- terra::project(eco_vec, pa_raster)
  
  occ <- readRDS(file.path(occ_path, paste0("occThin_", sp_code, ".Rdata")))
  if (!inherits(occ, "SpatVector")) occ <- terra::vect(occ)
  
  eco_mask_pts <- terra::intersect(eco_vec, occ)
  eco_occ <- unique(eco_mask_pts$NA_L2CODE)
  eco_vec_crop <- eco_vec[eco_vec$NA_L2CODE %in% eco_occ, ]
  
  # ---- Load Predictions and Thresholds ----
  preds <- list(
    hist      = readRDS(file.path(pred_path, paste0(sp_code, "_pred_hist_subs.Rdata"))),
    ssp245_30 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp245_30_subs.Rdata"))),
    ssp245_50 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp245_50_subs.Rdata"))),
    ssp245_70 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp245_70_subs.Rdata"))),
    ssp585_30 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp585_30_subs.Rdata"))),
    ssp585_50 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp585_50_subs.Rdata"))),
    ssp585_70 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp585_70_subs.Rdata")))
  )
  
  thresholds <- list(
    low  = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_1_subs.Rdata"))),
    mod  = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_10_subs.Rdata"))),
    high = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_50_subs.Rdata")))
  )
  
  # ---- Restrict to Ecoregions ----
  preds <- lapply(preds, function(x) terra::crop(x, eco_vec_crop, mask = TRUE))
  template <- preds[[1]]
  
  # ---- Historical Mask and Resampled PA ----
  hist_masked <- terra::ifel(preds[["hist"]] > thresholds[["high"]], 1, NA)
  pa_mask <- pa_raster == 1
  pa_mask_resampled <- terra::resample(pa_mask, template, method = "near")
  
  # ---- Metric Functions ----
  calculate_srsin <- function(pred, pa, occ, threshold) {
    r_masked <- pred > threshold
    names(r_masked) <- "binary"
    r_pa <- terra::mask(pred, pa_mask_resampled) > threshold
    names(r_pa) <- "binary"
    
    occ_pred <- terra::extract(r_masked, occ) %>% as.data.frame() %>% dplyr::filter(binary == 1)
    occ_pa   <- terra::extract(r_pa, occ)     %>% as.data.frame() %>% dplyr::filter(binary == 1)
    
    if (nrow(occ_pred) == 0) return(0)
    return(nrow(occ_pa) / nrow(occ_pred) * 100)
  }
  
  calculate_grsin <- function(pred, threshold, hist_masked) {
    r_masked <- (pred > threshold) & (hist_masked == 1)
    area_all <- terra::expanse(r_masked, byValue = TRUE, unit = "km") %>% 
      as.data.frame() %>% dplyr::filter(value == 1) %>% dplyr::pull(area)
    
    r_pa <- (terra::mask(pred, pa_mask_resampled) > threshold) & (hist_masked == 1)
    area_pa <- terra::expanse(r_pa, byValue = TRUE, unit = "km") %>% 
      as.data.frame() %>% dplyr::filter(value == 1) %>% dplyr::pull(area)
    
    if (length(area_all) == 0 || sum(area_all) == 0) return(0)
    return(sum(area_pa, na.rm = TRUE) / sum(area_all, na.rm = TRUE) * 100)
  }
  
  calculate_ersin <- function(pred, threshold, eco_vec, hist_masked) {
    r_masked <- (pred > threshold) & (hist_masked == 1)
    r_pa     <- (terra::mask(pred, pa_mask_resampled) > threshold) & (hist_masked == 1)
    
    eco_raster <- terra::rasterize(eco_vec, pred, field = "NA_L2CODE")
    eco_raster <- as.factor(eco_raster)
    
    eco_suit <- terra::mask(eco_raster, r_masked)
    eco_pa   <- terra::mask(eco_raster, r_pa)
    
    suit_ecos <- unique(na.omit(terra::values(eco_suit)))
    pa_ecos   <- unique(na.omit(terra::values(eco_pa)))
    
    if (length(suit_ecos) == 0) return(0)
    return(length(intersect(suit_ecos, pa_ecos)) / length(suit_ecos) * 100)
  }
  
  # ---- Run Gap Analysis ----
  out <- list()
  for (pname in names(preds)) {
    for (thresh_name in names(thresholds)) {
      th <- thresholds[[thresh_name]]
      pr <- preds[[pname]]
      
      cat("[INFO] Processing:", sp_code, pname, thresh_name, "\n")
      
      SRSin <- calculate_srsin(pr, pa_raster, occ, th)
      GRSin <- calculate_grsin(pr, th, hist_masked)
      ERSin <- calculate_ersin(pr, th, eco_vec, hist_masked)
      FCSin <- mean(c(SRSin, GRSin, ERSin))
      
      out[[paste(pname, thresh_name, sep = "_")]] <- tibble(
        species     = sp_name,
        sp_code     = sp_code,
        ssp         = ifelse(pname == 'hist', 'historical', ifelse(grepl('245', pname), '245', '585')),
        period      = case_when(
          pname == 'hist' ~ 2000,
          pname == 'ssp245_30' ~ 2030,
          pname == 'ssp245_50' ~ 2050,
          pname == 'ssp245_70' ~ 2070,
          pname == 'ssp585_30' ~ 2030,
          pname == 'ssp585_50' ~ 2050,
          pname == 'ssp585_70' ~ 2070
        ),
        suitability = thresh_name,
        SRSin       = SRSin,
        GRSin       = GRSin,
        ERSin       = ERSin,
        FCSin       = FCSin
      )
    }
  }
  
  # ---- Store per species ----
  result <- dplyr::bind_rows(out)
  full_result[[sp_code]] <- result
  
  End <- Sys.time()
  cat("Elapsed time for", sp_code, ":", End - Start, "\n")
}


# Combine Results ---------------------------------------------------------
End_all <- Sys.time()
cat("\nTotal elapsed time:", End_all - Start_all, "\n")

all_species_results <- dplyr::bind_rows(full_result)
readr::write_csv(all_species_results, "malus_gap_analysis_all_species.csv")
