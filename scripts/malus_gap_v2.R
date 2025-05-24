# Malus Gap Analysis — Future Projections Masked by Historical Suitability

library(terra)
library(tidyverse)
library(geodata)

# Species list
species_list <- c("cor", "fus", "ion", "ang", "chl")
species_names <- c("Malus coronaria", "Malus fusca", "Malus ioensis", "Malus angustifolia", "Sect. Chloromeles")

# Local selector
index <- 1
sp_code <- species_list[index]
sp_name <- species_names[index]
Start <- Sys.time()

# Threads and paths
terraOptions(threads = 15)
base_path   <- file.path("./sdm_output", sp_code)
pred_path   <- file.path(base_path, "subs")
thresh_path <- file.path(pred_path, "threshold")
occ_path    <- file.path("./occ_data", sp_code)

# Load NA boundaries and protected area raster
ca_bound <- geodata::gadm("CAN", level = 0, resolution = 1, path = "./maps/base_maps")
us_bound <- geodata::gadm("USA", level = 0, resolution = 1, path = "./maps/base_maps")
na_bound <- rbind(ca_bound, us_bound)
pa_raster <- terra::rast("./gap_analysis/pa_raster_us_can.tif")

# Load ecoregions and project
eco_vec <- readRDS(file.path("./maps/eco_regions", paste0("ecoNA_", sp_code, ".Rdata")))
eco_vec <- terra::project(eco_vec, pa_raster)

# Load occurrences
occ <- readRDS(file.path(occ_path, paste0("occThin_", sp_code, ".Rdata")))
if (!inherits(occ, "SpatVector")) occ <- terra::vect(occ)

# Restrict to ecoregions with occurrences
eco_mask_pts <- terra::intersect(eco_vec, occ)
eco_occ <- unique(eco_mask_pts$NA_L2CODE)
eco_vec_crop <- eco_vec[eco_vec$NA_L2CODE %in% eco_occ, ]

# Load predictions and thresholds
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

# Restrict predictions to known ecoregions
preds <- lapply(preds, function(x) terra::crop(x, eco_vec_crop, mask = TRUE))

# PA mask aligned to prediction raster
template <- preds[[1]]
pa_mask <- pa_raster == 1
pa_mask_resampled <- terra::resample(pa_mask, template, method = "near")

# Historical binary mask (used for masking future projections)
r_hist_masked <- preds[["hist"]] > thresholds[["low"]]

# ---- Metric Functions ----
calculate_srsin <- function(pred, occ, threshold) {
  r_masked <- pred > threshold
  names(r_masked) <- "binary"
  r_pa <- terra::mask(pred, pa_mask_resampled) > threshold
  names(r_pa) <- "binary"
  
  occ_pred <- terra::extract(r_masked, occ) %>% as.data.frame() %>% filter(binary == 1)
  occ_pa   <- terra::extract(r_pa, occ)     %>% as.data.frame() %>% filter(binary == 1)
  
  if (nrow(occ_pred) == 0) return(0)
  return(nrow(occ_pa) / nrow(occ_pred) * 100)
}

calculate_grsin <- function(pred, threshold) {
  r_masked <- pred > threshold
  r_masked <- terra::mask(r_masked, r_hist_masked)
  area_all <- terra::expanse(r_masked, byValue = TRUE, unit = "km") %>% as.data.frame() %>% filter(value == 1) %>% pull(area)
  
  r_pa <- terra::mask(pred, pa_mask_resampled) > threshold
  r_pa <- terra::mask(r_pa, r_hist_masked)
  area_pa <- terra::expanse(r_pa, byValue = TRUE, unit = "km") %>% as.data.frame() %>% filter(value == 1) %>% pull(area)
  
  if (length(area_all) == 0 || sum(area_all) == 0) return(0)
  return(sum(area_pa, na.rm = TRUE) / sum(area_all, na.rm = TRUE) * 100)
}

calculate_ersin <- function(pred, threshold, eco_vec) {
  r_masked <- pred > threshold
  r_masked <- terra::mask(r_masked, r_hist_masked)
  r_pa     <- terra::mask(pred, pa_mask_resampled) > threshold
  r_pa     <- terra::mask(r_pa, r_hist_masked)
  
  eco_raster <- terra::rasterize(eco_vec, pred, field = "NA_L2CODE")
  eco_raster <- as.factor(eco_raster)
  
  eco_suit <- terra::mask(eco_raster, r_masked)
  eco_pa   <- terra::mask(eco_raster, r_pa)
  
  suit_ecos <- unique(terra::values(eco_suit)) %>% na.omit()
  pa_ecos   <- unique(terra::values(eco_pa))   %>% na.omit()
  
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
    
    SRSin <- calculate_srsin(pr, occ, th)
    GRSin <- calculate_grsin(pr, th)
    ERSin <- calculate_ersin(pr, th, eco_vec)
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

result <- bind_rows(out)
End <- Sys.time()
cat("Total elapsed:", End - Start, "\n")
print(result)
