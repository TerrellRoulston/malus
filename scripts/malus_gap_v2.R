# Malus gap analysis V2
# I have been trying this anaylsis on the HPC cluster and its very quick, I think I can run it locally
# The analysis did complete initially, however the output was a bit suspicous
# So, I am going to debug locally
# And hopefully this runs...

library(terra)
library(tidyverse)
library(geodata)

# Species names and codes
species_list <- c("cor", "fus", "ion", "ang", "chl")
species_names <- c("Malus coronaria", "Malus fusca", "Malus ioensis", "Malus angustifolia", "Sect. Chloromeles")

# Accept SLURM_ARRAY_TASK_ID or local test arg
args <- commandArgs(trailingOnly = TRUE)
index <- if (length(args) > 0) as.integer(args[1]) + 1 else 1

# Debug Array Job - Fatal Halt
if (index < 1 || index > length(species_list)) {
  stop("[FATAL] Invalid index: ", index)
}

# Index species
sp_code <- species_list[index]
sp_name <- species_names[index]

cat("\n===== STARTING: ", sp_name, " =====\n")

# Set terra threading from allocated CPUs
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 4))
terraOptions(threads = n_threads)
cat("[INFO] Using", n_threads, "threads for terra\n")

# Paths
base_path   <- file.path("./sdm_output", sp_code)
pred_path   <- file.path(base_path, "subs")
thresh_path <- file.path(pred_path, "threshold")
occ_path    <- file.path("./occ_data", sp_code)

cat("[INFO] Loading North America boundaries...\n")
ca_bound <- gadm("CAN", level = 0, resolution = 1, path = "./maps/base_maps")
us_bound <- gadm("USA", level = 0, resolution = 1, path = "./maps/base_maps")
na_bound <- rbind(ca_bound, us_bound)

cat("[INFO] Loading pre-saved protected areas raster layer...\n")
pa_raster <- terra::rast("./gap_analysis/pa_raster_us_can.tif")

# ---- Analysis Function ----
process_species <- function(sp_code, sp_name, na_bound, pa_raster) {
  # Load predictions
  preds <- list(
    hist      = readRDS(file.path(pred_path, paste0(sp_code, "_pred_hist_subs.Rdata"))),
    ssp245_30 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp245_30_subs.Rdata"))),
    ssp245_50 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp245_50_subs.Rdata"))),
    ssp245_70 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp245_70_subs.Rdata"))),
    ssp585_30 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp585_30_subs.Rdata"))),
    ssp585_50 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp585_50_subs.Rdata"))),
    ssp585_70 = readRDS(file.path(pred_path, paste0(sp_code, "_pred_ssp585_70_subs.Rdata")))
  )
  
  # Load thresholds
  thresholds <- list(
    low  = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_1_subs.Rdata"))),
    mod  = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_10_subs.Rdata"))),
    high = readRDS(file.path(thresh_path, paste0(sp_code, "Pred_threshold_50_subs.Rdata")))
  )
  
  # Load occurrence data
  occ <- readRDS(file.path(occ_path, paste0("occThin_", sp_code, ".Rdata")))
  
  # Load species-specific ecoregion vector and rasterize
  cat("[INFO] Rasterizing species-specific ecoregion layer...\n")
  eco_vec <- readRDS(file.path("./maps/eco_regions", paste0("ecoNA_", sp_code, ".Rdata")))
  eco_codes <- sort(unique(eco_vec$NA_L2CODE))
  eco_lut <- setNames(seq_along(eco_codes), eco_codes)
  eco_vec$eco_id <- eco_lut[eco_vec$NA_L2CODE]
  
  # Ensure CRS match
  cat("[DEBUG] CRS before project:\n")
  print(terra::crs(eco_vec))
  print(terra::crs(pa_raster))
  eco_vec <- terra::project(eco_vec, pa_raster)
  cat("[DEBUG] CRS after project:\n")
  print(terra::crs(eco_vec))
  
  # Check raster status
  cat("[DEBUG] pa_raster summary:\n")
  print(summary(pa_raster))
  cat("[DEBUG] pa_raster extent:\n")
  print(terra::ext(pa_raster))
  cat("[DEBUG] eco_vec extent:\n")
  print(terra::ext(eco_vec))
  
  if (is.null(terra::ext(pa_raster)) || all(is.na(terra::values(pa_raster)))) {
    stop("[FATAL] pa_raster appears to be empty or has no non-NA values.")
  }
  
  # Create a raster template with same resolution/alignment as pa_raster, but cropped to eco_vec extent
  template_raster <- terra::crop(pa_raster, terra::ext(eco_vec))
  
  cat("[DEBUG] template_raster summary:\n")
  print(summary(template_raster))
  cat("[DEBUG] template_raster extent:\n")
  print(terra::ext(template_raster))
  
  # Rasterize ecoregion vector using aligned template
  eco_raster <- terra::rasterize(eco_vec, template_raster, field = "eco_id", background = NA)
  
  # Begin loop over predictions and thresholds
  out <- list()
  for (pname in names(preds)) {
    for (thresh_name in names(thresholds)) {
      cat("[INFO] ", sp_code, ":", pname, "-", thresh_name, "\n")
      
      th <- thresholds[[thresh_name]]
      r_sdm <- terra::crop(preds[[pname]], na_bound, mask = TRUE)
      r_masked <- r_sdm > th
      
      pa_raster_aligned <- terra::resample(pa_raster, r_masked, method = "near")
      r_pa <- terra::mask(r_masked, pa_raster_aligned, maskvalues = NA)
      
      occ_vals    <- terra::extract(r_masked, occ)[[1]]
      occ_vals_pa <- terra::extract(r_pa, occ)[[1]]
      SRSin <- ifelse(sum(!is.na(occ_vals)) == 0, 0,
                      sum(!is.na(occ_vals_pa)) / sum(!is.na(occ_vals)) * 100)
      
      area_total <- terra::expanse(r_masked, unit = "km", byValue = TRUE) %>%
        as.data.frame() %>%
        dplyr::filter(value == 1) %>%
        dplyr::pull(area)
      
      area_pa <- terra::expanse(r_pa, unit = "km", byValue = TRUE) %>%
        as.data.frame() %>%
        dplyr::filter(value == 1) %>%
        dplyr::pull(area)
      
      GRSin <- ifelse(length(area_total) == 0, 0, area_pa / area_total * 100)
      
      r_masked_aligned <- terra::resample(r_masked, eco_raster, method = "near")
      r_pa_aligned     <- terra::resample(r_pa, eco_raster, method = "near")
      
      eco_total <- unique(na.omit(terra::values(terra::mask(eco_raster, r_masked_aligned, maskvalues = 0))))
      eco_total <- eco_total[eco_total != 0.0]
      
      eco_pa <- unique(na.omit(terra::values(terra::mask(eco_raster, r_pa_aligned, maskvalues = 0))))
      eco_pa <- eco_pa[eco_pa != 0.0]
      
      ERSin <- ifelse(length(eco_total) == 0, 0, length(eco_pa) / length(eco_total) * 100)
      FCSin <- mean(c(SRSin, GRSin, ERSin))
      
      out[[paste(pname, thresh_name, sep = "_")]] <- tibble::tibble(
        species     = sp_name,
        sp_code     = sp_code,
        ssp         = ifelse(pname == 'hist', 'historical', ifelse(grepl('245', pname), '245', '585')),
        period      = dplyr::case_when(
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
  
  dplyr::bind_rows(out)
}


# Run
result <- process_species(sp_code, sp_name, na_bound, pa_raster)

# Save
outfile <- file.path("./gap_analysis", paste0("result_", sp_code, ".rds"))
cat("[INFO] Writing result to:", outfile, "\n")
saveRDS(result, outfile)

cat("===== DONE: ", sp_name, " =====\n")
cat("[INFO] ELASPED TIME:", Sys.time() - Start, "\n")
