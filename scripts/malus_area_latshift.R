# Top ---------------------------------------------------------------------
# Habitat Suitability Summarization
# Terrell Roulston
# Started May 5, 2025

library(tidyverse)
library(geodata)
library(terra)
library(dplyr)
library(tidyr)
library(purrr)


# Load Habitat Pred -------------------------------------------------------
# Only summarizing SSP585

# Malus fusca
fus_pred_hist <- readRDS(file = './sdm_output/fus/subs/fus_pred_hist_subs.Rdata')
fus_pred_ssp585_30 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_30_subs.Rdata')
fus_pred_ssp585_50 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_50_subs.Rdata')
fus_pred_ssp585_70 <- readRDS(file = './sdm_output/fus/subs/fus_pred_ssp585_70_subs.Rdata')
# Malus coronaria
cor_pred_hist <- readRDS(file = './sdm_output/cor/subs/cor_pred_hist_subs.Rdata')
cor_pred_ssp585_30 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_30_subs.Rdata')
cor_pred_ssp585_50 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_50_subs.Rdata')
cor_pred_ssp585_70 <- readRDS(file = './sdm_output/cor/subs/cor_pred_ssp585_70_subs.Rdata')
# Malus ioensis
ion_pred_hist <- readRDS(file = './sdm_output/ion/subs/ion_pred_hist_subs.Rdata')
ion_pred_ssp585_30 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_30_subs.Rdata')
ion_pred_ssp585_50 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_50_subs.Rdata')
ion_pred_ssp585_70 <- readRDS(file = './sdm_output/ion/subs/ion_pred_ssp585_70_subs.Rdata')
# Malus angustifolia
ang_pred_hist <- readRDS(file = './sdm_output/ang/subs/ang_pred_hist_subs.Rdata')
ang_pred_ssp585_30 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_30_subs.Rdata')
ang_pred_ssp585_50 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_50_subs.Rdata')
ang_pred_ssp585_70 <- readRDS(file = './sdm_output/ang/subs/ang_pred_ssp585_70_subs.Rdata')
# Sect. Chloromeles
chl_pred_hist <- readRDS(file = './sdm_output/chl/subs/chl_pred_hist_subs.Rdata')
chl_pred_ssp585_30 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_30_subs.Rdata')
chl_pred_ssp585_50 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_50_subs.Rdata')
chl_pred_ssp585_70 <- readRDS(file = './sdm_output/chl/subs/chl_pred_ssp585_70_subs.Rdata')

# Load Thresholds
# Only reporting high suitable area
# = 50th percentile suitability

# Malus fusca
fusPred_threshold_50 <- readRDS(file = './sdm_output/fus/subs/threshold/fusPred_threshold_50_subs.Rdata')
#Malus coronaria
corPred_threshold_50 <- readRDS(file = './sdm_output/cor/subs/threshold/corPred_threshold_50_subs.Rdata')
# Malus ioensis
ionPred_threshold_50 <- readRDS(file = './sdm_output/ion/subs/threshold/ionPred_threshold_50_subs.Rdata')
# Malus angustifolia
angPred_threshold_50 <- readRDS(file = './sdm_output/ang/subs/threshold/angPred_threshold_50_subs.Rdata')
# Sect Chloromeles
chlPred_threshold_50 <- readRDS(file = './sdm_output/chl/subs/threshold/chlPred_threshold_50_subs.Rdata')


# Area summarization ------------------------------------------------------
# Create function 
summarize_habitat <- function(pred_raster, threshold) {
  pred_bin <- pred_raster > threshold
  cell_area_km2 <- prod(res(pred_raster)) / 1e6  # convert from m² to km²
  
  coords_vals <- as.data.frame(terra::xyFromCell(pred_bin, which(values(pred_bin) == 1)))
  coords_vals$suitability <- values(pred_raster)[which(values(pred_bin) == 1)]
  
  total_area_km2 <- nrow(coords_vals) * cell_area_km2
  median_lat <- median(coords_vals$y)
  
  return(data.frame(total_area_km2 = total_area_km2, median_latitude = median_lat))
}

# Malus fusca
summarize_habitat(fus_pred_hist, fusPred_threshold_50)



# new function
compare_suitability <- function(t0_rast, t1_rast, threshold) {
  # 1. Threshold both rasters to binary
  t0_bin <- t0_rast >= threshold
  t1_bin <- t1_rast >= threshold
  
  # 2. Calculate true cell areas (km²) per pixel
  area_rast <- terra::cellSize(t0_rast, unit = "km")
  
  # 3. Identify valid (non-NA) cells in both layers
  valid_mask <- !is.na(t0_bin) & !is.na(t1_bin)
  valid_cells <- which(values(valid_mask) == 1)
  
  # 4. Extract values for valid cells
  t0_vals <- values(t0_bin)[valid_cells]
  t1_vals <- values(t1_bin)[valid_cells]
  area_vals <- values(area_rast)[valid_cells]
  
  # 5. Classify change status
  status <- dplyr::case_when(
    t0_vals == 1 & t1_vals == 1 ~ "stable",
    t0_vals == 1 & t1_vals == 0 ~ "contraction",
    t0_vals == 0 & t1_vals == 1 ~ "expansion",
    TRUE ~ "unsuitable"
  )
  
  # 6. Area per status
  area_by_change <- tibble(change = status, area_km2 = area_vals) |>
    group_by(change) |>
    summarize(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
  
  # 7. Total suitable areas
  total_area_t0 <- sum(area_vals[t0_vals == 1], na.rm = TRUE)
  total_area_t1 <- sum(area_vals[t1_vals == 1], na.rm = TRUE)
  
  total_area_summary <- tibble(
    time = c("t0", "t1"),
    total_area_km2 = c(total_area_t0, total_area_t1)
  )
  
  # 8. Latitude summaries (use coordinates from original raster)
  coords <- terra::xyFromCell(t0_rast, valid_cells)
  lat_t0 <- coords[t0_vals == 1, 2]
  lat_t1 <- coords[t1_vals == 1, 2]
  
  lat_summary <- tibble(
    time = c("t0", "t1"),
    median_latitude = c(median(lat_t0), median(lat_t1)),
    mean_latitude = c(mean(lat_t0), mean(lat_t1))
  )
  
  return(list(
    area_by_change = area_by_change,
    total_area_summary = total_area_summary,
    lat_summary = lat_summary
  ))
}


compare_suitability(cor_pred_hist, cor_pred_ssp585_30, fusPred_threshold_50)

# COMPLETE SUMMARY
# List of all inputs ------------------------------------------------------

# Updated taxa list including 2050 and 2070
taxa_list <- list(
  fus = list(hist = fus_pred_hist, ssp30 = fus_pred_ssp585_30, ssp50 = fus_pred_ssp585_50, ssp70 = fus_pred_ssp585_70, threshold = fusPred_threshold_50),
  cor = list(hist = cor_pred_hist, ssp30 = cor_pred_ssp585_30, ssp50 = cor_pred_ssp585_50, ssp70 = cor_pred_ssp585_70, threshold = corPred_threshold_50),
  ion = list(hist = ion_pred_hist, ssp30 = ion_pred_ssp585_30, ssp50 = ion_pred_ssp585_50, ssp70 = ion_pred_ssp585_70, threshold = ionPred_threshold_50),
  ang = list(hist = ang_pred_hist, ssp30 = ang_pred_ssp585_30, ssp50 = ang_pred_ssp585_50, ssp70 = ang_pred_ssp585_70, threshold = angPred_threshold_50),
  chl = list(hist = chl_pred_hist, ssp30 = chl_pred_ssp585_30, ssp50 = chl_pred_ssp585_50, ssp70 = chl_pred_ssp585_70, threshold = chlPred_threshold_50)
)

# Wrapper to summarize all time comparisons per taxon
# Long format summary function
summarize_long_format <- function(taxon_code, hist, ssp30, ssp50, ssp70, threshold) {
  
  # Step 1: t0 baseline row with total area and latitude summary
  t0_vals <- compare_suitability(hist, hist, threshold) # dummy comparison
  t0_total_area <- t0_vals$total_area_summary %>%
    filter(time == "t0") %>%
    select(total_area_km2)
  
  t0_lat <- t0_vals$lat_summary %>%
    filter(time == "t0") %>%
    select(median_latitude, mean_latitude)
  
  baseline_row <- tibble(
    taxon = taxon_code,
    timeseries = "t0",
    change = "total",
    area_km2 = t0_total_area$total_area_km2,
    median_latitude = t0_lat$median_latitude,
    mean_latitude = t0_lat$mean_latitude
  )
  
  # Step 2: compare each future timepoint to t0
  future_layers <- list(t1 = ssp30, t2 = ssp50, t3 = ssp70)
  
  future_rows <- imap_dfr(future_layers, function(layer, label) {
    res <- compare_suitability(hist, layer, threshold)
    
    lat <- res$lat_summary %>% filter(time == "t1")
    
    res$area_by_change %>%
      filter(change %in% c("contraction", "expansion", "stable")) %>%
      mutate(taxon = taxon_code,
             timeseries = label,
             median_latitude = lat$median_latitude,
             mean_latitude = lat$mean_latitude) %>%
      select(taxon, timeseries, change, area_km2, median_latitude, mean_latitude)
  })
  
  bind_rows(baseline_row, future_rows)
}

# Apply to all taxa
summary_long_df <- purrr::pmap_dfr(
  list(
    taxon_code = names(taxa_list),
    hist = map(taxa_list, "hist"),
    ssp30 = map(taxa_list, "ssp30"),
    ssp50 = map(taxa_list, "ssp50"),
    ssp70 = map(taxa_list, "ssp70"),
    threshold = map(taxa_list, "threshold")
  ),
  summarize_long_format
)

# View result
head(summary_long_df)

# Pivot wide for publication ----------------------------------------------



# Step 1: Extract t0 summaries (where change == "total")
df_t0 <- summary_long_df %>%
  filter(timeseries == "t0", change == "total") %>%
  select(taxon, area_km2, median_latitude, mean_latitude) %>%
  rename(
    total_area_t0 = area_km2,
    median_latitude_t0 = median_latitude,
    mean_latitude_t0 = mean_latitude
  )

# Step 2: Extract change values for t1–t3
df_changes <- summary_long_df %>%
  filter(timeseries != "t0") %>%
  select(taxon, timeseries, change, area_km2) %>%
  pivot_wider(
    names_from = c(timeseries, change),
    values_from = area_km2,
    values_fill = 0
  )

# Step 3: Extract lat summaries for t1–t3 (from any row per group since values are the same per group)
df_latitudes <- summary_long_df %>%
  filter(timeseries != "t0") %>%
  group_by(taxon, timeseries) %>%
  summarize(
    median_lat = median(median_latitude, na.rm = TRUE),
    mean_lat = median(mean_latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = timeseries,
    values_from = c(median_lat, mean_lat)
  ) %>%
  rename_with(~ paste0(.x, "_t", gsub(".*_", "", .x)), starts_with("median_lat")) %>%
  rename_with(~ paste0(.x, "_t", gsub(".*_", "", .x)), starts_with("mean_lat"))

# Step 4: Join all together
summary_wide_df <- df_t0 %>%
  left_join(df_changes, by = "taxon") %>%
  left_join(df_latitudes, by = "taxon")

# View result
print(summary_wide_df)


# Save for publication ----------------------------------------------------
# Save as CSV

write.csv(summary_wide_df, file = './sdm_output/summ_area/summary_wide_df.csv')


# Exploratory maps to make sure I am not going crazy ----------------------
cor_30_vect <- terra::as.polygons(cor_pred_ssp585_30 > corPred_threshold_50)
cor_50_vect <- terra::as.polygons(cor_pred_ssp585_50 > corPred_threshold_50)
cor_70_vect <- terra::as.polygons(cor_pred_ssp585_70 > corPred_threshold_50)



terra::plot(cor_pred_hist > corPred_threshold_50, col = c("#E8E8E8", 'black'), legend = F)
terra::plot(cor_30_vect, border = c("#FFFFFF00", 'red'), add = T, legend = T)
terra::plot(cor_50_vect, border = c("#FFFFFF00", 'blue'), add = T, legend = T)
terra::plot(cor_70_vect, border = c("#FFFFFF00", 'green'), add = T, legend = T)


