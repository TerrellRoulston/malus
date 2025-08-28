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


# Load ecoregions to constrain area calculations
ecoNA_fus <- readRDS(file = './maps/eco_regions/ecoNA_fus.Rdata')
ecoNA_cor <- readRDS(file = './maps/eco_regions/ecoNA_cor.Rdata')
ecoNA_ion <- readRDS(file = './maps/eco_regions/ecoNA_ion.Rdata')
ecoNA_ang <- readRDS(file = './maps/eco_regions/ecoNA_ang.Rdata')
ecoNA_chl <- readRDS(file = './maps/eco_regions/ecoNA_chl.Rdata')

# Load base ecoregion shape file and select ecoregions that are expanded into in the future!
# Cross refereing visually the predicted areas of future suitability that are north of the historical
# with the ecoregion map and selected ones that are expanded into northward and adjacent to recent ecoregion
ecoNA <- vect(x = "maps/eco_regions/na_cec_eco_l2/NA_CEC_Eco_Level2.shp")
ecoNA <- project(ecoNA, 'WGS84') # project ecoregion vector to same coords ref as basemap


# M. fusca
# Add 6.1?? a small area of suitability in southern Alaska is missed under 2070
# Historic: "7.1""6.2"  "10.1" "11.1" "10.2"
# Addition: 6.1, 2.2
eco_fus_code_add <- c("7.1", "6.2", "10.1", "11.1", "10.2", "6.1", "2.2")
ecoNA_fus_add <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_fus_code_add)

# M. coronaria
# Historic: "8.1" "8.2" "5.3" "8.4" "8.3" "8.5" "9.2" "9.4" "5.2"
# Addition: 5.1 and 3.4
eco_cor_code_add <- c("8.1", "8.2", "5.3", "8.4", "8.3", "8.5", "9.2", "9.4", "5.2", "5.1", '3.4')
ecoNA_cor_add <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_cor_code_add)

# M. ioensis
# Historic: "5.2" "8.1" "8.2" "8.3" "8.4" "8.5" "9.2" "9.4"
# Addition: 5.1, 4.1, 5.4 and 3.4
eco_ion_code_add <- c("5.2", "8.1", "8.2", "8.3", "8.4", "8.5", "9.2", "9.4", "5.1", "4.1", "5.4", "3.4")
ecoNA_ion_add <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_ion_code_add)

# M. angustifolia
# Historic: "5.3" "8.1" "8.2" "8.3" "8.4" "8.5" "9.5"
# Addition: "9.2", "9.4", "5.2", "5.1", '3.4', 1.1 (adding more than need i think but mirroring coronaria)
eco_ang_code_add <- c("5.3", "8.1", "8.2", "8.3", "8.4", "8.5", "9.5", "9.2", "9.4", "5.2", "5.1", '3.4', "1.1")
ecoNA_ang_add <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_ang_code_add)

# Sect. Chloromeles
# Historic: "5.3" "8.1" "8.2" "8.3" "8.4" "8.5" "9.5"
# Addition: "5.2" "9.2" "9.4" "5.1" "4.1" "5.4" "3.4" "1.1"
eco_chl_code_add <- c("5.2", "8.1", "8.2", "8.3", "8.4", "8.5", "9.2", "9.4", "5.1", "4.1", "5.4", "3.4", "5.3", "9.5", "1.1")
ecoNA_chl_add <- terra::subset(ecoNA, ecoNA$NA_L2CODE %in% eco_chl_code_add) 


# Build masked rasters for each species -----------------------------------

# fus: hist masked to ecoNA_fus; futures masked to ecoNA_fus_add
fus_hist_m <- mask(fus_pred_hist, ecoNA_fus)
fus_30_m <- mask(fus_pred_ssp585_30, ecoNA_fus_add)
fus_50_m <- mask(fus_pred_ssp585_50, ecoNA_fus_add)
fus_70_m  <- mask(fus_pred_ssp585_70, ecoNA_fus_add)

# cor
cor_hist_m <- mask(cor_pred_hist, ecoNA_cor)
cor_30_m <- mask(cor_pred_ssp585_30, ecoNA_cor_add)
cor_50_m <- mask(cor_pred_ssp585_50, ecoNA_cor_add)
cor_70_m <- mask(cor_pred_ssp585_70, ecoNA_cor_add)

# ion
ion_hist_m <- mask(ion_pred_hist, ecoNA_ion)
ion_30_m <- mask(ion_pred_ssp585_30, ecoNA_ion_add)
ion_50_m <- mask(ion_pred_ssp585_50, ecoNA_ion_add)
ion_70_m <- mask(ion_pred_ssp585_70, ecoNA_ion_add)

# ang
ang_hist_m <- mask(ang_pred_hist, ecoNA_ang)
ang_30_m <- mask(ang_pred_ssp585_30, ecoNA_ang_add)
ang_50_m <- mask(ang_pred_ssp585_50, ecoNA_ang_add)
ang_70_m <- mask(ang_pred_ssp585_70, ecoNA_ang_add)

# chl
chl_hist_m <- mask(chl_pred_hist, ecoNA_chl)
chl_30_m <- mask(chl_pred_ssp585_30, ecoNA_chl_add)
chl_50_m <- mask(chl_pred_ssp585_50, ecoNA_chl_add)
chl_70_m <- mask(chl_pred_ssp585_70, ecoNA_chl_add)

# Area summarization ------------------------------------------------------
# Create function 
summarize_habitat <- function(pred_raster, threshold) {
  pred_bin <- pred_raster > threshold
  cell_area_km2 <- terra::cellSize(pred_raster, unit = "km")  # Geodesic per-cell areas (km²) for a lon/lat raster
  total_area_km2 <- as.numeric(terra::global(cell_area_km2 * (pred_bin == 1), "sum", na.rm = TRUE)[1,1])
  idx <- which(pred_bin[] == 1)
  ys  <- terra::yFromCell(pred_raster, idx)
  data.frame(total_area_km2 = total_area_km2,
             median_latitude = median(ys),
             stdev_latitude  = sd(ys))
}
# Malus fusca
summarize_habitat(fus_hist_m, fusPred_threshold_50)


compare_suitability <- function(t0_rast, t1_rast, threshold) {
  # Align rasters to the union extent (same grid; same res/alignment assumed)
  ext_u <- terra::union(terra::ext(t0_rast), terra::ext(t1_rast))
  r0 <- terra::extend(t0_rast, ext_u)
  r1 <- terra::extend(t1_rast, ext_u)
  
  # Threshold
  b0 <- r0 >= threshold
  b1 <- r1 >= threshold
  
  # Union of valid cells
  valid <- !is.na(b0) | !is.na(b1)
  idx <- which(valid[])
  
  # Treat NA as FALSE
  t0v <- as.logical(b0[])[idx]; t0v[is.na(t0v)] <- FALSE
  t1v <- as.logical(b1[])[idx]; t1v[is.na(t1v)] <- FALSE
  
  # Geodesic per-cell areas (prefer r0; fall back to r1)
  a0 <- terra::cellSize(r0, unit = "km")[]
  a1 <- terra::cellSize(r1, unit = "km")[]
  area <- ifelse(!is.na(a0), a0, a1)[idx]
  
  # Change classes
  status <- dplyr::case_when(
    t0v &  t1v ~ "stable",
    t0v & !t1v ~ "contraction",
    !t0v & t1v ~ "expansion",
    TRUE       ~ "unsuitable"
  )
  
  area_by_change <- tibble::tibble(change = status, area_km2 = area) |>
    dplyr::group_by(change) |>
    dplyr::summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
  
  total_area_summary <- tibble::tibble(
    time = c("t0","t1"),
    total_area_km2 = c(sum(area[t0v], na.rm = TRUE),
                       sum(area[t1v], na.rm = TRUE))
  )
  
  # Latitude summaries
  coords <- terra::xyFromCell(r0, idx)
  lat_summary <- tibble::tibble(
    time = c("t0","t1"),
    median_latitude = c(if (any(t0v)) median(coords[t0v,2]) else NA_real_,
                        if (any(t1v)) median(coords[t1v,2]) else NA_real_),
    mean_latitude   = c(if (any(t0v)) mean(coords[t0v,2])   else NA_real_,
                        if (any(t1v)) mean(coords[t1v,2])   else NA_real_),
    sd_latitude     = c(if (any(t0v)) sd(coords[t0v,2])     else NA_real_,
                        if (any(t1v)) sd(coords[t1v,2])     else NA_real_)
  )
  
  list(area_by_change = area_by_change,
       total_area_summary = total_area_summary,
       lat_summary = lat_summary)
}

compare_suitability(fus_hist_m, fus_70_m, fusPred_threshold_50)

# COMPLETE SUMMARY
# List of all inputs ------------------------------------------------------

# Use masked rasters (hist → ecoNA_*, futures → ecoNA_*_add)
taxa_list <- list(
  fus = list(hist = fus_hist_m, ssp30 = fus_30_m, ssp50 = fus_50_m, ssp70 = fus_70_m, threshold = fusPred_threshold_50),
  cor = list(hist = cor_hist_m, ssp30 = cor_30_m, ssp50 = cor_50_m, ssp70 = cor_70_m, threshold = corPred_threshold_50),
  ion = list(hist = ion_hist_m, ssp30 = ion_30_m, ssp50 = ion_50_m, ssp70 = ion_70_m, threshold = ionPred_threshold_50),
  ang = list(hist = ang_hist_m, ssp30 = ang_30_m, ssp50 = ang_50_m, ssp70 = ang_70_m, threshold = angPred_threshold_50),
  chl = list(hist = chl_hist_m, ssp30 = chl_30_m, ssp50 = chl_50_m, ssp70 = chl_70_m, threshold = chlPred_threshold_50)
)

# Wrapper to summarize all time comparisons per taxon
# Long format summary function
# Wrapper to summarize all time comparisons per taxon
summarize_long_format <- function(taxon_code, hist, ssp30, ssp50, ssp70, threshold) {
  
  # --- Baseline (hist vs hist): total + lat summaries
  t0_vals <- compare_suitability(hist, hist, threshold)  # masked hist already
  t0_total_area <- t0_vals$total_area_summary %>%
    dplyr::filter(time == "t0") %>%
    dplyr::pull(total_area_km2) %>% .[1]
  
  t0_lat <- t0_vals$lat_summary %>%
    dplyr::filter(time == "t0") %>%
    dplyr::slice(1)
  
  baseline_row <- tibble::tibble(
    taxon = taxon_code,
    timeseries = "t0",
    change = "total",
    area_km2 = t0_total_area,
    median_latitude = t0_lat$median_latitude,
    mean_latitude = t0_lat$mean_latitude,
    sd_latitude = t0_lat$sd_latitude
  )
  
  # --- Futures relative to t0
  future_layers <- list(ssp30 = ssp30, ssp50 = ssp50, ssp70 = ssp70)
  
  future_rows <- purrr::imap_dfr(future_layers, function(layer, label) {
    res <- compare_suitability(hist, layer, threshold)
    
    lat <- res$lat_summary %>% dplyr::filter(time == "t1") %>% dplyr::slice(1)
    
    changes <- res$area_by_change %>%
      dplyr::filter(change %in% c("contraction", "expansion", "stable")) %>%
      dplyr::mutate(
        taxon = taxon_code,
        timeseries = label,
        median_latitude = lat$median_latitude,
        mean_latitude = lat$mean_latitude,
        sd_latitude = lat$sd_latitude
      ) %>%
      dplyr::select(taxon, timeseries, change, area_km2,
                    median_latitude, mean_latitude, sd_latitude)
    
    # Optional: include total suitable area at the future timepoint
    total_future <- res$total_area_summary %>%
      dplyr::filter(time == "t1") %>%
      dplyr::transmute(
        taxon = taxon_code,
        timeseries = label,
        change = "total",
        area_km2 = total_area_km2,
        median_latitude = lat$median_latitude,
        mean_latitude = lat$mean_latitude,
        sd_latitude = lat$sd_latitude
      )
    
    dplyr::bind_rows(changes, total_future)
  })
  
  dplyr::bind_rows(baseline_row, future_rows)
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
  dplyr::filter(timeseries == "t0", change == "total") %>%
  dplyr::select(taxon, area_km2, median_latitude, mean_latitude, sd_latitude) %>%
  dplyr::rename(
    total_area_t0 = area_km2,
    median_latitude_t0 = median_latitude,
    mean_latitude_t0 = mean_latitude,
    sd_latitude_t0 = sd_latitude
  )

# Step 2: Extract change values for t1–t3
df_changes <- summary_long_df %>%
  dplyr::filter(timeseries != "t0") %>%
  dplyr::select(taxon, timeseries, change, area_km2) %>%
  pivot_wider(
    names_from = c(timeseries, change),
    values_from = area_km2,
    values_fill = 0
  )

# Step 3: Extract lat summaries for t1–t3 (from any row per group since values are the same per group)
df_latitudes <- summary_long_df %>%
  dplyr::filter(timeseries != "t0") %>%
  dplyr::group_by(taxon, timeseries) %>%
  dplyr::slice(1) %>%  # each row per taxon/time has same lat summary, just grab the first
  dplyr::ungroup() %>%
  dplyr::select(taxon, timeseries, median_latitude, mean_latitude, sd_latitude) %>%
  pivot_wider(
    names_from = timeseries,
    values_from = c(median_latitude, mean_latitude, sd_latitude)
  ) %>%
  dplyr::rename_with(~ paste0(.x, "_t", gsub(".*_", "", .x)), starts_with("median_latitude")) %>%
  dplyr::rename_with(~ paste0(.x, "_t", gsub(".*_", "", .x)), starts_with("mean_latitude")) %>%
  dplyr::rename_with(~ paste0(.x, "_t", gsub(".*_", "", .x)), starts_with("sd_latitude"))

# Step 4: Join all together
summary_wide_df <- df_t0 %>%
  left_join(df_changes, by = "taxon") %>%
  left_join(df_latitudes, by = "taxon")

# View result
print(summary_wide_df)


# Save for publication ----------------------------------------------------
# Save as CSV

write.csv(summary_wide_df, file = './sdm_output/summ_area/summary_wide_df_v5.csv')


