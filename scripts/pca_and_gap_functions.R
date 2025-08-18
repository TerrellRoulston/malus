# This script contains the functions for running the PCA analysis
# This functions were adapted from the ade4 and ecospat packages!

# Niche Equivalency Test (PCA-based)
# Permutational Test of Importance (pseudo p-value)
# Species occurrences are randomized in repeteted runs
# Build on ecospat.niche.overlap function from ecospat

niche_equivalency_test <- function(sp1_scores, sp2_scores,
                                   sp1_bg_scores, sp2_bg_scores,
                                   bg_scores, R = 100, reps = 999,
                                   quant = 0.1,
                                   parallel = TRUE, ncores = 2,
                                   verbose = TRUE) {
  require(ecospat)
  require(foreach)
  start_time <- Sys.time()  # Start timer
  
  if (parallel) {
    require(doParallel)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
  }
  
  # Observed overlap grid
  grid1 <- ecospat.grid.clim.dyn(glob = bg_scores, glob1 = sp1_bg_scores, sp = sp1_scores, R = R)
  grid2 <- ecospat.grid.clim.dyn(glob = bg_scores, glob1 = sp2_bg_scores, sp = sp2_scores, R = R)
  
  # Try calculating observed overlap safely
  obs_vals <- tryCatch(
    ecospat.niche.overlap(grid1, grid2, cor = FALSE),
    error = function(e) return(NULL)
  )
  
  if (is.null(obs_vals) || !is.finite(obs_vals$D) || !is.finite(obs_vals$I)) {
    stop("Observed overlap metrics could not be calculated. Check if grid1/grid2 are valid.")
  }
  
  obs_D <- obs_vals$D
  obs_I <- obs_vals$I
  
  # Combine all occurrences
  all_occ <- rbind(sp1_scores, sp2_scores)
  n1 <- nrow(sp1_scores)
  n2 <- nrow(sp2_scores)
  
  permute_fun <- function(i) {
    set.seed(i)
    sampled <- all_occ[sample(nrow(all_occ)), ]
    sp1_null <- sampled[1:n1, ]
    sp2_null <- sampled[(n1 + 1):(n1 + n2), ]
    
    g1 <- tryCatch(ecospat.grid.clim.dyn(glob = bg_scores, glob1 = sp1_bg_scores, sp = sp1_null, R = R), error = function(e) return(NULL))
    g2 <- tryCatch(ecospat.grid.clim.dyn(glob = bg_scores, glob1 = sp2_bg_scores, sp = sp2_null, R = R), error = function(e) return(NULL))
    
    if (is.null(g1) || is.null(g2)) return(c(D = NA, I = NA))
    
    # Use uncorrected occupancy
    vals <- tryCatch(ecospat.niche.overlap(g1, g2, cor = FALSE), error = function(e) return(c(D = NA, I = NA)))
    D <- vals$D
    I <- vals$I
    c(D = D, I = I)
  }
  
  if (parallel) {
    null_vals <- foreach(i = 1:reps, .combine = rbind, .packages = "ecospat") %dopar% permute_fun(i)
  } else {
    null_vals <- foreach(i = 1:reps, .combine = rbind, .packages = "ecospat") %do% permute_fun(i)
  }
  
  D_null <- null_vals[, "D"]
  I_null <- null_vals[, "I"]
  
  p_D <- mean(D_null <= obs_D, na.rm = TRUE)
  p_I <- mean(I_null <= obs_I, na.rm = TRUE)
  
  end_time <- Sys.time()  # End timer
  duration <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("Observed Schoener's D:", round(obs_D, 3), "\n")
    cat("Observed Warren's I:", round(obs_I, 3), "\n")
    cat("P(D <= observed):", round(p_D, 4), "\n")
    cat("P(I <= observed):", round(p_I, 4), "\n")
    cat("Elapsed time:", round(duration, 2), "seconds\n")
  }
  
  return(tibble::tibble(
    D_obs = obs_D,
    I_obs = obs_I,
    D_null_mean = mean(D_null, na.rm = TRUE),
    I_null_mean = mean(I_null, na.rm = TRUE),
    p_D = p_D,
    p_I = p_I,
    run_time_sec = as.numeric(duration)
  ))
}

# This is a modification of Ecospat's ecospat.plot.niche.dyn() function
# Their function is focused on invasion ecology and uses terms and colouring that is not relevant to comparing two species niches
# My function similfies and makes the plotting more relevant to my research questions
# Comparing niches of two species with questionanle taxonomy and shared backgrounds


ecospat.plot.niche.pair <- function(z1, z2,
                                    use.zcor = FALSE,
                                    quant = 0.1,
                                    quant_occ = 0.1,
                                    quant_bg = 0.1,
                                    drawOccLayer = TRUE,
                                    show.legendOcc = TRUE,
                                    legend.pos = "topleft",
                                    legend.xpd = NA,
                                    legend.inset = 0.02,
                                    legend.cex = 1,
                                    legend.bty = "o",
                                    legend.bg = "white",
                                    legend.box.lwd = 1,
                                    legend.box.col = "black",
                                    drawKD = FALSE,
                                    drawBGextent = TRUE,
                                    drawOccOutline = TRUE,
                                    col.unique_i = "magenta",
                                    col.unique_j = "#E88E00",
                                    col.shared = "#009E73",
                                    col.bg_i = "magenta",
                                    col.bg_j = "#E88E00",
                                    name.axis1 = "PC 1",
                                    name.axis2 = "PC 2",
                                    title = "",
                                    cex.main = 1.2,
                                    cex.lab = 1,
                                    cex.axis = 1,
                                    max.alpha = 0.5) {
  
  if (drawKD) drawOccLayer <- FALSE
  
  # Shared axis and extent
  common_x <- z1$x
  common_y <- z1$y
  ext_obj <- terra::ext(min(common_x), max(common_x),
                        min(common_y), max(common_y))
  
  plot(1, type = "n",
       xlim = range(common_x), ylim = range(common_y),
       xlab = name.axis1, ylab = name.axis2,
       main = title, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
  
  # Choose z layer
  zname <- if (use.zcor) "z.cor" else "z.uncor"
  z1_vals <- terra::values(z1[[zname]])
  z2_vals <- terra::values(z2[[zname]])
  Z1_vals <- terra::values(z1$Z)
  Z2_vals <- terra::values(z2$Z)
  
  # Threshold low background values before computing preference
  Z1_thresh_mask <- Z1_vals > quantile(Z1_vals, probs = quant_bg, na.rm = TRUE)
  Z2_thresh_mask <- Z2_vals > quantile(Z2_vals, probs = quant_bg, na.rm = TRUE)
  
  # Apply mask and compute preference
  if (use.zcor) {
    norm1 <- ifelse(Z1_thresh_mask, z1_vals / Z1_vals, NA)
    norm2 <- ifelse(Z2_thresh_mask, z2_vals / Z2_vals, NA)
  } else {
    norm1 <- z1_vals
    norm2 <- z2_vals
  }
  norm1_mat <- matrix(norm1, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
  norm2_mat <- matrix(norm2, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
  
  # Binary occupancy
  threshold1 <- quantile(norm1, probs = quant, na.rm = TRUE)
  threshold2 <- quantile(norm2, probs = quant, na.rm = TRUE)
  occ1 <- norm1_mat > threshold1
  occ2 <- norm2_mat > threshold2
  
  occ_cat <- matrix(NA, nrow = length(common_y), ncol = length(common_x))
  occ_cat[occ1 & !occ2] <- 1
  occ_cat[occ2 & !occ1] <- 2
  occ_cat[occ1 & occ2]  <- 3
  
  # Draw occupancy raster BEFORE other overlays
  if (drawOccLayer) {
    r_occ_cat <- terra::rast(occ_cat, extent = ext_obj)
    terra::plot(r_occ_cat, add = TRUE,
                col = c(col.unique_i, col.unique_j, col.shared),
                legend = FALSE)
  }
  
  # Kernel shading
  if (drawKD) {
    z1_mat <- matrix(z1_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    z2_mat <- matrix(z2_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    
    plot_density_layer <- function(density_mat, col_base) {
      q_low <- quantile(density_mat, probs = 0.01, na.rm = TRUE)
      q_high <- quantile(density_mat, probs = 0.99, na.rm = TRUE)
      scaled_vals <- pmin((density_mat - q_low) / (q_high - q_low), 1)
      scaled_vals[scaled_vals < 0] <- 0
      scaled_vals[is.na(scaled_vals)] <- 0
      rgb_base <- grDevices::col2rgb(col_base) / 255
      rgba_vec <- rgb(rgb_base[1], rgb_base[2], rgb_base[3], alpha = scaled_vals * max.alpha)
      color_mat <- matrix(rgba_vec, nrow = nrow(density_mat), ncol = ncol(density_mat))
      rasterImage(as.raster(color_mat),
                  xleft = min(common_x), xright = max(common_x),
                  ybottom = min(common_y), ytop = max(common_y))
    }
    
    plot_density_layer(z1_mat, col.unique_i)
    plot_density_layer(z2_mat, col.unique_j)
  }
  
  # Background extent outlines (Z)
  if (drawBGextent) {
    Z1_mat <- matrix(Z1_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    Z2_mat <- matrix(Z2_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    Z1_thresh <- quantile(Z1_vals, probs = quant_bg, na.rm = TRUE)
    Z2_thresh <- quantile(Z2_vals, probs = quant_bg, na.rm = TRUE)
    Z1_mask <- ifelse(Z1_mat > Z1_thresh, 1, NA)
    Z2_mask <- ifelse(Z2_mat > Z2_thresh, 1, NA)
    
    rZ1_outline <- terra::as.polygons(terra::rast(Z1_mask, extent = ext_obj), dissolve = TRUE)
    rZ2_outline <- terra::as.polygons(terra::rast(Z2_mask, extent = ext_obj), dissolve = TRUE)
    
    terra::plot(rZ1_outline, add = TRUE, border = col.bg_i, lty = 3, lwd = 5)
    terra::plot(rZ2_outline, add = TRUE, border = col.bg_j, lty = 3, lwd = 5)
  }
  
  # Occurrence extent outlines (z)
  if (drawOccOutline) {
    z1_mat <- matrix(z1_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    z2_mat <- matrix(z2_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    
    occ_thresh1 <- quantile(z1_vals, probs = quant_occ, na.rm = TRUE)
    occ_thresh2 <- quantile(z2_vals, probs = quant_occ, na.rm = TRUE)
    
    occ_extent1 <- ifelse(z1_mat > occ_thresh1, 1, NA)
    occ_extent2 <- ifelse(z2_mat > occ_thresh2, 1, NA)
    
    rZ1_occ_outline <- terra::as.polygons(terra::rast(occ_extent1, extent = ext_obj), dissolve = TRUE)
    rZ2_occ_outline <- terra::as.polygons(terra::rast(occ_extent2, extent = ext_obj), dissolve = TRUE)
    
    terra::plot(rZ1_occ_outline, add = TRUE, border = col.bg_i, lty = 1, lwd = 1)
    terra::plot(rZ2_occ_outline, add = TRUE, border = col.bg_j, lty = 1, lwd = 1)
  }
  
  if (show.legendOcc && drawOccLayer) {
    legend(legend.pos,
           legend = c("Unique (i)", "Unique (j)", "Shared"),
           fill = c(col.unique_i, col.unique_j, col.shared),
           title = "Occurrence",
           bty = legend.bty,
           bg = legend.bg,
           box.lwd = legend.box.lwd,
           box.col = legend.box.col,
           xpd = legend.xpd,
           inset = legend.inset,
           cex = legend.cex)
  }
  
  invisible(NULL)
}


# Gap analysis functions --------------------------------------------------

#### Gap Analysis functions ####
calculate_srsin <- function(pred, pa, occ, threshold) {
  r_masked <- pred > threshold
  names(r_masked) <- "binary"
  r_pa <- terra::mask(pred, pa) > threshold
  names(r_pa) <- "binary"
  
  occ_pred <- terra::extract(r_masked, occ) %>% as.data.frame() %>% dplyr::filter(binary == 1)
  occ_pa   <- terra::extract(r_pa, occ)     %>% as.data.frame() %>% dplyr::filter(binary == 1)
  
  if (nrow(occ_pred) == 0) return(0)
  return(nrow(occ_pa) / nrow(occ_pred) * 100)
}

calculate_grsin <- function(pred, threshold, hist_masked, pa_mask) {
  r_pred <- pred > threshold
  r_masked <- r_pred & (hist_masked == 1)
  area_all <- terra::expanse(r_masked, byValue = TRUE, unit = "km") %>%
    as.data.frame() %>% filter(value == 1) %>% pull(area)
  
  r_pa_masked <- r_pred & (hist_masked == 1) & (pa_mask == 1)
  area_pa <- terra::expanse(r_pa_masked, byValue = TRUE, unit = "km") %>%
    as.data.frame() %>% filter(value == 1) %>% pull(area)
  
  if (length(area_all) == 0 || sum(area_all) == 0) return(0)
  return(sum(area_pa, na.rm = TRUE) / sum(area_all, na.rm = TRUE) * 100)
}

calculate_ersin <- function(pred, threshold, eco_vec, hist_masked, pa_mask_resampled) {
  vals <- as.data.frame(terra::values(c(pred, hist_masked, pa_mask_resampled, eco_vec)))
  colnames(vals) <- c("pred", "hist", "pa", "eco")
  
  # Remove rows with NA or 0 eco values (assumes 0 = outside eco region)
  vals <- vals[!is.na(vals$eco) & vals$eco != 0, ]
  
  # Historical suitable ecoregions (regardless of PA)
  hist_filter <- vals$hist == 1
  eco_hist_codes <- sort(unique(vals$eco[hist_filter]))
  
  # Future suitable + protected *within historical range*
  future_filter <- vals$pred > threshold & vals$hist == 1 & vals$pa == 1
  eco_future_codes <- sort(unique(vals$eco[future_filter]))
  
  # Coerce to numeric to avoid factor/character mismatch issues
  eco_hist_codes <- as.numeric(as.character(eco_hist_codes))
  eco_future_codes <- as.numeric(as.character(eco_future_codes))
  
  # Filter out invalid or NA codes
  eco_hist_codes <- eco_hist_codes[!is.na(eco_hist_codes)]
  eco_future_codes <- eco_future_codes[!is.na(eco_future_codes)]
  
  # Compute overlap
  overlap_codes <- intersect(eco_hist_codes, eco_future_codes)
  overlap_codes <- overlap_codes[!is.na(overlap_codes)]
  
  # Debug output
  message("[ERSin] Unique ecoregions (historical suitable): ", paste(eco_hist_codes, collapse = " "))
  message("[ERSin] Unique ecoregions (future suitable *in hist range* + protected): ", paste(eco_future_codes, collapse = " "))
  message("[ERSin] Overlap: ", paste(overlap_codes, collapse = " "))
  
  # Return 0 if no overlap
  if (length(overlap_codes) == 0) {
    message("[ERSin] No overlapping ecoregions — returning 0.")
    return(0)
  }
  
  # Final ERSin: proportion of historical suitable ecoregions still suitable+protected in future
  return((length(overlap_codes) / length(eco_hist_codes)) * 100)
}

