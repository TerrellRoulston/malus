# Top ---------------------------------------------------------------------
# Performing PCA of ecological niches for each Malus species
# Terrell Roulston
# Started March 2nd 2025

library(tidyverse)
library(terra)
library(geodata)
library(ade4)
library(ecospat)

#source("scripts/malus_bg.R")
# Load Occurrence data ----------------------------------------------------
occThin_cor <- readRDS(file = './occ_data/cor/occThin_cor.Rdata')
occThin_ion <- readRDS(file = './occ_data/ion/occThin_ion.Rdata')
occThin_ang <- readRDS(file = './occ_data/ang/occThin_ang.Rdata')
occThin_chl <- readRDS(file = './occ_data/chl/occThin_chl.Rdata')

# Load Cropped WClim Data -------------------------------------------------
# These WorldClim rasters are cropped to the background extend used for the SDMs
# See script malus_bg for more details
wclim_cor <- readRDS(file = './wclim_data/wclim_cor.Rdata')
wclim_ion <- readRDS(file = './wclim_data/wclim_ion.Rdata')
wclim_ang <- readRDS(file = './wclim_data/wclim_ang.Rdata')
wclim_chl <- readRDS(file = './wclim_data/wclim_chl.Rdata')


## # Subset the variables included in modeling
## # See WorldClim 2.1 documents for more details on vars
wclim_cor_subs <- wclim_cor %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_ion_subs <- wclim_ion %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_ang_subs <- wclim_ang %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_chl_subs <- wclim_chl %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))


# Convert BG to Matrix for Analysis ---------------------------------------
# Get the full background area by combining Sect. Chloromeles BG together (already did in SDM workflow)
# This is the BG that is used for PCA - models the entire niche space of ALL Sect. Chloromeles species

# Convert BG Rasters to Matrices for PCA ----------------------------------
# Note: Sect. Chloromeles is the full bg of the eastern species
bg_mat_full <- values(wclim_chl_subs) %>% na.omit() # remove NA Values
bg_mat_cor <- values(wclim_cor_subs) %>% na.omit()
bg_mat_ion <- values(wclim_ion_subs) %>% na.omit()
bg_mat_ang <- values(wclim_ang_subs) %>% na.omit()

# Extract Climate Vars from Points ----------------------------------------
# Extract the wclim raster values from occurrence points then bind them with the SpatVector points
# Remove NA values = True
## NB: na.rm = T is not doing anything here!

wclim_cor_occ <- cbind(crds(occThin_cor), extract(wclim_cor_subs, occThin_cor, ID = FALSE))
wclim_ion_occ <- cbind(crds(occThin_ion), extract(wclim_ion_subs, occThin_ion, ID = FALSE))
wclim_ang_occ <- cbind(crds(occThin_ang), extract(wclim_ang_subs, occThin_ang, ID = FALSE))

wclim_cor_occ <- wclim_cor_occ[complete.cases(wclim_cor_occ), ]
wclim_ion_occ <- wclim_ion_occ[complete.cases(wclim_ion_occ), ]
wclim_ang_occ <- wclim_ang_occ[complete.cases(wclim_ang_occ), ]

# PCA ---------------------------------------------------------------------
# Now create a PCA using the FULL BG EXTENT matrix
# We want to generate components of the entire enviromental variability
# Make sure to center (subtract from mean) and scale (divide by SD to make it 0-1)
# Set scannf to false to skip selecting the number of axes and set it equal to 2 with nf

pca_full <- dudi.pca(bg_mat_full, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
pca_score <- pca_full$li # extract PCA scores

# Extract Eigenvalues and compute statistics
eig <- pca_full$eig
prop_var <- eig / sum(eig) * 100    # Proportion of variance in percentage
cum_var <- cumsum(prop_var)         # Cumulative variance

# Build a summary table for principal components
pca_summary <- data.frame(
  PC = paste0("PC", 1:length(eig)),
  Eigenvalue = round(eig, 2),
  Proportion = round(prop_var, 2),
  Cumulative = round(cum_var, 2)
)

# Display the PCA summary table
print(pca_summary)

# Extract the variable loadings for the retained PCs 
loadings <- pca_full$c1
print(loadings)

# Occurrence point PCA scores
# Not sure why but some are producing NA values
# The occurrences should match the BG fine and NA values should have been omited above
cor_occ_score <- suprow(pca_full, data.frame(wclim_cor_occ)[, colnames(bg_mat_full)])$li
ion_occ_score <- suprow(pca_full, data.frame(wclim_ion_occ)[, colnames(bg_mat_full)])$li
ang_occ_score <- suprow(pca_full, data.frame(wclim_ang_occ)[, colnames(bg_mat_full)])$li

# BG Point PCA scores
cor_bg_score <- suprow(pca_full, bg_mat_cor)$li
ion_bg_score <- suprow(pca_full, bg_mat_ion)$li
ang_bg_score <- suprow(pca_full, bg_mat_ang)$li


# Plotting Pairwise Matrix ------------------------------------------------
# Use PCA ordination scores to map the environmental BG and ecological niche in a grid Z of dimension RxR pixels
# Max 1-2 principal components input
# Remove marginal PCA scores in the 95% percentile with th.env = 0.05 (default = 0)
# Can be useful to compare if the shape of the 95 percentile matches the remaining portion of grid
# Create list of <ecospat> grids Z space (kernel density).
# List of SpatRaster Objects <terra>

grids <- list(
  cor = ecospat.grid.clim.dyn(glob = pca_score, glob1 = cor_bg_score, sp = cor_occ_score, R = 100),
  ion = ecospat.grid.clim.dyn(glob = pca_score, glob1 = ion_bg_score, sp = ion_occ_score, R = 100),
  ang = ecospat.grid.clim.dyn(glob = pca_score, glob1 = ang_bg_score, sp = ang_occ_score, R = 100))


## TYLERS FUNCTION
## 
# dev.new(height = 10, width = 7)
# par(mfrow = c(3, 2))
# for(i in c(1, 4, 10, 11, 15, 16)){
#   VAR <- paste("wc2.1_2.5m_bio", i, sep = "_")
#   corbio1 <- ecospat.grid.clim.dyn(glob = bg_mat_full[, VAR],
#                                    glob1 = bg_mat_cor[, VAR],  
#                                    data.frame(wclim_cor_occ)[, VAR],
#                                    R = 100) 
#   ionbio1 <- ecospat.grid.clim.dyn(glob = bg_mat_full[, VAR],
#                                    glob1 = bg_mat_ion[, VAR], 
#                                    data.frame(wclim_ion_occ)[, VAR], R = 100)
#   tws.plot.niche.dyn(corbio1, ionbio1)
#   title(VAR)
# }


# Niche Overlap Metrics and Permutational Signif. -------------------------
# Niche Equivalency Test (PCA-based)
# Permutational test of importance
# Species occurrences are randomized in repeteted runs


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

# Test Species Pairs
# Run on HPC instead

# cor_ion_test <- niche_equivalency_test_broennimann(
#                     sp1_scores = cor_occ_score,
#                     sp2_scores = ion_occ_score,
#                     sp1_bg_scores = cor_bg_score,
#                     sp2_bg_scores = ion_bg_score,
#                     bg_scores = pca_score,
#                     reps = 999,
#                     R = 100,
#                     parallel = TRUE,
#                     ncores = 15,
#                     verbose = TRUE
#                   )
# 
# cor_ang_test <- niche_equivalency_test_broennimann(
#                     sp1_scores = cor_occ_score,
#                     sp2_scores = ion_occ_score,
#                     sp1_bg_scores = cor_bg_score,
#                     sp2_bg_scores = ion_bg_score,
#                     bg_scores = pca_score,
#                     reps = 999,
#                     R = 100,
#                     parallel = TRUE,
#                     ncores = 15,
#                     verbose = TRUE
#                   )
# 
# ion_ang_test <- niche_equivalency_test_broennimann(
#                   sp1_scores = cor_occ_score,
#                   sp2_scores = ion_occ_score,
#                   sp1_bg_scores = cor_bg_score,
#                   sp2_bg_scores = ion_bg_score,
#                   bg_scores = pca_score,
#                   reps = 999,
#                   R = 100,
#                   parallel = TRUE,
#                   ncores = 15,
#                   verbose = TRUE
#                 )


# Load results from HPC
files <- list.files("./", pattern = "niche_equiv_.*\\.rds", full.names = TRUE)

results <- map_dfr(files, readRDS, .id = "pair")
results$pair <- gsub("niche_equiv_|\\.rds", "", basename(results$pair))

print(results)


# Modify Ecospat Function -------------------------------------------------
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




# PCA Overlap Plot w/ Custom Func. ----------------------------------------
# Assuming 'grids' is a list of grid objects generated via ecospat.grid.clim.dyn.
# NOTE: Change variables, save file name and legend colours each time before plotting

# There is a bug when plotting the Chloromeles where it takes the aestic of layer i instead of j

#Legend colour
# Coronaria = #F763C2
# Ioensis = #E88E00
# Angustifolia = #007CBE
# Shared niche = #009E73

# Open a new high-resolution graphics device
jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/ion_ang_pair.jpeg',
     width = 3333, height = 3333, res = 300)

par(mar = c(7, 8, 4, 2), mgp = c(5, 2, 0))

ecospat.plot.niche.pair(
  z1 = grids[["ion"]], 
  z2 = grids[["ang"]],
  use.zcor = F,        # Toggle between using z.cor (TRUE) and z.uncor (FALSE). Z.cor = z/Z = niche preference, z.uncor = z = realized niche
  quant = 0.1,         # occupancy threshold
  quant_occ = 0.5,     # niche (z.cor) outline threshold
  quant_bg = 0.01,     # background (Z) extent threshold
  drawKD = F,          # Draw kernel density shading
  drawOccOutline = F,
  drawBGextent = T,
  col.unique_i = "#E88E00",  # Color for niche unique to species i (z1)
  col.unique_j = "#007CBE",# Color for niche unique to species j (z2)
  col.bg_i = "#E88E00",      # Color for background outline from z1's Z component
  col.bg_j = "#007CBE",    # Color for background outline from z2's Z component
  name.axis1 = "PC 1",    # X-axis label
  name.axis2 = "PC 2",    # Y-axis label
  title = NULL,  # Plot title
  show.legendOcc = F,
  max.alpha = 0.6,
  cex.lab = 3.5,
  cex.axis = 3.5,
)


  ####
# Compute overlap metrics for displaying
# Note excluded D metric, only reporting Warren's I in the MS
### 
metrics <- ecospat.niche.overlap(grids[["ion"]], grids[["ang"]], cor = F)
legend("topright", legend = paste0(
                                    #"D = ", round(metrics$D, 2), "***",
                                   "I = ", round(metrics$I, 2), "***"),
       bty = "n", cex = 3.5, inset = -0.01)


### SAVE
dev.off() 



# Save Legend for Plot
legend_labs <- c(expression(italic('Malus coronaria')), expression(italic('Malus ioensis')), expression(italic('Malus angustifolia')), 'Shared Niche')
fill_cols <- c("magenta", "#E88E00", "#007CBE", "#009E73")

jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/pca_legend.jpeg", width = 3333, height = 3333, res = 300)
par(mar = c(0,0,0,0))
# Plot a legend that can be saved on its own
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#legend('center', xpd = NA, title = c(as.expression(bquote(bold('Habitat Suitability')))), legend = legend_labs, fill = fill_cols, cex = 3)
legend('center', xpd = NA, box.lwd = 2, legend = legend_labs, fill = fill_cols, cex = 5, horiz = F, bty = "o", title = expression(underline("Realized Niche")))
dev.off()

# Biplot with each species ------------------------------------------------
occ_scores_list <- list(
  "Malus coronaria"      = cor_occ_score,
  "Malus ioensis"        = ion_occ_score,
  "Malus angustifolia"   = ang_occ_score,
  "Sect. Chloromeles"    = chl_occ_score
)


# Compute global axis limits using "Axis1" and "Axis2" across all species
all_scores <- do.call(rbind, occ_scores_list)
x_min <- min(all_scores$Axis1, na.rm = TRUE)
x_max <- max(all_scores$Axis1, na.rm = TRUE)
y_min <- min(all_scores$Axis2, na.rm = TRUE)
y_max <- max(all_scores$Axis2, na.rm = TRUE)

# Axis labels
xlab_text <- paste0("PC 1 (", signif(pca_summary$Proportion[1], 3), "%)")
ylab_text <- paste0("PC 2 (", signif(pca_summary$Proportion[2], 3), "%)")

# Set up the plotting area (adjust mfrow as needed)

# Scale factor for loadings (optional tweak to improve visibility)
scale_factor <- 1.75

jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/bi_plot.jpeg',
     width = 6666, height = 6666, res = 300)

par(mar = c(7, 8, 6, 2), mgp = c(5, 2, 0))
# Plot occurrence points using Axis1 and Axis2
# Plot M. coronaria
plot(cor_occ_score$Axis1, cor_occ_score$Axis2,
       xlim = c(x_min, x_max), ylim = c(y_min, y_max),
       xlab = xlab_text, ylab = ylab_text,
       main = NULL,
       pch = 19, col = adjustcolor("magenta", alpha.f = 0.6),
       cex.axis = 3.5, cex.lab = 3.5, cex.main = 5, cex = 2.5)


# Plot M. ioensis
points(ion_occ_score$Axis1, ion_occ_score$Axis2,
     pch = 19, col = adjustcolor("#E88E00", alpha.f = 0.6), cex = 2.5)

# Plot M. angustifolia
points(ang_occ_score$Axis1, ang_occ_score$Axis2,
     pch = 19, col = adjustcolor("#007CBE", alpha.f = 0.6), cex = 2.5)

# Draw reference lines at 0
abline(h = 0, v = 0, col = "gray60", lty = 2, lwd = 2)
  
# Overlay PCA loadings (from your global PCA stored in pca_full$c1)
# Multiply by scale_factor to enhance visibility if needed
arrows(0, 0,
         pca_full$c1[,1] * scale_factor,
         pca_full$c1[,2] * scale_factor,
         length = 0.1, col = "red", lwd = 5)
  
# Create simplified labels from the original rownames
simplified_labels <- sub("wc2.1_2.5m_bio_", "Bio ", rownames(pca_full$c1))
  
# Add text labels for each variable loading
text(pca_full$c1[,1] * scale_factor,
       pca_full$c1[,2] * scale_factor,
       labels = simplified_labels,
       col = "black", pos = 3, cex = 3.5, font = 2)

# Add legend
# legend("topright",
#   legend = c(expression(italic("M. coronaria")), 
#              expression(italic("M. ioensis")),
#              expression(italic("M. angustifolia")),
#              "Sect. Chloromeles"),
#   fill = adjustcolor(c("#6A3D9A", "#E88E00", "#007CBE", "#009E73"), alpha.f = 0.6),
#   box.col = "black",  
#   bg = "white",
#   text.col = 'black',
#   bty = "n",
#   cex = 1.4
# )

dev.off() # Save