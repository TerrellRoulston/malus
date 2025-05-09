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
bg_mat_chl <- values(wclim_chl_subs) %>% na.omit()

# Extract Climate Vars from Points ----------------------------------------
# Extract the wclim raster values from occurrence points then bind them with the SpatVector points
# Remove NA values = True
## NB: na.rm = T is not doing anything here!

wclim_cor_occ <- cbind(crds(occThin_cor), extract(wclim_cor_subs, occThin_cor, ID = FALSE))
wclim_ion_occ <- cbind(crds(occThin_ion), extract(wclim_ion_subs, occThin_ion, ID = FALSE))
wclim_ang_occ <- cbind(crds(occThin_ang), extract(wclim_ang_subs, occThin_ang, ID = FALSE))
wclim_chl_occ <- cbind(crds(occThin_chl), extract(wclim_chl_subs, occThin_chl, ID = FALSE))

wclim_cor_occ <- wclim_cor_occ[complete.cases(wclim_cor_occ), ]
wclim_ion_occ <- wclim_ion_occ[complete.cases(wclim_ion_occ), ]
wclim_ang_occ <- wclim_ang_occ[complete.cases(wclim_ang_occ), ]
wclim_chl_occ <- wclim_chl_occ[complete.cases(wclim_chl_occ), ]

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
chl_occ_score <- suprow(pca_full, data.frame(wclim_chl_occ)[, colnames(bg_mat_full)])$li

# BG Point PCA scores
cor_bg_score <- suprow(pca_full, bg_mat_cor)$li
ion_bg_score <- suprow(pca_full, bg_mat_ion)$li
ang_bg_score <- suprow(pca_full, bg_mat_ang)$li
chl_bg_score <- suprow(pca_full, bg_mat_chl)$li


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
  ang = ecospat.grid.clim.dyn(glob = pca_score, glob1 = ang_bg_score, sp = ang_occ_score, R = 100),
  chl = ecospat.grid.clim.dyn(glob = pca_score, glob1 = chl_bg_score, sp = chl_occ_score, R = 100)
)


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


# Return A list with observed D/I, null D/I vectors, and p-values

niche_similarity_shift_test <- function(z1, z2,
                                        quant = 0.05,
                                        reps = 999,
                                        parallel = TRUE,
                                        ncores = 2,
                                        verbose = TRUE) {
  require(terra)
  require(foreach)
  require(tibble)
  require(dplyr)
  if (parallel) {
    require(doParallel)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  # Extract rasters
  r1 <- z1$z.cor / z1$Z
  r2 <- z2$z.cor / z2$Z
  r1[!is.finite(r1)] <- 0
  r2[!is.finite(r2)] <- 0
  
  # Binary occupancy
  thresh1 <- quantile(values(r1), probs = quant, na.rm = TRUE)
  thresh2 <- quantile(values(r2), probs = quant, na.rm = TRUE)
  occ1 <- r1 > thresh1
  occ2 <- r2 > thresh2
  
  # Matrix conversions
  mat1 <- as.matrix(occ1)
  mat2 <- as.matrix(occ2)
  bg_mask <- as.matrix((z1$Z > 0) & (z2$Z > 0))
  
  get_DI <- function(m1, m2, mask) {
    idx <- which(mask & (m1 | m2))
    v1 <- m1[idx] / sum(m1[idx])
    v2 <- m2[idx] / sum(m2[idx])
    D <- 1 - sum(abs(v1 - v2)) / 2
    I <- 1 - sum((sqrt(v1) - sqrt(v2))^2) / 2
    c(D = D, I = I)
  }
  
  obs_metrics <- get_DI(mat1, mat2, bg_mask)
  
  shift_fun <- function(seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    dx <- sample(-(ncol(mat1) - 1):(ncol(mat1) - 1), 1)
    dy <- sample(-(nrow(mat1) - 1):(nrow(mat1) - 1), 1)
    roll_matrix <- function(mat, dx, dy) {
      out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
      x1 <- max(1, 1 + dx); x2 <- min(ncol(mat), ncol(mat) + dx)
      y1 <- max(1, 1 + dy); y2 <- min(nrow(mat), nrow(mat) + dy)
      i1 <- max(1, 1 - dx); i2 <- min(ncol(mat), ncol(mat) - dx)
      j1 <- max(1, 1 - dy); j2 <- min(nrow(mat), nrow(mat) - dy)
      out[j1:j2, i1:i2] <- mat[y1:y2, x1:x2]
      return(out)
    }
    roll1 <- roll_matrix(mat1, dx, dy)
    roll2 <- roll_matrix(mat2, -dx, -dy)
    tryCatch(get_DI(roll1, roll2, bg_mask), error = function(e) c(D = NA, I = NA))
  }
  
  null_vals <- foreach(i = 1:reps, .combine = rbind, .packages = "terra") %dopar% {
    shift_fun()
  }
  
  D_null <- null_vals[, "D"]
  I_null <- null_vals[, "I"]
  
  p_D <- mean(D_null >= obs_metrics["D"], na.rm = TRUE)
  p_I <- mean(I_null >= obs_metrics["I"], na.rm = TRUE)
  
  if (verbose) {
    cat("Observed Schoener's D:", round(obs_metrics["D"], 3), "\n")
    cat("Observed Warren's I:", round(obs_metrics["I"], 3), "\n")
    cat("P(D >= observed):", round(p_D, 4), "\n")
    cat("P(I >= observed):", round(p_I, 4), "\n")
  }
  
  return(tibble(
    metric = c("Schoener's D", "Warren's I"),
    observed = c(obs_metrics["D"], obs_metrics["I"]),
    p_value = c(p_D, p_I),
    null_mean = c(mean(D_null, na.rm = TRUE), mean(I_null, na.rm = TRUE)),
    null_sd   = c(sd(D_null, na.rm = TRUE), sd(I_null, na.rm = TRUE)),
    null_min  = c(min(D_null, na.rm = TRUE), min(I_null, na.rm = TRUE)),
    null_max  = c(max(D_null, na.rm = TRUE), max(I_null, na.rm = TRUE))
  ))
}


niche_similarity_shift_test(grids[["cor"]], grids[["ion"]],
                            quant = 0.1,
                            reps = 1000,
                            parallel = TRUE,
                            ncores = 4,
                            verbose = TRUE)



species_pairs <- list(
  c("cor", "ion"),
  c("cor", "ang"),
  c("ion", "ang")
)

results <- purrr::map_dfr(species_pairs, ~ run_similarity_test(
  grids[[.x[1]]], grids[[.x[2]]],
  name_i = .x[1], name_j = .x[2],
  reps = 1000, cores = 4
))


# Modify Ecospat Function -------------------------------------------------
# Work i PROGRESS. Still not producing what I want 
# Under the hood the dyn. fuction plots the kernel density of the occurrence points of species j
# To prevent this I am going to define a new function, and only select certain arguments from
# the ecospat.plot.niche.dyn function

# col.stab = "#9467bd",   # shared niche region (stability)
# col.unf = "#1f77b4",    # niche exclusive to species 1 (unfilling)
# col.exp = "#ff7f0e",

# new new new new ---------------------------------------------------------

ecospat.plot.niche.pair <- function(z1, z2,
                                    use.zcor = TRUE,         # Logical, T use the z.cor value, F use z.uncor value
                                    quant = 0.1,             # Set what percentile to threshold occupancy (z and Z)
                                    drawOccLayer = TRUE,     # Draw categorical occurrence layer
                                    show.legendOcc = TRUE,   # Toggle legend for occurrence layer
                                    legend.pos = "topleft",  # Legend position
                                    legend.xpd = NA,         # Allow legend outside plot region (NA = anywhere on device)
                                    legend.inset = 0.02,     # Legend inset (can be negative to push outside)
                                    legend.cex = 1,          # Legend text scaling
                                    legend.bty = "o",        # Legend border type ("o" for box, "n" for none)
                                    legend.bg = "white",     # Legend background fill
                                    legend.box.lwd = 1,      # Legend box line width
                                    legend.box.col = "black",# Legend box border color
                                    drawKD = TRUE,           # Add occupancy Kernel density 
                                    transparency = 50,       # Set base transparency
                                    col.shared = "forestgreen", # Colour shared niche
                                    col.unique_i = "blue",   # Colour unique i niche
                                    col.unique_j = "orange", # Colour unique j niche
                                    col.bg_i = "blue",       # Colour background boundary i
                                    col.bg_j = "orange",     # Colour background boundary j
                                    name.axis1 = "PC 1",     # X axis title
                                    name.axis2 = "PC 2",     # Y axis title
                                    title = "",              # Main title
                                    cex.main = 1.2,          # Main title font size
                                    cex.lab = 1,             # Axis label font size
                                    cex.axis = 1,            # Axis tick font size
                                    compute.metrics = TRUE,  # Toggle metric calculation and display
                                    cex.metrics = 1.2,       # Text size for printed metrics
                                    metrics.pos = "topright", # Where to print metrics
                                    metrics.inset = 0.02     # Fine control of inset
) {
  # Get shared grid dimensions
  common_x <- z1$x
  common_y <- z1$y
  ext_obj <- terra::ext(min(common_x), max(common_x),
                        min(common_y), max(common_y))
  
  # Plot setup (axes span the full PCA range, including negatives)
  plot(1, type = "n",
       xlim = range(common_x), ylim = range(common_y),
       xlab = name.axis1, ylab = name.axis2,
       main = title, 
       cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
  
  # Choose z layer
  zname <- if (use.zcor) "z.cor" else "z.uncor"
  
  # Extract values from SpatRasters
  z1_vals <- terra::values(z1[[zname]])
  Z1_vals <- terra::values(z1$Z)
  z2_vals <- terra::values(z2[[zname]])
  Z2_vals <- terra::values(z2$Z)
  
  # Normalize by background
  norm1 <- ifelse(Z1_vals != 0, z1_vals / Z1_vals, 0)
  norm2 <- ifelse(Z2_vals != 0, z2_vals / Z2_vals, 0)
  
  # Convert normalized values to matrices
  norm1_mat <- matrix(norm1, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
  norm2_mat <- matrix(norm2, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
  
  # Compute thresholds
  threshold1 <- quantile(norm1, probs = quant, na.rm = TRUE)
  threshold2 <- quantile(norm2, probs = quant, na.rm = TRUE)
  
  # Binary occupancy
  occ1 <- norm1_mat > threshold1
  occ2 <- norm2_mat > threshold2
  
  # Create categorical occupancy layer
  if (drawOccLayer) {
    occ_cat <- matrix(NA, nrow = length(common_y), ncol = length(common_x))
    occ_cat[occ1 & !occ2] <- 1 # Unique z species i
    occ_cat[occ2 & !occ1] <- 2 # Unique z species j
    occ_cat[occ1 & occ2]  <- 3 # Shared z species i & j
    
    r_occ_cat <- terra::rast(occ_cat, extent = ext_obj)
    
    # Plot the occupancy raster without an automatic legend.
    terra::plot(r_occ_cat, add = TRUE,
                col = c(col.unique_i, col.unique_j, col.shared),
                legend = FALSE)
  }
  
  # Optional density shading (drawn on top of the categorical occupancy layer)
  if (drawKD) {
    z1_mat <- matrix(z1_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    z2_mat <- matrix(z2_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
    
    r_z1 <- terra::rast(z1_mat, extent = ext_obj)
    r_z2 <- terra::rast(z2_mat, extent = ext_obj)
    
    terra::plot(r_z1, add = TRUE,
                col = colorRampPalette(c("white", col.unique_i))(transparency + 10),
                legend = FALSE)
    terra::plot(r_z2, add = TRUE,
                col = colorRampPalette(c("white", col.unique_j))(transparency + 10),
                legend = FALSE)
  }
  
  # Draw background extent outlines
  Z1_mat <- matrix(Z1_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
  Z2_mat <- matrix(Z2_vals, nrow = length(common_y), ncol = length(common_x), byrow = TRUE)
  
  Z1_thresh <- quantile(Z1_vals, probs = quant, na.rm = TRUE)
  Z2_thresh <- quantile(Z2_vals, probs = quant, na.rm = TRUE)
  
  Z1_mask <- ifelse(Z1_mat > Z1_thresh, 1, NA)
  Z2_mask <- ifelse(Z2_mat > Z2_thresh, 1, NA)
  
  rZ1_mask <- terra::rast(Z1_mask, extent = ext_obj)
  rZ2_mask <- terra::rast(Z2_mask, extent = ext_obj)
  
  rZ1_outline <- terra::as.polygons(rZ1_mask, dissolve = TRUE)
  rZ2_outline <- terra::as.polygons(rZ2_mask, dissolve = TRUE)
  
  terra::plot(rZ1_outline, add = TRUE, border = col.bg_i, lty = 2, lwd = 1.5)
  terra::plot(rZ2_outline, add = TRUE, border = col.bg_j, lty = 2, lwd = 1.5)
  
  # Add a manual legend for the categorical occupancy layer if requested
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
  
  # Compute and optionally display niche overlap metrics
  if (compute.metrics) {
    # Flatten binary occupancy layers and restrict to overlapping cells
    v1 <- as.vector(occ1)
    v2 <- as.vector(occ2)
    valid_idx <- which(!is.na(v1) & !is.na(v2) & (v1 | v2))
    v1 <- v1[valid_idx]
    v2 <- v2[valid_idx]
    
    # Normalize binary occupancy to sum to 1
    v1 <- v1 / sum(v1)
    v2 <- v2 / sum(v2)
    
    # Overlap metrics
    D_val <- 1 - sum(abs(v1 - v2)) / 2
    I_val <- 1 - sum((sqrt(v1) - sqrt(v2))^2) / 2
    
    # Display metrics on the plot
    metrics_text <- paste0("D = ", round(D_val, 2),
                           "\nI = ", round(I_val, 2))
    legend(metrics.pos, legend = metrics_text,
           bty = "n", cex = cex.metrics, inset = metrics.inset,
           text.col = "black")
    
    return(invisible(list(D = D_val, I = I_val)))
  } else {
    return(invisible(NULL))
  }
}



# PCA Overlap Plot w/ Custom Func. ----------------------------------------
# Assuming 'grids' is a list of grid objects generated via ecospat.grid.clim.dyn.
# NOTE: Change variables, save file name and legend colours each time before plotting

# There is a bug when plotting the Chloromeles where it takes the aestic of layer i instead of j

#Legend colour
# Coronaria = magenta
# Ioensis = #E88E00
# Angustifolia = #007CBE
# Shared niche = #009E73

# Open a new high-resolution graphics device
jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/cor_ion_pair.jpeg',
     width = 3333, height = 3333, res = 300)

par(mar = c(7, 8, 4, 2), mgp = c(5, 2, 0))

ecospat.plot.niche.pair(
  z1 = grids[["cor"]], 
  z2 = grids[["ion"]],
  use.zcor = TRUE,        # Toggle between using z.cor (TRUE) and z.uncor (FALSE)
  quant = 0.1,           # Quantile threshold for determining occupied cells
  drawKD = F,          # Draw kernel density shading
  transparency = 50,      # Transparency setting for the density color ramp
  col.shared = "#009E73",  # Color for shared niche region
  col.unique_i = "magenta",  # Color for niche unique to species i (z1)
  col.unique_j = "#E88E00",# Color for niche unique to species j (z2)
  col.bg_i = "magenta",      # Color for background outline from z1's Z component
  col.bg_j = "#E88E00",    # Color for background outline from z2's Z component
  name.axis1 = "PC 1",    # X-axis label
  name.axis2 = "PC 2",    # Y-axis label
  title = "Ang. vs. Chl.",  # Plot title
  show.legendOcc = F,
  cex.axis = 2,
  cex.lab = 2,
  cex.main = 2,
  cex.metrics = 2,
  metrics.inset = 0.05
)

dev.off() # Return "Null" to save


# Save Legend for Plot
legend_labs <- c('M. coronaria', 'M. ioensis', 'M. angustifolia', 'Sect. Chloromeles', 'Shared Niche')
fill_cols <- c("#6A3D9A", "#E88E00", "#007CBE", "#333f07", "#009E73")

jpeg(filename = "C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/pca_legend.jpeg", width = 9999, height = 6666, res = 300)
# Plot a legend that can be saved on its own
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#legend('center', xpd = NA, title = c(as.expression(bquote(bold('Habitat Suitability')))), legend = legend_labs, fill = fill_cols, cex = 3)
legend('center', xpd = NA, box.lwd = 2, legend = legend_labs, fill = fill_cols, cex = 9, horiz = F, bty = "o", title = 'Taxons Niche')
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
#par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Scale factor for loadings (optional tweak to improve visibility)
scale_factor <- 1.75

jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/bi_plot.jpeg',
     width = 6666, height = 6666, res = 300)

par(mar = c(7, 8, 4, 2), mgp = c(5, 2, 0))
# Plot occurrence points using Axis1 and Axis2
# Plot M. coronaria
plot(cor_occ_score$Axis1, cor_occ_score$Axis2,
       xlim = c(x_min, x_max), ylim = c(y_min, y_max),
       xlab = xlab_text, ylab = ylab_text,
       main = "Sect. Chloromeles PCA Biplot",
       pch = 19, col = adjustcolor("magenta", alpha.f = 0.6),
       cex.axis = 3, cex.lab = 3, cex.main = 4.5, cex = 2.5)


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
         length = 0.1, col = "red", lwd = 3)
  
# Create simplified labels from the original rownames
simplified_labels <- sub("wc2.1_2.5m_bio_", "Bio ", rownames(pca_full$c1))
  
# Add text labels for each variable loading
text(pca_full$c1[,1] * scale_factor,
       pca_full$c1[,2] * scale_factor,
       labels = simplified_labels,
       col = "black", pos = 3, cex = 3, font = 2)

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


####
# PCA Biplot faceted by species
####

# occ_scores_list <- list(
#   "Malus coronaria"      = cor_occ_score,
#   "Malus ioensis"        = ion_occ_score,
#   "Malus angustifolia"   = ang_occ_score,
#   "Sect. Chloromeles"    = chl_occ_score
# )
# 
# 
# # Compute global axis limits using "Axis1" and "Axis2" across all species
# all_scores <- do.call(rbind, occ_scores_list)
# x_min <- min(all_scores$Axis1, na.rm = TRUE)
# x_max <- max(all_scores$Axis1, na.rm = TRUE)
# y_min <- min(all_scores$Axis2, na.rm = TRUE)
# y_max <- max(all_scores$Axis2, na.rm = TRUE)
# 
# # Axis labels
# xlab_text <- paste0("PC 1 (", signif(pca_summary$Proportion[1], 3), "%)")
# ylab_text <- paste0("PC 2 (", signif(pca_summary$Proportion[2], 3), "%)")
# 
# # Set up the plotting area (adjust mfrow as needed)
# par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# 
# # Scale factor for loadings (optional tweak to improve visibility)
# scale_factor <- 1.75
# 
# # Loop over each species in occ_scores_list
# for (sp in names(occ_scores_list)) {
#   occ_data <- occ_scores_list[[sp]]
#   
#   # Plot occurrence points using Axis1 and Axis2
#   plot(occ_data$Axis1, occ_data$Axis2,
#        xlim = c(x_min, x_max), ylim = c(y_min, y_max),
#        xlab = xlab_text, ylab = ylab_text,
#        main = paste("PCA Biplot for", sp),
#        pch = 16, col = adjustcolor("darkgray", alpha.f = 0.6))
#   
#   # Draw reference lines at 0
#   abline(h = 0, v = 0, col = "gray60", lty = 2)
#   
#   # Overlay PCA loadings (from your global PCA stored in pca_full$c1)
#   # Multiply by scale_factor to enhance visibility if needed
#   arrows(0, 0,
#          pca_full$c1[,1] * scale_factor,
#          pca_full$c1[,2] * scale_factor,
#          length = 0.1, col = "red", lwd = 2)
#   
#   # Create simplified labels from the original rownames
#   simplified_labels <- sub("wc2.1_2.5m_bio_", "Bio ", rownames(pca_full$c1))
#   
#   # Add text labels for each variable loading
#   text(pca_full$c1[,1] * scale_factor,
#        pca_full$c1[,2] * scale_factor,
#        labels = simplified_labels,
#        col = "blue", pos = 3, cex = 1.2)
# }

