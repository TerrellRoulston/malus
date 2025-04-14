# Top ---------------------------------------------------------------------
# Performing PCA of ecological niches for each Malus species
# Terrell Roulston
# Started March 2nd 2025

library(tidyverse)
library(ecospat)
library(terra)
library(geodata)
library(ade4)
library(grid)


# Load Occurrence data ----------------------------------------------------
occThin_cor <- readRDS(file = './occ_data/cor/occThin_cor.Rdata')
occThin_fus <- readRDS(file = './occ_data/fus/occThin_fus.Rdata')
occThin_ion <- readRDS(file = './occ_data/ion/occThin_ion.Rdata')
occThin_ang <- readRDS(file = './occ_data/ang/occThin_ang.Rdata')
occThin_chl <- readRDS(file = './occ_data/chl/occThin_chl.Rdata')

# Load Cropped WClim Data -------------------------------------------------
# These WorldClim rasters are cropped to the background extend used for the SDMs
# See script malus_bg for more details
wclim_cor <- readRDS(file = './wclim_data/wclim_cor.Rdata') 
wclim_fus <- readRDS(file = './wclim_data/wclim_fus.Rdata')
wclim_ion <- readRDS(file = './wclim_data/wclim_ion.Rdata')
wclim_ang <- readRDS(file = './wclim_data/wclim_ang.Rdata')
wclim_chl <- readRDS(file = './wclim_data/wclim_chl.Rdata')


# Subset the variables included in modeling
# See WorldClim 2.1 documents for more details on vars
wclim_cor_subs <- wclim_cor %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_fus_subs <- wclim_fus %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_ion_subs <- wclim_ion %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_ang_subs <- wclim_ang %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))
wclim_chl_subs <- wclim_chl %>% terra::subset(c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16'))


# Convery BG to Matrix for Analysis ---------------------------------------
# Get the full background area by combing the M. fusca and Sect. Chloromeles BG together
# This is the BG that is used for PCA - models the entire niche space of ALL species
wclim_full <- terra::mosaic(wclim_fus_subs, wclim_chl_subs)

# Convert BG Rasters to Matrices for PCA ----------------------------------
bg_mat_full <- values(wclim_full) %>% na.omit() # remove NA Values
bg_mat_cor <- values(wclim_cor_subs) %>% na.omit()
bg_mat_fus <- values(wclim_fus_subs) %>% na.omit()
bg_mat_ion <- values(wclim_ion_subs) %>% na.omit()
bg_mat_ang <- values(wclim_ang_subs) %>% na.omit()
bg_mat_chl <- values(wclim_chl_subs) %>% na.omit()

# Extract Climate Vars from Points ----------------------------------------
# Extract the wclim raster values from occurrence points then bind them with the SpatVector points
# Remove NA values = True
wclim_cor_occ <- cbind(occThin_cor, extract(wclim_cor, occThin_cor, na.rm = T))
wclim_fus_occ <- cbind(occThin_fus, extract(wclim_fus, occThin_fus, na.rm = T))
wclim_ion_occ <- cbind(occThin_ion, extract(wclim_ion, occThin_ion, na.rm = T))
wclim_ang_occ <- cbind(occThin_ang, extract(wclim_ang, occThin_ang, na.rm = T))
wclim_chl_occ <- cbind(occThin_chl, extract(wclim_chl, occThin_chl, na.rm = T))

# PCA ---------------------------------------------------------------------
# Now create a PCA using the FULL BG EXTENT matrix
# We want to generate components of the entire enviromental variability
# Make sure to center (subtract from mean) and scale (divide by SD to make it 0-1)
# Set scannf to false to skip selecting the number of axes and set it equal to 2 with nf

pca_full <- dudi.pca(bg_mat_full, center = TRUE,scale = TRUE, scannf = FALSE, nf = 2)
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
cor_occ_score <- suprow(pca_full, data.frame(wclim_cor_occ)[, colnames(bg_mat_full)])$li %>% na.omit() 
fus_occ_score <- suprow(pca_full, data.frame(wclim_fus_occ)[, colnames(bg_mat_full)])$li %>% na.omit()
ion_occ_score <- suprow(pca_full, data.frame(wclim_ion_occ)[, colnames(bg_mat_full)])$li %>% na.omit()
ang_occ_score <- suprow(pca_full, data.frame(wclim_ang_occ)[, colnames(bg_mat_full)])$li %>% na.omit()
chl_occ_score <- suprow(pca_full, data.frame(wclim_chl_occ)[, colnames(bg_mat_full)])$li %>% na.omit()

# BG Point PCA scores
cor_bg_score <- suprow(pca_full, bg_mat_cor)$li
fus_bg_score <- suprow(pca_full, bg_mat_fus)$li
ion_bg_score <- suprow(pca_full, bg_mat_ion)$li
ang_bg_score <- suprow(pca_full, bg_mat_ang)$li
chl_bg_score <- suprow(pca_full, bg_mat_chl)$li


# Pseudo-code for plotting ------------------------------------------------
# The plotting function for the ecospat plot is a bit limited and I am trying to replicate it from source
# Below is an example of the data structure that I want the plot density to show.

# Assume D1 and D2 are matrices of density values for species 1 and species 2, respectively.
# And T1 and T2 are thresholds determined by quantiles of D1 and D2.

# For each cell (i,j) in the grid:
#   if (D1[i,j] > T1 and D2[i,j] > T2)
#     then assign cell to category "stability" (shared niche)
# else if (D1[i,j] <= T1 and D2[i,j] > T2)
#   then assign cell to category "expansion" (species 2 occupies unique space)
# else if (D1[i,j] > T1 and D2[i,j] <= T2)
#   then assign cell to category "unfilling" (species 1 occupies
# else
#   assign cell as background or zero

# Optionally, define an "abandoned" niche category for species 2.



# Plotting Pairwise Matrix ------------------------------------------------
# Use PCA ordination scores to map the environmental BG and ecological niche in a grid Z of dimension RxR pixels
# Max 1-2 principal components input
# Remove marginal PCA scores in the 95% percentile with th.env = 0.05 (default = 0)
# Can be useful to compare if the shape of the 95 percentile matches the remaining portion of grid
# Create list of <ecospat> grids Z space (kernel density).
# List of SpatRaster Objects <terra>

grids <- list(
  fus = ecospat.grid.clim.dyn(glob = pca_score, glob1 = fus_bg_score, sp = fus_occ_score, R = 100),
  cor = ecospat.grid.clim.dyn(glob = pca_score, glob1 = cor_bg_score, sp = cor_occ_score, R = 100),
  ion = ecospat.grid.clim.dyn(glob = pca_score, glob1 = ion_bg_score, sp = ion_occ_score, R = 100),
  ang = ecospat.grid.clim.dyn(glob = pca_score, glob1 = ang_bg_score, sp = ang_occ_score, R = 100),
  chl = ecospat.grid.clim.dyn(glob = pca_score, glob1 = chl_bg_score, sp = chl_occ_score, R = 100)
)


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
                                    quant = 0.05,            # Set what percentile to threshold occupancy (z and Z)
                                    drawOccLayer = TRUE,     # Draw categorical occurrence layer
                                    show.legendOcc = TRUE,   # Toggle legend for occurrence layer
                                    legend.pos = "topright", # Legend position
                                    legend.xpd = NA,         # Allow legend outside plot region (NA = anywhere on device)
                                    legend.inset = 0.02,     # Legend inset (can be negative to push outside)
                                    legend.cex = 1,          # Legend text scaling
                                    legend.bty = "o",        # Legend border type ("o" for box, "n" for none)
                                    legend.bg = "white",     # Legend background fill
                                    legend.box.lwd = 1,      # Legend box line width
                                    legend.box.col = "black",# Legend box border color
                                    drawKD = TRUE,           # Add occupancy Kernel density 
                                    transparency = 50,       # Set base transparency
                                    col.shared = "purple",   # Colour shared niche
                                    col.unique_i = "blue",   # Colour unique i niche
                                    col.unique_j = "orange", # Colour unique j niche
                                    col.bg_i = "blue",       # Colour background boundary i
                                    col.bg_j = "orange",     # Colour background boundary j
                                    name.axis1 = "PC 1",     # X axis title
                                    name.axis2 = "PC 2",     # Y axis title
                                    title = "")              # Main title
{
  # Get shared grid dimensions
  common_x <- z1$x
  common_y <- z1$y
  ext_obj <- terra::ext(min(common_x), max(common_x),
                        min(common_y), max(common_y))
  
  # Plot setup (axes span the full PCA range, including negatives)
  plot(1, type = "n",
       xlim = range(common_x), ylim = range(common_y),
       xlab = name.axis1, ylab = name.axis2,
       main = title)
  
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
  
  invisible(NULL)
}


# EXAMPLE PLOT
# Assuming 'grids' is a list of grid objects generated via ecospat.grid.clim.dyn.
# Here we test the function with species "fus" and "chl".
ecospat.plot.niche.pair(
  z1 = grids[["fus"]], 
  z2 = grids[["chl"]],
  use.zcor = TRUE,        # Toggle between using z.cor (TRUE) and z.uncor (FALSE)
  quant = 0.05,           # Quantile threshold for determining occupied cells
  drawKD = F,          # Draw kernel density shading
  transparency = 50,      # Transparency setting for the density color ramp
  col.shared = "purple",  # Color for shared niche region
  col.unique_i = "blue",  # Color for niche unique to species i (z1)
  col.unique_j = "orange",# Color for niche unique to species j (z2)
  col.bg_i = "blue",      # Color for background outline from z1's Z component
  col.bg_j = "orange",    # Color for background outline from z2's Z component
  name.axis1 = "PC 1",    # X-axis label
  name.axis2 = "PC 2",    # Y-axis label
  title = "Niche Plot: fus vs chl",  # Plot title
  show.legendOcc = F
)

# Plot Niche Overlap, Schoener's D and Warren et al. I metrics ------------

olaps <- function(sp.i, sp.j, sample_size = 5000, bandwidth = "nrd0", 
                   xlab = "Envir. Grad.", ylab = "Occupancy",
                   show.legend = FALSE) {
  # Scale occupancy values so they sum to 1
  sp.i <- sp.i / sum(sp.i)
  sp.j <- sp.j / sum(sp.j)
  
  # Define the environmental gradient (assumes sp.i and sp.j have equal length)
  env <- seq_along(sp.i)
  
  # Calculate overlap metrics:
  D_val <- 1 - sum(abs(sp.i - sp.j)) / 2
  I_val <- 1 - sum((sqrt(sp.i) - sqrt(sp.j))^2) / 2
  cor_val <- cor(sp.i, sp.j, method = "spearman")
  
  # Compute kernel density for sp.i
  sp.i_samples <- sample(env, size = sample_size, replace = TRUE, prob = sp.i)
  d.i <- density(sp.i_samples, bw = bandwidth, from = min(env), to = max(env))
  scale_factor.i <- max(sp.i) / max(d.i$y)
  d.i$y <- d.i$y * scale_factor.i
  
  # Compute kernel density for sp.j
  sp.j_samples <- sample(env, size = sample_size, replace = TRUE, prob = sp.j)
  d.j <- density(sp.j_samples, bw = bandwidth, from = min(env), to = max(env))
  scale_factor.j <- max(sp.j) / max(d.j$y)
  d.j$y <- d.j$y * scale_factor.j
  
  # Set up the plot with proper limits and the provided axis labels
  y_max <- max(max(sp.i), max(sp.j))
  plot(env, sp.i, type = "n",
       xlab = xlab, ylab = ylab,
       ylim = c(0, y_max * 1.1), main = "")
  
  # Draw the smoothed kernel density lines using the desired colours:
  # sp.i will use "#1f77b4" and sp.j will use "#ff7f0e"
  lines(d.i, col = "#1f77b4", lwd = 3)
  lines(d.j, col = "#ff7f0e", lwd = 3)
  
  # Optionally add a legend
  if (show.legend) {
    legend("topright", legend = c("sp.i", "sp.j"), col = c("#1f77b4", "#ff7f0e"), lwd = 3)
  }
  
  # Annotate the plot with overlap metrics
  usr <- par("usr")  # returns c(x_min, x_max, y_min, y_max)
  x_pos <- usr[1] + 0.05 * diff(usr[1:2])
  y_pos <- usr[4] - 0.05 * diff(usr[3:4])
  text(x = x_pos, y = y_pos,
       labels = paste("D =", round(D_val, 2),
                      "\nI =", round(I_val, 2),
                      "\nCor =", round(cor_val, 2)),
       adj = c(0, 1), cex = 1.2)
  
  invisible(list(D = D_val, I = I_val, cor = cor_val, d.i = d.i, d.j = d.j))
}


# Test olaps function
sp.i_vals <- as.vector(grids[['fus']]$z.cor)
sp.j_vals <- as.vector(grids[['chl']]$z.cor)
res <- olaps(sp.i_vals, sp.j_vals)

# Start Plotting Matrix ---------------------------------------------------
# Species codes = SpatRaster obj names
# Species labs = Full name for plot labels
# Define your species codes and labels
species_codes <- names(grids)
species_labs <- c('Malus fusca', 'Malus coronaria', 'Malus ioensis', 'Malus angustifolia', 'Sect. Chloromeles')

# Close any existing graphics device
while (!is.null(dev.list())) dev.off()

# Open a new high-resolution graphics device
jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/pca_grid_pairs.jpeg',
     width = 3333, height = 3333, res = 300)

# Set up a 5×5 panel with global settings
par(mfrow = c(5, 5), mar = c(3, 3, 2, 1), oma = c(0, 0, 0, 0),
    cex.axis = 1.4, cex.lab = 1.6, cex.main = 1.5)

# Loop over species pairs
for (i in seq_along(species_codes)) {
  for (j in seq_along(species_codes)) {
    
    # Diagonal: species names
    if (i == j) {
      plot.new()
      text(0.5, 0.5, species_labs[i], cex = 1.4, font = 2)
      
      # Lower triangle: PCA niche plots (ecospat.plot.niche.dyn)
    } else if (i > j) {
      ecospat.plot.niche.dyn(
        z1 = grids[[species_codes[j]]],
        z2 = grids[[species_codes[i]]],
        name.axis1 = 'PC 1',
        name.axis2 = 'PC 2',
        col.stab = "#9467bd", 
        col.unf = "#1f77b4", 
        col.abn = "#1f77b4", 
        col.exp = "#ff7f0e", 
        colZ1 = "#1f77b4", 
        colZ2 = "#ff7f0e",
        title = '',
        th.env = 0.05,
        sp.env = 0.05,
        interest = 0  # transparency = 0 for the kernel
      )
      
      # Upper triangle: occupancy kernel plots using the global par settings
    } else {
      # Retrieve bias-corrected occupancy (z.cor)
      sp.i_vals <- as.vector(grids[[species_codes[i]]]$z.cor)
      sp.j_vals <- as.vector(grids[[species_codes[j]]]$z.cor)
      
      olaps(sp.i_vals, sp.j_vals,
            sample_size = 5000,
            bandwidth = "nrd0",
            xlab = "Envir. Grad.",  # shorter label to fit
            ylab = "Occupancy",
            show.legend = TRUE)
    }
  }
}

dev.off()


# Plot with new function --------------------------------------------------
species_codes <- names(grids)
species_labs <- c('Malus fusca', 'Malus coronaria', 'Malus ioensis', 'Malus angustifolia', 'Sect. Chloromeles')

# Close any existing graphics device
while (!is.null(dev.list())) dev.off()

# Open a new high-resolution graphics device
jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/pca_plots/pca_grid_pairs_new.jpeg',
     width = 3333, height = 3333, res = 300)

# Set up a 5×5 panel with global settings
par(mfrow = c(5, 5), mar = c(3, 3, 2, 1), oma = c(0, 0, 0, 0),
    cex.axis = 1.4, cex.lab = 1.6, cex.main = 1.5)

# Loop over species pairs
for (i in seq_along(species_codes)) {
  for (j in seq_along(species_codes)) {
    
    # Diagonal: species names
    if (i == j) {
      plot.new()
      text(0.5, 0.5, species_labs[i], cex = 1.4, font = 2)
      
      # Lower triangle: PCA niche plots (ecospat.plot.niche.dyn)
    } else if (i > j) {
      ecospat.plot.niche.pair(
        z1 = grids[["fus"]], 
        z2 = grids[["chl"]],
        use.zcor = TRUE,        # Toggle between using z.cor (TRUE) and z.uncor (FALSE)
        quant = 0.05,           # Quantile threshold for determining occupied cells
        drawKD = F,          # Draw kernel density shading
        transparency = 50,      # Transparency setting for the density color ramp
        col.shared = "#009E73",  # Color for shared niche region
        col.unique_i = "blue",  # Color for niche unique to species i (z1)
        col.unique_j = "orange",# Color for niche unique to species j (z2)
        col.bg_i = "blue",      # Color for background outline from z1's Z component
        col.bg_j = "orange",    # Color for background outline from z2's Z component
        name.axis1 = "PC 1",    # X-axis label
        name.axis2 = "PC 2",    # Y-axis label
        title = "",  # Plot title
        show.legendOcc = F
      )
      
      # Upper triangle: occupancy kernel plots using the global par settings
    } else {
      # Retrieve bias-corrected occupancy (z.cor)
      sp.i_vals <- as.vector(grids[[species_codes[i]]]$z.cor)
      sp.j_vals <- as.vector(grids[[species_codes[j]]]$z.cor)
      
      olaps(sp.i_vals, sp.j_vals,
            sample_size = 5000,
            bandwidth = "nrd0",
            xlab = "Envir. Grad.",  # shorter label to fit
            ylab = "Occupancy",
            show.legend = TRUE)
    }
  }
}

dev.off()


# PCA Biplot (variable loading) by Species --------------------------------

# Compute global axis limits using "Axis1" and "Axis2" across all species
all_scores <- do.call(rbind, occ_scores_list)
x_min <- min(all_scores$Axis1, na.rm = TRUE)
x_max <- max(all_scores$Axis1, na.rm = TRUE)
y_min <- min(all_scores$Axis2, na.rm = TRUE)
y_max <- max(all_scores$Axis2, na.rm = TRUE)

# Axis labels
xlab_text <- "PC 1"
ylab_text <- "PC 2"

# Set up the plotting area (adjust mfrow as needed)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), asp = 1)

# Scale factor for loadings (optional tweak to improve visibility)
scale_factor <- 1.75

# Loop over each species in occ_scores_list
for (sp in names(occ_scores_list)) {
  occ_data <- occ_scores_list[[sp]]
  
  # Plot occurrence points using Axis1 and Axis2
  plot(occ_data$Axis1, occ_data$Axis2,
       xlim = c(x_min, x_max), ylim = c(y_min, y_max),
       xlab = xlab_text, ylab = ylab_text,
       main = paste("PCA Biplot for", sp),
       pch = 16, col = adjustcolor("darkgray", alpha.f = 0.6))
  
  # Draw reference lines at 0
  abline(h = 0, v = 0, col = "gray60", lty = 2)
  
  # Overlay PCA loadings (from your global PCA stored in pca_full$c1)
  # Multiply by scale_factor to enhance visibility if needed
  arrows(0, 0,
         pca_full$c1[,1] * scale_factor,
         pca_full$c1[,2] * scale_factor,
         length = 0.1, col = "red", lwd = 2)
  
  # Create simplified labels from the original rownames
  simplified_labels <- sub("wc2.1_2.5m_bio_", "Bio ", rownames(pca_full$c1))
  
  # Add text labels for each variable loading
  text(pca_full$c1[,1] * scale_factor,
       pca_full$c1[,2] * scale_factor,
       labels = simplified_labels,
       col = "blue", pos = 3, cex = 1.2)
}


# PCA Bioplot for Sect Chlormeles -----------------------------------------

# --- Subset Occurrence Scores for Eastern Species (Excluding fusca) ---
# Assuming occ_scores_list has names like "Malus fusca", "Malus coronaria", etc.
occ_scores_list_eastern <- occ_scores_list[!(names(occ_scores_list) %in% "Malus fusca")]

# --- Global Axis Limits from Eastern Species Only ---
all_scores_east <- do.call(rbind, occ_scores_list_eastern)
x_min <- min(all_scores_east$Axis1, na.rm = TRUE)
x_max <- max(all_scores_east$Axis1, na.rm = TRUE)
y_min <- min(all_scores_east$Axis2, na.rm = TRUE)
y_max <- max(all_scores_east$Axis2, na.rm = TRUE)

# --- Create Simplified Labels for the Variable Loadings ---
simplified_labels <- sub("wc2.1_2.5m_bio_", "Bio ", rownames(pca_full$c1))

# Optionally, if you computed the proportion of variance:
# (Assuming prop_var exists from your PCA summary code)
var_PC1 <- round(prop_var[1], 1)
var_PC2 <- round(prop_var[2], 1)
xlab_text <- paste0("Axis1 (", var_PC1, "%)")
ylab_text <- paste0("Axis2 (", var_PC2, "%)")

# --- Set up the plotting area (adjust layout as needed) ---
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), asp = 1)

# Loop over only the eastern species
for(sp in names(occ_scores_list_eastern)) {
  occ_data <- occ_scores_list_eastern[[sp]]
  
  # Plot occurrence points using Axis1 and Axis2 (from occ_data)
  plot(occ_data$Axis1, occ_data$Axis2,
       xlim = c(x_min, x_max), ylim = c(y_min, y_max),
       xlab = xlab_text, ylab = ylab_text,
       main = paste("PCA Biplot for", sp),
       pch = 16, col = adjustcolor("darkgray", alpha.f = 0.6))
  
  # Draw dashed reference lines at 0 (helpful for interpretation)
  abline(h = 0, v = 0, col = "gray60", lty = 2)
  
  # Scale the loadings so they are clearly visible relative to occurrence points
  scale_factor <- 1.2
  
  # Overlay arrows for variable loadings from the global PCA object
  arrows(0, 0,
         pca_full$c1[,1] * scale_factor,
         pca_full$c1[,2] * scale_factor,
         length = 0.1, col = "red", lwd = 2)
  
  # Add simplified text labels for the loadings
  text(pca_full$c1[,1] * scale_factor,
       pca_full$c1[,2] * scale_factor,
       labels = simplified_labels,
       col = "blue", pos = 3, cex = 0.8)
}

