# Top ---------------------------------------------------------------------
# Script for running Niche Equivalency test remotely on HPC
# run_malus_pca_equiv.R
# Schedule with RUN_PCA.sh
# Terrell Roulston
# May 14 2025 
# Note this script has be modified to run on HPC and differs from malus_pca slighlty

args <- commandArgs(trailingOnly = TRUE)
pair_id <- as.numeric(args[1])
cat("Starting run for pair ID:", pair_id, "\n")

setwd("/project/6074193/mig_lab/malus_pca/")

# Load packages
cat("Loading packages...\n")
library(tidyverse)
library(terra)
library(geodata)
library(ade4)
library(ecospat)

cat("Loading inputs...\n")
# Global variables
vars <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_4', 'wc2.1_2.5m_bio_10',
          'wc2.1_2.5m_bio_11', 'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_16')

# Species input file map
species_map <- list(
  cor = list(occ = "occThin_cor.Rdata", clim = "wclim_cor.Rdata"),
  ion = list(occ = "occThin_ion.Rdata", clim = "wclim_ion.Rdata"),
  ang = list(occ = "occThin_ang.Rdata", clim = "wclim_ang.Rdata")
)

# Pairwise tests (species codes)
pairwise_tests <- list(
  cor_ion = c("cor", "ion"),
  cor_ang = c("cor", "ang"),
  ion_ang = c("ion", "ang")
)

# Resolve species pair
pair_label <- names(pairwise_tests)[pair_id + 1]
species_pair <- pairwise_tests[[pair_id + 1]]
sp1 <- species_pair[1]
sp2 <- species_pair[2]

cat("Running pair:", pair_label, "→", sp1, "vs", sp2, "\n")

# Load full background for PCA (Sect. Chloromeles)
cat("Loading background for PCA...\n")
wclim_chl <- readRDS("./inputs/wclim_data/wclim_chl.Rdata") %>%
  terra::subset(vars)
bg_mat_full <- values(wclim_chl) %>% na.omit()

# Run PCA
cat("Performing PCA...\n")
pca_full <- dudi.pca(bg_mat_full, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
pca_score <- pca_full$li

# Function to load and process species data
load_species_data <- function(code) {
  occ_path <- paste0("./inputs/occ_data/", species_map[[code]]$occ)
  clim_path <- paste0("./inputs/wclim_data/", species_map[[code]]$clim)
  
  cat("  Loading data for:", code, "\n")
  occ <- readRDS(occ_path)
  clim <- readRDS(clim_path) %>% terra::subset(vars)
  
  bg_mat <- values(clim) %>% na.omit()
  occ_vals <- cbind(crds(occ), extract(clim, occ, ID = FALSE)) %>% na.omit()
  
  list(
    occ_score = suprow(pca_full, data.frame(occ_vals)[, colnames(bg_mat_full)])$li,
    bg_score = suprow(pca_full, bg_mat)$li
  )
}

# Load + project species 1 and 2
sp1_data <- load_species_data(sp1)
sp2_data <- load_species_data(sp2)

# Define test function
niche_equivalency_test <- function(sp1_scores, sp2_scores,
                                               sp1_bg_scores, sp2_bg_scores,
                                               bg_scores, R = 100, reps = 999,
                                               quant = 0.1,
                                               parallel = TRUE, ncores = 2,
                                               verbose = TRUE) {
  require(ecospat)
  require(foreach)
  start_time <- Sys.time()
  
  if (parallel) {
    require(doParallel)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
  }
  
  grid1 <- ecospat.grid.clim.dyn(glob = bg_scores, glob1 = sp1_bg_scores, sp = sp1_scores, R = R)
  grid2 <- ecospat.grid.clim.dyn(glob = bg_scores, glob1 = sp2_bg_scores, sp = sp2_scores, R = R)
  
  obs_vals <- tryCatch(ecospat.niche.overlap(grid1, grid2, cor = FALSE), error = function(e) return(NULL))
  if (is.null(obs_vals) || !is.finite(obs_vals$D) || !is.finite(obs_vals$I)) {
    stop("Observed overlap metrics could not be calculated.")
  }
  
  obs_D <- obs_vals$D
  obs_I <- obs_vals$I
  
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
    
    vals <- tryCatch(ecospat.niche.overlap(g1, g2, cor = FALSE), error = function(e) return(c(D = NA, I = NA)))
    c(D = vals$D, I = vals$I)
  }
  
  null_vals <- if (parallel) {
    foreach(i = 1:reps, .combine = rbind, .packages = "ecospat") %dopar% permute_fun(i)
  } else {
    foreach(i = 1:reps, .combine = rbind, .packages = "ecospat") %do% permute_fun(i)
  }
  
  D_null <- null_vals[, "D"]
  I_null <- null_vals[, "I"]
  p_D <- mean(D_null <= obs_D, na.rm = TRUE)
  p_I <- mean(I_null <= obs_I, na.rm = TRUE)
  duration <- difftime(Sys.time(), start_time, units = "secs")
  
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

# Run the test
cat("Running niche equivalency test for:", pair_label, "\n")
result <- niche_equivalency_test(
  sp1_scores = sp1_data$occ_score,
  sp2_scores = sp2_data$occ_score,
  sp1_bg_scores = sp1_data$bg_score,
  sp2_bg_scores = sp2_data$bg_score,
  bg_scores = pca_score,
  reps = 999,
  R = 100,
  parallel = TRUE,
  ncores = 15,
  verbose = TRUE
)

# Save results
output_path <- paste0("results/niche_equiv_", pair_label, ".rds")
cat("Saving result to", output_path, "\n")
saveRDS(result, file = output_path)

cat("Job complete for pair:", pair_label, "\n")