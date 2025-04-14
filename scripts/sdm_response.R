library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package
library(predicts) # SDM package
library(geodata) # basemaps
library(rJava) # MaxEnt models are dependant on JDK
library(ENMeval) # Another modeling package, useful for data partitioning (Checkerboarding)
library(dismo) # For plotting response curves (accepts Maxent objects)


# Load Maxent Models ------------------------------------------------------
# Note these are the models made with the subsetted variables described in the SDM script
cor_maxent <- readRDS(file = './sdm_output/cor/subs/cor_maxent_subs.Rdata')
fus_maxent <- readRDS(file = './sdm_output/fus/subs/fus_maxent_subs.Rdata')
ion_maxent <- readRDS(file = './sdm_output/ion/subs/ion_maxent_subs.Rdata')
ang_maxent <- readRDS(file = './sdm_output/ang/subs/ang_maxent_subs.Rdata')
chl_maxent <- readRDS(file = './sdm_output/chl/subs/chl_maxent_subs.Rdata')


# Extract best model ------------------------------------------------------
best_cor_maxent <- subset(cor_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_cor_maxent <- eval.models(cor_maxent)[[best_cor_maxent$tune.args]] # extracts the best model - returns MaxEnt object

best_fus_maxent <- subset(fus_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_fus_maxent <- eval.models(fus_maxent)[[best_fus_maxent$tune.args]] # extracts the best model - returns MaxEnt object

best_ion_maxent <- subset(ion_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_ion_maxent <- eval.models(ion_maxent)[[best_ion_maxent$tune.args]] # extracts the best model - returns MaxEnt object

best_ang_maxent <- subset(ang_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_ang_maxent <- eval.models(ang_maxent)[[best_ang_maxent$tune.args]] # extracts the best model - returns MaxEnt object

best_chl_maxent <- subset(chl_maxent@results, delta.AICc == 0) # selects the best performing model based on delta AICc - returns data frame object
mod.best_chl_maxent <- eval.models(chl_maxent)[[best_chl_maxent$tune.args]] # extracts the best model - returns MaxEnt object


# Variable importance plots -----------------------------------------------

plot(mod.best_cor_maxent)

# SDM Response Curve ------------------------------------------------------
# Species info: each element is a list with
#   $model : the fitted SDM
#   $env   : the environmental predictors (SpatRaster or RasterStack)
species_info <- list(
  "Malus fusca"        = list(model = mod.best_fus_maxent, env = wclim_fus_subs),
  "Malus coronaria"    = list(model = mod.best_cor_maxent, env = wclim_cor_subs),
  "Malus ioensis"      = list(model = mod.best_ion_maxent, env = wclim_ion_subs),
  "Malus angustifolia" = list(model = mod.best_ang_maxent, env = wclim_ang_subs),
  "Sect. Chloromeles"  = list(model = mod.best_chl_maxent, env = wclim_chl_subs)
)

# Define the order you want the species to appear in (each row)
species_order <- c("Malus fusca", "Malus coronaria", "Malus ioensis",
                   "Malus angustifolia", "Sect. Chloromeles")

# Get the names of the predictors (assume they are identical across species)
predictor_names <- names(species_info[[1]]$env)

# Define custom labels for the predictors (to display on the top of each column)
custom_names <- c("Bio 1 (Tmean)", "Bio 4 (Tvar)", "Bio 10 (Twarmq)", 
                  "Bio 11 (Tcoldq)", "Bio 15 (Pvar)", "Bio 16 (Pwetq)")


# Define the desired species order (rows, from top to bottom)
species_order <- names(species_info)

preds0 <- names(species_info[[1]]$env)

global_ranges <- list()
for (p in preds0) {
  global_min <- min(sapply(species_info, function(sp) min(values(sp$env[[p]]), na.rm = TRUE)))
  global_max <- max(sapply(species_info, function(sp) max(values(sp$env[[p]]), na.rm = TRUE)))
  global_ranges[[p]] <- c(global_min, global_max)
}
names(global_ranges) <- preds0

# Save the file 
# Close any existing graphics device
while (!is.null(dev.list())) dev.off()

# Open a new high-resolution graphics device
jpeg(file = 'C:/Users/terre/Documents/Acadia/Malus Project/sdm_plot/sdm_response_curves.jpeg',
     width = 3333, height = 3333, res = 300)

par(mfrow = c(5, 6), mar = c(3, 1, 4, 1), oma = c(0, 5, 0, 0))


# Loop Over Species & Predictors to Create Response Curves
for (sp in species_order) {
  
  # Extract the model and environmental predictors for the current species.
  mod <- species_info[[sp]]$model
  env <- species_info[[sp]]$env
  preds <- names(env)
  
  for (j in seq_along(preds)) {
    
    varname <- preds[j]
    
    # Use the global range for this predictor (consistent x-axis across species)
    xlim <- global_ranges[[varname]]
    xseq <- seq(xlim[1], xlim[2], length.out = 100)
    
    # For all predictors, use the species-specific mean values.
    const_vals <- sapply(preds, function(p) mean(values(env[[p]]), na.rm = TRUE))
    
    # Construct a data frame for prediction:
    # All predictors are set to their mean except varname, which is varied across xseq.
    newdata <- as.data.frame(matrix(rep(const_vals, each = 100),
                                    nrow = 100, byrow = FALSE))
    colnames(newdata) <- preds
    newdata[[varname]] <- xseq
    
    # Obtain predictions from your model.
    predsuit <- predict(mod, newdata)
    
    # Draw the plot: set up axes manually.
    # The y-axis is fixed from 0 to 1.
    plot(xseq, predsuit, type = "n", ylim = c(0, 1),
         xlab = "", ylab = "", main = "", axes = FALSE)
    lines(xseq, predsuit, col = "red", lwd = 2)
    box()
    
    # Draw x-axis with default tick marks.
    axis(1)
    # Draw a y-axis with ticks at 0.0, 0.6, and 1.0.
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0.0", "0.25", "0.5", "0.75", "1.0"))
    
    # --- Add Blue Rug Ticks (Deciles) ---
    # Compute deciles from the actual training values for varname.
    var_values <- values(env[[varname]])
    var_values <- var_values[!is.na(var_values)]
    deciles <- quantile(var_values, probs = seq(0.1, 0.9, by = 0.1))
    rug(deciles, col = "blue", ticksize = 0.05)
    
    # --- Labeling ---
    # For the left-most plot in each species row (first predictor), label with the species name.
    if (j == 1) {
      # Increase the line number (e.g., line = 5) so that the species name is further to the right,
      # leaving room for the global y-axis label.
      mtext(sp, side = 2, line = 5, cex = 0.9)
    }
    # For the top row (first species), add the custom predictor label at the top.
    if (sp == species_order[1]) {
      mtext(custom_names[j], side = 3, line = 1.5, cex = 0.9)
    }
  } # end for each predictor
} # end for each species


# Add Global Y-Axis Label in Outer Margin
mtext("Predicted Suitability", side = 2, outer = TRUE, line = 2.5, cex = 1.1)

dev.off()
