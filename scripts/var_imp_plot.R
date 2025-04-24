# Top ---------------------------------------------------------------------
# Plotting SDM Variable Permutation Importance
# Converting table of data into a figure to increase ease of comparison
# Terrell Roulston
# Started April 15 2025

library(tidyverse)
library(readxl)
library(tidytext) # for custom labeling of facets with italics
library(cowplot)
library(ggh4x)

# Import data -------------------------------------------------------------
importance_df <- read_xlsx("C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/variable importance/sdm_vars_R_import.xlsx")

# Pivoting from wide to long format
importance_long <- importance_df %>%
  pivot_longer(
    cols = starts_with("Bio"), 
    names_to = "Variable", 
    values_to = "Importance"
  ) %>%
  # Reverse the factor order for a top-down ordering of variables
  mutate(Variable = factor(Variable, levels = rev(unique(Variable))))

# Facet by Variables ------------------------------------------------------

# Function for species that italicizes the entire name unless it's "Sect. Chloromeles"
italicize_species <- function(x) {
  # Remove the suffix that reorder_within() attaches (e.g., "___Bio_1")
  x_clean <- gsub("___.*$", "", x)
  
  sapply(x_clean, function(spec) {
    if (spec == "Sect. Chloromeles") {
      spec  # Leave unchanged
    } else {
      # Add an extra trailing space inside the quotes so that the italicized text has some padding.
      parse(text = paste0('italic("', spec, ' ")'))
    }
  }, USE.NAMES = FALSE)
}



# Factor order or variables for plotting
importance_long_reorder <- importance_long %>%
  mutate(Variable = factor(Variable, levels = c("BIO1", "BIO4", "BIO10", "BIO11", "BIO15", "BIO16")))

importance_ranked <- importance_long_reorder %>%
  group_by(Species) %>%
  # "min_rank(desc())" sets the biggest Importance to rank=1, next biggest=2, etc.
  mutate(Rank = min_rank(desc(Importance))) %>%
  ungroup() %>% 
  group_by(Variable) %>%
  mutate(VarRank = min_rank(desc(Importance))) %>%
  ungroup()

# Reorder species within each variable facet using tidytext's reorder_within():
plot_data_by_var <- importance_ranked %>%
  mutate(Species_reordered = reorder_within(Species, Importance, Variable))

p <- ggplot(plot_data_by_var, aes(x = Importance, y = Species_reordered, color = Species)) +
  geom_segment(aes(x = 0, xend = Importance, 
                   y = Species_reordered, yend = Species_reordered),
               linetype = "dashed", linewidth = 1.25) +
  geom_point(size = 5) +
  geom_text(aes(label = Rank), 
            nudge_x = 2,
            size = 5, 
            color = "black",
            fontface = "bold") +
  # Facet by Variable in a single column, free y to allow separate ordering per facet
  facet_wrap2(~ Variable, ncol = 1, scales = "free_y") +
  # Remove reorder_within's suffix from the y-axis labels AND italicize
  scale_y_reordered(labels = italicize_species) +
  # Also apply the same style to the legend labels
  scale_color_discrete(labels = italicize_species) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Variable's Permutational Importance (%)",
    y = "Taxon"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  theme(text = element_text(family = "Arial")) +
  theme(axis.text = element_text(size = 16, colour = 'black'),
        legend.text = element_text(size = 16),
        axis.title  = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 20), size = 20),
        strip.text.x = element_text(size = 18))


# Export to PNG -----------------------------------------------------------

# Open png device with your settings
png(filename = "C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/variable importance/var_important_draft_plot.png", 
     width = 4000, height = 3333, res = 300)

# Print your plot (if using ggplot, make sure the plot object is printed)
print(p)

# Close the device
dev.off()

# Plot with labels
p_tagged <- ggdraw(p) +
  draw_plot_label(label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                  x = c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95),  # adjust x-coordinates as needed
                  y = c(0.99, 0.85, 0.69, 0.54, 0.39, 0.235),  # adjust y-coordinates for each facet
                  size = 22,  # adjust the size
                  fontfamily = "Arial")

png(filename = "C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/variable importance/var_important_draft_plot_v2.png", 
   width = 4000, height = 3333, res = 300)

print(p_tagged)

dev.off()


# Var Importance by Species -----------------------------------------------
# Now lets plot faceting by species
# Reorder the variables from most to least important by species

# Colour blind friendly pallette
cbPalette <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00"   # vermilion
)

# Parse text for labeller of Facetes
species_exprs <- c(
  "Malus fusca"         = "italic(Malus~fusca)",
  "Malus coronaria"     = "italic(Malus~coronaria)",
  "Malus ioensis"       = "italic(Malus~ioensis)",
  "Malus angustifolia"  = "italic(Malus~angustifolia)",
  "Sect. Chloromeles"   = "\"Sect. Chloromeles\""       # plain text
)

# Prepare data
importance_ranked <- importance_ranked %>%
  # set Species factor levels in the desired facet order
  mutate(Species = factor(Species, levels = c(
    "Malus fusca",
    "Malus coronaria",
    "Malus ioensis",
    "Malus angustifolia",
    "Sect. Chloromeles"
  )))

plot_data_by_sp <- importance_ranked %>%
  mutate(Variable_reordered = reorder_within(Variable, Importance, Species)) 

plot_data_by_sp$Species <- factor(plot_data_by_sp$Species,
                            levels = names(species_exprs),
                            labels = species_exprs)

# PLOT
# Plot species var importance
p_sp <- ggplot(plot_data_by_sp, aes(x = Importance, y = Variable_reordered, colour = Variable)) +
  # dashed guides
  geom_segment(aes(x = 0, xend = Importance,
                   y = Variable_reordered, yend = Variable_reordered),
               linetype = "dashed", linewidth = 1.25) +
  # points colored by Variable
  geom_point(aes(color = Variable), size = 5) +
  # rank labels
  #geom_text(aes(label = VarRank), nudge_x = 1, fontface = "bold", size = 5, color = "black") +
  geom_text(aes(x = 0, label = VarRank),
            hjust = 4, vjust = 0.5,
            fontface = "bold", size = 4,
            color = "black") +
  # facet by Species
  facet_wrap(~ Species, 
             ncol   = 1, 
             scales = "free_y",
             labeller = label_parsed) +
  # clean up the reorder_within labels
  scale_y_reordered() +
  # a colorblind‐friendly palette for variables
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = cbPalette, name = "Variable") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  labs(
    x = "Variable's Permutational Importance (%)",
    y = "Bioclimatic Variable",
    color = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  theme(text = element_text(family = "Arial")) +
  theme(axis.text = element_text(size = 16, colour = 'black'),
        legend.text = element_text(size = 16),
        axis.title  = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 20), size = 20),
        strip.text.x = element_text(size = 18))

# Print or save
print(p_sp)

png(filename = "C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/variable importance/var_important_species.png", 
    width = 4000, height = 3333, res = 300)

# Print your plot (if using ggplot, make sure the plot object is printed)
print(p_sp)

# Close the device
dev.off()
