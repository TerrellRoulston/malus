library(tidyverse)
library(ggplot2)

# Color palette
fill_cols <- c('#C7E9B4', "#7FCDBB", "#41AE76", '#0868AC')

# Load and prepare data
in_situ <- read_csv("malus_gap_analysis_all_species.csv") %>%
  filter(suitability == "high", ssp %in% c("historical", "585"), period %in% c(2000, 2030, 2050, 2070)) %>%
  mutate(species = case_when(
    species == "Malus coronaria"      ~ "italic('Malus coronaria')",
    species == "Malus fusca"          ~ "italic('Malus fusca')",
    species == "Malus ioensis"        ~ "italic('Malus ioensis')",
    species == "Malus angustifolia"   ~ "italic('Malus angustifolia')",
    species == "Sect. Chloromeles"    ~ "'Sect. Chloromeles'"
  )) %>%
  mutate(species = factor(species, levels = c(
    "italic('Malus fusca')",
    "italic('Malus coronaria')",
    "italic('Malus ioensis')",
    "italic('Malus angustifolia')",
    "'Sect. Chloromeles'"
  ))) %>%
  pivot_longer(cols = c(SRSin, GRSin, ERSin, FCSin), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("SRSin", "GRSin", "ERSin", "FCSin")),
         period = factor(period, levels = c(2000, 2030, 2050, 2070),
                         labels = c("1970–2000", "2030", "2050", "2070")))


# # ggplot figure -----------------------------------------------------------
# # Plot
# facet_plot <- ggplot(in_situ, aes(x = period, y = value, fill = metric)) +
#   geom_col(position = position_dodge()) +
#   facet_wrap(~ species, ncol = 2, scales = 'free', labeller = label_parsed) +
#   scale_fill_manual(values = fill_cols,
#                     breaks = c("SRSin", "GRSin", "ERSin", "FCSin")) +
#   scale_y_continuous(limits = c(0, 105),
#                      breaks = seq(0, 100, 25),
#                      expand = c(0, 0)) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 22, colour = "black"),
#     axis.text = element_text(colour = "black"),
#     axis.title.x = element_blank(),
#     legend.position = c(0.55, 0.02),
#     legend.justification = c("left", "bottom"),
#     legend.direction = "vertical",
#     legend.background = element_rect(fill = "white", color = NULL),
#     legend.key.size = unit(1.25, "cm"),
#     legend.text = element_text(size = 20),
#     strip.text = element_text(size = 20),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   ylab(bquote(atop("SSP585", "Conservation Score"))) +
#   geom_hline(yintercept = 25, linetype = "dashed", color = "black", linewidth = 1.1) +
#   guides(fill = guide_legend(title = expression(italic("In Situ") ~ "Indicator")))
# 
# facet_plot
# 
# ggsave(plot = facet_plot, filename = "C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/gap_analysis/gap_analysis_plot.png", width = 5000, height = 3333, dpi = 300, units = 'px')
# 


# Making plot in base R ---------------------------------------------------
png("C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/gap_analysis/gap_analysis_baseR_plot.png", width = 5000, height = 3333, res = 300)

# Base R version
fill_cols <- c("SRSin" = '#C7E9B4', "GRSin" = "#7FCDBB", "ERSin" = "#41AE76", "FCSin" = '#0868AC')
period_levels <- c("1970–2000", "2030", "2050", "2070")
metric_levels <- c("SRSin", "GRSin", "ERSin", "FCSin")

# Create layout: 3 rows x 2 columns
par(mfrow = c(3, 2), mar = c(4, 4.5, 2, 1), oma = c(1, 6.5, 1, 0))

# Species values as they exist in the data
species_codes <- c(
  "italic('Malus fusca')",
  "italic('Malus coronaria')",
  "italic('Malus ioensis')",
  "italic('Malus angustifolia')",
  "'Sect. Chloromeles'"
)

# Plotting labels
species_labels <- c(
  "italic('Malus fusca')"        = "italic('Malus fusca')",
  "italic('Malus coronaria')"    = "italic('Malus coronaria')",
  "italic('Malus ioensis')"      = "italic('Malus ioensis')",
  "italic('Malus angustifolia')" = "italic('Malus angustifolia')",
  "'Sect. Chloromeles'"          = expression("Sect. *italic(' Chloromeles')")
)

for (sp in species_codes) {
  dat <- subset(in_situ, species == sp)
  
  # Create empty plot
  plot(NA, xlim = c(0.5, 4.5), ylim = c(0, 105),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = NULL, bty = "n")
  
  # Add y-axis and x-axis
  axis(2, at = seq(0, 100, 25), las = 1, cex.axis = 2.25)
  axis(1, at = 1:4, labels = period_levels, cex.axis = 2.25)
  
  # Plot bars
  for (i in 1:4) {
    for (j in 1:4) {
      value <- dat$value[dat$period == period_levels[i] & dat$metric == metric_levels[j]]
      if (length(value) == 1) {
        rect(xleft = i - 0.4 + (j - 1) * 0.2,
             xright = i - 0.4 + j * 0.2,
             ybottom = 0, ytop = value,
             col = fill_cols[metric_levels[j]], border = NA)
      }
    }
  }
  
  abline(h = 25, lty = 2, lwd = 2)
  
  # Add species label
  mtext(parse(text = species_labels[sp]), side = 3, line = 0.5, cex = 2)
}

# Add common y-axis title
mtext("SSP585\nConservation Score", side = 2, outer = TRUE, line = 0.5, cex = 2)

# Legend
# Leave previous plotting area as is
par(xpd = NA)  # Allow drawing outside plot region

# Place legend in user coordinates (outside plot layout)
legend(x = 5, y = 85,  # Adjust as needed
       legend = metric_levels, 
       fill = fill_cols[metric_levels],
       title = expression(italic("In Situ") ~ "Indicator"),
       horiz = FALSE, bty = "n", cex = 2.5)
dev.off()