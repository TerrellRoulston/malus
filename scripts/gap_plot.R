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

# Plot
facet_plot <- ggplot(in_situ, aes(x = period, y = value, fill = metric)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ species, ncol = 2, scales = 'free', labeller = label_parsed) +
  scale_fill_manual(values = fill_cols,
                    breaks = c("SRSin", "GRSin", "ERSin", "FCSin")) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(size = 22, colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.02),
    legend.justification = c("left", "bottom"),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = NULL),
    legend.key.size = unit(1.25, "cm"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5)
  ) +
  ylab(bquote(atop("SSP585", "Conservation Score"))) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "black", linewidth = 1.1) +
  guides(fill = guide_legend(title = expression(italic("In Situ") ~ "Indicator")))

facet_plot

ggsave(plot = facet_plot, filename = "C:/Users/terre/Documents/Acadia/Malus Project/statistical analysis/gap_analysis/gap_analysis_plot.png", width = 5000, height = 3333, dpi = 300, units = 'px')
