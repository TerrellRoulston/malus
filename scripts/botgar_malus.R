# Top ---------------------------------------------------------------------
# Terrell Roulston
# Started Oct. 10th 2024

library(tidyverse)
library(readxl)
library(geodata)
library(terra)
library(tidyterra)

# Load data
malus_data <- read_excel("C:/Users/terre/Documents/Acadia/Malus Project/ex situ data/botanical_garden_data.xlsx")

# calculate some statistics about the data

malus_data %>% filter(Class == 'Native, Wild') %>% nrow()

# Summarize data for visualization

malus_summary <- malus_data %>% 
  filter(Class != 'Unknown') %>% # remove 26 records that have insufficenit description - listed as only 'Malus' or 'Malus sp.'
  group_by(Province, Class) %>% 
  summarize(Total = n()) %>% 
  ungroup()


# Load Basemap data -------------------------------------------------------
# Load Lambert projected data from Stats Canada
can_lambert <- vect('C:/Users/terre/Documents/Acadia/Malus Project/maps/canada_lambert')


# Prep spatial data --------------------------------------------------------
# Rename "Newfoundland" to "Newfoundland and Labrador" to be compatiable with the spatial data

total_y <- c(280, 198, 127, 69, 38, 33, 17) # Set height for text of overall totals (I use a +15 offset above the greatest class per province)
province_order <- c('Ontario', 'British Columbia', 'Quebec', 'Alberta', 'Saskatchewan', 'Nova Scotia', 'Newfoundland and Labrador')

malus_summary <- malus_summary %>%
  mutate(Province = ifelse(Province == "Newfoundland", "Newfoundland and Labrador", Province)) %>%
  mutate(Province = factor(Province, levels = province_order))

total_species_per_province <- malus_summary %>%
  group_by(Province) %>%
  summarize(Total_Species = sum(Total)) %>%
  mutate(Province_numeric = as.numeric(factor(Province, levels = province_order)))

malus_summary <- malus_summary %>% left_join(total_species_per_province)

can_map_sp <- can_lambert %>% 
  left_join(total_species_per_province, by = c('PRENAME' = 'Province'))

cbind(total_species_per_province, total_y)

# Choropleth map ----------------------------------------------------------

dev.new() # open graphics window

ggplot() +
  geom_spatvector(data = can_map_sp, aes(fill = Total_Species), color = "white") +
  scale_fill_gradient(low = "#ccece6", high = "#00441b", na.value = "lightgray") + 
  theme_minimal() +
  labs(title = NULL,
       fill = "Total\nAccessions") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    legend.key.size = unit(1.5, "cm"),         
    legend.spacing = unit(0.5, "cm")           
  )


# Accession boxplot -------------------------------------------------------

ggplot() +
  geom_col(data = malus_summary, aes(x = Province, y = Total, fill = Class), position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(data = malus_summary, aes(label = Total, group = Class, x = Province, y = Total), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5,  # Adjust vertical position above the bars
            size = 8) +
  geom_text(data = total_species_per_province, aes(x = Province, y = total_y, label = Total_Species), vjust = -0.5, size = 8) +
  geom_segment(data = total_species_per_province, aes(x = Province_numeric - 0.4, xend = Province_numeric + 0.4, y = total_y - 5, yend = total_y - 5), color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = NULL,
       x = NULL, y = "Number of Accessions",
       fill = 'Accession Type') +
  scale_fill_manual(values = c('#e6f5d0', '#a1d99b', '#41ab5d', '#006d2c'))+
  theme(
    legend.key.size = unit(1.5, "cm"),         
    legend.text = element_text(size = 24, color = 'black'),     
    legend.title = element_text(size = 26, color = 'black'),
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(size = 24, color = 'black'),     
    axis.text.y = element_text(size = 24, color = 'black'),     
    axis.title.y = element_text(size = 24, color = 'black')
  ) 

  
