# Top ---------------------------------------------------------------------
# Conservation Gap Analysis
# Started September 26th, 2024

library(tidyverse) # Grammar and data management
library(terra) # Spatial Data package

# The following gap analysis is referencing Carver et al. (2021)
# https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.05430

# For the purposes of the analysis we are restricting the study area to just Canada

data <- terra::vect("C:/Users/terre/Documents/Acadia/Malus Project/maps/canada_pa/ProtectedConservedArea_2023/ProtectedConservedArea_2023/ProtectedConservedArea_2023.gdb")

plot(data)
