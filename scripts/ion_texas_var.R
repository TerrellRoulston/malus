# Exploratory script
# Get counties from Points of Malus ionesis isolated in Texas

# Libs
library(terra)
library(tigris)  # For US county shapefiles
options(tigris_use_cache = TRUE)

# Your coordinates
coords <- matrix(c(
  -98.72972, 29.78556,
  -98.603547, 29.97294,
  -98.880074, 30.42365,
  -99.273229, 29.760225,
  -98.542258, 30.109597,
  -98.94185, 30.325093,
  -98.399217, 30.266456,
  -98.753121, 30.377799,
  -98.568931, 30.041149,
  -98.603547, 29.97294,
  -98.55565, 30.097918,
  -98.386405, 29.816091,
  -98.4105, 29.840551,
  -98.387154, 29.818214,
  -98.344116, 30.10593
), ncol = 2, byrow = TRUE)

# Create a SpatVector from the coordinates
pts <- vect(coords, type = "points", crs = "EPSG:4326")

# Download US counties
counties <- tigris::counties(cb = TRUE, year = 2022)

# Convert to SpatVector
counties_vect <- vect(counties)

# Intersect to get county names
county_names <- extract(counties_vect, pts)

# Print county names
print(unique(county_names$NAME))
