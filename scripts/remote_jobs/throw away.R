library(tidyverse)

files <- list.files("./", pattern = "niche_equiv_.*\\.rds", full.names = TRUE)

results <- map_dfr(files, readRDS, .id = "pair")
results$pair <- gsub("niche_equiv_|\\.rds", "", basename(results$pair))

print(results)
