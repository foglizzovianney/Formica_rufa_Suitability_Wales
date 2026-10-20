# Script: 07_distance_forest_edge.R
# Description: Calculate distance to forest edge and forest mask buffer

setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)

# LOAD DATA 
wales_27700 <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")

cover_vom <- rast("spatial/vom/cover_VOM_30m_27700.tif")

# CREATE FOREST MASK
# Forest = canopy cover above 30% following Robinson et al. 2025
forest_mask <- cover_vom
forest_mask[forest_mask < 0.3] <- NA
forest_mask[!is.na(forest_mask)] <- 1

# CALCULATE DISTANCE TO FOREST EDGE
distance_forest <- terra::distance(forest_mask, unit = "m") %>%
  mask(wales_27700)

# CREATE BUFFERED FOREST MASK
# 70m buffer beyond forest edge - zone where 95% of nests were found
forest_mask_buff <- ifel(distance_forest <= 70, 1, NA)

# SAVE
writeRaster(distance_forest,
            "spatial/vom/distance_forest_30m.tif",
            overwrite = TRUE)

writeRaster(forest_mask_buff,
            "spatial/vom/forest_mask_70m_buff.tif",
            overwrite = TRUE)

# PLOT
par(mfrow = c(1, 2))

plot(distance_forest,
     main = "Distance to forest edge (m)",
     col = hcl.colors(20, "YlOrRd"),
     axes = FALSE)

plot(forest_mask_buff,
     main = "Forest mask + 70m buffer",
     col = "darkgreen",
     axes = FALSE)
plot(wales_27700, add = TRUE, border = "grey40")

