# Script: 03_topographic_variables.R
# Description: Derive topographic variables from NASA SRTM 30m DEM

setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)

# Load Wales boundary
wales_27700 <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")

# Load SRTM VOM
srtm <- rast("spatial/vom/srtm_wales_30m_27700.tif")

# Derive topographic variables
aspect_30m <- terrain(srtm, v = "aspect")
northness_30m <- cos(aspect_30m * pi / 180)
eastness_30m <- sin(aspect_30m * pi / 180)
slope_30m <- terrain(srtm, v = "slope")
hillshade_30m <- shade(slope_30m, aspect_30m, angle = 45)

# Stack all together
topo_stack <- c(northness_30m, eastness_30m, hillshade_30m, slope_30m)
names(topo_stack) <- c("northness", "eastness", "hillshade", "slope")

# Save
writeRaster(topo_stack, 
            "spatial/vom/topo_stack_30m_27700.tif", 
            overwrite = TRUE)

