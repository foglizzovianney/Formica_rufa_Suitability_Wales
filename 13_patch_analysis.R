# Script: 13_patch_analysis.R
# Description: Identify ON, SNO and translocation patches for F. rufa in Wales


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(sf)

# PARAMETERS 
suit_threshold <- 0.5
point_buff     <- 50      # metres
dispersal_dist <- 100     # metres
trans_dist_1   <- 2000    # metres
trans_dist_2   <- 5000    # metres

# LOAD DATA 
suitability  <- rast("outputs/suitability_wales_current.tif")
forest_mask  <- rast("spatial/vom/forest_mask_70m_buff.tif")
all_records  <- st_read("data/processed/Formica_rufa_Wales_all_records.gpkg") %>%
  vect()
wales_27700  <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")

# BINARY SUITABILITY 
suit_binary <- ifel(suitability >= suit_threshold, 1, NA) %>%
  mask(forest_mask)

# OCCURRENCE BUFFER 
records_buff_rast <- all_records %>%
  buffer(width = point_buff) %>%
  rasterize(suitability, field = 1)

# ON PATCHES - suitable + within occurrence buffer 
ON_rast <- ifel(!is.na(suit_binary) & !is.na(records_buff_rast), 1, NA)

# SNO PATCHES - suitable + outside occurrence buffer 
SNO_rast <- ifel(!is.na(suit_binary) & is.na(records_buff_rast), 1, NA)

# DISPERSAL GRADIENT - 100m buffer around ON patches 
ON_poly <- as.polygons(ifel(!is.na(ON_rast), 1, NA))
dispersal_rast <- buffer(ON_poly, width = dispersal_dist) %>%
  rasterize(suitability, field = 1)
SNO_dispersal <- ifel(!is.na(SNO_rast) & !is.na(dispersal_rast), 1, NA)

# TRANSLOCATION PATCHES 
trans_2km_rast <- buffer(ON_poly, width = trans_dist_1) %>%
  rasterize(suitability, field = 1)
trans_5km_rast <- buffer(ON_poly, width = trans_dist_2) %>%
  rasterize(suitability, field = 1)

trans_2km <- ifel(!is.na(SNO_rast) & !is.na(trans_2km_rast) &
                    is.na(dispersal_rast), 1, NA)
trans_5km <- ifel(!is.na(SNO_rast) & !is.na(trans_5km_rast) &
                    is.na(trans_2km_rast) & is.na(dispersal_rast), 1, NA)

# SAVE 
writeRaster(ON_rast,       "outputs/ON_patches.tif",        overwrite = TRUE)
writeRaster(SNO_rast,      "outputs/SNO_patches.tif",       overwrite = TRUE)
writeRaster(SNO_dispersal, "outputs/dispersal_patches.tif", overwrite = TRUE)
writeRaster(trans_2km,     "outputs/trans_2km_patches.tif", overwrite = TRUE)
writeRaster(trans_5km,     "outputs/trans_5km_patches.tif", overwrite = TRUE)

# SUMMARY STATISTICS 
ha <- function(r) global(r, fun = "sum", na.rm = TRUE)$sum * 900 / 10000

forest_area <- ha(forest_mask)
suit_area   <- ha(suit_binary)
ON_area     <- ha(ON_rast)
SNO_area    <- ha(SNO_rast)
disp_area   <- ha(SNO_dispersal)
t2km_area   <- ha(trans_2km)
t5km_area   <- ha(trans_5km)

cat("\n=== PATCH ANALYSIS SUMMARY ===\n")
cat("Total forest area:        ", round(forest_area), "ha\n")
cat("Suitable area (>0.5):     ", round(suit_area), "ha (", round(suit_area/forest_area*100, 1), "%)\n")
cat("ON patches:               ", round(ON_area), "ha (", round(ON_area/forest_area*100, 1), "%)\n")
cat("SNO patches:              ", round(SNO_area), "ha (", round(SNO_area/forest_area*100, 1), "%)\n")
cat("Dispersal patches (100m): ", round(disp_area), "ha\n")
cat("Translocation 2km:        ", round(t2km_area), "ha\n")
cat("Translocation 5km:        ", round(t5km_area), "ha\n")

# PLOT 
plot(wales_27700, border = "grey40", 
     main = "F. rufa patch analysis - Wales")
plot(SNO_rast,      col = "orange",    add = TRUE, legend = FALSE)
plot(SNO_dispersal, col = "lightblue", add = TRUE, legend = FALSE)
plot(trans_2km,     col = "red",       add = TRUE, legend = FALSE)
plot(trans_5km,     col = "darkred",   add = TRUE, legend = FALSE)
plot(ON_rast,       col = "darkgreen", add = TRUE, legend = FALSE)
plot(all_records,   col = "blue", pch = 16, cex = 0.8, add = TRUE)

legend("bottomright",
       legend = c("ON - likely occupied",
                  "Dispersal (100m)",
                  "Translocation 2km",
                  "Translocation 5km",
                  "SNO - suitable unoccupied",
                  "Occurrence records"),
       fill = c("darkgreen", "lightblue", "red", "darkred", "orange", NA),
       pch  = c(NA, NA, NA, NA, NA, 16),
       col  = c(NA, NA, NA, NA, NA, "blue"),
       cex  = 0.8)

