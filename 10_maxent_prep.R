# Script: 10_maxent_prep.R
# Description: Prepare occurrence and background points for MaxEnt


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(sf)
library(dismo)
library(raster)

# PARAMETERS 
n_background <- 10000 
seed <- 42
thin_dist <- 100  # 100m thinning

# LOAD DATA
predictors <- rast("spatial/vom/predictor_stack_wales.tif")
effort_rast <- raster("spatial/vom/effort_masked_30m.tif")
forest_mask_buff <- rast("spatial/vom/forest_mask_70m_buff.tif")

# Load occurrence records
all_records <- st_read("data/processed/Formica_rufa_Wales_all_records.gpkg")

# Convert to dataframe with coordinates
occs <- st_coordinates(all_records) %>%
  as.data.frame() %>%
  rename(x = X, y = Y)

cat("Number of occurrence records:", nrow(occs), "\n")

# BACKGROUND POINTS weighted by effort
set.seed(seed)
bg_bias <- dismo::randomPoints(mask = effort_rast, 
                               p = as.matrix(occs),
                               n = n_background, 
                               prob = TRUE)

bg_bias_df <- as.data.frame(bg_bias)
names(bg_bias_df) <- c("x", "y")

cat("Background points selected:", nrow(bg_bias_df), "\n")

# VALIDATION BACKGROUND POINTS
set.seed(seed + 100)
val_bg_bias <- dismo::randomPoints(mask = effort_rast,
                                   p = as.matrix(occs),
                                   n = n_background,
                                   prob = TRUE)

val_bg_bias_df <- as.data.frame(val_bg_bias)
names(val_bg_bias_df) <- c("x", "y")

# SAVE 
write.csv(occs, "data/processed/occs_maxent.csv", row.names = FALSE)
write.csv(bg_bias_df, "data/processed/bg_bias_n10000.csv", row.names = FALSE)
write.csv(val_bg_bias_df, "data/processed/val_bg_bias_n10000.csv", row.names = FALSE)

# PLOT to verify
plot(forest_mask_buff, col = "lightgreen", axes = FALSE,
     main = "Occurrence records and background points")
points(bg_bias_df$x, bg_bias_df$y, col = "grey50", pch = 16, cex = 0.3)
points(occs$x, occs$y, col = "red", pch = 16, cex = 1)
legend("bottomright", legend = c("Occurrences", "Background"),
       col = c("red", "grey50"), pch = 16)

