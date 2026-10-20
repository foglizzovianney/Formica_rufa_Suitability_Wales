# Script: 08_forest_PCA.R
# Description: Apply PCA to forest structure variables


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(tidyterra)

# LOAD DATA 
cover_vom  <- rast("spatial/vom/cover_VOM_30m_27700.tif")
perc09_vom <- rast("spatial/vom/perc09_VOM_30m_27700.tif")
sd_vom     <- rast("spatial/vom/sd_VOM_30m_27700.tif")
forest_mask_buff <- rast("spatial/vom/forest_mask_70m_buff.tif")

# PREPARE FOREST STACK 
# Square root transform following Joris
forest_stack <- c(cover_vom, 
                  sqrt(perc09_vom), 
                  sqrt(sd_vom))

names(forest_stack) <- c("cover_VOM", 
                         "perc09_height_VOM_sqrt", 
                         "sd_height_VOM_sqrt")

# Mask to forest buffer
forest_stack <- mask(forest_stack, forest_mask_buff)

# RUN PCA
rpc_forest <- terra::prcomp(forest_stack, center = TRUE, scale. = TRUE)
summary(rpc_forest)

# Print loadings
forest_loadings_df <- data.frame(rpc_forest$rotation) %>%
  mutate(PCA_comp = rownames(.)) %>%
  dplyr::select(PCA_comp, everything())

print(forest_loadings_df)

# PREDICT PCA SCORES 
forest_PC_stack <- predict(forest_stack, rpc_forest) %>%
  subset(1:2) %>%
  c(., forest_mask_buff)

names(forest_PC_stack) <- c("forest_PC1", "forest_PC2", "forest_mask_buff")

# SAVE
write.csv(forest_loadings_df, "outputs/tables/forest_PCA_loadings.csv")

writeRaster(forest_PC_stack,
            "spatial/vom/forest_PCA_30m.tif",
            overwrite = TRUE)

# PLOT 
par(mfrow = c(1, 2))
plot(forest_PC_stack$forest_PC1,
     main = "Forest PC1",
     col = hcl.colors(20, "Spectral"),
     axes = FALSE)
plot(forest_PC_stack$forest_PC2,
     main = "Forest PC2", 
     col = hcl.colors(20, "Spectral"),
     axes = FALSE)

