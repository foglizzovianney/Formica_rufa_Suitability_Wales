# Script:      09_covariate_prep.R
# Description: Prepare covariate stack, check collinearity, select final variables


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(tidyterra)
library(data.table)
library(usdm)
library(corrplot)


# LOAD COVARIATES

forest_covariates <- rast("spatial/vom/forest_PCA_30m.tif")
forest_mask_buff  <- forest_covariates$forest_mask_buff

climate_current <- rast("spatial/vom/clim_temp_stack_rcp45.tif") %>%
  resample(forest_covariates)

rain_coldest <- rast("spatial/vom/rain_coldest_quarter_rcp45.tif") %>%
  resample(forest_covariates)

rain_hottest <- rast("spatial/vom/rain_hottest_quarter_rcp45.tif") %>%
  resample(forest_covariates)

dry_spell <- rast("spatial/vom/dry_spell_rcp45.tif") %>%
  resample(forest_covariates)

rain_coldest <- log(rain_coldest + 1)
rain_hottest <- log(rain_hottest + 1)

names(rain_coldest) <- "perc09_total_rain_coldest"
names(rain_hottest) <- "perc09_total_rain_hottest"
names(dry_spell)    <- "dry_duration_09perc"

topo_covariates <- rast("spatial/vom/topo_stack_30m_27700.tif") %>%
  resample(forest_covariates)

topo_covariates$slope <- log(topo_covariates$slope + 1)

distance_forest <- rast("spatial/vom/distance_forest_30m.tif") %>%
  resample(forest_covariates)

distance_forest        <- log(distance_forest + 0.001)
names(distance_forest) <- "distance_forest"

sbio_files <- list.files("spatial/vom",
                         pattern    = "SBIO.*wales_27700.tif",
                         full.names = TRUE)

soil_temp <- rast(sbio_files) %>%
  resample(forest_covariates)


# STACK ALL COVARIATES 

covariates_stack <- c(forest_covariates,
                      climate_current,
                      rain_coldest,
                      rain_hottest,
                      dry_spell,
                      topo_covariates,
                      distance_forest,
                      soil_temp) %>%
  terra::mask(forest_mask_buff)

writeRaster(covariates_stack,
            "spatial/vom/full_covariates_stack.tif",
            overwrite = TRUE)

effort_masked_30m <- rast("spatial/vom/effort_rast_10km_wales.tif") %>%
  resample(forest_mask_buff) %>%
  mask(forest_mask_buff)

writeRaster(effort_masked_30m,
            "spatial/vom/effort_masked_30m.tif",
            overwrite = TRUE)


# CHECK COLLINEARITY 

keep <- which(values(forest_mask_buff) == 1)

set.seed(1)
cells <- sample(keep, min(1000000, length(keep)))

sample_df <- terra::extract(covariates_stack, cells) %>%
  dplyr::select(-ID, -forest_mask_buff) %>%
  as.data.table() %>%
  drop_na()

usdm::vif(sample_df)
usdm::vifstep(sample_df, th = 5)
usdm::vifstep(sample_df, th = 10)

cor_matrix <- cor(sample_df, use = "complete.obs")

png("outputs/figures/covariate_correlation.png",
    width  = 2000,
    height = 2000,
    res    = 150)

corrplot::corrplot(cor_matrix,
                   method      = "circle",
                   type        = "upper",
                   order       = "original",
                   addCoef.col = "black",
                   number.cex  = 0.65,
                   tl.cex      = 1,
                   tl.col      = "black",
                   tl.srt      = 45,
                   diag        = FALSE,
                   col         = colorRampPalette(c("red", "white", "blue"))(200),
                   addgrid.col = "grey")

dev.off()

