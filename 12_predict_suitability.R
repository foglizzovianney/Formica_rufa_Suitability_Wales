# Script: 12_predict_suitability.R
# Description: Predict F. rufa habitat suitability across Wales


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(dismo)
library(rJava)

# LOAD DATA 
predictors <- rast("spatial/vom/predictor_stack_wales.tif")
best_model <- readRDS("outputs/best_model_full.RDS")
wales_27700 <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")


# PREDICT CURRENT SUITABILITY 

suitability_current <- predict(predictors, best_model,
                               type = "cloglog",
                               wopt = list(steps = 250),
                               na.rm = TRUE) %>%
  mask(wales_27700)

writeRaster(suitability_current,
            "outputs/suitability_wales_current.tif",
            overwrite = TRUE)


# PREDICT FUTURE SUITABILITY RCP8.5 

# Load future climate
clim_future <- rast("spatial/vom/clim_temp_stack_rcp85.tif") %>%
  resample(predictors)

rain_coldest_85 <- rast("spatial/vom/rain_coldest_quarter_rcp85.tif") %>%
  resample(predictors)
dry_spell_85 <- rast("spatial/vom/dry_spell_rcp85.tif") %>%
  resample(predictors)

# Log transform rainfall
rain_coldest_85 <- log(rain_coldest_85 + 1)

# Build future predictor stack - matching current variable names
predictors_future <- predictors
predictors_future$current_perc09_temp_hottest <- clim_future$rcp85_perc09_temp_hottest
predictors_future$current_temp_seasonality <- clim_future$rcp85_temp_seasonality
predictors_future$perc09_total_rain_coldest <- rain_coldest_85
predictors_future$dry_duration_09perc <- dry_spell_85

suitability_future <- predict(predictors_future, best_model,
                              type = "cloglog",
                              wopt = list(steps = 250),
                              na.rm = TRUE) %>%
  mask(wales_27700)

writeRaster(suitability_future,
            "outputs/suitability_wales_future_rcp85.tif",
            overwrite = TRUE)

# PLOT
par(mfrow = c(1, 2))

plot(suitability_current,
     main = "Current suitability (2010-2013)",
     col = hcl.colors(100, "RdYlGn"),
     axes = FALSE)
plot(wales_27700, add = TRUE, border = "grey40")

plot(suitability_future,
     main = "Future suitability RCP8.5 (2070-2073)",
     col = hcl.colors(100, "RdYlGn"),
     axes = FALSE)
plot(wales_27700, add = TRUE, border = "grey40")

# SAVE FIGURES 
png("outputs/figures/suitability_current.png",
    width = 800, height = 1000, res = 150)
plot(suitability_current,
     main = "F. rufa habitat suitability - Wales (2010-2013)",
     col = hcl.colors(100, "RdYlGn"),
     axes = FALSE)
plot(wales_27700, add = TRUE, border = "grey40")
dev.off()

png("outputs/figures/suitability_future_rcp85.png",
    width = 800, height = 1000, res = 150)
plot(suitability_future,
     main = "F. rufa habitat suitability - Wales (RCP8.5 2070-2073)",
     col = hcl.colors(100, "RdYlGn"),
     axes = FALSE)
plot(wales_27700, add = TRUE, border = "grey40")
dev.off()


