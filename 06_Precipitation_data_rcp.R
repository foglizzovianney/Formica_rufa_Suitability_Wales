# Script: 06_climate_variables.R
# Description: Process CHESS-SCAPE climate data into seasonal rainfall variables


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)

# SETUP 
calc_quantiles <- function(x, probs = 0.9) {
  quantile(x, probs = probs, na.rm = TRUE)
}

# Load Wales boundary
wales_27700 <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")

# Parameters
n_quarters <- 4

# Rainfall files
rain_files_45 <- list.files("spatial/precipitation/rcp45/pr", 
                            pattern = ".nc", full.names = TRUE)
rain_files_85 <- list.files("spatial/precipitation/rcp85/pr", 
                            pattern = ".nc", full.names = TRUE)

# FUNCTIONS 
# Load and prepare temperature stack for a given period
load_temp_stack <- function(climate_file, start_year, end_year) {
  climate <- rast(climate_file)
  dates <- as.Date(time(climate))
  indices <- which(dates >= as.Date(paste0(start_year, "-01-01")) & 
                     dates <= as.Date(paste0(end_year, "-12-31")))
  stack <- climate[[indices]]
  values(stack) <- values(stack) - 273.15
  names(stack) <- sprintf("%d%02d", rep(start_year:end_year, each = 12), 1:12)
  wales_crs <- project(wales_27700, crs(stack))
  stack <- crop(stack, wales_crs) %>% mask(wales_crs)
  return(list(stack = stack, wales_crs = wales_crs))
}

# Identify coldest and hottest quarters
get_quarters <- function(temp_stack, start_year, n_quarters) {
  quarters <- ceiling(1:nlyr(temp_stack) / 3)
  temp_stack_qtr <- tapp(temp_stack, quarters, fun = median)
  avg_qtr_temp <- global(temp_stack_qtr, stat = 'median', na.rm = TRUE)
  
  coldest_indices <- order(avg_qtr_temp$mean)[1:n_quarters]
  hottest_indices <- order(avg_qtr_temp$mean, decreasing = TRUE)[1:n_quarters]
  
  to_dates <- function(indices) {
    lapply(indices, function(q) {
      year <- start_year + ((q - 1) %/% 4)
      qtr <- (q - 1) %% 4 + 1
      paste(year, qtr, sep = "-Q")
    })
  }
  
  list(coldest = to_dates(coldest_indices),
       hottest = to_dates(hottest_indices))
}

# Calculate quarterly rainfall
get_quarterly_rainfall <- function(qtrs_dates, rain_files) {
  quarterly_rainfall_rasters <- list()
  
  for (i in seq_along(qtrs_dates)) {
    q <- qtrs_dates[[i]]
    year_qtr <- strsplit(q, "-Q")[[1]]
    year <- as.numeric(year_qtr[1])
    qtr <- as.numeric(year_qtr[2])
    q_start <- (qtr - 1) * 3 + 1
    q_end <- q_start + 2
    q_rainfall_sum <- NULL
    
    for (m in q_start:q_end) {
      month_pattern <- sprintf("%04d%02d", year, m)
      month_rain_files <- grep(month_pattern, rain_files, value = TRUE)
      
      if (length(month_rain_files) > 0) {
        month_rain_stack <- rast(month_rain_files)
        values(month_rain_stack) <- values(month_rain_stack) * 86400
        month_rain_sum <- sum(month_rain_stack, na.rm = TRUE)
        q_rainfall_sum <- if (is.null(q_rainfall_sum)) month_rain_sum else q_rainfall_sum + month_rain_sum
      }
    }
    quarterly_rainfall_rasters[[i]] <- q_rainfall_sum
  }
  app(rast(quarterly_rainfall_rasters), fun = calc_quantiles, probs = 0.9)
}

# PROCESS RCP4.5 (future 2070-2073)
result_45 <- load_temp_stack(
  "spatial/vom/chess-scape_rcp45_bias-corrected_01_tas_uk_1km_monthly_19801201-20801130.nc",
  2070, 2073)

quarters_45 <- get_quarters(result_45$stack, 2070, n_quarters)

rain_coldest_45 <- get_quarterly_rainfall(quarters_45$coldest, rain_files_45) %>%
  crop(result_45$wales_crs) %>% mask(result_45$wales_crs)

rain_hottest_45 <- get_quarterly_rainfall(quarters_45$hottest, rain_files_45) %>%
  crop(result_45$wales_crs) %>% mask(result_45$wales_crs)

# PROCESS RCP8.5 (future 2070-2073) 
result_85 <- load_temp_stack(
  "spatial/vom/chess-scape_rcp85_bias-corrected_01_tas_uk_1km_monthly_19801201-20801130.nc",
  2070, 2073)

quarters_85 <- get_quarters(result_85$stack, 2070, n_quarters)

rain_coldest_85 <- get_quarterly_rainfall(quarters_85$coldest, rain_files_85) %>%
  crop(result_85$wales_crs) %>% mask(result_85$wales_crs)

rain_hottest_85 <- get_quarterly_rainfall(quarters_85$hottest, rain_files_85) %>%
  crop(result_85$wales_crs) %>% mask(result_85$wales_crs)

# SAVE 
writeRaster(rain_coldest_45, "spatial/vom/rain_coldest_quarter_rcp45.tif", overwrite = TRUE)
writeRaster(rain_hottest_45, "spatial/vom/rain_hottest_quarter_rcp45.tif", overwrite = TRUE)
writeRaster(rain_coldest_85, "spatial/vom/rain_coldest_quarter_rcp85.tif", overwrite = TRUE)
writeRaster(rain_hottest_85, "spatial/vom/rain_hottest_quarter_rcp85.tif", overwrite = TRUE)


# PLOT 
par(mfrow = c(2, 2))
plot(rain_coldest_45, main = "Coldest quarter rainfall RCP4.5", col = hcl.colors(20, "Blues"))
plot(rain_hottest_45, main = "Hottest quarter rainfall RCP4.5", col = hcl.colors(20, "Blues"))
plot(rain_coldest_85, main = "Coldest quarter rainfall RCP8.5", col = hcl.colors(20, "Blues"))
plot(rain_hottest_85, main = "Hottest quarter rainfall RCP8.5", col = hcl.colors(20, "Blues"))




# Dry spell function - following Joris
calculate_longest_dryspell <- function(year_stack, threshold = 1) {
  dry_day_count <- terra::init(year_stack[[1]], 0) %>%
    mask(year_stack[[1]])
  
  for (i in 1:nlyr(year_stack)) {
    daily_rain <- year_stack[[i]]
    is_dry <- daily_rain < threshold
    dry_day_count <- ifel(is_dry, dry_day_count + 1, 0)
    
    if (i == 1) {
      max_dry_spell <- dry_day_count
    } else {
      max_dry_spell <- max(max_dry_spell, dry_day_count)
    }
  }
  return(max_dry_spell)
}

# Calculate dry spell for RCP4.5

dryspell_rasters_45 <- list()

for (i in seq_along(quarters_45$hottest)) {
  q <- quarters_45$hottest[[i]]
  year_qtr <- strsplit(q, "-Q")[[1]]
  year <- as.numeric(year_qtr[1])
  qtr <- as.numeric(year_qtr[2])
  
  q_start <- (qtr - 1) * 3 + 1
  q_end <- q_start + 2
  
  # Get all daily files for this quarter
  month_files <- c()
  for (m in q_start:q_end) {
    month_pattern <- sprintf("%04d%02d", year, m)
    month_files <- c(month_files, grep(month_pattern, rain_files_45, value = TRUE))
  }
  
  if (length(month_files) > 0) {
    quarter_stack <- rast(month_files)
    values(quarter_stack) <- values(quarter_stack) * 86400
    dryspell_rasters_45[[i]] <- calculate_longest_dryspell(quarter_stack)
  }
}

dry_duration_45 <- app(rast(dryspell_rasters_45), fun = calc_quantiles, probs = 0.9) %>%
  crop(result_45$wales_crs) %>% mask(result_45$wales_crs)


# Calculate dry spell for RCP8.5

dryspell_rasters_85 <- list()

for (i in seq_along(quarters_85$hottest)) {
  q <- quarters_85$hottest[[i]]
  year_qtr <- strsplit(q, "-Q")[[1]]
  year <- as.numeric(year_qtr[1])
  qtr <- as.numeric(year_qtr[2])
  
  q_start <- (qtr - 1) * 3 + 1
  q_end <- q_start + 2
  
  month_files <- c()
  for (m in q_start:q_end) {
    month_pattern <- sprintf("%04d%02d", year, m)
    month_files <- c(month_files, grep(month_pattern, rain_files_85, value = TRUE))
  }
  
  if (length(month_files) > 0) {
    quarter_stack <- rast(month_files)
    values(quarter_stack) <- values(quarter_stack) * 86400
    dryspell_rasters_85[[i]] <- calculate_longest_dryspell(quarter_stack)
  }
}

dry_duration_85 <- app(rast(dryspell_rasters_85), fun = calc_quantiles, probs = 0.9) %>%
  crop(result_85$wales_crs) %>% mask(result_85$wales_crs)


# Save
writeRaster(dry_duration_45, "spatial/vom/dry_spell_rcp45.tif", overwrite = TRUE)
writeRaster(dry_duration_85, "spatial/vom/dry_spell_rcp85.tif", overwrite = TRUE)

# Plot
par(mfrow = c(1, 2))
plot(dry_duration_45, main = "Dry spell duration RCP4.5", col = hcl.colors(20, "YlOrRd"))
plot(dry_duration_85, main = "Dry spell duration RCP8.5", col = hcl.colors(20, "YlOrRd"))

