# Script: 05_climate_variables.R
# Description: Process CHESS-SCAPE climate and precipitation data into seasonal variables for current (2010-2013) and future (RCP8.5 2070-2073)


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)

# FUNCTIONS 
calc_quantiles <- function(x, probs = 0.9) {
  quantile(x, probs = probs, na.rm = TRUE)
}

load_temp_stack <- function(climate_file, start_year, end_year, wales_crs) {
  climate <- rast(climate_file)
  dates <- as.Date(time(climate))
  indices <- which(dates >= as.Date(paste0(start_year, "-01-01")) &
                     dates <= as.Date(paste0(end_year, "-12-31")))
  stack <- climate[[indices]]
  values(stack) <- values(stack) - 273.15
  names(stack) <- sprintf("%d%02d", rep(start_year:end_year, each = 12), 1:12)
  stack <- crop(stack, wales_crs) %>% mask(wales_crs)
  return(stack)
}

get_quarters <- function(temp_stack, start_year, n_quarters = 4) {
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

get_quarterly_temperature <- function(qtrs_dates, temp_stack, quantile = 0.9) {
  quarterly_temp_rasters <- list()
  for (i in seq_along(qtrs_dates)) {
    q <- qtrs_dates[[i]]
    year_qtr <- strsplit(q, "-Q")[[1]]
    year <- as.numeric(year_qtr[1])
    qtr <- as.numeric(year_qtr[2])
    q_start <- (qtr - 1) * 3 + 1
    q_end <- q_start + 2
    month_temp_rasters <- list()
    for (k in seq_along(q_start:q_end)) {
      month <- (q_start:q_end)[k]
      month_pattern <- sprintf("%d%02d", year, month)
      month_temp_sub <- temp_stack[[grep(month_pattern, names(temp_stack))]]
      month_temp_rasters[[k]] <- month_temp_sub
    }
    quarterly_temp_rasters[[i]] <- rast(month_temp_rasters) %>% median()
  }
  app(rast(quarterly_temp_rasters), fun = calc_quantiles, probs = quantile)
}

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

calculate_longest_dryspell <- function(year_stack, threshold = 1) {
  dry_day_count <- terra::init(year_stack[[1]], 0) %>% mask(year_stack[[1]])
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

get_quarterly_dryspell <- function(qtrs_dates, rain_files) {
  dryspell_rasters <- list()
  for (i in seq_along(qtrs_dates)) {
    q <- qtrs_dates[[i]]
    year_qtr <- strsplit(q, "-Q")[[1]]
    year <- as.numeric(year_qtr[1])
    qtr <- as.numeric(year_qtr[2])
    q_start <- (qtr - 1) * 3 + 1
    q_end <- q_start + 2
    month_files <- c()
    for (m in q_start:q_end) {
      month_pattern <- sprintf("%04d%02d", year, m)
      month_files <- c(month_files, grep(month_pattern, rain_files, value = TRUE))
    }
    if (length(month_files) > 0) {
      quarter_stack <- rast(month_files)
      values(quarter_stack) <- values(quarter_stack) * 86400
      dryspell_rasters[[i]] <- calculate_longest_dryspell(quarter_stack)
    }
  }
  app(rast(dryspell_rasters), fun = calc_quantiles, probs = 0.9)
}

# SETUP 
wales_27700 <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")

tas_file <- "spatial/vom/chess-scape_rcp45_bias-corrected_01_tas_uk_1km_monthly_19801201-20801130.nc"
tas_file_85 <- "spatial/vom/chess-scape_rcp85_bias-corrected_01_tas_uk_1km_monthly_19801201-20801130.nc"

rain_files_current <- list.files("spatial/precipitation/rcp45_current/pr",
                                 pattern = ".nc", full.names = TRUE)
rain_files_85 <- list.files("spatial/precipitation/rcp85/pr",
                            pattern = ".nc", full.names = TRUE)

# CURRENT CLIMATE (2010-2013) 

wales_crs_current <- project(wales_27700, crs(rast(tas_file)))
temp_current <- load_temp_stack(tas_file, 2010, 2013, wales_crs_current)
quarters_current <- get_quarters(temp_current, 2010)

# Temperature variables
perc01_temp_coldest_current <- get_quarterly_temperature(quarters_current$coldest, temp_current, quantile = 0.1)
perc09_temp_hottest_current <- get_quarterly_temperature(quarters_current$hottest, temp_current, quantile = 0.9)

yearly_temp <- tapp(temp_current, index = rep(2010:2013, each = 12), fun = mean, na.rm = TRUE)
annual_mean_temp_current <- mean(yearly_temp, na.rm = TRUE)

year_sd_list <- list()
year_range_list <- list()
for (year in 2010:2013) {
  year_layers <- temp_current[[grep(paste0("^", year), names(temp_current))]]
  year_sd_list[[as.character(year)]] <- app(year_layers, fun = sd, na.rm = TRUE)
  year_range_list[[as.character(year)]] <- max(year_layers) - min(year_layers)
}
annual_sd_temp_current <- mean(rast(year_sd_list))
temp_seasonality_current <- mean(rast(year_range_list))

# Rainfall variables
wales_crs_rain <- project(wales_27700, crs(rast(rain_files_current[1])))

rain_coldest_current <- get_quarterly_rainfall(quarters_current$coldest, rain_files_current) %>%
  crop(wales_crs_rain) %>% mask(wales_crs_rain)

rain_hottest_current <- get_quarterly_rainfall(quarters_current$hottest, rain_files_current) %>%
  crop(wales_crs_rain) %>% mask(wales_crs_rain)

dry_spell_current <- get_quarterly_dryspell(quarters_current$hottest, rain_files_current) %>%
  crop(wales_crs_rain) %>% mask(wales_crs_rain)


# FUTURE CLIMATE RCP8.5 (2070-2073) 

wales_crs_85 <- project(wales_27700, crs(rast(tas_file_85)))
temp_future_85 <- load_temp_stack(tas_file_85, 2070, 2073, wales_crs_85)
quarters_85 <- get_quarters(temp_future_85, 2070)


# Temperature variables
perc01_temp_coldest_85 <- get_quarterly_temperature(quarters_85$coldest, temp_future_85, quantile = 0.1)
perc09_temp_hottest_85 <- get_quarterly_temperature(quarters_85$hottest, temp_future_85, quantile = 0.9)

yearly_temp_85 <- tapp(temp_future_85, index = rep(2070:2073, each = 12), fun = mean, na.rm = TRUE)
annual_mean_temp_85 <- mean(yearly_temp_85, na.rm = TRUE)

year_sd_list_85 <- list()
year_range_list_85 <- list()
for (year in 2070:2073) {
  year_layers <- temp_future_85[[grep(paste0("^", year), names(temp_future_85))]]
  year_sd_list_85[[as.character(year)]] <- app(year_layers, fun = sd, na.rm = TRUE)
  year_range_list_85[[as.character(year)]] <- max(year_layers) - min(year_layers)
}
annual_sd_temp_85 <- mean(rast(year_sd_list_85))
temp_seasonality_85 <- mean(rast(year_range_list_85))

# Rainfall variables
wales_crs_rain_85 <- project(wales_27700, crs(rast(rain_files_85[1])))

rain_coldest_85 <- get_quarterly_rainfall(quarters_85$coldest, rain_files_85) %>%
  crop(wales_crs_rain_85) %>% mask(wales_crs_rain_85)

rain_hottest_85 <- get_quarterly_rainfall(quarters_85$hottest, rain_files_85) %>%
  crop(wales_crs_rain_85) %>% mask(wales_crs_rain_85)

dry_spell_85 <- get_quarterly_dryspell(quarters_85$hottest, rain_files_85) %>%
  crop(wales_crs_rain_85) %>% mask(wales_crs_rain_85)

# SAVE

# Current
writeRaster(rain_coldest_current, "spatial/vom/rain_coldest_quarter_current.tif", overwrite = TRUE)
writeRaster(rain_hottest_current, "spatial/vom/rain_hottest_quarter_current.tif", overwrite = TRUE)
writeRaster(dry_spell_current,    "spatial/vom/dry_spell_current.tif",            overwrite = TRUE)

clim_stack_current <- c(perc01_temp_coldest_current, perc09_temp_hottest_current,
                        annual_mean_temp_current, annual_sd_temp_current,
                        temp_seasonality_current)
names(clim_stack_current) <- c("current_perc01_temp_coldest", "current_perc09_temp_hottest",
                               "current_annual_mean_temp", "current_annual_sd_temp",
                               "current_temp_seasonality")
writeRaster(clim_stack_current, "spatial/vom/clim_temp_stack_current.tif", overwrite = TRUE)

# Future RCP8.5
writeRaster(rain_coldest_85, "spatial/vom/rain_coldest_quarter_rcp85.tif", overwrite = TRUE)
writeRaster(rain_hottest_85, "spatial/vom/rain_hottest_quarter_rcp85.tif", overwrite = TRUE)
writeRaster(dry_spell_85,    "spatial/vom/dry_spell_rcp85.tif",            overwrite = TRUE)

clim_stack_85 <- c(perc01_temp_coldest_85, perc09_temp_hottest_85,
                   annual_mean_temp_85, annual_sd_temp_85,
                   temp_seasonality_85)
names(clim_stack_85) <- c("rcp85_perc01_temp_coldest", "rcp85_perc09_temp_hottest",
                          "rcp85_annual_mean_temp", "rcp85_annual_sd_temp",
                          "rcp85_temp_seasonality")
writeRaster(clim_stack_85, "spatial/vom/clim_temp_stack_rcp85.tif", overwrite = TRUE)


