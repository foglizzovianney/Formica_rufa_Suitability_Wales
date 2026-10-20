# Script: 04_sampling_effort.R
# Description: Create sampling effort raster for Wales from full BWARS dataset


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(sf)

# LOAD DATA 
wales_27700 <- vect("spatial/boundaries/wales_boundary.shp") %>%
  project("EPSG:27700")

wales_sf <- st_read("spatial/boundaries/wales_boundary.shp") %>%
  st_transform(27700)

bwars_all <- read.csv("data/raw/BWARS_all_species.csv")

# FILTER TO WALES AND 2010+
bwars_wales_effort <- bwars_all %>%
  filter(easting >= 146000 & easting <= 356000,
         northing >= 165000 & northing <= 395000) %>%
  mutate(upper_date = as.Date(upper_date, format = "%d/%m/%Y"),
         year = as.numeric(format(upper_date, "%Y"))) %>%
  filter(year >= 2010, !is.na(upper_date)) %>%
  mutate(grid_10km = case_when(
    nchar(osgr_to_gridref) >= 6 ~ paste0(substr(osgr_to_gridref, 1, 2),
                                         substr(osgr_to_gridref, 3, 3),
                                         substr(osgr_to_gridref, 5, 5)),
    nchar(osgr_to_gridref) == 4 ~ osgr_to_gridref,
    TRUE ~ NA_character_))

# COUNT UNIQUE DAYS PER 10KM SQUARE 
effort_summary <- bwars_wales_effort %>%
  group_by(grid_10km) %>%
  summarise(days_sampled = n_distinct(upper_date)) %>%
  filter(!is.na(grid_10km))

# RASTERISE 
bwars_effort_sf <- bwars_wales_effort %>%
  filter(!is.na(easting), !is.na(northing)) %>%
  left_join(effort_summary, by = "grid_10km") %>%
  st_as_sf(coords = c("easting", "northing"), crs = 27700)

effort_rast <- rast(ext = ext(wales_27700),
                    crs = "EPSG:27700",
                    resolution = 10000)

effort_rast_wales <- rasterize(vect(bwars_effort_sf), effort_rast,
                               field = "days_sampled", fun = "max") %>%
  crop(wales_27700) %>%
  mask(wales_27700)

effort_rast_wales[is.na(effort_rast_wales)] <- 0
effort_rast_log <- log(effort_rast_wales + 1)

# SAVE RASTER
writeRaster(effort_rast_log,
            "spatial/vom/effort_rast_10km_wales.tif",
            overwrite = TRUE)

# PLOT 
bbox <- st_as_sfc(st_bbox(st_buffer(wales_sf, 50000)))
outside_wales <- st_difference(bbox, wales_sf)

effort_df <- as.data.frame(effort_rast_log, xy = TRUE)
names(effort_df) <- c("x", "y", "effort")
effort_df <- effort_df %>% filter(!is.na(effort))

ggplot() +
  geom_raster(data = effort_df %>% mutate(effort = replace_na(effort, 0)),
              aes(x = x, y = y, fill = effort)) +
  geom_sf(data = outside_wales, fill = "white", colour = NA) +
  geom_sf(data = wales_sf, fill = NA, colour = "grey40", linewidth = 0.8) +
  scale_fill_viridis_c(name = "log(Visits)", option = "viridis") +
  coord_sf(xlim = c(146000, 356000), ylim = c(165000, 395000), crs = 27700) +
  labs(title = "Sampling effort - Number of volunteer visits in Wales", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.ontop = TRUE,
        panel.background = element_rect(fill = NA, colour = NA),
        panel.grid = element_line(colour = "lightgrey", linewidth = 0.3))

ggsave("outputs/figures/sampling_effort_wales.png", width = 8, height = 10, dpi = 300)


