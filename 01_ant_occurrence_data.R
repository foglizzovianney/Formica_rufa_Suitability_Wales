# Script:      02_occurrence_data.R
# Description: Clean and filter Formica rufa occurrence records from NBN Atlas, BWARS, and Charlotte's data

setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(sf)
library(tidyverse)


# NBN DATA 

nbn_raw <- read.csv("data/raw/Formica_rufa_NBNatlas.csv")

nbn_clean <- nbn_raw %>%
  filter(State.Province == "Wales") %>%
  filter(Year >= 2010) %>%
  filter(Coordinate.Uncertainty..m. <= 100)

write.csv(nbn_clean, "data/processed/Formica_rufa_Wales_clean.csv", row.names = FALSE)

st_write(wales, "spatial/boundaries/wales_boundary.shp")


# BWARS DATA 

bwars_raw <- read.csv("data/raw/BWARS_Welsh_data.csv")

bwars_sf <- st_as_sf(bwars_raw,
                     coords = c("easting", "northing"),
                     crs = 27700)

wales_27700 <- st_transform(wales, crs = 27700)

bwars_wales <- st_intersection(bwars_sf, wales_27700)

bwars_clean <- bwars_wales %>%
  filter(accuracy <= 100) %>%
  filter(upper_year >= 2010)

bwars_final <- bwars_clean %>%
  mutate(source = "BWARS",
         year = upper_year) %>%
  dplyr::select(source, year, geometry)

nbn_final <- nbn_clean %>%
  st_as_sf(coords = c("Longitude..WGS84.", "Latitude..WGS84."), crs = 4326) %>%
  st_transform(crs = 27700) %>%
  mutate(source = "NBN",
         year = Year) %>%
  dplyr::select(source, year, geometry)


# COMBINE NBN + BWARS 

combined_records <- rbind(bwars_final, nbn_final)

combined_unique <- combined_records %>%
  distinct(geometry, .keep_all = TRUE)

coords <- st_coordinates(combined_unique)

thinned <- thin(
  loc.data = data.frame(
    Longitude = coords[,1],
    Latitude  = coords[,2],
    Species   = "Formica rufa"),
  lat.col  = "Latitude",
  long.col = "Longitude",
  spec.col = "Species",
  thin.par = 0.1,
  reps     = 10,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  verbose     = FALSE
)

thinned_df <- thinned[[1]]

final_records <- st_as_sf(thinned_df,
                          coords = c("Longitude", "Latitude"),
                          crs = 27700)

st_write(final_records,
         "data/processed/Formica_rufa_Wales_final.gpkg",
         delete_dsn = TRUE)

plot(st_geometry(wales_27700), col = "lightgreen", border = "grey40",
     main = "Final Formica rufa records - Wales")
plot(st_geometry(final_records), add = TRUE, col = "red", pch = 16, cex = 0.8)


# CHARLOTTE DATA 

final_records <- st_read("data/processed/Formica_rufa_Wales_final.gpkg") %>%
  rename(geometry = geom) %>%
  st_set_geometry("geometry") %>%
  mutate(source = "NBN_BWARS", year = NA)

charlotte_raw <- read.csv("data/raw/Rufa_data_Charlotte.csv")

charlotte_raw$Coordinates <- trimws(gsub("\n", "", charlotte_raw$Coordinates))

charlotte_coords <- charlotte_raw %>%
  separate(Coordinates, into = c("Latitude", "Longitude"), sep = ",") %>%
  mutate(Latitude  = as.numeric(trimws(Latitude)),
         Longitude = as.numeric(trimws(Longitude)))

charlotte_sf <- st_as_sf(charlotte_coords,
                         coords = c("Longitude", "Latitude"),
                         crs = 4326)

charlotte_27700 <- st_transform(charlotte_sf, crs = 27700)

plot(st_geometry(wales_27700), col = "lightgreen", border = "grey40",
     main = "Charlotte's F. rufa nest locations (pre-thinning)")
plot(st_geometry(charlotte_27700), add = TRUE, col = "red", pch = 16, cex = 0.8)

charlotte_27700_clean <- charlotte_27700 %>%
  mutate(source = "Charlotte", year = 2025) %>%
  dplyr::select(source, year, geometry)

keep <- rep(TRUE, nrow(charlotte_27700_clean))

for (i in 1:(nrow(charlotte_27700_clean) - 1)) {
  if (keep[i]) {
    dists     <- st_distance(charlotte_27700_clean[i,], charlotte_27700_clean[keep,])
    too_close <- which(dists < units::set_units(250, "m") & dists > units::set_units(0, "m"))
    if (length(too_close) > 0) {
      keep[which(keep)[too_close]] <- FALSE
    }
  }
}

charlotte_thinned_manual <- charlotte_27700_clean[keep,]


# COMBINE ALL RECORDS

all_records <- rbind(final_records, charlotte_thinned_manual)

plot(st_geometry(wales_27700), col = "lightgreen", border = "grey40",
     main = "All F. rufa records - Wales (combined)")
plot(st_geometry(all_records[all_records$source == "NBN_BWARS",]),
     add = TRUE, col = "red", pch = 16, cex = 0.8)
plot(st_geometry(all_records[all_records$source == "Charlotte",]),
     add = TRUE, col = "blue", pch = 16, cex = 0.8)
legend("bottomright", legend = c("NBN/BWARS", "Charlotte"),
       col = c("red", "blue"), pch = 16)

st_write(all_records,
         "data/processed/Formica_rufa_Wales_all_records.gpkg",
         delete_dsn = TRUE)

