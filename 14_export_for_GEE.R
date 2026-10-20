# Script:      14_export_for_GEE.R
# Description: Prepare raster files for Google Earth Engine


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)


# CREATE OUTPUT DIRECTORIES 

dir.create("outputs/for_GEE",    showWarnings = FALSE)
dir.create("outputs/gee_upload", showWarnings = FALSE, recursive = TRUE)


# SUITABILITY MAPS 

suit_current <- rast("outputs/suitability_wales_current.tif")
writeRaster(project(suit_current, "EPSG:4326"),
            "outputs/for_GEE/Formica_rufa_Wales_suitability_current_4326.tif",
            overwrite = TRUE)

suit_future <- rast("outputs/suitability_wales_future_rcp85.tif")
writeRaster(project(suit_future, "EPSG:4326"),
            "outputs/for_GEE/Formica_rufa_Wales_suitability_future_rcp85_4326.tif",
            overwrite = TRUE)


# PATCH RASTERS

patches <- list(
  ON_patches                  = rast("outputs/ON_patches.tif"),
  SNO_patches                 = rast("outputs/SNO_patches.tif"),
  dispersal_patches           = rast("outputs/dispersal_patches.tif"),
  reintroduction_patches_50ha = rast("outputs/reintroduction_patches_50ha.tif")
)

for (nm in names(patches)) {
  r <- project(patches[[nm]], "EPSG:4326", method = "near")
  writeRaster(r,
              paste0("outputs/gee_upload/", nm, ".tif"),
              overwrite = TRUE)
}
