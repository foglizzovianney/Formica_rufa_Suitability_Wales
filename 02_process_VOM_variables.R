# Script: 02_process_VOM_variables.R
# Description: Derive forest structure variables from VOM at 30m resolution


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)

# Load VOM
vom <- rast("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales/wales_ndsm/wales_ndsm_75m_VOM.tif")

# Define tile size - process in 10x10km chunks
tile_size <- 10000  # 10km in metres

# Get extent of VOM
e <- ext(vom)

# Create tile extents
x_starts <- seq(e$xmin, e$xmax, by = tile_size)
y_starts <- seq(e$ymin, e$ymax, by = tile_size)

# Output files
sd_file <- "spatial/vom/sd_VOM_30m_27700.tif"
perc09_file <- "spatial/vom/perc09_VOM_30m_27700.tif"

# Process each tile
sd_tiles <- list()
perc09_tiles <- list()

for(i in seq_along(x_starts)){
  for(j in seq_along(y_starts)){
    
    # Define tile extent
    tile_ext <- ext(x_starts[i], 
                    min(x_starts[i] + tile_size, e$xmax),
                    y_starts[j], 
                    min(y_starts[j] + tile_size, e$ymax))
    
    # Crop VOM to tile
    tile <- crop(vom, tile_ext)
    
    # Skip if tile is empty
    if(all(is.na(values(tile)))) next
    
    # Aggregate to 30m using SD
    sd_tile <- aggregate(tile, fact = 30, fun = "sd")
    
    # Aggregate to 30m using 90th percentile
    perc09_tile <- aggregate(tile, fact = 30, 
                             fun = function(x) quantile(x, 0.9, na.rm = TRUE))
    
    sd_tiles[[length(sd_tiles) + 1]] <- sd_tile
    perc09_tiles[[length(perc09_tiles) + 1]] <- perc09_tile
    
  }
}


# Merge all tiles
sd_merged <- do.call(merge, sd_tiles)
perc09_merged <- do.call(merge, perc09_tiles)

# Save
writeRaster(sd_merged, sd_file, overwrite = TRUE)
writeRaster(perc09_merged, perc09_file, overwrite = TRUE)
