## Function to get information from target raster (dem) and store in a list

#-------------------------------------------------------------------------------

library(terra)
# dem <- rast("E:/MSc/LiDAR/ChipmunkCreek/ChipmunkCreek_DEM_Clip.tif")

#-------------------------------------------------------------------------------

target_rast <- function(dem) {
  
  # check that the class is SpatRaster
  if (!inherits(dem, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }
  
  # get extent
  e_vect <- c(
    xmin = terra::xmin(dem),
    xmax = terra::xmax(dem),
    ymin = terra::ymin(dem),
    ymax = terra::ymax(dem)
  )
  
  list(
    crs = terra::crs(dem, proj = TRUE),
    extent = e_vect,
    res = terra::res(dem),
    nrow = terra::nrow(dem),
    ncol = terra::ncol(dem)
  )
  
}