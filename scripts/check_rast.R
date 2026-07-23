## function to check rasters against target and reproject/resample to match

#-------------------------------------------------------------------------------

library(terra)

#-------------------------------------------------------------------------------

check_rast <- function(r, target, resample_method = "bilinear") {
  
  if (!inherits(r, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }
  print(target$extent)
  
  template <- terra::rast(
    nrows = target$nrow,
    ncols = target$ncol,
    xmin = target$extent["xmin"],
    xmax = target$extent["xmax"],
    ymin = target$extent["ymin"],
    ymax = target$extent["ymax"],
    crs = target$crs
  )
  
  # reproject if needed
  if (!identical(terra::crs(r, proj=TRUE), target$crs)) {
    r <- terra::project(r, target$crs)
  }
  
  ext_r <- terra::ext(r)
  ext_r_vect <- c(xmin=ext_r[1], xmax=ext_r[2], ymin=ext_r[3], ymax=ext_r[4])
  
  ext_t <- terra::ext(template)
  ext_t_vect <- c(xmin=ext_t[1], xmax=ext_t[2], ymin=ext_t[3], ymax=ext_t[4])
  
  
  needs_resample <- (
    !all(abs(terra::res(r) - target$res) < 1e-6) | 
      !all(abs(ext_r_vect - ext_t_vect) < 1e-6)
  )
  
  if (needs_resample) {
    message("Resampling raster to match target resolution and extent.")
    r <- terra::resample(r, template, method = resample_method)
  }
  
  return(r)
}