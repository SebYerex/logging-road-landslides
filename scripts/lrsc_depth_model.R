## LRSC soil depth model translated to R from REGOLITH (fortran), based on: https://doi.org/10.5066/P9U2RDWJ

# ============================================================
# Libraries
# ============================================================
library(terra)


## LRSC PARAMETER NOTES ##

# - C0 represents background soil thickness, the y-intercept in the linear regression

# - C1 represents the sensitivity of soil thickness to slope curvature

# - C2 represents the control of slope angle on soil thickness, larger C2 = more thinning due to increasing slope angle


lrsc_depth_model <- function(
    dem,
    ca,
    C0, C1, C2,
    theta_c,
    depth_min, depth_max,
    chan_thresh, chan_depth,
    smooth_topo = FALSE,
    smooth_soil = FALSE
) {
  stopifnot(inherits(dem, "SpatRaster"), inherits(ca, "SpatRaster"))
  if (!compareGeom(dem, ca)) stop("DEM and CA rasters must match in extent and resolution.")
  
  message("Computing slope (radians) and Laplacian curvature ...")
  
  if (isTRUE(smooth_topo)) {
    message("Applying topographic smoothing (3x3 mean kernel)...")
    kernel <- matrix(1, nrow = 3, ncol = 3)
    kernel <- kernel / sum(kernel)
    dem <- focal(dem, w = kernel, fun = mean, na.policy = "omit", pad = TRUE, padValue = NA)
  }
  
  slope_rad <- terrain(dem, v = "slope", unit = "radians", neighbors = 8)
  theta_c_rad <- theta_c * pi / 180
  
  res_m <- res(dem)[1]
  laplace_kernel <- matrix(c(0,1,0,1,-4,1,0,1,0), nrow = 3, byrow = TRUE)
  laplacian <- focal(dem, w = laplace_kernel, fun = sum, na.policy = "omit", pad = TRUE, padValue = NA) / (res_m^2)
  
  message("Derived slope and Laplacian curvature successfully.")
  
  soil_depth <- app(c(slope_rad, laplacian, ca), fun = function(x) {
    slope_angle <- x[1]; lap <- x[2]; area <- x[3]
    if (is.na(slope_angle) || is.na(lap) || is.na(area)) return(NA)
    
    mag_del_z <- tan(slope_angle)
    sc <- tan(theta_c_rad)
    
    if (mag_del_z > sc) {
      depth <- 0
    } else {
      h1 <- C2 * (sc - mag_del_z)
      h2 <- C0 + C1 * lap
      depth <- h1 + h2
      depth <- max(min(depth, depth_max), depth_min)
    }
    
    if (area > chan_thresh && slope_angle > 0.1 * theta_c_rad) {
      depth <- min(depth, chan_depth)
    }
    
    return(depth)
  })
  
  if (isTRUE(smooth_soil)) {
    message("Applying soil depth smoothing (3x3 mean kernel)...")
    kernel <- matrix(1, nrow = 3, ncol = 3)
    kernel <- kernel / sum(kernel)
    soil_depth <- focal(soil_depth, w = kernel, fun = mean, na.policy = "omit", pad = TRUE, padValue = NA)
  }
  
  names(soil_depth) <- "soil_depth"
  return(soil_depth)
}