## function to classify Ksat and
# ROSETTA based on:https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-research-unit/docs/model/rosetta-class-average-hydraulic-parameters/
# eu based on: https://www.futurewater.nl/wp-content/uploads/2020/10/HiHydroSoil-v2.0-High-Resolution-Soil-Maps-of-Global-Hydraulic-Properties.pdf
# & https://doi.org/10.1111/ejss.12192

#-------------------------------------------------------------------------------

library(terra)

#-------------------------------------------------------------------------------

# ---- having multiple methods ---

transmissivity <- function(texture = NULL, 
                           thickness = NULL, 
                           return = c("transmissivity", "ksat"), 
                           method = c("rosetta", "eu"), 
                           clay = NULL,
                           sand = NULL,
                           silt = NULL, 
                           ph = NULL,
                           cec = NULL,
                           is_topsoil = FALSE) {
  return <- match.arg(return)
  method <- match.arg(method)
  
  if (method == "rosetta") {
    if (is.null(texture) || !inherits(texture, "SpatRaster")) {
      stop("For method = 'rosetta', 'texture' must be provided as a SpatRaster.")
    }
    
    # Assign Ksat (m/hr) based on USDA texture codes
    ksat <- terra::lapp(texture, fun = function(texture) {
      k <- rep(NA_real_, length(texture))
      k[texture == 1]  <- 0.2678 # sand
      k[texture == 2]  <- 0.0438 # loamy sand
      k[texture == 3]  <- 0.0159 # sandy loam
      k[texture == 4]  <- 0.0050 # loam
      k[texture == 5]  <- 0.0076 # silt loam
      k[texture == 6]  <- 0.0183 # silt
      k[texture == 7]  <- 0.0055 # sandy clay loam
      k[texture == 8]  <- 0.0034 # clay loam
      k[texture == 9]  <- 0.0046 # silty clay loam
      k[texture == 10] <- 0.0048 # sandy clay
      k[texture == 11] <- 0.0040 # silty clay
      k[texture == 12] <- 0.0062 # clay
      return(k)
    })
    names(ksat) <- "ksat_m_per_hr"
    
  } else if (method == "eu") {
    
    # check that all of the rasters are valid
    if (any(sapply(list(clay, sand, ph, cec, silt), is.null))) {
      stop("For method = 'eu', clay, sand, silt, ph, and cec must be provided as SpatRaster objects.")
    }
    if (!all(sapply(list(clay, sand, ph, cec, silt), inherits, what = "SpatRaster"))) {
      stop("All inputs for method = 'eu' must be SpatRasetr objects.")
    }
    if (!all(terra::compareGeom(clay, sand, ph, cec, silt, stopOnError = FALSE))) {
      stop("All rasters (clay, sand, silt, ph, cec) must have matching extent, resolution, and CRS.")
    }
    
    # apply EU-SoilHydroGrids PTF
    ts_val <- ifelse(is_topsoil, 1, 0)
    
    log_ksat <- 0.40220 + 0.26122 * ph + 0.44565 * ts_val - 0.02329 * clay - 0.01265 * silt - 0.01038 * cec
    
    ksat <- (10^log_ksat) / 100 # cm/day -> m/day
    ksat <- ksat / 24 # m/day -> m/hr
    names(ksat) <- "ksat_m_per_hr"
  }
  
  if (return == "ksat") {
    return(ksat)
  }
  
  # transmissivity 
  if (is.null(thickness) || !inherits(thickness, "SpatRaster")) {
    stop("thickness must be provided as a SpatRaster to calculate transmissivity.")
  }
  
  if (!terra::compareGeom(ksat, thickness, stopOnError = FALSE)) {
    stop("ksat and thickness rasters must have matching extent, resolution, and CRS.")
  }
  
  trans <- ksat * thickness
  names(trans) <- "transmissivity_m2_per_hr"
  
  return(trans)
}