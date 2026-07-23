## Fix soilgrids layers for soil texture classification

#-------------------------------------------------------------------------------

library(terra)

#------------------------------------------------------------------------------

fix_soil <- function(soil) {
  if (!inherits(soil, "SpatRaster")) {
    stop("Input must be a SpatRaster.")
  }
  
  # store original name(s)
  layer_names <- names(soil)
  
  # determine scale factor based on variable name
  scale <- sapply(layer_names, function(name) {
    name <- tolower(name)
    if (grepl("sand|silt|clay", name)) {
      return(1000) # convert g/1000g to fraction
    } else if (grepl("ph", name)) {
      return(10) # convert to pH units
    } else if (grepl("cec", name)) {
      return(10) # convert to cmol(c)/kg
    } else if (grepl("bdod", name)) {
      return(100) # convert to kg/dm^3
    } else {
      warning(paste("Unknown variable in layer name:", name, "- no conversion applied."))
      return(1) # no conversion
    }
  })
  
  # apply conversion to each layer
  if (length(scale) == 1) {
    soil_fixed <- soil / scale
  } else {
    # if multiple layers with different scaling factors, use 'lapp'
    soil_fixed <- terra::lapp(soil, fun = function(...) {
      vals <- list(...)
      sapply(seq_along(vals), function(i) vals[[i]] / scale [i])
    })
  }
  
  return(soil_fixed)
}
