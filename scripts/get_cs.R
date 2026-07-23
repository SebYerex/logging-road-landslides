library(terra)

#-------------------------------------------------------------------------------
# coarse sand function
## based on: http://refhub.elsevier.com/S2095-6339(19)30004-8/sref159
## same as ESDAC model from the function above

get_cs <- function(sand) {
  
  #--- check the data type ---
  if (!inherits(sand, "SpatRaster")) {
    stop("Inputs must be SpatRaster objects.")
  }
  
  #--- apply the model ---
  1/5 * sand
}
