library(terra)

#-------------------------------------------------------------------------------
# friction angle function
## based on: http://refhub.elsevier.com/S0167-1987(18)31213-3/sbref0020
int_friction <- function(sand = NULL, 
                         silt = NULL, 
                         clay = NULL,
                         fineSand = NULL, 
                         veryFineSand = NULL,
                         method = c("subfraction", "GMD"),
                         D_sand = 0.5, 
                         D_silt = 0.026,
                         D_clay = 0.001) {
  method <- match.arg(method)
  
  #--- check data ---
  if (method == "subfraction") {
    if (any(sapply(list(fineSand, veryFineSand), is.null)) || 
        !all(sapply(list(fineSand, veryFineSand), inherits, what = "SpatRaster"))) {
      stop("For method = 'subfraction', 'fineSand' and 'veryFineSand' must be provided as SpatRaster objects.")
    }
    # ---- subfraction equation ----
    FA <- 1.40 + 0.0001 * (fineSand^2) + 0.0001 * (veryFineSand^2)
    
    # convert to linear scaling and round
    FA <- round(10^FA, 3)
    return(FA)
  }
  
  else if (method == "GMD") {
    if (any(sapply(list(sand, silt, clay), is.null)) || 
        !all(sapply(list(sand, silt, clay), inherits, what = "SpatRaster"))) {
      stop("For method = 'GMD', 'sand', 'silt', and 'clay' must be provided as SpatRaster objects.")
    }
    #---- GMD equation ----
    # based on: https://doi.org/10.1155/2022/2122554
    # & http://refhub.elsevier.com/S0167-8809(20)30084-0/sbref0445
    GMD <- sand * log10(D_sand) + silt * log10(D_silt) + clay * log10(D_clay)
    
    # convert to linear
    GMD <- 10^GMD
    
    # calculate friction angle
    FA <- 1.43 + 1.23 * GMD
    
    # convert to linear and round
    FA <- round(10^FA, 3)
    return(FA)
  }
}