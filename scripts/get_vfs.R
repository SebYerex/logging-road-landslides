#-------------------------------------------------------------------------------

library(terra)

#-------------------------------------------------------------------------------
# very fine sand function
## based on: http://refhub.elsevier.com/S0048-9697(20)35138-X/rf0035

get_vfs <- function(texture, sand, silt, clay) {
  
  # --- check inputs ---
  if (!inherits(texture, "SpatRaster") ||
      !inherits(sand, "SpatRaster") ||
      !inherits(silt, "SpatRaster") ||
      !inherits(clay, "SpatRaster")) {
    stop("All inputs must be SpatRaster objects.")
  }
  
  # --- model functions ---
  ESDAC_model <- function(sand) (1/5)*sand
  RUSLE2_model <- function(sand) (0.74 - 0.62 * sand) * sand
  SB_model <- function(sand) {
    sand <- pmin(pmax(sand, 0.0001), 0.9999)
    phi_inv <- qnorm(1 - sand)
    phi_arg <- 0.698810 + 0.812098 * phi_inv
    pnorm(phi_arg) - 1 + sand
  }
  
  # --- core apply function (vectorized over entire rasters) ---
  model_vfs <- function(texture, sand, silt, clay) {
    n <- length(texture)
    vfs <- rep(NA_real_, n)
    
    na_idx <- is.na(texture) | is.na(sand) | is.na(silt) | is.na(clay)
    
    idx_rusle2 <- texture %in% c(2, 3) & !na_idx
    vfs[idx_rusle2] <- RUSLE2_model(sand[idx_rusle2])
    
    idx_esdac <- texture %in% c(1, 4, 6, 7, 8, 10, 11, 12) & !na_idx
    vfs[idx_esdac] <- ESDAC_model(sand[idx_esdac])
    
    idx_sb <- texture %in% c(5, 9) & !na_idx
    vfs[idx_sb] <- SB_model(sand[idx_sb])
    
    vfs[na_idx] <- NA_real_
    return(vfs)
  }
  
  # Apply over rasters (remove na.rm)
  rast_stack <- c(texture, sand, silt, clay)
  predicted_vfs <- terra::lapp(rast_stack, fun = model_vfs)
  names(predicted_vfs) <- "predicted_vfs"
  
  return(predicted_vfs)
}
