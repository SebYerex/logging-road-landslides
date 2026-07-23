## SURFACE DETRENDING TRANSLATED TO PYTHON ##
## translated for MATLAB script: https://github.com/HydrogeomorphologyTools/DTM-Inpainting-surface-roughness-restitution/tree/master
## which applies methods used in:  https://doi.org/10.1002/esp.4739

############################
### INSTALL DEPENDANCIES ###
############################

# install.packages("reticulate") # if required
# install.packages("terra") # if required

library(reticulate)
library(terra)

###########################
### DEFINE THE FUNCTION ###
###########################
# V2
# Python function
py_run_string("
import numpy as np
from scipy.ndimage import convolve, generic_filter

def detrend_surface(Z, ker_size=9, big_window=41, n_samples=1, mask=None):
    '''
    Compute residual topography from Z and sample local variability.

    Parameters:
      Z : 2D numpy array (DEM)
      ker_size : int, size of smoothing kernel
      big_window : int, size of moving window for sampling residuals
      n_samples : int, number of samples to take in each window
      mask : 2D array, 1 for void areas to add residuals, 0 elsewhere (optional)
    '''
    Z = np.array(Z, dtype=np.float64)

    # Step 1: Smooth DEM with mean filter
    kernel = np.ones((ker_size, ker_size)) / (ker_size ** 2)
    Z_m = convolve(Z, kernel, mode='reflect')

    # Step 2: Residual topography
    Z_res = Z - Z_m

    # Step 3: Sampling function
    def sample_non_nan(values):
        valid = values[~np.isnan(values)]
        if len(valid) > 0:
            return np.median(np.random.choice(valid, size=min(len(valid), n_samples), replace=False))
        else:
            return np.nan

    # Step 4: Apply moving window sampling
    Z_window_random = generic_filter(Z_res, sample_non_nan, size=(big_window, big_window), mode='reflect')

    # Step 5: Apply mask if provided
    if mask is not None:
        mask = np.array(mask, dtype=float)
        Z_window_random = Z_window_random * mask

    return Z_window_random
")

# R wrapper
detrend_surface <- function(Z, ker_size = 9L, big_window = 41L, n_samples = 1L, mask = NULL) {
  
  m <- terra::as.matrix(Z, wide = TRUE)
  
  mask_m <- if (!is.null(mask)) terra::as.matrix(mask, wide = TRUE) else NULL
  
  results <- py$detrend_surface(m, ker_size = as.integer(ker_size),
                                big_window = as.integer(big_window),
                                n_samples = as.integer(n_samples),
                                mask = mask_m)
  
  r <- terra::rast(results, ext = terra::ext(Z), crs = terra::crs(Z))
  
  return(r)
}



##################
### SMALL TEST ###
##################

# library(terra)
# # Create raster
# r <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin=0, ymax=100)
# values(r) <- runif(ncell(r), min=1, max=100)
# 
# # Introduce a missing data region
# r[40:60, 40:60] <- NA
# #plot(r, main = "Raster with NA's")
# 
# # Make a mask layer for later
# A <- ifel(is.na(r), 0, r)
# # Fill in NA's with simple mean
# r <- ifel(is.na(r), mean(na.omit(values(r))), r)
# #plot(r, main = "Raster with filled NA's")
# 
# # Apply the detrend to get a noise layer
# noise <- detrend_surface(r)
# 
# # Where there were originally NA's add noise to the filled values
# d <- ifel(A < 1, (r + noise), r)
# #plot(d, main = "Raster with detrending")
