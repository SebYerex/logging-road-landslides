# --------------------------------------------------------------------
# terra-based approximation of REGOLITH's gauss_approx smoothing, based on: https://doi.org/10.5066/P9U2RDWJ

library(terra)

gauss_filter <- function(r, npoints = 3, n_iter = 4, naflag = NA) {
  
  # Make sure NA values are recognized
  if (!is.na(naflag)) {
    r[r == naflag] <- NA
  }
  
  # Define moving-average kernel
  w <- matrix(1, npoints, npoints)
  
  # Build 1D kernels
  w_y <- matrix(1, npoints, 1) / npoints
  w_x <- matrix(1, 1, npoints) / npoints
  
  # Run iterative smoothing
  r_temp <- r
  
  for (i in 1:n_iter) {
    # Moving average in Y (rows)
    r_temp <- focal(r_temp, w = w_y, fun = mean, na.policy = "omit", fillvalue = NA, na.rm = TRUE, expand = TRUE, filename = "", wopt = list())
    # Moving average in X (cols)
    r_temp <- focal(r_temp, w = w_x, fun = mean, na.policy = "omit", fillvalue = NA, na.rm = TRUE, expand = TRUE, filename = "", wopt = list())
  }
  
  return(r_temp)
}