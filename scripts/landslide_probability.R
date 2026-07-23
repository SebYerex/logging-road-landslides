
library(terra)

# ===============================
# Helper functions
# ===============================

# helper functions

## for converting degrees to radians
deg2rad <- function(x) {
  x * (pi / 180)
}

## for computing sin_slope
compute_sin_slope <- function(slope_rad) {
  sin_slope <- sin(slope_rad)
  sin_slope[abs(sin_slope) < 1e-3] <- 1e-3 # to prevent zero denom
  return(sin_slope)
}

## for computing wetness
compute_wetness <- function(R, transmissivity, sca, sin_slope) {
  terra::clamp((R / transmissivity) * (sca / sin_slope), lower = 0, upper = 1)
}


# main FS function
SINMAP_FS <- function(
    slope_rad,
    sin_slope,
    friction_angle,
    bulk_density = 2000,
    g = 9.81,
    density_w = 1000,
    R = 0.0001,
    transmissivity,
    sca,
    cohesion_s,
    cohesion_r,
    soil_depth,
    wetness_override = NULL,
    cohesion_star_override = NULL,
    tan_fa_override = NULL,
    cos_slope_override = NULL,
    clamp_FS = FALSE,
    intermediates = FALSE
) {
  # ---- check inputs ----
  
  # check that all non-default inputs are provided
  if (any(sapply(list(friction_angle,
                      bulk_density,
                      g,
                      density_w,
                      R,
                      transmissivity,
                      sca,
                      cohesion_s,
                      cohesion_r,
                      soil_depth), is.null))) {
    stop("All inputs must be provided.")
  }
  
  # check data type
  if (!all(sapply(list(friction_angle,
                       transmissivity,
                       sca, 
                       cohesion_s, 
                       cohesion_r, 
                       soil_depth), inherits, what = "SpatRaster"))) {
    stop("Slope, friction angle, bulk density, transmissivity, specific catchment area, cohesion, and soil depth must be SpatRaster objects.")
  }
  
  # check resolution, crs, and extent
  if (!all(terra::compareGeom(friction_angle,
                              transmissivity,
                              sca,
                              cohesion_s,
                              cohesion_r,
                              soil_depth, 
                              stopOnError = FALSE))) {
    stop("All rasters must have matching extent, resolution, and CRS.")
  }
  
  # ---- convert units ----
  
  # degrees --> radians
  ## slope
  # slope_rad <- deg2rad(slope)
  ## friction angle
  # friction_angle <- deg2rad(friction_angle)
  
  # kPa --> Pa
  ## soil cohesion
  # cohesion_s <- cohesion_s * 1000
  ## root cohesion
  # cohesion_r <- cohesion_r * 1000
  
  # ---- calculate intermediates ----
  
  # sin slope
  # sin_slope <- compute_sin_slope(slope_rad)
  
  # compute wetness or use precomputed
  wetness <- if (!is.null(wetness_override)) {
    wetness_override
  } else {
    compute_wetness(R, transmissivity, sca, sin_slope)
  }
  
  # compute cohesion star or use precomputed
  ## if computing converts from kPa to Pa
  cohesion_star <- if (!is.null(cohesion_star_override)) {
    cohesion_star_override
  } else {
    (cohesion_s * 1000 + cohesion_r * 1000) / (bulk_density * g * soil_depth)
  }
  
  tan_fa <- if (!is.null(tan_fa_override)) tan_fa_override else tan(deg2rad(friction_angle))
  cos_slope <- if (!is.null(cos_slope_override)) cos_slope_override else cos(slope_rad)
  
  
  # ---- main FS calculation ----
  
  # calculate factor of safety
  FS <- (cohesion_star + cos_slope * (1 - wetness * (density_w / bulk_density)) * tan_fa) / sin_slope
  
  
  if (clamp_FS) {
    FS <- terra::clamp(FS, lower = 0, upper = 5)
  }
  
  # ---- post-processing ----
  # smoothing option?
  
  # ---- return ----
  
  if (intermediates) {
    return(list(
      FS = FS,
      wetness = wetness,
      cohesion_dim = cohesion_star
    ))
  } else {
    return(FS)
  }
}

# ===============================
# Optimized Green-Ampt Solver
# ===============================
# Optimized Green-Ampt Solver
compute_transient_wetness <- function(K_sat, psi, d_theta, rain_intensity, duration, soil_depth, F_template) {
  F_final <- F_template 
  ponding_possible <- rain_intensity > K_sat
  
  if (global(ponding_possible, "sum", na.rm = TRUE)[1,1] > 0) {
    denom <- rain_intensity - K_sat
    denom[denom <= 0] <- 1e-6 
    SM <- max(psi * d_theta, 1e-6)
    F_p <- (K_sat * SM) / denom
    t_p <- F_p / rain_intensity
    is_ponded <- ponding_possible & (duration > t_p)
    
    if (global(is_ponded, "sum", na.rm = TRUE)[1,1] > 0) {
      F_guess <- F_p + K_sat * (duration - t_p)
      F_iter <- mask(F_guess, is_ponded, maskvalues = FALSE)
      target_Kt <- K_sat * duration
      
      for(k in 1:3) {  
        log_term <- log(1 + F_iter / SM)
        f_val <- F_iter - (SM * log_term) - target_Kt
        f_prime <- F_iter / (SM + F_iter)
        
        # Update guess
        F_iter <- F_iter - (f_val / f_prime)
      }
      F_final <- cover(F_iter, F_final)
    }
  }
  return(clamp(F_final / (d_theta * soil_depth), lower = 0, upper = 1.0))
}

landslide_probability <- function(
    slope, friction_angle, bulk_density, ksat, sca, cohesion_s, cohesion_r, soil_depth,
    R_steady, rain_intensity, duration, psi,
    g = 9.81, density_w = 1000, alpha = 10, n_bins = 10,
    perturb_var = c("friction_angle", "cohesion_s"),
    perturb_settings = list(
      friction_angle = list(sd = 5, min = 5, max = 50),
      ksat = list(sd = 0.01, min = 0.0001, max = 1),
      cohesion_s = list(sd = 2, min = 0, max = 20)
    ),
    output = NULL,
    progress = TRUE
) {
  
  # This affects all new rasters created during the function
  old_opt <- terraOptions(datatype = "FLT4S")
  on.exit(terraOptions(datatype = old_opt$datatype)) # Revert when done
  
  
  rain_intensity <- max(as.numeric(rain_intensity[1]), 0)
  duration <- max(as.numeric(duration[1]), 0)
  R_steady <- max(as.numeric(R_steady[1]), 0)
  n_bins <- as.integer(n_bins[1])
  
  to_rast <- function(x) if (is.character(x)) rast(x) else x
  slope <- to_rast(slope)
  friction_angle <- to_rast(friction_angle)
  ksat <- to_rast(ksat)
  sca <- to_rast(sca)
  cohesion_s <- to_rast(cohesion_s)
  cohesion_r <- to_rast(cohesion_r)
  soil_depth <- to_rast(soil_depth)
  bulk_density <- to_rast(bulk_density)
  
  # --- PRE-COMPUTE SECTION ---
  slope_rad <- deg2rad(slope)
  sin_slope <- compute_sin_slope(slope_rad)
  cos_slope <- cos(slope_rad)
  soil_depth_normal <- soil_depth * cos_slope
  porosity_base <- 1 - (bulk_density / 2650)
  F_total_template <- (slope * 0) + (rain_intensity * duration)
  inv_n_bins <- 1 / n_bins
  
  # --- LHS SETUP ---
  lhs_matrix <- matrix(NA, nrow = n_bins, ncol = length(perturb_var))
  colnames(lhs_matrix) <- perturb_var
  for (j in seq_along(perturb_var)) {

    p_values <- (seq_len(n_bins) - 0.5) / n_bins 
    lhs_matrix[, j] <- qnorm(p_values, mean = 0, sd = perturb_settings[[perturb_var[j]]]$sd)
  }
  
  prob_landslide <- rast(slope)
  values(prob_landslide) <- 0
  
  
  if (progress) {
    pb <- txtProgressBar(min = 0, max = n_bins, style = 3)
  }
  
  # ---- SIMULATION LOOP ----
  for (i in 1:n_bins) {
    params <- list(friction_angle = friction_angle, cohesion_s = cohesion_s,
                   ksat = ksat, soil_depth = soil_depth_normal,
                   cohesion_r = cohesion_r, bulk_density = bulk_density, 
                   R_steady = R_steady)
    
    for (var in perturb_var) {
      shift <- lhs_matrix[i, var]
      params[[var]] <- clamp(params[[var]] + shift, 
                             lower = perturb_settings[[var]]$min, 
                             upper = perturb_settings[[var]]$max)
    }
    
    current_trans <- params$ksat * params$soil_depth
    m1 <- compute_wetness(params$R_steady, current_trans, sca, sin_slope)
    d_theta_dynamic <- clamp(porosity_base * (1 - m1), lower = 1e-4)
    
    # Updated transient wetness
    m2 <- compute_transient_wetness(params$ksat, psi, d_theta_dynamic, 
                                    rain_intensity, duration, params$soil_depth, 
                                    F_total_template)
    
    wetness_total <- clamp(m1 + m2, lower = 0, upper = 1.0)
    
    # Calculate saturated density using porosity and density of water (1000 kg/m3)
    # porosity_base is already 1 - (bulk_density / 2650) 
    rho_sat <- params$bulk_density + (porosity_base * density_w)
    dynamic_bulk_density <- (wetness_total * rho_sat) + ((1 - wetness_total) * params$bulk_density)
    
    fs <- SINMAP_FS(
      slope_rad = slope_rad, sin_slope = sin_slope,
      friction_angle = params$friction_angle,
      bulk_density = dynamic_bulk_density,
      R = R_steady, transmissivity = current_trans, sca = sca,
      cohesion_s = params$cohesion_s, cohesion_r = params$cohesion_r,
      soil_depth = params$soil_depth,
      wetness_override = wetness_total,
      cos_slope_override = cos_slope
    )
    
    prob_landslide <- prob_landslide + (inv_n_bins / (1 + exp(alpha * (fs - 1))))
    
    if (progress) setTxtProgressBar(pb, i)
    
    # Force memory release at the end of each bin
    rm(params, m1, m2, fs, wetness_total) 
    gc()
    
    # terra::tmpFiles(remove = TRUE)
  }
  
  if (progress) close(pb)
  if (!is.null(output)) writeRaster(prob_landslide, output, overwrite = TRUE)
  
  return(prob_landslide)
}