## Get bulk density using ROSETTA volumetric water content
# based on : https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-research-unit/docs/model/rosetta-model/

library(terra)

# construct the ROSETTA dataframe
rosetta <- data.frame(class_id = 1:12)

rosetta$texture_class <- c("sand", "l_sand", "s_loam", "loam", "si_loam", "silt",
                           "s_c_l", "c_loam", "si_c_l", "s_clay", "si_clay", "clay")

rosetta$N <- c(308, 201, 476, 242, 330, 6, 87, 140, 172, 11, 28, 84)

rosetta$theta_r <- c(0.053, 0.049, 0.039, 0.061, 0.065, 0.050, 0.063, 0.079, 
                     0.090, 0.117, 0.111, 0.098)

rosetta$theta_s <- c(0.375, 0.390, 0.387, 0.399, 0.439, 0.489, 0.384, 0.442, 
                     0.482, 0.385, 0.481, 0.459)

rosetta$log_alpha <- c(-1.453, -1.459, -1.574, -1.954, -2.296, -2.182, -1.676, 
                       -1.801, -2.076, -1.476, -1.790, -1.825)

rosetta$log_n <- c(0.502, 0.242, 0.161, 0.168, 0.221, 0.225, 0.124, 0.151, 
                   0.182, 0.082, 0.121, 0.098)

rosetta$k_s <- c(2.808, 2.022, 1.583, 1.081, 1.261, 1.641, 1.120, 0.913, 1.046, 
                 1.055, 0.983, 1.169)

rosetta$k_o <- c(1.389, 1.386, 1.190, 0.568, 0.243, 0.524, 0.841, 0.699, 0.349, 
                 0.637, 0.501, 0.472)

rosetta$L <- c(-0.930, -0.874, -0.861, -0.371, 0.365, 0.624, -1.280, -0.763, 
               -0.156, -3.665, -1.287, -1.561)

# build the function to convert texture class into bulk density
bulk_density <- function(texture, assumed_density = 2650) {
  
  # make a vector to hold the texture raster
  v <- c(texture)
  
  # reclassify to the theta_s value for the corresponding texture class
  v[texture == 1] <- rosetta$theta_s[rosetta$class_id == 1]
  v[texture == 2] <- rosetta$theta_s[rosetta$class_id == 2]
  v[texture == 3] <- rosetta$theta_s[rosetta$class_id == 3]
  v[texture == 4] <- rosetta$theta_s[rosetta$class_id == 4]
  v[texture == 5] <- rosetta$theta_s[rosetta$class_id == 5]
  v[texture == 6] <- rosetta$theta_s[rosetta$class_id == 6]
  v[texture == 7] <- rosetta$theta_s[rosetta$class_id == 7]
  v[texture == 8] <- rosetta$theta_s[rosetta$class_id == 8]
  v[texture == 9] <- rosetta$theta_s[rosetta$class_id == 9]
  v[texture == 10] <- rosetta$theta_s[rosetta$class_id == 10]
  v[texture == 11] <- rosetta$theta_s[rosetta$class_id == 11]
  v[texture == 12] <- rosetta$theta_s[rosetta$class_id == 12]
  
  # calculate bulk density with constant density (kg/m^3), default = quartz
  bd <- assumed_density * (1-v)
  
  return(bd)
}