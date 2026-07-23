## function to get ground water well depth to bedrock data
# ---- required libraries ----

library(terra)
library(dplyr)

# ---- define the function ----
well_depth <- function(well_points, target, buff = 100) {
  
  # check that the projects match
  if (crs(well_points) != target$crs) {
    well_points <- project(well_points, target$crs)
  }
  
  # make a bounding box from the extent
  bbox <- as.polygons(ext(target$extent[1], 
                          target$extent[2], 
                          target$extent[3], 
                          target$extent[4]), 
                      crs = target$crs)
  # add the buffer distance
  bbox <- buffer(bbox, width = buff)
  
  # subset the well points
  selected_points <- mask(well_points, bbox)
  
  return(selected_points)
}