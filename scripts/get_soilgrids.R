## function to get SoilGrids data

# adapted from https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/commit/23fe857b81fea0149526fbdee2115d1480b1568c

#-------------------------------------------------------------------------------

library(terra) # required

#-------------------------------------------------------------------------------
get_soilgrids <- function(
    resolution = "250m",
    voi = "sand",
    depth = "5-15cm",
    quantile = "mean",
    target,
    load = TRUE,
    write_to = NULL
) {
  
  # enable progress bar
  old_progress <- terraOptions()$progress
  terraOptions(progress = 1)
  on.exit(terraOptions(progress = old_progress), add = TRUE)
  
  # get extent from target
  aoi_extent <- terra::ext(
    target$extent["xmin"],
    target$extent["xmax"],
    target$extent["ymin"],
    target$extent["ymax"]
  )
  
  # get crs from target
  aoi_crs <- target$crs
  
  # set SoilGrids crs (Goode Homolosine)
  soilgrids_crs <- "ESRI:54052"
  
  # standard depths (cm)
  available_depths <- list(
    "0-5cm" = 5,
    "5-15cm" = 10,
    "15-30cm" = 15,
    "30-60cm" = 30,
    "60-100cm" = 40,
    '100-200cm' = 100
  )
  
  # expand to multiple layers if needed
  if (depth %in% names(available_depths)) {
    depth_layers <- list(depth)
    thickness <- available_depths[[depth]]
  } else {
    # parse composite depth (e.g., "0-100cm")
    depth_parts <- as.numeric(unlist(regmatches(depth, gregexpr("[0-9]+", depth))))
    if (length(depth_parts) !=2) stop("Invalid depth format.")
    
    lower <- depth_parts[1]
    upper <- depth_parts[2]
    
    # filter layers that fall within composite depth
    layer_names <- names(available_depths)
    layer_bounds <- t(sapply(strsplit(layer_names, "-|cm"), function(x) as.numeric(x[1:2])))
    selected <- which(layer_bounds[,2] > lower & layer_bounds[,1] < upper)
    
    if (length(selected) == 0) stop("No matching layers for this depth range.")
    
    depth_layers <- layer_names[selected]
    thickness <- sapply(depth_layers, function(x) available_depths[[x]])
  }
  
  # project target extent to SoilGrids crs
  projected_bb <- terra::project(aoi_extent, from = aoi_crs, to = soilgrids_crs)
  
  # load each raster and crop
  rasters <- list()
  for (i in seq_along(depth_layers)) {
    voi_layer <- paste(voi, depth_layers[[i]], quantile, sep = "_")
    
    # build SoilGrids URL
    rast_url <- switch(
      resolution,
      "250m" = paste0("/vsicurl/https://files.isric.org/soilgrids/latest/data/", voi, "/", voi_layer, ".vrt"),
      "1000m" = paste0("/vsicurl/https://files.isric.org/soilgrids/latest/data_aggregated/1000m/", voi, "/", voi_layer, "_1000.tif"),
      "5000m" = paste0("/vsicurl/https://files.isric.org/soilgrids/latest/data_aggregated/5000m/", voi, "/", voi_layer, "_5000.tif"),
      stop("Invalid resolution. Choose from '250m', '1000m', or '5000m'.")
    )
    
    r <- terra::rast(rast_url)
    terra::crs(r) <- soilgrids_crs
    r <- terra::crop(r, projected_bb)
    rasters[[i]] <- r * thickness[[i]] # weighted by thickness  
  }
  
  # combine weighted layers and divide by the total thickness
  total_thickness <- sum(thickness)
  combined_rast <- Reduce('+', rasters) / total_thickness
  
  # reproject back to target crs
  combined_rast <- terra::project(combined_rast, aoi_crs)
  
  # write to disk if required
  if (!is.null(write_to)) {
    terra::writeRaster(combined_rast, write_to, overwrite = TRUE)
    if (!load) return(write_to)
  }
  
  
  return(combined_rast)
  
}
