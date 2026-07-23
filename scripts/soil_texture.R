## classifying SoilGrid to USDA classification
## based on: https://www.mathworks.com/matlabcentral/fileexchange/45468-soil_classification-sand-clay-t-varargin
## &: https://code.usgs.gov/ghsc/lhp/regiongrow3d/-/blob/main/lib/functions/soil_classification_NM.m?ref_type=heads

## 1 = sand
## 2 = loamy sand
## 3 = sandy loam
## 4 = loam 
## 5 = silt loam
## 6 = silt
## 7 = sandy clay loam
## 8 = clay loam 
## 9 = silty clay loam
## 10 = sandy clay
## 11 = silty clay
## 12 = clay

#-------------------------------------------------------------------------------
library(terra) # required
#-------------------------------------------------------------------------------

soil_texture <- function(sand, clay){
  stopifnot(terra::compareGeom(sand, clay)) # check inputs
  
  # start main function process
  rast.out <- terra::lapp(c(sand, clay), fun = function(sand, clay) {
    silt <- 1 - (sand + clay) # make silt layer
    SC <- rep(NA_integer_, length(sand)) # make a vector the size of the inputs
    
    SC[(silt + 1.5 * clay) < 0.15] <- 1
    SC[(silt + 1.5 * clay) >= 0.15 & (silt + 2 * clay) < 0.3] <- 2
    cond3 <- (clay >= 0.07 & clay <= 0.2 & sand > 0.52 & (silt + 2 * clay) >= 0.3) |
      (clay <0.07 & silt < 0.5 & (silt + 2 * clay) >= 0.3) # require because of 'OR' statement
    SC[cond3] <- 3
    SC[clay >=0.07 & clay <= 0.27 & silt >= 0.28 & silt < 0.5 & sand <= 0.52] <- 4
    cond5 <- (silt >= 0.5 & clay >= 0.12 & clay < 0.27) | (silt >=0.5 & silt < 0.8 & clay < 0.12) # required becasue of 'OR" statement
    SC[cond5] <- 5
    SC[silt >= 0.8 & clay < 0.12] <- 6
    SC[clay >= 0.2 & clay < 0.35 & silt < 0.28 & sand > 0.45] <- 7
    SC[clay >= 0.27 & clay < 0.4 & sand > 0.2 & sand <=0.45] <- 8
    SC[clay >= 0.27 & clay < 0.4 & sand <= 0.2] <- 9
    SC[clay >= 0.35 & sand >= 0.45] <- 10
    SC[clay >= 0.4 & silt >= 0.4] <- 11
    SC[clay >= 0.4 & sand <= 0.45 & silt <0.4] <- 12
    
    return(SC)
  })
  names(rast.out) <- "soil_class"
  return(rast.out)
}