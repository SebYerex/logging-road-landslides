library(terra)

#-------------------------------------------------------------------------------
# unsat cohesion function
## based on: https://doi.org/10.1016/j.still.2018.07.006

unsat_cohesion <- function(clay, coarseSand, veryFineSand) {
  
  coh = -0.75 + 2.07 * clay^0.5 - 5.87 * log10(coarseSand) - 0.035 * veryFineSand^2
  
  coh <- round(coh, 3)
  
  coh <- ifel(coh < 0, 0, coh)
  
  return(coh)
}
