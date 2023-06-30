#
# Settings tree
#

construct_flipping_point_vector <- function(nAlt, nCS){
  flipping_point_vector <- matrix(0, ncol = nCS)
  for (iCS in 1:nCS){
    flipping_point_vector[iCS] = nAlt^(nCS - iCS)
  }
  return(flipping_point_vector)
}

# lookup point
setting_to_index <- function(setting, flipping_point_vector){
  
  return( as.double(1 + flipping_point_vector %*% (setting-1)) )

}

# convert an index to a setting
index_to_setting <- function(index, flipping_point_vector, nAlt){
  
  nCS <- length(flipping_point_vector)
  
  setting <- as.vector(floor( (index - 1) / flipping_point_vector) %% nAlt + 1)
  
  return(setting)
  
}


# convert a setting to a response vector
setting_to_response <- function(setting, nAlt){
  
  Y <- matrix(0, nrow = nAlt * length(setting))
  
  for (iCS in 1:length(setting)){
    
    Y[ nAlt*(iCS-1) + setting[iCS] ] <- 1
    
  }
  
  return(Y)
  
}



