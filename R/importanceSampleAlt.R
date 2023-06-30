importanceSampleAlt <- function(X, mu, sig, nChoiceSet,
                                nU){
  
  nAlternative <- nrow(X) / nChoiceSet
  nRE <- sum(sig > 0)
  
  fpv <- construct_flipping_point_vector(nAlternative, nChoiceSet)
  
  
  
  conditional_probs <- rep(0.0, nU)
  
  CJ <- nChoiceSet * nAlternative
  
  Delta <- matrix(0, nrow = CJ, ncol = CJ)
  PPt   <- matrix(0, nrow = CJ, ncol = CJ)
  

  # Eu(p [u^2 / sigma^3, ...])
  PRatiosT <- matrix(0, nrow = CJ, ncol = nRE)
  
  responseProbs <- matrix(0, nrow = CJ, ncol = nU)
  
  U22 <- matrix(0, nrow = nRE, ncol = nU)
  
  # Monte Carlo Sampling for the terms that only depend on a single quantity (u)
  impsamp_update_easy_terms_inplace(Delta, PPt, PRatiosT,
                                    responseProbs, U22,
                                    X,
                                    mu,
                                    sig,
                                    nChoiceSet,
                                    nU = nU)
  
  # Now the parts that are conditional on the outcome vector
  
  Epypyt <- matrix(0, nrow = CJ, ncol = CJ)
  EpyRatiot <- matrix(0, nrow = CJ, ncol = nRE)
  EReatioRatiot <- matrix(0, nrow = nRE, ncol = nRE)
  
  sigma_sigmaprime <- (1 / sqrt(sig)) %*% t(1/sqrt(sig))
  
  marginal_p_setting <- rep(NA, 2^nChoiceSet)
  
  # loop starts here
  print("Todo calculate number of settings correctly!")
  for (iY in 1:2^nChoiceSet){
    
    setting_i <- index_to_setting(iY, flipping_point_vector = fpv, nAlt = nAlternative)
    Y_i        <- setting_to_response(setting_i, nAlt = nAlternative)
    impsamp_calc_conditional_prob_vector_inplace_R(conditional_probs, Y_i, responseProbs, nChoiceSet)
    
    # Marginal probability of this combination of outcomes (P Y=y)
    Py_eq_y <- sum(conditional_probs) / nU; 
    marginal_p_setting[iY] <- Py_eq_y
    
    # Term 1
    Py     <- responseProbs %*% conditional_probs / nU / Py_eq_y
    Epypyt <- Epypyt + Py %*% t(Py) * Py_eq_y
    
    # Term 2
    RatioY <- U22 %*% conditional_probs / nU / Py_eq_y
    PyRatiot  <- Py %*% t(RatioY)
    EpyRatiot <-  EpyRatiot + PyRatiot * Py_eq_y
    
    # Term 3
    RatioRatiot <- RatioY %*% t(RatioY)
    EReatioRatiot <-  EReatioRatiot + RatioRatiot  * Py_eq_y
    
  }
  
  I11 <- t(X) %*% (Delta - PPt + Epypyt) %*% X
  
  I12 <- t(X) %*% (PRatiosT - EpyRatiot)
  
  I22 <- EReatioRatiot - sigma_sigmaprime
  #I22 <- EReatioRatiot
  
  info_matrix_appr <- rbind( cbind(I11, I12), cbind(t(I12), I22) )
  
  
}