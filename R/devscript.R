#
# Looking at alternative version of importance sampling ----
#

if (0 == 1){
  
  library(PanelMixedDesign)
  library(DoE.base)
  
  nattr <- 3
  nChoiceSet <- 6
  nAlternative <- 2 
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.3, 0.4)
  st = NULL
  fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
  for(i in 1:nChoiceSet){
    st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=FALSE) )
  }
  designm=(fullfac)[st,]
  contr=rep('contr.sum',nattr)
  names(contr)=names(designm)
  contr=as.list(contr)
  M = model.matrix(~., designm, contrasts = contr)[,-1] 
  Y <- gen_all_choice_seq(nChoiceSet, nAlternative)
  infoAppr <- PMLInfoApprox(M, method = "importance",
                            nChoiceSet = nChoiceSet, effectMean = mu,
                            effectVar = sig,
                            Y = Y,
                            opts = list(nU = 50000))
  
  temp <- PMLInfoApprox(M, method = "importance",
                nChoiceSet = nChoiceSet, effectMean = mu,
                effectVar = sig,
                opts = list(nY = 10000, nU = 1000))
  
  # temp2 <- PMLInfoApprox(M, method = "MSM",
  #               nChoiceSet = nChoiceSet, effectMean = mu,
  #               effectVar = sig)
  
  # 
  # temp2 <- PMLInfoApprox(M, method = "importance",
  #                       nChoiceSet = nChoiceSet, effectMean = mu,
  #                       effectVar = sig)
  # 
  plot(c(infoAppr), c(temp))
  
  plot(c(temp2), c(infoAppr))
  
  # 
  # New ----
  #
  
  # Finding the "most likely" outcome
  XB <- M %*% mu
  expXB <- exp(XB)
  index_set <- 1:nAlternative
  myY <- matrix(0, nrow = nChoiceSet*nAlternative, ncol = ncol(Y))
  
  # Finding the probabilities if each parameter were excluded, in turn
  nBeta <- 3
  XBparamdrop <- matrix(0, nrow = nrow(M), ncol = 3)
  expXBp <- matrix(0, nrow = nrow(M), ncol = 3)
  for (p in 1:3){
    mupdrop    <- mu
    mupdrop[p] <- 0
    XBp <- M %*% mupdrop
    expXBp[,p] <- exp(XBp)
  }
  
  rankings <- matrix(0, nrow =  nChoiceSet*nAlternative, ncol = 4)
  
  delta_prob_param_drop <- matrix(0, nrow = nrow(M), ncol = 3)
  
  info_parameter <- matrix(0, nrow = nChoiceSet*nAlternative, ncol = length(mu))
  for (iCS in 1:nChoiceSet){
    
    # Convert to probs
    probs <- expXB[index_set] / sum(expXB[index_set])
    rankings[index_set,1] <- probs
    
    # Get the probs if each variable were to be excluded:
    for (p in 1:nBeta){
      probs_param_p_dropped <- expXBp[index_set,p] / sum(expXBp[index_set,p])
      delta_prob_param_drop[index_set,p] <- probs - probs_param_p_dropped
    }
    
    # cols are which is selected, rows are parameter it provides information about. Might need abs?
    XPX <- t(M[index_set, ]) %*% (diag(1, nAlternative) - cbind(probs, probs)) %*% diag(probs)
    info_parameter[index_set, ] <- t(XPX)
    
    varhat <- diag(M[index_set,] %*% diag(sig) %*% t(M[index_set, ]))
    rankings[index_set,2] <- as.numeric(varhat)
    
    maxind <- which(expXB[index_set] == max(expXB[index_set] ))
    ranks_within_cs <- rank(-probs, ties.method = "random")
    rankings[index_set,3] <- ranks_within_cs
    
    myY[index_set[maxind], ] <- 1
    
    index_set <- index_set + nAlternative
  }
  rankings[,4] <- rank(-rankings[,1])
  
  
  
  unweighted_value = FALSE
  mstore <- matrix(0, nrow = 5)
  total_sum_marginal_probs <- 0.0
  setting    <- c(1,1,1,1,1,1)
  Y_i        <- get_setting(setting, nAlternative, rankings[,3])
  infoAppr_i <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000, return_unweighted = unweighted_value, mstore)
  total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
  
  store_margin <- mstore[1]
  
  max_change <- max(abs(infoAppr_i[1:3, 1:3]))
  n_draw <- 1
  
  diag_mat_store           <- matrix(diag(infoAppr_i))
  predicted_diag_mat_store <- matrix(diag(infoAppr_i) / total_sum_marginal_probs)
  
  
  initial <- infoAppr_i
  
  prob_store <- rep(0, nCS)
  mat_store <- list()
  
  # second level
  for (iCS in 1:nCS){
    setting    <- c(1,1,1,1,1,1)
    setting[iCS] <- 2
    Y_i        <- get_setting(setting, nAlternative, rankings[,3])
    prob_store[iCS] <- prod(rankings[Y_i == 1, 1])
    info_setting    <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
    total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
    mat_store[[iCS]] <- info_setting
    infoAppr_i       <- infoAppr_i + info_setting
    max_change <- c(max_change, max(abs(info_setting[1:3, 1:3])))
    diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
    predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
    n_draw <- n_draw + 1
    store_margin <- c(store_margin, mstore[1])
  }
  
  #testdf <- tibble(p = prob_store, m = mat_store)
  
  
  # # compare each to each other (upper third)
  # any_gt_0 <- matrix(NA, nrow = nCS, ncol = nCS)
  # for (i in 1:nCS){
  #   for (j in 1:nCS){
  #     if (i != j){
  #       mat1 <- testdf$m[[i]]
  #       mat2 <- testdf$m[[j]]
  #       diffmat <- mat1 - mat2
  #       diffmat <- diffmat[1:3, 1:3]
  #       ev <- eigen(diffmat)$values
  #       any_gt_0[i, j] <- as.numeric(any(ev < 0))
  #     }
  #   }
  # }
  
  # sets of P ?
  
  
  
  
  
  
  # third level
  for (iCS in 1:nCS){
    for (jCS in iCS:nCS){
      if (iCS != jCS){
        setting    <- c(1,1,1,1,1,1)
        setting[iCS] <- 2
        setting[jCS] <- 2
        Y_i        <- get_setting(setting, nAlternative, rankings[,3])
        info_setting <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
        total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
        infoAppr_i <- infoAppr_i + info_setting
        max_change <- c(max_change, max(abs(info_setting[1:3, 1:3])))
        diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
        predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
        
        n_draw <- n_draw + 1
        store_margin <- c(store_margin, mstore[1])
        }
    }
  }
  
  # fourth level
  for (iCS in 1:nCS){
    for (jCS in iCS:nCS){
      for (kCS in jCS:nCS){
        if (iCS != jCS){
          if (iCS != kCS){
            if (jCS != kCS){
              setting    <- c(1,1,1,1,1,1)
              setting[iCS] <- 2
              setting[jCS] <- 2
              setting[kCS] <- 2
              print(setting)
              Y_i        <- get_setting(setting, nAlternative, rankings[,3])
              info_setting <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
              total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
              infoAppr_i <- infoAppr_i + info_setting
              max_change <- c(max_change, max(abs(info_setting[1:3, 1:3])))
              diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
              predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
              
              n_draw <- n_draw + 1
              store_margin <- c(store_margin, mstore[1])
            }
          }
        }
      }
    }
  }
  
  # fifth level
  for (iCS in 1:nCS){
    for (jCS in iCS:nCS){
      for (kCS in jCS:nCS){
        for (lCS in kCS:nCS){
          if (iCS != jCS){
            if (iCS != kCS){
              if (jCS != kCS){
                if (lCS != kCS){
                  setting    <- c(1,1,1,1,1,1)
                  setting[iCS] <- 2
                  setting[jCS] <- 2
                  setting[kCS] <- 2
                  setting[lCS] <- 2
                  print(setting)
                  Y_i        <- get_setting(setting, nAlternative, rankings[,3])
                  info_setting <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
                  total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
                  infoAppr_i <- infoAppr_i + info_setting
                  max_change <- c(max_change, max(abs(info_setting[1:3, 1:3])))
                  diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
                  predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
                  
                  n_draw <- n_draw + 1
                  store_margin <- c(store_margin, mstore[1])
                }
              }
            }
          }
        }
      }
    }
  }
  
  # sixth level
  for (iCS in 1:nCS){
    for (jCS in iCS:nCS){
      for (kCS in jCS:nCS){
        for (lCS in kCS:nCS){
          for (mCS in lCS:nCS){
          if (iCS != jCS){
            if (iCS != kCS){
              if (jCS != kCS){
                if (lCS != kCS){
                  if (mCS != lCS){
                      setting    <- c(1,1,1,1,1,1)
                      setting[iCS] <- 2
                      setting[jCS] <- 2
                      setting[kCS] <- 2
                      setting[lCS] <- 2
                      setting[mCS] <- 2
                      print(setting)
                      Y_i        <- get_setting(setting, nAlternative, rankings[,3])
                      info_setting <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
                      total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
                      infoAppr_i <- infoAppr_i + info_setting
                      max_change <- c(max_change, max(abs(info_setting[1:3, 1:3])))
                      diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
                      predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
                      
                      n_draw <- n_draw + 1
                      store_margin <- c(store_margin, mstore[1])
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  setting    <- c(2,2,2,2,2,2)
  Y_i        <- get_setting(setting, nAlternative, rankings[,3])
  infoAppr_i <- infoAppr_i + R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
  total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
  #mat_11_store <- c(mat_11_store, infoAppr_i[1,1])
  diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
  predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
  
  store_margin <- c(store_margin, mstore[1])
  n_draw <- n_draw + 1
  
  # Looking at biggest drivers for each parameter
  param_changes <- t(diag_mat_store - cbind(rep(0, 6),  diag_mat_store[,1:63]))
  combined_info <- cbind(store_margin, param_changes)
  
  plot(combined_info[,1], combined_info[,2])
  plot(combined_info[,1], combined_info[,3])
  plot(combined_info[,1], combined_info[,4], main = "Beta 3")
  plot(log(combined_info[,1]), combined_info[,5], main = "Sigma 1")
  plot(log(combined_info[,1]), combined_info[,6], main = "Sigma 2")
  plot(log(combined_info[,1]), combined_info[,7], main = "Sigma 3")
  
  
  
  
  infoAppr_i <- PMLInfoApprox(M, method = "importance",
                             nChoiceSet = nChoiceSet, effectMean = mu,
                             effectVar = sig,
                             Y = Y_i)

  
  
  get_setting <- function(setting, J, within_cs_rankings){
    
    nCS <- length(setting)
    
    index_set <- 1:J
    Y <- matrix(0, nrow = J * nCS)
    
    for (i in 1:nCS){
      
      flip_index      <- which(within_cs_rankings[index_set] == setting[i])
      Y[index_set[flip_index], 1] <- 1
      
      index_set <- index_set + nAlternative
    }
    
    return(Y)
    
  }

  
  
  
  
  
  # Completely Random Version
  nrep <- 5000
  diag_mat_store           <- NULL
  predicted_diag_mat_store <- NULL
  infoAppr_i <- matrix(0, nrow = 6, ncol = 6)
  total_sum_marginal_probs <- 0.0
  
  # second level
  for (ii in 1:nrep){
    setting    <- sample(1:2, size = 6, replace = TRUE)
    Y_i        <- get_setting(setting, nAlternative, rankings[,3])
    info_setting    <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=10000,  return_unweighted = unweighted_value, mstore)
    total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]

    infoAppr_i       <- infoAppr_i + info_setting
    
    diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
    predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
    
    n_draw <- n_draw + 1
  }
  
  infoAppr_i / total_sum_marginal_probs
  
  plot(predicted_diag_mat_store[1,])
  plot(predicted_diag_mat_store[2,])
  plot(predicted_diag_mat_store[3,])
  plot(predicted_diag_mat_store[4,])
  plot(predicted_diag_mat_store[5,])
  plot(predicted_diag_mat_store[6,])
  
  
  
  
  
  
  
  #
  # New Random Order ----
  #
  
  nattr <- 3
  nChoiceSet <- 10
  nAlternative <- 3 
  nAlt <- nAlternative
  nCS <- nChoiceSet
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.3, 0.4)
  st = NULL
  fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
  for(i in 1:nChoiceSet){
    st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=FALSE) )
  }
  designm=(fullfac)[st,]
  contr=rep('contr.sum',nattr)
  names(contr)=names(designm)
  contr=as.list(contr)
  M = model.matrix(~., designm, contrasts = contr)[,-1] 
  Y <- gen_all_choice_seq(nChoiceSet, nAlternative)
  infoAppr <- PMLInfoApprox(M, method = "importance",
                            nChoiceSet = nChoiceSet, effectMean = mu,
                            effectVar = sig,
                            Y = Y,
                            opts = list(nU = 1000))
  
  
  # Random Version without repetition
  nrep <- 60000
  diag_mat_store           <- NULL
  predicted_diag_mat_store <- NULL
  infoAppr_i <- matrix(0, nrow = 6, ncol = 6)
  total_sum_marginal_probs <- 0.0
  unweighted_value = FALSE
  mstore <- matrix(0, nrow = 5)
  n_draw <- 0
  
  fpv <- construct_flipping_point_vector(nAlt, nCS)
  
  visit_order <- sample(1:(nAlt)^(nCS))
  
  # second level
  for (ii in 1:min(nrep, (nAlt)^(nCS)) ){
    setting    <- index_to_setting(visit_order[ii], fpv, nAlt)
    Y_i        <- setting_to_response(setting, nAlt)
    info_setting    <- R_impsamp_calc_info_given_Y(Y_i, M, mu, sig, nChoiceSet, nU=5000,  return_unweighted = unweighted_value, mstore)
    total_sum_marginal_probs = total_sum_marginal_probs + mstore[1]
    
    infoAppr_i       <- infoAppr_i + info_setting
    
    diag_mat_store <- cbind(diag_mat_store, matrix(diag(infoAppr_i)))
    predicted_diag_mat_store <- cbind(predicted_diag_mat_store, matrix(diag(infoAppr_i)) / total_sum_marginal_probs)
    
    n_draw <- n_draw + 1
  }
  
  infoAppr_i / total_sum_marginal_probs
  
  par(mfrow = c(3, 2))
  
  plot(predicted_diag_mat_store[1,], main = "VarCov Element (1,1)", ylab = "value", xlab = "Sampled outcome")
  abline(h = infoAppr[1,1], col = "red")
  
  plot(predicted_diag_mat_store[2,], main = "VarCov Element (2,2)", ylab = "value", xlab = "Sampled outcome")
  abline(h = infoAppr[2,2], col = "red")
  
  plot(predicted_diag_mat_store[3,], main = "VarCov Element (3,3)", ylab = "value", xlab = "Sampled outcome")
  abline(h = infoAppr[3,3], col = "red")
  
  plot(predicted_diag_mat_store[4,], main = "VarCov Element (4,4)", ylab = "value", xlab = "Sampled outcome")
  abline(h = infoAppr[4,4], col = "red")
  
  plot(predicted_diag_mat_store[5,], main = "VarCov Element (5,5)", ylab = "value", xlab = "Sampled outcome")
  abline(h = infoAppr[5,5], col = "red")
  
  plot(predicted_diag_mat_store[6,], main = "VarCov Element (6,6)", ylab = "value", xlab = "Sampled outcome")
  abline(h = infoAppr[6,6], col = "red")
  
  dev.off()
  
  
  
}



#
# New (term specific) importance sampling ----
#

if (0 == 1){
  
  
  library(PanelMixedDesign)
  library(DoE.base)
  
  nUsel <- 150000
  
  
  nattr <- 3
  nChoiceSet <- 6
  nAlternative <- 2 
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.3, 0.4)
  nRE <- sum(sig > 0)
  st = NULL
  fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
  for(i in 1:nChoiceSet){
    st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=FALSE) )
  }
  designm=(fullfac)[st,]
  contr=rep('contr.sum',nattr)
  names(contr)=names(designm)
  contr=as.list(contr)
  M = model.matrix(~., designm, contrasts = contr)[,-1] 
  
  Y <- gen_all_choice_seq(nChoiceSet, nAlternative)
  infoAppr <- PMLInfoApprox(M, method = "importance",
                            nChoiceSet = nChoiceSet, effectMean = mu,
                            effectVar = sig,
                            Y = Y,
                            opts = list(nU = 100000))
  infoAppr2 <- PMLInfoApprox(M, method = "importanceAlt",
                            nChoiceSet = nChoiceSet, effectMean = mu,
                            effectVar = sig,
                            opts = list(nU = 1000000))
  
  det(infoAppr)
  det(infoAppr2)
  sum(diag(solve(infoAppr)))
  sum(diag(solve(infoAppr2)))
  
  # For manual checking
  if (0==1){
    X = M
  }
  
  fpv <- construct_flipping_point_vector(nAlternative, nChoiceSet)
  
  conditional_probs <- rep(0.0, nUsel)
  
  CJ <- nChoiceSet * nAlternative
  
  Delta <- matrix(0, nrow = CJ, ncol = CJ)
  PPt   <- matrix(0, nrow = CJ, ncol = CJ)
  
  # Eu(p [u^2 / sigma^3, ...])
  PRatiosT <- matrix(0, nrow = CJ, ncol = nRE)
  
  responseProbs <- matrix(0, nrow = CJ, ncol = nUsel)
  
  U22 <- matrix(0, nrow = nRE, ncol = nUsel)
  
  # Monte Carlo Sampling for the terms that only depend on a single quantity (u)
  impsamp_update_easy_terms_inplace(Delta, PPt, PRatiosT,
                                    responseProbs, U22,
                                    M,
                                    mu,
                                    sig,
                                    nChoiceSet,
                                    nU = nUsel)
  
  #uprob <- dnorm(U)
  PUsq <- (dchisq(U^2, df = 1))
  
  # Now the parts that are conditional on the outcome vector
  
  Epypyt <- matrix(0, nrow = CJ, ncol = CJ)
  EpyRatiot <- matrix(0, nrow = CJ, ncol = nRE)
  EReatioRatiot <- matrix(0, nrow = nRE, ncol = nRE)
  
  sigma_sigmaprime <- (1 / sqrt(sig)) %*% t(1/sqrt(sig))
  
  marginal_p_setting <- rep(NA, 2^nChoiceSet)
  
  # loop starts here
  for (iY in 1:2^nChoiceSet){
    
    setting_i <- index_to_setting(iY, flipping_point_vector = fpv, nAlt = nAlternative)
    Y_i        <- setting_to_response(setting_i, nAlt = nAlternative)
    impsamp_calc_conditional_prob_vector_inplace_R(conditional_probs, Y_i, responseProbs, nChoiceSet)
    
    # Marginal probability of this combination of outcomes (P Y=y)
    Py_eq_y <- sum(conditional_probs) / nUsel; 
    marginal_p_setting[iY] <- Py_eq_y
    
    # Term 1
    Py     <- responseProbs %*% conditional_probs / nUsel / Py_eq_y
    Epypyt <- Epypyt + Py %*% t(Py) * Py_eq_y
    
    # Term 2
    RatioY <- U22 %*% conditional_probs / nUsel / Py_eq_y
    PyRatiot  <- Py %*% t(RatioY)
    EpyRatiot <-  EpyRatiot + PyRatiot * Py_eq_y

    # Term 3
    RatioRatiot <- RatioY %*% t(RatioY)
    EReatioRatiot <-  EReatioRatiot + RatioRatiot  * Py_eq_y
    
  }
  
  I11 <- t(M) %*% (Delta - PPt + Epypyt) %*% M
  
  I12 <- t(M) %*% (PRatiosT - EpyRatiot)
    
  I22 <- -sigma_sigmaprime + EReatioRatiot
  
  info_matrix_appr <- rbind( cbind(I11, I12), cbind(t(I12), I22) )

  det(info_matrix_appr)
  sum(diag(solve(info_matrix_appr)))
  
  plot(c(info_matrix_appr), c(infoAppr))
  
}



