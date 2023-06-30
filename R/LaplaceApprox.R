#
# Version of Laplace approach that uses R as an outer wrapper
#

test_laplace <- function(ny,
                         X,
                         b_mean,
                         var_vec,
                         n_choice_set){

  
  ## Actual Work
  
  ## Pre-Loop
  
  # Get the number of regression coefficients
  n_beta <- length(b_mean)
  
  # Detect random versus fixed effects
  if (any(var_vec < 0.0)){
    stop("Negative variance provided. Please ensure all variances are >= 0.0")
  }
  n_random_effect <- sum(var_vec > 0)
  n_fixed_effect  <- sum(var_vec == 0)
  
  # Storage for the information matrix, we will build onto this at each iteration
  information_matrix_approximation <- matrix(0,
                                             nrow = n_random_effect + n_beta,
                                             ncol =  n_random_effect + n_beta)
  
  
  # Some bookkeeping quantities
  n_alt_total <- nrow(X)
  choice_set_size <- n_alt_total / n_choice_set;
  
  # Get the inverse of the prior covariance
  omega = diag(1 / (var_vec))
  omega_re <- omega[1:n_random_effect, 1:n_random_effect]
  sigma_R <- var_vec # rename for elsewhere
  
  
  constant_term <- 100
  
  # Obs values for each sample (beta)
  sd_vec    <- sqrt(var_vec)
  u         <- matrix(rnorm(n_beta*ny), nrow = n_beta, ncol = ny)
  beta      <- diag(sd_vec) %*% u + b_mean
  
  # Calculate the epsilon terms (get added to linear predictor later)
  epsilon <- matrix(runif(n_alt_total*ny), nrow = n_alt_total)
  epsilon <- -log(-1.0 * log(epsilon))
  
  # Get the "most likely" Y for each draw of beta
  expXB <- exp(X %*% beta + epsilon)
  Y     <- matrix(0, nrow = n_alt_total, ncol = ny)
  
  # Get Xr, which just contains the columns with a random effect
  Xr <- X[,1:n_random_effect]
  delta_matrix <- matrix(0, nrow = n_alt_total, ncol= n_alt_total)
  
  # Update Y to have 1s in the "most likely" locations
  laplace_set_maxprob_Y_inplace(Y, expXB, choice_set_size)
  
  # Storage
  prob      <- matrix(0, nrow = n_alt_total)
  numerator <- rep(0, n_alt_total)
  re_ratio_condexpect_numerator <- rep(0, n_random_effect)
  
  # Information matrix
  I11 = matrix(0, nrow = n_beta, ncol = n_beta)
  I12 = matrix(0, nrow = n_beta, ncol = n_random_effect)
  I22 = matrix(0, nrow = n_beta, ncol = n_beta)
  
  n_converged = 0
  
  # Loop over draws of Y
  for (iy in 1:ny){
    
    u_init <- u[,iy]
    beta_i <- beta[,iy]
    Y_i    <- Y[,iy]
    
    opt_result <- optim(u_init,
                        fn = laplace_u1_fn_minimize,
                        b_mean=b_mean,
                        OmegaRE=omega_re,
                        X=X,
                        beta = beta_i,
                        prob = prob,
                        nRandEffect  = n_random_effect,
                        nFixedEffect = n_fixed_effect,
                        nChoiceSet   = nChoiceSet)
                  
    converged = FALSE
    if (opt_result$convergence == 0){converged = TRUE}
    
    # if converged, calculate information matrix
    if (converged == TRUE){
      
      beta_current <- beta_i
      beta_current[1:n_random_effect] <- beta_current[1:n_random_effect]  + opt_result$par
      
      # Get the probabilites for each selection under the current u
      update_response_probabilities_inplace(prob, X, beta_current, nChoiceSet)
      
      # Update delta, which is the matrix containing XXX
      update_delta_inplace(delta_matrix, prob, nChoiceSet)
      
      # Need to get H
      H = -t(Xr) %*% delta_matrix %*% Xr - omega_re;
      Hd <- -solve(H)
      H_inv_det <- det(Hd)
      
      # this is denom from c++ version
      denom <- sqrt(det(Hd)) *
        prod(prob^Y_i) *
        exp(-0.5 * sum( (1.0 / sd_vec) * opt_result$par^2))
        
      # Second layer of optimization
      u0_start <- opt_result$par
      
      for (row_j in 1:n_alt_total){
        
        u0j_mle <- u0_start
        cs_membership <- ceiling(row_j / choice_set_size)
        #i1        = (j1 / choice_set_size); // indexes which choice set this alternative is a member of
        
        cs_indices <- (choice_set_size * (cs_membership - 1) + 1):(choice_set_size * cs_membership)
        
        # Get the corresponding choice set
        X_cs     <- X[cs_indices, ]
        X_cs_rdm <- X[cs_indices, 1:n_random_effect]
        
        # Use optim to find the mle for u0j
        opt_result_rowj <- optim(u0j_mle,
                            fn = laplace_usj_fn_minimize,
                            b_mean=b_mean,
                            OmegaRE=omega_re,
                            X=X,
                            beta = beta_i,
                            prob = prob,
                            nRandEffect  = n_random_effect,
                            nFixedEffect = n_fixed_effect,
                            nChoiceSet   = nChoiceSet,
                            irow  = (row_j-1),
                            iCS = (cs_membership-1),
                            choice_set_size = choice_set_size,
                            X_cs_rdm = X_cs_rdm)
        
        converged = FALSE
        if (opt_result_rowj$convergence == 0){converged = TRUE}
        
        # TODO add quit out here if converged is false
        
        # Evaluate H at the MLE
        H_sj <- laplace_usj_evaluate_H(u0j_mle = opt_result_rowj$par,
                                     b_mean=b_mean,
                                     OmegaRE=omega_re,
                                     X=X,
                                     XR=Xr,
                                     beta = beta_i,
                                     prob = prob,
                                     nRandEffect  = n_random_effect,
                                     nFixedEffect = n_fixed_effect,
                                     nChoiceSet   = nChoiceSet,
                                     irow  = (row_j-1),
                                     iCS = (cs_membership-1),
                                     choice_set_size = choice_set_size,
                                     X_cs_rdm = X_cs_rdm)
    
        # Construct the portion of the numerator corresponding to this row
        H_sj_inv     <- -solve(H_sj);
        H_sj_inv_det <- det(H_sj_inv);
        
        # This is going to typically run into numeric problems if exponential term is small
        # TODO - log form of ratio
        f_theta_u0jmle <- exp(-(opt_result_rowj$par %*% omega_re %*% opt_result_rowj$par)/2)

        # Probability for the current row
        beta_current <- beta_i
        beta_current[1:n_random_effect] <- beta_current[1:n_random_effect]  + opt_result_rowj$par
        update_response_probabilities_inplace(prob, X, beta_current, nChoiceSet)
        prob_sj_cond_u0jmle <- prob[row_j]
        
        # Probability for the "observed" Y conditioned on u (maximized for this row)
        prob_Y_cond_u0jmle <- prod(prob^Y_i) 

        # This is the numerator of E_{u1} (p_{1sj} | y), pg 112 of Wei Zhang Thesis
        numerator[row_j] <- sqrt(H_sj_inv_det) * prob_sj_cond_u0jmle * prob_Y_cond_u0jmle * f_theta_u0jmle
      }
      
      #
      #---- E(uj^2 / sigmaj^3 | y_montecarlo)
      #
      # This is the second main piece that we need to evaluate the 
      # Laplace approximation to the information matrix
      
      # Now we loop over the individual random effects to achieve effect-specific maximizers
      for (i_re in 1:n_random_effect){
        
        u_init <- u0_start
        
        opt_result_effect <- optim(u_init,
                                 fn = laplace_uj_fn_minimize,
                                 b_mean=b_mean,
                                 sigma_re=sd_vec,
                                 OmegaRE = omega_re,
                                 X=X,
                                 beta = beta_i,
                                 prob = prob,
                                 constant_term = constant_term,
                                 nRandEffect  = n_random_effect,
                                 nFixedEffect = n_fixed_effect,
                                 nChoiceSet   = nChoiceSet,
                                 iRE = (i_re - 1),
                                 choice_set_size = choice_set_size,
                                 X_cs_rdm = X_cs_rdm)
        
        # Check convergence
        converged = FALSE
        if (opt_result_effect$convergence == 0){converged = TRUE}
        
        # Evaluate H at the MLE
        H_j <- laplace_uj_evaluate_H(u = opt_result_effect$par,
                                       b_mean=b_mean,
                                       OmegaRE=omega_re,
                                       sigma_re=sd_vec,
                                       X=X,
                                       XR=Xr,
                                       beta = beta_i,
                                       prob = prob,
                                       nRandEffect  = n_random_effect,
                                       nFixedEffect = n_fixed_effect,
                                       nChoiceSet   = nChoiceSet,
                                       iRE  = (i_re - 1),
                                       choice_set_size = choice_set_size,
                                       X_cs_rdm = X_cs_rdm,
                                      constant_term = constant_term)
        
        # Construct the portion of the numerator corresponding to this row
        H_j_inv     <- -solve(H_j);
        H_j_inv_det <- det(H_j_inv);
        
        # This is going to typically run into numeric problems if exponential term is small
        # TODO - log form of ratio
        f_theta_ujmle <- exp(-(opt_result_effect$par %*% omega_re %*% opt_result_effect$par)/2)
        
        # Probability for the current row
        beta_current <- beta_i
        beta_current[1:n_random_effect] <- beta_current[1:n_random_effect]  + opt_result_effect$par
        update_response_probabilities_inplace(prob, X, beta_current, nChoiceSet)

        # Probability for the "observed" Y conditioned on u (maximized for this row)
        prob_Y_cond_ujmle <- prod(prob^Y_i) 
        
        ratio_term <- (opt_result_effect$par[i_re]^2 + constant_term * sd_vec[i_re]^3) / sd_vec[i_re]^3
        re_ratio_condexpect_numerator[i_re] <- sqrt(H_j_inv_det) * ratio_term * prob_Y_cond_ujmle * f_theta_ujmle
        
      }
      
      #
      #---- E_u(p | y)
      #
      
      Ep_given_y_unnormalized <- numerator / denom
      Ep_given_y <- normalizeProbabilities(Ep_given_y_unnormalized, n_choice_set)
      
      #
      #---- E_u( f(u, sig) | y)
      #
     
      Eratio_given_y_uncorrected <- re_ratio_condexpect_numerator / denom - constant_term
      Eratio_given_y             <- abs(Eratio_given_y_uncorrected)
      
      #
      #---- Variance-Covariance
      #
      
      score1 <- matrix( t(X) %*% (Y_i - Ep_given_y) )
      #score2 <- Eratio_given_y + expm::sqrtm(omega_re)
      score2 <- matrix( Eratio_given_y - sqrt(diag(omega_re)) )
      
      I11i = score1 %*% t(score1)
      I12i = score1 %*% t(score2)
      I22i = score2 %*% t(score2);
      
      # Add to current state of information matrix
      I11 = I11 + I11i
      I12 = I12 + I12i
      I22 = I22 + I22i
      
      # increment counter for successful runs
      n_converged <- n_converged + 1
      
    }
    
  }
  
  variance_covariance <- rbind( cbind(I11, I12), cbind(t(I12), I22))
  variance_covariance <- variance_covariance / n_converged
  
  return(variance_covariance)
}



if (0 == 1){
  
  library(PanelMixedDesign)
  library(DoE.base)
  
  ## PRELIM
  nattr <- 3
  nChoiceSet <- 6
  nAlternative <- 2 
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.3, 0.4)
  st = NULL
  fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
  for(i in 1:nChoiceSet){
    st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=TRUE) )
  }
  designm=(fullfac)[st,]
  contr=rep('contr.sum',nattr)
  names(contr)=names(designm)
  contr=as.list(contr)
  M = model.matrix(~., designm, contrasts = contr)[,-1] 
  
  # ## FUNCTION ARGUMENTS
  # ny = 100
  # X = M
  # b_mean = mu
  # var_vec = sig
  # n_choice_set = nChoiceSet
  
  
  test_laplace(ny = 1000, X = M, b_mean = mu, var_vec = sig, n_choice_set = nChoiceSet)
  
  
  x2 <- PMLInfoApprox(X = M, method = "Laplace", nChoiceSet = nChoiceSet,
                            effectMean = mu,
                            effectVar = sig,
                            opts = list(nY=10000,
                                        nU=1000) )
  det(x2)
  sum(diag(solve(x2)))
  
}
