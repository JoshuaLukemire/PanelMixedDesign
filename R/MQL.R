#' Calculate choice probabilities for each alternative in each choice set
#' 
#' @param b Vector of effects
#' @param model_matrix The model matrix. Dimension should be 
#' @param alternative_ordering A vector of alternative ids (eg 1, 2, 3, 1, 2, 3, 1, 2, 3)
#' @return The probability of each alternative being selected in each choice set
#' @keywords internal
evaluate_choice_probabilities <- function(b, model_matrix, alternative_ordering){
  rep=exp(model_matrix %*% b)                               
  sumrep=apply(rep,2,function(x) stats::ave(x, alternative_ordering, FUN=sum)) 
  p=rep/sumrep                              
  return(p)                                          
}

#' Marginal Quasi-Likelihood approach to approximate the variance covariance matrix for the
#' model parameters
#' 
#' @param model_matrix The model matrix (effects coded), should be (n_alternative * n_choice_set) x n_beta
#' @param effect_mean Vector of means for the effects coded attribute effects
#' @param effect_var Vector of variances for the effects coded attribute effects
#' @param n_choice_set Number of choice sets
#' @return The approximation to the variance covariance matrix for the model parameters
#' @keywords internal
MQL_approx_vcov <- function(model_matrix,
                      effect_mean,
                      effect_var,
                      n_choice_set){  
  
  # Total number of effects
  dimb   <- length(effect_mean)
  
  # Number of random effects
  dimbr  <- sum(effect_var > 0)
  
  # Corresponding variances
  sigmar <- effect_var[1:dimbr]
  
  # Number of alternatives per choice set
  n_alternative <- nrow(model_matrix) / n_choice_set
  
  # vector containing number of alternatives per choice set, only needed
  # to simplify some sequence/vector commands later
  alternative_counts = rep(n_alternative, n_choice_set)   
  
  # index set used when later doing setdiff
  setdiff_indices = seq(n_alternative, by=n_alternative, length=n_choice_set)
  
  # a vector of ids of the alternatives in each choice set, e.g. 1, 2, 3, 1, 2, 3,...
  alternative_ordering <- rep(factor(1:n_choice_set), times=alternative_counts)
  
  # "differenced" version of the model matrix
  # This is obtained by taking each alternative EXCEPT the last and subtracting
  # the settings for the last alternative within each choice set
  model_matrix_differenced = matrix(0, n_choice_set*(n_alternative-1), dimb)
  
  # Difference within each choice set
  for(i in 1:n_choice_set){
    # Final alternative of this choice set
    reference_row <- matrix(model_matrix[i*n_alternative,], nrow = 1)
    # Subtract final alternative from all other alternatives in this choice set
    if (n_alternative == 2){
      model_matrix_differenced[((i-1)*(n_alternative-1)+1):(i*(n_alternative-1)),] = 
        model_matrix[((i-1)*n_alternative+1):((i-1)*n_alternative+n_alternative-1),] - reference_row
    } else {
      model_matrix_differenced[((i-1)*(n_alternative-1)+1):(i*(n_alternative-1)),] = 
        sweep(model_matrix[((i-1)*n_alternative+1):((i-1)*n_alternative+n_alternative-1),],
              2, reference_row)
    }
  }
  
  # Inverse of the variance-covariance matrix
  I = matrix(0,dimb+dimbr,dimb+dimbr)
  
  # This is Delta_n = diag(p) - pp'
  delta = matrix(0, nrow = n_choice_set*n_alternative, ncol = n_choice_set*n_alternative)
  # p evaluated at tilde(u), which is 0 for MQL
  p = evaluate_choice_probabilities(effect_mean, model_matrix, alternative_ordering)
  for (j in 1:n_choice_set){
    # Indices for this choice set
    cs_indices = ((j-1)*n_alternative+1):(j*n_alternative)
    delta[cs_indices,cs_indices] <- diag(p[cs_indices,1]) - p[cs_indices,1] %*% t(p[cs_indices,1])
  }
  
  # V = Delta^{-1} | u=0 + X Sigma X'
  delta_inverse <- solve(delta[setdiff(1:(n_choice_set*n_alternative), setdiff_indices),
                               setdiff(1:(n_choice_set*n_alternative), setdiff_indices)])
  Vn            <- delta_inverse + model_matrix_differenced %*% diag(effect_var) %*% t(model_matrix_differenced)
  Vn_inverse    <- solve(Vn)
  Vn_inverse_X  <- Vn_inverse %*% model_matrix_differenced
  
  # Upper left block of information matrix
  I11 <- t(model_matrix_differenced) %*% Vn_inverse_X 
  
  # Lower right block of information matrix
  I22 <- 0.5 * diag(sqrt(effect_var)) %*% I11^2 %*% diag(sqrt(effect_var))

  # # Variance Covariance Matrix via inversion of information matrix
  # vcov_mql                        <- matrix(0,dimb+dimbr,dimb+dimbr)
  # vcov_mql[1:dimb,1:dimb]         <- solve(I11)
  # vcov_mql[(dimb+1):(dimb+dimbr),
  #          (dimb+1):(dimb+dimbr)] <- solve(I22[1:dimbr,1:dimbr])
  
  vcov_mql <- tryCatch(
    expr = {
      vcov_mql                        <- matrix(0,dimb+dimbr,dimb+dimbr)
      vcov_mql[1:dimb,1:dimb]         <- solve(I11)
      vcov_mql[(dimb+1):(dimb+dimbr),
               (dimb+1):(dimb+dimbr)] <- solve(I22[1:dimbr,1:dimbr])
      vcov_mql
    },
    error = function(e){
      vcov_mql = matrix(NA, nrow = dimb + dimbr, ncol = dimb + dimbr)
    }
  )
  
  return(vcov_mql)
}