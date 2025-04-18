#' Approximate the variance covariance matrix for the model parameters
#' 
#' @param X The model matrix for the choice experiment. Each row is an effects coded alternative.
#' @param method Choice of variance-covariance approximation method. Options are MQL, PQL, MSM, Laplace, or Importance. MQL is the default.
#' @param n_choice_set Number of choice sets
#' @param effect_means Vector of means for the effects coded attribute effects
#' @param effect_vars Vector of variances for the effects coded attribute effects. Must be equal in length to the effect_means argument. Use a 0 for any non-random terms.
#' @param PQL_n_samples Number of random samples to use in PQL evaluation
#' @param MSM_n_samples Number of random samples to use in MSM evaluation
#' @param Laplace_n_samples Number of random samples to use in the Laplace evaluation
#' @param IS_n_samples_u Number of random samples of random effects to use in importance sampling
#' @param IS_n_samples_Y Number of random samples of the response vector to use in importance sampling
#' @param Y an optional enumeration of all possible responses. This is only used for importance sampling. If provided, importance sampling will use the full enumeration of outcomes instead of randomly sampling responses.
#' @return The approximation to the variance covariance matrix for the model parameters.
#' @examples
#' # Libraries
#' library(DoE.base)
#' 
#' # 3 attributes
#' nattr <- 3
#' 
#' # 6 choice sets
#' n_choice_set <- 6
#' 
#' # 2 alternatives per choice set
#' n_alternative <- 2 
#' 
#' # 3 factors, each at 2 levels
#' n_level_attribute <- c(2, 2, 2)
#' 
#' # effect means
#' mu <- c(1.0, -0.4, -0.8)
#' 
#' # effect variances
#' sig <- c(0.5, 0.3, 0.4)
#' 
#' # indices for alternatives
#' st = NULL
#' 
#' # all possible design points
#' fullfac <- purrr::quietly(fac.design)(nlevels=n_level_attribute, random=FALSE)$result
#' 
#' # Generate a random design
#' for(i in 1:n_choice_set){
#'     st = c( st, sample(prod(n_level_attribute), n_alternative, replace=TRUE) )
#' }
#' designm=(fullfac)[st,]
#' 
#' # Effects coding
#' contr=rep('contr.sum',nattr)
#' names(contr)=names(designm)
#' contr=as.list(contr)
#' M = model.matrix(~., designm, contrasts = contr)[,-1] 
#' 
#' varcovAppr <- varcov_approx_PML(M, method = "PQL",
#'     n_choice_set = n_choice_set, effect_means = mu, effect_vars = sig) 
#' 
#' @export
varcov_approx_PML <- function(X,
                          method = "MQL",
                          n_choice_set,
                          effect_means,
                          effect_vars,
                          PQL_n_samples = 1000,
                          MSM_n_samples = 1000,
                          Laplace_n_samples = 1000,
                          IS_n_samples_u = 1000,
                          IS_n_samples_Y = 1000,
                          Y = NULL)
{
  
  # Verify that dimensions match
  if (length(effect_means) != length(effect_vars)){
    stop(paste("Length of means:", length(effect_means),
               "does not match length of effect_vars:",
               length(effect_vars)))
  }
  
  # Check if the model matrix is (clearly) the wrong dimension
  if (nrow(X) / n_choice_set != round(nrow(X) / n_choice_set)){
    error_string <- paste0("Number of rows in the model matrix (",
                           nrow(X), ") must be a multiple of the number of choice sets (",
                           n_choice_set, "). Please check input.")
    stop(error_string)
  }
  
  # Check and warn the user if there are any duplicated alternatives within
  # each choice set
  n_alternative <- nrow(X) / n_choice_set
  for (i_choice_set in 1:n_choice_set){
    indices <- ((i_choice_set-1)*n_alternative+1):(i_choice_set*n_alternative)
    if (any(duplicated(X[indices,]))){
      print(X)
      warning(paste0("Choice set ", i_choice_set, " has duplicated alternatives."))
      stop(paste0("Choice set ", i_choice_set, " has duplicated alternatives."))
    }
  }
  
  # Re-arrange effect means and vars so that the random effects come first.
  # this is just for convenience; the sorting will be reversed when
  # the approximation is complete during construction of variance_covariance_final
  random_effect_inds  <- which(effect_vars != 0)
  n_random_effect     <- length(random_effect_inds)
  fixed_effect_inds   <- (1:length(effect_means))[-random_effect_inds]
  b_mean              <- c(effect_means[random_effect_inds], effect_means[fixed_effect_inds])
  b_var               <- c(effect_vars[random_effect_inds],  effect_vars[fixed_effect_inds])
  X                   <- X[,c(random_effect_inds, fixed_effect_inds)]
  
  
  variance_covariance_approx <- NULL
  
  if (pracma::strcmpi(method, "Laplace")){
    
    variance_covariance_approx <- Laplace_approx_vcov(Laplace_n_samples, X, b_mean, b_var, n_choice_set)
    
  } else if (pracma::strcmpi(method, "MQL")){
    
    variance_covariance_approx <- MQL_approx_vcov(X, b_mean, b_var, n_choice_set)
    
  } else if (pracma::strcmpi(method, "PQL")){
    
    variance_covariance_approx <- PQL_approx_vcov(PQL_n_samples, X, b_mean, b_var, n_choice_set)
    
  } else if (pracma::strcmpi(method, "Importance")){
    
    if (is.null(Y)){
      variance_covariance_approx <- importance_sample_random_Y(X, b_mean, b_var, n_choice_set, IS_n_samples_u, IS_n_samples_Y)
    } else {
      variance_covariance_approx <- importance_sample_all_Y(Y, X, b_mean, b_var, n_choice_set, IS_n_samples_u)
    }
    
  } else if (pracma::strcmpi(method, "MSM")){
    
    variance_covariance_approxTemp <- MSM_approx_vcov(MSM_n_samples, X, b_mean, b_var, n_choice_set)
    
    # Finalize quantities from MSM
    f <- variance_covariance_approxTemp[1:(nrow(variance_covariance_approxTemp) / 2), ]
    d <- variance_covariance_approxTemp[(nrow(variance_covariance_approxTemp) / 2 +1):nrow(variance_covariance_approxTemp),]
    
    # Return matrix of NA if numeric issues occur
    variance_covariance_approx <- tryCatch(
      expr = {
        solve(t(f) %*% MASS::ginv(d) %*% f)
      },
      error = function(e){
        matrix(NA, nrow = ncol(f), ncol = ncol(f))
      }
    )
    
    
  } else {
    
    stop(paste("Approximation Method:", method, "not recognized."))
    
  }
  
  # Undo the re-arrangement by random effects
  variance_covariance_final <- matrix(0, nrow=length(b_mean) + n_random_effect, ncol=length(b_mean) + n_random_effect)

  variance_covariance_final[c(random_effect_inds, fixed_effect_inds, (length(b_mean)+1):(length(b_mean) + n_random_effect) ),
                  c(random_effect_inds, fixed_effect_inds, (length(b_mean)+1):(length(b_mean) + n_random_effect) )] <- 
    variance_covariance_approx[ c(1:length(b_mean), 1:n_random_effect + length(b_mean)), c(1:length(b_mean), 1:n_random_effect + length(b_mean))]
  
  return( variance_covariance_final )
  
}