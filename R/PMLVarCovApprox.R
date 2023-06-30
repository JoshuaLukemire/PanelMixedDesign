#' Approximate the variance covariance matrix for the model parameters
#' 
#' @param X model matrix (a 3d array).
#' @param method Choice of approximation method. Options are MQL, PQL, MSM, Laplace, or Importance
#' @param nChoiceSet Number of choice sets
#' @param effectMean Vector of means for the effects coded attribute effects
#' @param effectVar Vector of variances for the effects coded attribute effects
#' @param opts Additional options controlling the number of samples used in evaluation
#' @param Y an enumeration of all possible responses. This is only required for importance sampling
#' @return The approximation to the variance covariance matrix for the model parameters
#' @examples
#' # Libraries
#' library(DoE.base)
#' 
#' # 3 attributes
#' nattr <- 3
#' 
#' # 6 choice sets
#' nChoiceSet <- 6
#' 
#' # 2 alternatives per choice set
#' nAlternative <- 2 
#' 
#' # 3 factors, each at 2 levels
#' nLevelAttribute <- c(2, 2, 2)
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
#' fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
#' 
#' # Generate a random design
#' for(i in 1:nChoiceSet){
#'     st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=TRUE) )
#' }
#' designm=(fullfac)[st,]
#' 
#' # Effects coding
#' contr=rep('contr.sum',nattr)
#' names(contr)=names(designm)
#' contr=as.list(contr)
#' M = model.matrix(~., designm, contrasts = contr)[,-1] 
#' 
#' varcovAppr <- PMLVarCovApprox(M, method = "PQL",
#'     nChoiceSet = nChoiceSet, effectMean = mu, effectVar = sig) 
#' 

PMLVarCovApprox <- function(X,
                          method = "MQL",
                          nChoiceSet,
                          effectMean,
                          effectVar,
                          opts = list(nY=1000,
                                      nU=1000),
                          Y = NULL)
{
  
  # Verify that dimensions match
  if (length(effectMean) != length(effectVar)){
    stop(paste("Length of means:", length(effectMean),
               "does not match length of effectVar:",
               length(effectVar)))
  }
  
  # Calculate number of alternatives per choice set
  nAltPerChoiceSet <- nrow(X) / nChoiceSet
  
  # Re-arrange effect means and vars so that the random effects come first.
  # this is just for convenience; the sorting will be reversed when
  # the approximation is complete.
  randomEffectIndices <- which(effectVar != 0)
  nRandEffect         <- length(randomEffectIndices)
  fixedEffectIndices  <- (1:length(effectMean))[-randomEffectIndices]
  b_mean              <- c(effectMean[randomEffectIndices], effectMean[fixedEffectIndices])
  b_var               <- c(effectVar[randomEffectIndices],  effectVar[fixedEffectIndices])
  X                   <- X[,c(randomEffectIndices, fixedEffectIndices)]
  
  VarCovApprox <- NULL
  
  if (pracma::strcmpi(method, "Laplace")){
    
    VarCovApprox <- laplaceApproximation(opts$nY, X, b_mean, b_var, nChoiceSet)
    
  } else if (pracma::strcmpi(method, "MQL")){
    
    VarCovApprox <- MQLApprox(X, b_mean, b_var, nChoiceSet)
    
  } else if (pracma::strcmpi(method, "PQL")){
    
    VarCovApprox <- PQLApprox(opts$nU, X, b_mean, b_var, nChoiceSet)
    
  } else if (pracma::strcmpi(method, "Importance")){
    
    if (is.null(Y)){
      #stop("Argument: Y must be provided for importance sampling. Y can be generated
      #     usign the gen_all_choice_seq function.")
      VarCovApprox <- importanceSample(X, b_mean, b_var, nChoiceSet, opts$nU, opts$nY)
    } else {
      VarCovApprox <- importanceSampleFixedY(Y, X, b_mean, b_var, nChoiceSet, opts$nU)
    }
    
  } else if (pracma::strcmpi(method, "ImportanceAlt")){
    
    VarCovApprox <- importanceSampleAlt(X, b_mean, b_var, nChoiceSet, opts$nU)

  } else if (pracma::strcmpi(method, "MSM")){
    
    VarCovApproxTemp <- MSMApprox(opts$nU, X, b_mean, b_var, nChoiceSet)
    f <- VarCovApproxTemp[1:(nrow(VarCovApproxTemp) / 2), ]
    d <- VarCovApproxTemp[(nrow(VarCovApproxTemp) / 2 +1):nrow(VarCovApproxTemp),]
    VarCovApprox = t(f) %*% MASS::ginv(d) %*% f
    
  } else {
    
    stop(paste("Approximation Method:", method, "not recognized."))
    
  }
  
  # Undo the re-arrangement by random effects
  VarCovApproxFinal <- matrix(0, nrow=length(b_mean) + nRandEffect, ncol=length(b_mean) + nRandEffect)

  VarCovApproxFinal[c(randomEffectIndices, fixedEffectIndices, (length(b_mean)+1):(length(b_mean) + nRandEffect) ),
                  c(randomEffectIndices, fixedEffectIndices, (length(b_mean)+1):(length(b_mean) + nRandEffect) )] <- 
    VarCovApprox[ c(1:length(b_mean), 1:nRandEffect + length(b_mean)), c(1:length(b_mean), 1:nRandEffect + length(b_mean))]
  
  return( solve(VarCovApproxFinal) )
  
}