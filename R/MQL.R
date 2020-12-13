#' Calculate choice probabilities
#' 
#' @param b Vector of effects
#' @param modmat The model matrix
#' @param qes A vector of alternative ids (eg 1, 2, 3, 1, 2, 3, 1, 2, 3)
#' @return The probability of each alternative being selected in each choice set
cprob<-function(b,modmat, qes)                        # each column of b is a vector of parameter values(nrow = dimb), modmat is the coded design matrix (sum(cs)*dimb)
{
  rep=exp(modmat%*%b)                                # each column of rep is a vector of representative utilities of alternatives, nrow = sum(cs)
  sumrep=apply(rep,2,function(x) stats::ave(x,qes,FUN=sum)) # each column of sumrep is a vector of sum of utilities in each choice set, nrow = ns
  p=rep/sumrep                                       # each column of p is a vector of choice probabilities of alternatives, nrow = ns
  return(p)                                          
}



#' Marginal Quasi Likelihood approach to approximate the information matrix for the
#' model parameters
#' 
#' @param modmat The model matrix
#' @param effectMean Vector of means for the effects coded attribute effects
#' @param effectVar Vector of variances for the effects coded attribute effects
#' @param nChoiceSet Number of choice sets
#' @return The approximation to the information matrix for the model parameters
MQLApprox <- function(modmat,
                      effectMean,
                      effectVar,
                      nChoiceSet){  
  
  dimb <- length(effectMean)
  dimbr <- sum(effectVar > 0)
  sigmar <- effectVar[1:dimbr]
  
  n_alternative <- nrow(modmat) / nChoiceSet
  
  modmatd = matrix(0, nChoiceSet*(n_alternative-1), dimb)
  S = nChoiceSet
  J = n_alternative
  delta = matrix(0,S*J,S*J)
  cs=rep(n_alternative, nChoiceSet)           
  ns=length(cs) 
  ind2 = seq(J,by=J,length=S)
  qes<- rep(factor(1:ns),times=cs) # a vector of ids of the alternatives in each choice set. 
  for(i in 1:nChoiceSet){
    
    modmatd[((i-1)*(J-1)+1):(i*(J-1)),] = modmat[((i-1)*J+1):((i-1)*J+J-1),]-modmat[i*J,]
    
  }
  I = matrix(0,dimb+dimbr,dimb+dimbr)
  I11 = 0
  I22 = 0
  u = 0
  p = cprob(effectMean,modmat, qes) # p evaluated at tilde(u)
  for( j in 1:nChoiceSet){
    st = (j-1)*n_alternative+1
    en = j*n_alternative
    pos = st:en
    delta[pos,pos]=diag(p[pos,1])-p[pos,1]%*%t(p[pos,1])
  }
  v = solve(delta[setdiff(1:(nChoiceSet*n_alternative),ind2),setdiff(1:(nChoiceSet*n_alternative),ind2)])+modmatd%*%diag(effectVar)%*%t(modmatd)  # vn in the formula
  vi = solve(v)
  vix = vi%*%modmatd
  I11 = I11 + t(modmatd)%*%vix
  I22 = 	I22 + 2*diag(sqrt(effectVar))%*%(t(modmatd)%*%vi%*%modmatd)^2%*%diag(sqrt(effectVar))
  I[1:dimb,1:dimb] = I11
  I[(dimb+1):(dimb+dimbr),(dimb+1):(dimb+dimbr)] = I22[1:dimbr,1:dimbr]
  return(I)
}