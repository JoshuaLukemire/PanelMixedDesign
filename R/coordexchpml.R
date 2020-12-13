#' Search for an optimal design under the panel mixed logit model
#' 
#' @param nChoiceSet Number of choice sets
#' @param nAlternative The number of alternatives in each choice set
#' @param nLevelAttribute Vector where each element is the number of levels for a factor
#' @param effectMean Vector of means for the effects coded attribute effects
#' @param effectVar Vector of variances for the effects coded attribute effects
#' @param approx Type of approximation to use (MQL, PQL, Importance, Laplace, or MSM)
#' @param nStartDesign Number of times to run the coordinate exchange algorithm
#' @param optcrit Optimality criterion. Either "A" or "D"
#' @return Object containing the best design found by the coordinate exchange algorithm
coordexchpml <- function(nChoiceSet,
                  nAlternative,
                  nLevelAttribute,
                  effectMean, 
                  effectVar,
                  approx = "pql",
                  nStartDesign = 10,
                  optcrit = "D"){
  
  cs=rep(nAlternative, nChoiceSet)           # a vector for numbers of alternatives in choice sets, used 
  ns=length(cs)         # number of choice sets reading from the length of cs
  nattr = length(nLevelAttribute)
  dimb <- length(effectMean)
  dimbr <- sum(effectVar > 0)
  
  design  <- array(0, c(nStartDesign, sum(cs), nattr)) # optimal designs from nStartDesign starting designs
  designc <-array(0,c(nStartDesign, sum(cs), dimb))    # coded values of design 
  derr <- rep(0, nStartDesign)                         # the determinant of the nStartDesign optimal designs 
  count=0                                
  
  do = Inf
  
  # Y is only used if using importance sampling (not rec.)
  Y <- NULL
  if (pracma::strcmpi(approx, "importance")){
    Y <- gen_all_choice_seq(nChoiceSet, nAlternative)
  }
  
  searchTimes <- rep(0, nStartDesign)
  
  for(r in 1:nStartDesign){
    # starting time
    startTime = tictoc::tic()
    
    st=NULL
    #fullfac=fac.design(nlevels=nLevelAttribute,random=F)
    fullfac <- purrr::quietly(DoE.base::fac.design)(nlevels=nLevelAttribute, random=F)$result
    for(i in 1:nChoiceSet){
      st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=T) )
    }
    designm=(fullfac)[st,]
    
    contr=rep('contr.sum',nattr)
    names(contr)=names(designm)
    contr=as.list(contr)
    
    modmat=stats::model.matrix(~.,designm,contrasts = contr)[,-1]  #contr is used to get effects type coding, current coded design
    
    Ic = PMLInfoApprox(modmat,
                       method = approx,
                       nChoiceSet,
                       effectMean = effectMean,
                       effectVar = effectVar, 
                       Y = Y)
    
    Ic=round(Ic,digits=10)
    
    if(all(is.na(Ic))==T){
      dc = Inf
    } else{
      # Check if less than full rank design (failure)
      if(Matrix::rankMatrix(Ic)<nrow(Ic) || isSymmetric(Ic)==F){
        dc = Inf
      } else{
        if (optcrit == "A"){
          dc=sum(eigen(solve(Ic))$values) / nrow(Ic)
          if(dc < 0){
            dc = Inf
          }
        } else {
          if(det(Ic)<0){
            dc=Inf
          }else{
            dc = abs(det(Ic))^{-1/nrow(Ic)}
          }
        }
      } # end of else statement for valid (non-nan) Ic value
    } # end of check that not all na
    
    new=matrix()              # new design  
    
    m=1
    while(m!=0){               # if no exchange is made, then m=0
      n=0                                    # number of exchange
      for(i in 1:(sum(cs)*nattr)){            # i goes through all elements in the uncoded design matrix
        
        j=(i%%nattr)                       # column in uncoded design matrix, i.e., jth attribute
        if(j==0) {j=nattr}                  
        k=(i-j)/nattr+1                    # row in uncoded design matrix, i.e., kth row
        ch=ceiling(k/nAlternative)                    # the 'ch'th choice set 
        diff=setdiff(1:nLevelAttribute[j],designm[k,j]) # possible levels for exchange
        for(l in diff){
          new=designm
          new[k,j]=l                                          # uncoded design matrix after exchange
          modmatnew=stats::model.matrix(~.,new,contrasts=contr)[,-1] # coded matrix of new        result1 = fi(modmat)
          
          I1 = PMLInfoApprox(modmatnew,
                             method = approx,
                             nChoiceSet,
                             effectMean = effectMean,
                             effectVar = effectVar, 
                             Y = Y)
          I1=round(I1,digits=10)
          
          if(all(is.na(I1))==T){
            d1 = Inf
          } else{
            if(Matrix::rankMatrix(I1)<nrow(I1) || isSymmetric(I1)==F){
              d1 = Inf
            } else{
              if (optcrit == "A"){
                d1=sum(eigen(solve(I1))$values)/nrow(I1)
                if(d1 < 0){
                  d1 = Inf
                }
              } else {
                if(det(I1)<0){
                  d1=Inf
                }else{
                  d1 = abs(det(I1))^{-1/(nrow(I1))}
                }
              }
            }						     	 
          }
          
          if (d1<dc){
            designm = new 
            modmat=modmatnew
            Ic=I1
            dc=d1
            n=n+1       # exchange is kept, add 1 to number of change 
          }
          # Handle case where new design is equally as good as old design
          if(d1 == dc){ 
            u = stats::runif(1,0,1)
            if(u < 0.5){
              designm=new 
              modmat=modmatnew
              Ic=I1
              dc=d1
              n=n+1       # exchange is kept, add 1 to number of change 
            }
          }       		                            
          count=count+1
        } #l loop end
        #print(dc)
      } # end of loop over i
      m=n
    } # end of while loop
    design[r,,]=as.matrix(designm)
    designc[r,,]=modmat
    derr[r]=dc
    if(dc<do){
      Io=Ic        # information matrix of the optimal design
      opt=designm  # optimal design
      mo=modmat    # coded optimal design
      do=dc        # determinant of optimal design
    }
    endTime = purrr::quietly(tictoc::toc)()$result
    searchTimes[r] = as.numeric(endTime$toc)
  }# r
  
  return(list(
    designs      = design,
    modmats      = designc,
    criteriavals = derr,
    optdesign    = opt,
    optInfo      = Io,
    optModmat    = mo,
    optCriteria  = do,
    searchTimes  = searchTimes
  ))
  
} # end of cepml function
