#' Search for an optimal design under the panel mixed logit model
#' 
#' @param n_choice_set Number of choice sets
#' @param n_alternative The number of alternatives in each choice set
#' @param n_level_attribute Vector where each element is the number of levels for a factor. For example, for two 2-level factors and one 3-level factor this would contain c(2, 2, 3).
#' @param effect_mean Vector of means for the effects coded attribute effects
#' @param effect_vars Vector of variances for the effects coded attribute effects
#' @param approx Type of approximation to use (MQL, PQL, Importance, Laplace, or MSM)
#' @param n_run Number of times to run the coordinate exchange algorithm
#' @param optcrit Optimality criterion. Either "A" or "D"
#' @return List object containing the following results of the coordinate exchange algorithm: \tabular{ll}{
#'    \code{design_list} \tab A list object with the final design from each coordinate exchange run \cr
#'    \tab \cr
#'    \code{modmat_list} \tab  A list object with the final model matrix from each coordinate exchange run \cr
#'    \tab \cr
#'    \code{criteria_list} \tab  A list object with the optimality criteria value for the final design from each coordinate exchange run \cr
#'    \tab \cr
#'    \code{best_design} \tab  The best design found by the coordinate exchange algorithm \cr
#'    \tab \cr
#'    \code{best_design_varcov} \tab  The variance-covariance matrix for the model parameters for the best design found by the coordinate exchange algorithm. The form will depend on which approximation method was used. \cr
#'    \tab \cr
#'    \code{best_model_matrix} \tab  The model matrix corresponding to the best design found by the coordinate exchange algorithm \cr
#'    \tab \cr
#'    \code{best_criteria_value} \tab  The value of the optimality criteria for the best design found by the coordinate exchange algorithm \cr
#'    \tab \cr
#'    \code{computation_times} \tab  The computation time for each run of the coordinate exchange algorithm \cr
#'    \tab \cr
#'    \code{n_accepted_exchange_list} \tab  The number of exchanges accepted during each run of the coordinate exchange algorithm. Can be used to verify that the algorithm is successfully moving away from the initial design \cr
#'    \tab \cr
#'    \code{n_proposal_evaled_list} \tab  The number of exchanges proposed during each run of the coordinate exchange algorithm \cr
#'    \tab \cr
#'    \code{overall_search_time} \tab  Overall time required to search for a design across all coordinate exchange algorithm runs \cr
#' }
#' @export
coord_exchange_PML <- function(n_choice_set,
                  n_alternative,
                  n_level_attribute,
                  effect_mean, 
                  effect_vars,
                  approx = "MQL",
                  n_run = 10,
                  optcrit = "D"){
  
  # Log the overall starting time
  overall_start_time <- Sys.time()
  
  # Number of alternatives per choice set
  cs=rep(n_alternative, n_choice_set)           # a vector for numbers of alternatives in choice sets, used 
  
  # number of choice sets reading from the length of cs
  n_choice_set <- length(cs)
  
  # Number of attributes
  nattr = length(n_level_attribute)
  
  # Total Number of coefficients (fixed and/or random)
  dimb <- length(effect_mean)
  
  # Number of coefficients that are random
  dimbr <- sum(effect_vars > 0)
  
  # Array storing designs found by each run of the coordinate exchange algorithm
  design  <- array(0, c(n_run, sum(cs), nattr)) # optimal designs from n_run starting designs
  
  # Array storing the model matrix for each run of the coordinate exchange algorithm
  designc <-array(0,c(n_run, sum(cs), dimb))    # coded values of design 
  
  # TODO remove this
  derr <- rep(0, n_run)                         # the determinant of the n_run optimal designs 
  count=0                                
  
  # Best criterion value across all CE runs (starting designs)
  best_criterion   <- Inf
  # Criterion-values specific to a current CE run (n_run of these)
  stored_criterion <- Inf
  criterion        <- Inf
  
  # Y is only needed if using importance sampling (not rec.)
  Y <- NULL
  if (pracma::strcmpi(approx, "importance")){
    Y <- gen_all_choice_seq(n_choice_set, n_alternative)
  }
  
  # Track how long each search took
  searchTimes <- NULL
  
  # How many exchanges were accepted for each run
  n_accepted_exchange_list <- rep(0, n_run)
  
  # How many proposals had their var-cov evaluated for each run
  n_proposal_evaled_list <- rep(0, n_run)
  
  # Initialize storage for optimal design
  opt <- NULL
  mo  <- NULL
  Io  <- NULL
  
  for(r in 1:n_run){
    
    # starting time
    startTime = tictoc::tic()
    
    # Reset criteria
    stored_criterion <- Inf
    criterion        <- Inf
    
    # Running tracker of how many exchanges were accepted for this run
    n_accepted_exchange_run <- 0
    
    # Running tracker of how many proposals were actually evaluated
    # (does not include invalid proposals)
    n_proposal_evaled <- 0
    
    # Obtain a starting design
    st=NULL
    fullfac <- purrr::quietly(DoE.base::fac.design)(nlevels=n_level_attribute, random=F)$result
    # If there are enough unique alternatives, use a starting design with no
    # repeated alternatives. Otherwise randomly select one choice set at a time
    if (nrow(fullfac) >= n_choice_set*n_alternative){
      st <- sample(1:nrow(fullfac), n_choice_set*n_alternative)
    } else {
      for(i in 1:n_choice_set){
        st = c( st, sample(prod(n_level_attribute), n_alternative, replace=FALSE) )
      }
    }
    designm=(fullfac)[st,]
    
    # Setup effects coding
    contr        <- rep('contr.sum', nattr)
    names(contr) <- names(designm)
    contr        <- as.list(contr)

    # Get the model matrix (effects coded)
    modmat=stats::model.matrix(~., designm, contrasts = contr)[,-1]
    
    # Calculate the variance covariance matrix
    VarCovMat = varcov_approx_PML(modmat,
                       method = approx,
                       n_choice_set,
                       effect_mean = effect_mean,
                       effect_vars = effect_vars, 
                       Y = Y)
    
    # Evaluate optimality criteria, with checks for failed designs
    if(any(is.na(VarCovMat)) == TRUE){
      stored_criterion = Inf
    } else{
      # Check if less than full rank design (failure)
      if(Matrix::rankMatrix(VarCovMat) < nrow(VarCovMat) || isSymmetric(VarCovMat) == F){
        stored_criterion = Inf
      } else{
        if (optcrit == "A"){
          stored_criterion <- sum(eigen(VarCovMat)$values) / nrow(VarCovMat)
          if(stored_criterion < 0){
            stored_criterion = Inf
          }
        } else {
          if(det(VarCovMat) < .Machine$double.eps){
            stored_criterion = Inf
          }else{
            # Determinant of the variance covariance matrix
            VarCov_det   <- det(VarCovMat)
            stored_criterion    <- VarCov_det^(1 / nrow(VarCovMat))
          }
        }
      } # end of else statement for valid (non-nan) VarCovMat value
    } # end of check that not all na
    
    # Now that the starting design has been constructed and evaluated, start 
    # proposing single coordinate exchanges
    
    new_design = matrix()
    
    # Number of full cycles through the design
    n_cycle = 0
    
    # convergence takes place when no additional coordinate exchanges
    # result in an improvement in the design
    converged = FALSE
    
    while(!converged){               
      
      n_accepted_exchange = 0  # number of exchange
      n_cycle = n_cycle + 1
      
      # Loop over each possible attribute setting
      for(i in 1:(sum(cs)*nattr)){
        
        # Determine which column we are currently considering for an exchange
        j <- i %% nattr               
        if(j==0) {j <- nattr} # case where is final column   
        
        # Determine which row we are currently considering for an exchange
        k <- (i - j) / nattr + 1                    
        
        # Figure out which choice set this corresponds to
        ch <- ceiling(k / n_alternative)                   
        
        # Find the possible exchanges for this location
        possible_exchanges <- setdiff(1:n_level_attribute[j], designm[k,j])
        
        # Loop over available settings for this location, assess improvement 
        # at each possible value
        for(l in possible_exchanges){
          
          # Make a copy of the current design
          new_design <- designm
          
          # Plug in the new setting
          new_design[k,j] <- l      
          
          # Verify that the corresponding choice set does not include
          # duplicate choices. If it does, move on to next
          current_choice_set_index <- ceiling(k / n_alternative)
          current_choice_set       <- new_design[
            ((current_choice_set_index-1)*n_alternative+1):(current_choice_set_index*n_alternative),]
          if (any(duplicated(current_choice_set))){next}
          
          # Check that flipping to this setting does not cause a loss of ability to
          # estimate the corresponding parameter (due to it being the same wn each 
          # choice set). First, only need to check this if it is a potential
          # problem within this choice set
          uniqueness_problem <- TRUE
          if (all(current_choice_set[,j] == l)){
            # If it is a potential problem, run a comparison against all
            # other choice sets
            for (iCS in 1:n_choice_set){
              current_choice_set_vals <- new_design[
                ((iCS-1)*n_alternative+1):(iCS*n_alternative), j] 
              if (length(unique(current_choice_set_vals)) > 1){
                uniqueness_problem <- FALSE
              }
            }
          }
          if (uniqueness_problem == TRUE){next}
          
          
          # Get the corresponding model matrix conversion
          modmatnew=stats::model.matrix(~.,new_design,contrasts=contr)[,-1] # coded matrix of new        result1 = fi(modmat)
          
          # Get the variance-covariance matrix
          var_cov_matrix_proposal = varcov_approx_PML(modmatnew,
                             method = approx,
                             n_choice_set,
                             effect_mean = effect_mean,
                             effect_vars = effect_vars, 
                             Y = Y)
          
          # Number of evaluated proposals
          n_proposal_evaled <- n_proposal_evaled + 1

          if(any(is.na(var_cov_matrix_proposal))==T){
            new_criterion = Inf
          } else{
            if(Matrix::rankMatrix(var_cov_matrix_proposal) < nrow(var_cov_matrix_proposal) || isSymmetric(var_cov_matrix_proposal)==F){
              new_criterion = Inf
            } else{
              if (optcrit == "A"){
                new_criterion=sum(eigen(solve(var_cov_matrix_proposal))$values)/nrow(var_cov_matrix_proposal)
                if(new_criterion < 0){
                  new_criterion = Inf
                }
              } else {
                if(det(var_cov_matrix_proposal) < .Machine$double.eps){
                  new_criterion=Inf
                }else{
                  new_criterion = abs(det(var_cov_matrix_proposal))^{1/(nrow(var_cov_matrix_proposal))}
                }
              }
            }						     	 
          }
          
          if (new_criterion < stored_criterion){
            designm   <- new_design 
            modmat    <- modmatnew
            VarCovMat <- var_cov_matrix_proposal
            stored_criterion    <- new_criterion
            n_accepted_exchange <- n_accepted_exchange + 1      
          }
          # Handle case where new design is equally as good as old design
          if(new_criterion == stored_criterion){ 
            u = stats::runif(1,0,1)
            if(u < 0.5){
              designm=new_design 
              modmat=modmatnew
              VarCovMat=var_cov_matrix_proposal
              stored_criterion    <- new_criterion
              n_accepted_exchange <- n_accepted_exchange + 1      
            }
          }       		                            
          count=count+1
        } #l loop end
      } # end of loop over i
      
      # Check if we looped through the entire design without accepting a single
      # proposed change
      if (n_accepted_exchange == 0){
        converged = TRUE
      }
      
      n_accepted_exchange_run <- n_accepted_exchange_run + n_accepted_exchange
      
    } # end of while loop
    
    # Store 
    design[r,,]=as.matrix(designm)
    designc[r,,]=modmat
    derr[r]=stored_criterion
    n_accepted_exchange_list[r] <- n_accepted_exchange_run
    n_proposal_evaled_list[r]   <- n_proposal_evaled
    if(stored_criterion < best_criterion){
      Io=VarCovMat        # variance covariance matrix of the optimal design
      opt=designm  # optimal design
      mo=modmat    # coded optimal design
      best_criterion=stored_criterion        # determinant of optimal design
    }
    endTime = purrr::quietly(tictoc::toc)()$result
    searchTimes[r] = endTime$callback_msg
  }# r
  
  # Log the overall starting ending
  overall_end_time <- Sys.time()
  overall_time <- difftime(overall_end_time, overall_start_time, units = "hours")
  
  return(list(
    design_list         = design,
    modmat_list         = designc,
    criteria_list       = derr,
    best_design         = opt,
    best_design_varcov  = Io,
    best_model_matrix   = mo,
    best_criteria_value = best_criterion,
    computation_times   = searchTimes,
    n_accepted_exchange_list = n_accepted_exchange_list,
    n_proposal_evaled_list = n_proposal_evaled_list,
    overall_search_time = overall_time
  ))
  
} # end of cepml function
