#' Generate all possible response sequences for a DCE. The resulting matrix can be used to conduct importance sampling to estimate the variance-covariance matrix for the parameters using a full enumeration of possible responses instead of sampling them.
#' 
#' @param nChoiceSet Number of choice sets
#' @param nAlternative Number of alternatives per choice set
#' @return An array with all possible response sequences to the experiment
#' @export
#' @examples
#' 
#' # DCE with 10 choice sets, and each choice set has 2 alternatives
#' n_choice_set <- 10
#' n_alternative <- 2
#' 
#' # Return a matrix where each column is a possible outcome of the experiment
#' Y <- gen_all_choice_seq(n_choice_set, n_alternative)
#' 
#' # The code below sets up a design matrix and converts it to an effects coded model matrix.
#' library(DoE.base)
#' 
#' # 3 attributes
#' nattr <- 3
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
#' # For this example, we do not already have a design, so we are going to
#' # construct a random design by first enumerating all possible alternatives
#' # and then randomly selecting a set of them for the experiment. 
#' # When inputting this quantity yourself, this should be a matrix of dimension
#' # (n_alternative * n_choice_set) x (n_attribute), where element i,j is the 
#' # integer-valued attribute setting for row i, attribute j. 
#' 
#' # Full enumeration of possible alternatives
#' all_possible_points <- purrr::quietly(fac.design)(nlevels=n_level_attribute, random=FALSE)$result
#' 
#' # Randomly select n_choice_set * n_alternative rows (this is the DCE)
#' for(i in 1:n_choice_set){
#'     st = c( st, sample(prod(n_level_attribute), n_alternative, replace=FALSE) )
#' }
#' design_matrix = (all_possible_points)[st,]
#' rownames(design_matrix) <- NULL # remove unneeded row labels
#' 
#' # Effects coding of the design matrix
#' contr=rep('contr.sum',nattr)
#' names(contr)=names(design_matrix)
#' contr=as.list(contr)
#' model_matrix = model.matrix(~., design_matrix, contrasts = contr)[,-1] 
#' 
#' varcovAppr <- varcov_approx_PML(model_matrix, method = "Importance",
#'     n_choice_set = n_choice_set, effect_means = mu, effect_vars = sig,
#'     Y = Y)
#' 
gen_all_choice_seq <- function(nChoiceSet, nAlternative){
  
  choice_sets <- rep(nAlternative, nChoiceSet)     
  suppressMessages(  full <- DoE.base::fac.design(nlevels = choice_sets))
  full = DoE.base::qua.design(full, quantitative = "all")
  full=as.matrix(full)
  
  y<-array(0, c(sum(choice_sets), prod(choice_sets)))
  
  for(i in 1:nChoiceSet)
  {
    for(j in 1:nAlternative) 
    { y[((i-1)*nAlternative+j),][full[,i]==j]=1 }
  }
  
  return(y)
}

.onUnload <- function (libpath) {
  library.dynam.unload("PanelMixedDesign", libpath)
}