#' Generate all possible response sequences for a DCE
#' 
#' @param nChoiceSet Number of choice sets
#' @param nAlternative Number of alternatives per choice set
#' @return An array with all possible response sequences to the experiment
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