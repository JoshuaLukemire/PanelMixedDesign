library(PanelMixedDesign)
library(DoE.base)
context("Laplace Tests")

test_that("Laplace Runs for random effects only", {
  nattr <- 3
  n_choice_set <- 6
  nAlternative <- 2 
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.3, 0.4)
  st = NULL
  fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
  for(i in 1:n_choice_set){
    st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=TRUE) )
  }
  designm=(fullfac)[st,]
  contr=rep('contr.sum',nattr)
  names(contr)=names(designm)
  contr=as.list(contr)
  M = model.matrix(~., designm, contrasts = contr)[,-1] 
  infoAppr <- varcov_approx_PML(M, method = "laplace",
                            n_choice_set = n_choice_set, effect_mean = mu, effect_vars = sig)
  expect_true(is.matrix(infoAppr))
})


test_that("Laplace Runs for random effects + fixed effects", {
  nattr <- 3
  n_choice_set <- 6
  nAlternative <- 2 
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.0, 0.4)
  st = NULL
  fullfac <- purrr::quietly(fac.design)(nlevels=nLevelAttribute, random=FALSE)$result
  for(i in 1:n_choice_set){
    st = c( st, sample(prod(nLevelAttribute), nAlternative, replace=TRUE) )
  }
  designm=(fullfac)[st,]
  contr=rep('contr.sum',nattr)
  names(contr)=names(designm)
  contr=as.list(contr)
  M = model.matrix(~., designm, contrasts = contr)[,-1] 
  infoAppr <- varcov_approx_PML(M, method = "laplace",
                            n_choice_set = n_choice_set, effect_mean = mu, effect_vars = sig)
  expect_true(is.matrix(infoAppr))
})