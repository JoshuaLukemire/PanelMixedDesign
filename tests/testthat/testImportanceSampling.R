library(PanelMixedDesign)
library(DoE.base)
context("Importance Sampling Tests")

test_that("Importance Sampling Runs for random effects only", {
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
  Y <- gen_all_choice_seq(nChoiceSet, nAlternative)
  infoAppr <- PMLInfoApprox(M, method = "importance",
                            nChoiceSet = nChoiceSet, effectMean = mu,
                            effectVar = sig,
                            Y = Y)
  expect_true(is.matrix(infoAppr))
})


test_that("Importance Sampling Runs for random effects + fixed effects", {
  nattr <- 3
  nChoiceSet <- 6
  nAlternative <- 2 
  nLevelAttribute <- c(2, 2, 2)
  mu <- c(1.0, -0.4, -0.8)
  sig <- c(0.5, 0.0, 0.4)
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
  Y <- gen_all_choice_seq(nChoiceSet, nAlternative)
  infoAppr <- PMLInfoApprox(M, method = "Importance",
                            nChoiceSet = nChoiceSet, effectMean = mu,
                            effectVar = sig,
                            Y = Y)
  expect_true(is.matrix(infoAppr))
})