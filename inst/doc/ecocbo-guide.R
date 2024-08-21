## ----setup, include = FALSE---------------------------------------------------
library(ecocbo)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.retina=2,
  fig.align='center',
  fig.width = 7, 
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----step0, eval=FALSE--------------------------------------------------------
#  # Load data and adjust it.
#  data(epiDat)
#  
#  simResults <- prep_data(data = epiDat, type = "counts", Sest.method = "average",
#                          cases = 5, N = 100, sites = 10,
#                          n = 5, m = 5, k = 30,
#                          transformation = "none", method = "bray",
#                          dummy = FALSE, useParallel = TRUE)

## ----step1--------------------------------------------------------------------
compVar <- scompvar(data = simResults)
compVar


## ----step21-------------------------------------------------------------------
cboCost <- sim_cbo(comp.var = compVar, ct = 20000, ck = 100, cj = 2500)
cboCost

## ----step22-------------------------------------------------------------------
cboPrecision <- sim_cbo(comp.var = compVar, multSE = 0.10, ck = 100, cj = 2500)
cboPrecision


## ----step3--------------------------------------------------------------------
betaResult <- sim_beta(simResults, alpha = 0.05)
betaResult


## ----step4--------------------------------------------------------------------
plot_power(data = betaResult, n = NULL, m = 3, method = "both")


