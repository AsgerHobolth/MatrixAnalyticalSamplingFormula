##---------------------------------
## Dependencies
##---------------------------------
## library(partitions)
source("BetaCoal.R")
##
##----------------------------------
##
## Author: Asger Hobolth
## Date: May 2024
## 
## Purpose
## Define the Beta-coalescent as an evolutionary model
## for the site frequency spectrum
##
## Input
## nSmpl: Sample size n
## alphaPrm: Parameter alpha
## mutaRate: Mutation rate 
##
## Output
## list giving the Beta-coalescent evolutionary model: 
## initVec: Initial vector
## SMat: Subintensity transition rate matrix
## RMat: Reward matrix 
## mutaRate: Mutation rate
##
## 
##--------------------------------------------------
BetaModel <- function( nSmpl,alphaPrm,mutaRate ){
  BetaRMASS <- BetaRateMatAndStateSpace( nSmpl,alphaPrm )
  m <- dim(BetaRMASS$RateM)[1]
  ## Obtain sub-intensity matrix
  SubIntMat <- BetaRMASS$RateM[-m,-m]
  ## Reward matrix is the state space matrix of the block counting process
  RewMat <- BetaRMASS$StSpM[-m,1:(nSmpl-1)]
  ## Initial probability vector
  InitPrb <- c(1, rep(0, m-2))
  ## Output the list
  out <- list()
  out$initVec <- InitPrb
  out$SMat <- SubIntMat
  out$RMat <- RewMat
  out$mutaRate <- mutaRate
  return(out)
}

