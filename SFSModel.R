##--------------------------------------------------------------
## Specification of evolutionary model for the SFS  
## in the case of the Kingman coalescent
##--------------------------------------------------------------
## 
## Name:
## SFSModel( nSmpl,mutaRate )
## 
## Input: 
## nSmpl: Number of samples
## mutaRate: Mutation rate
##
## Output:
## Object of type evolutionary model, i.e.
## list of the four components
## 
## Dependencies: 
## Function SFSRateMAndStateSpace()
source("SFSCoal.R")
##------------------------------------------------ 
SFSModel <- function( nSmpl,mutaRate ){
  BlockCP <- SFSRateMAndStateSpace( nSmpl )
  nSt <- nrow( BlockCP$RateM )
  out <- list()
  out$initVec <- c( 1,rep(0,(nSt-2)) )
  out$SMat <- BlockCP$RateM[-nSt,-nSt]
  nTypes <- ncol( BlockCP$StSpM )
  out$RMat <- BlockCP$StSpM[-nSt,-nTypes]
  out$mutaRate <- mutaRate
  return(out)
}
