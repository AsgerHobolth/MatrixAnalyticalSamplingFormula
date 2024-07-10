##
## Author: Hobolth
## Date: April 2024
## 
## Purpose: 
## Converts an evolutionary model to the desired
## initial probability vector, 
## exit (sub)probability vector, and
## (sub)transition probability matrices.
##
## Name: EvoMdlToPrbMdl
##
## Input:
## List giving the four parameters of the time-homogeneous 
## evolutionary model
## $initVec
## $SMat
## $RMat
## $mutaRate
##
## Output:
## List consisting of
## $initVec: Initial probability vector
## $PrbArray: Array of (sub)transition probability matrices
## $exitVec: Exit (sub)probability vector
##
EvoMdlToPrbMdl <- function(EvoMdl){
  SMat <- EvoMdl$SMat
  RMat <- EvoMdl$RMat
  lam <- EvoMdl$mutaRate
  ## Ptotal matrix
  IdMat <- diag( nrow=nrow(SMat) )
  rtotalVec <- rowSums( RMat ) 
  PtotalMat <- solve( IdMat - ( diag(1/rtotalVec)%*% SMat/lam) )
  ## exit vector
  nSt <- nrow( SMat )
  exitVec <- ( IdMat-PtotalMat ) %*% rep(1,nSt)                 
  ## Sub-transition matrices
  PrbArray <- array( 0, dim = c( ncol(RMat), rep(nSt, 2) ) )
  for (i in 1:ncol(RMat)){
    PrbArray[i,,] <- PtotalMat %*% ( diag( RMat[,i]/rtotalVec )) 
  }
  out <- list()
  out$initVec <- EvoMdl$initVec
  out$PrbArray <- PrbArray
  out$exitVec <- exitVec
  return(out)
}