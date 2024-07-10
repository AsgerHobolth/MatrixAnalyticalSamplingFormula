##-------------------------------------------------------------
##
## Purpose: Probability of a sample
##
## Name: SmplPrb
##
## Input: 
## smpl: Sample (mutation count vector)
## PrbMdl: Probability model (object) 
##
## Output:
## Probability of the sample
##
## Dependencies: Uses the R package 'partitions' 
## library('partitions')
##
##------------------------------------------------------------
SmplPrb <- function( smpl,PrbMdl ){
  initVec <- PrbMdl$initVec
  PrbArray <- PrbMdl$PrbArray
  exitVec <- PrbMdl$exitVec
  ## Zero mutations
  if (sum(smpl)==0) return( as.double( initVec %*% exitVec) )
  ## Non-zero number of mutations
  nSt <- length(initVec)
  mSet <- multiset( c( rep( 1:length(smpl), smpl ) ) )  
  sum.prb <- 0
  for (k in 1:ncol(mSet)){
    ProdMat <- diag(nSt)
    for (m in 1:nrow(mSet)){
      ProdMat <-  ProdMat %*% PrbArray[mSet[m,k],,]
    }
    sum.prb <- sum.prb + initVec %*% ProdMat %*% exitVec
  }
  return( as.double( sum.prb ) )
}