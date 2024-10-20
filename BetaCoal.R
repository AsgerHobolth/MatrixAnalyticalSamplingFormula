##-------------------------------------------------------------
## Dependencies: library("partitions")
##-------------------------------------------------------------
##
## Purpose:
## This function determines the state space and corresponding 
## rate matrix for the block counting process for a number of 
## samples n in the beta coalescent with parameter alpha.
##
## Name: BetaRateMatAndStateSpace
## Author: Hobolth
## Date: Spring 2024
##
## Input:
## n: Number of samples
## alpha: alpha parameter in the beta coalescent (1<alpha<2)
##
## Output:
## List consisting of
## RateM: Rate matrix of beta coalescent
## StSpM: Matrix with rows corresponding to the states
##        A state is a n-dimensional row vector (a1,...,an)
##        similar to Ewens sampling formula.
##        We always have a1+2*a2+...+n*an=n
##        The beginning state is (n,0,...,0),
##        the next state is (n-2,1,0,...,0) 
##        with rate (n choose 2), 
##        and the ending state (MRCA) is (0,...,0,1).
##
##----------------------------------------------------------------
BetaRateMatAndStateSpace <- function(n,alpha){
  ##----------------------------------------------------
  ## Possible states
  ##----------------------------------------------------
  ## Size of the state space (number of states)
  nSt <- P(n)
  ## Definition of the state space
  StSpM <- matrix( ncol=n,nrow=nSt )
  ## Set of partitions of [n]
  x <- parts(n)
  ## Rewriting the partitions as (a1,...,an)
  for (i in 1:nSt) {
    st <- x[,i]
    StSpM[i,] <- tabulate(x[,i],nbins=n)
  }
  ## Reordering
  StSpM <- StSpM[order(rowSums(StSpM),decreasing=TRUE),]
  ## Because of this ordering the process is feed-forward, i.e.
  ## below the diagonal the entries are always zero
  ##----------------------------------------------------
  ## Intensity matrix
  ##----------------------------------------------------
  RateM <- matrix(0,ncol=nSt,nrow=nSt)
  ## Algorithm for finding rates between states
  for (i in 1:(nSt-1)){
    for (j in (i+1):nSt){
      #cat("\n")
      #cat(i," state i",StSpM[i,])
      #cat(" ",j," state j",StSpM[j,])
      cvec <- StSpM[i,]-StSpM[j,]
      #cat(" cvec",cvec)
      cneg <- - cvec*ifelse(cvec<0,1,0)
      #cat(" cneg",cneg)
      ## Check that exactly one new branch is created
      check1 <- ( sum(cvec[cvec<0])==-1 ) ## TRUE/FALSE statement
      #cat(" check1",check1)
      if (check1){
        ## Number of branches before coalescent event
        nBranchBefore <- sum( StSpM[i,] )
        ## Number of branches that coalesce 
        nBranchCoal <- sum(cvec)+1 
        prelimRate <- 
          beta(nBranchCoal-alpha,nBranchBefore-nBranchCoal+alpha)/
          beta(alpha,2-alpha)
        prelimRate <- 
          prelimRate*prod( choose( StSpM[i,],cvec*ifelse(cvec>0,1,0) ) )
        RateM[i,j] <- prelimRate
      }
    }
  }
  ## Diagonal part of the rate matrix
  for (i in 1:nSt){
    RateM[i,i] <- -sum(RateM[i,])
  }
  return( list( RateM=RateM,StSpM=StSpM ) )
}
##--------------------------------------------------------
