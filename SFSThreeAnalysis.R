##-------------------------------------------------------
## Simple example: SFS for sample size n=3
##-------------------------------------------------------
library("partitions")
source("EvoModelToPrbModel.R")
source("PrbOfSample.R")
##-------------------------------------------------------
## Analyze SFS for coalescent model with sample size n=3
##-------------------------------------------------------
## Specify evolutionary model: Coalescent with n=3
SFSThree <- function(mutaRate){
  initVec <- c(1, 0)
  SMat <- matrix( c(-3,3, 0,-1), nrow=2, ncol=2, byrow=TRUE)
  RMat <- matrix( c(3,0, 1,1), nrow=2, ncol=2, byrow=TRUE)
  out <- list()
  out$initVec <- initVec
  out$SMat <- SMat
  out$RMat <- RMat
  out$mutaRate <- mutaRate
  return(out)
}
mutaRate <- 5
SFSMdl <- SFSThree( mutaRate=mutaRate )
## Convert the evolutionary model to the desired 
## discrete probability vectors and transition matrices 
PrbMdl <- EvoMdlToPrbMdl( SFSMdl )
##---------------------------------------------
## Simple examples of probabilities of samples:
## Compare analytical expressions and phase-type implementation
##---------------------------------------------
## SFS sample (0,0)
SmplPrb( c(0,0), PrbMdl=PrbMdl )
( 3/(3+3*mutaRate) )*( 1/(1+2*mutaRate) )
## SFS sample (1,0)
SmplPrb(c(1,0),PrbMdl=PrbMdl)
( 3/(3+3*mutaRate) )*( mutaRate/(2*mutaRate+1))*( 1/(1+2*mutaRate) )+
  ( 3*mutaRate/(3+3*mutaRate) )*( 3/(3*mutaRate+3))*( 1/(1+2*mutaRate) )
## SFS sample (0,1)
SmplPrb(c(0,1),PrbMdl=PrbMdl)
( 3/(3+3*mutaRate) )*( mutaRate/(1+2*mutaRate) )*( 1/(1+2*mutaRate) )
##---------------------------------------------------------------------
## Enumerate all compositions with less than a specified number of
## mutations and find their probabilities.
## Summarize in a table.
##---------------------------------------------------------------------
mutaRate <- 0.3
SFSMdl <- SFSThree( mutaRate=mutaRate )
PrbMdl <- EvoMdlToPrbMdl( SFSMdl )
## Maximum number of mutations
maxMuta <- 5
## All possible compositions with less than max number of mutations
SmplMat <- matrix( c(0,0), nrow=1,ncol=2 )
for (nMuta in 1:maxMuta){
  SmplMat <- rbind(SmplMat,t(as.matrix(compositions(nMuta,2))))
}
nSmpl <- nrow(SmplMat)
cat("Number of configurations:",nSmpl,"\n")
## Probability of each composition
prbSmpl <- rep(0,nSmpl)
for (i in 1:nSmpl){
  prbSmpl[i] <- SmplPrb(SmplMat[i,],PrbMdl=PrbMdl)
}
cat("Total probability:",sum(prbSmpl),"\n")
## Summarize the results in a matrix with four columns: 
## (number of singletons, number of doubletons,
##  total number of mutations, probability of SFS)
SmplMat <- cbind( SmplMat,SmplMat[,1]+SmplMat[,2] )
colnames(SmplMat) <- c("singletons","doubletons","total","probability")
ResMat <- cbind( SmplMat,prbSmpl )
print(ResMat)
##------------------------------------------------------------
## Plot the likelihood for 10 singletons and two doubletons
##-----------------------------------------------------------
lenMuta <- 40
mutaRateV <- c(seq(0.1,4.5,len=lenMuta/2),seq(5,15,len=lenMuta/2))
ResPrb <- rep(0,len=lenMuta)
for (i in 1:lenMuta){
  SFSMdl <- SFSThree(mutaRate=mutaRateV[i])
  PrbMdl <- ConvertEvoMdlToPrbMdl( SFSMdl )
  ResPrb[i] <- SmplPrb( c(10,2),PrbMdl=PrbMdl )
}
#pdf(file="ProbabilityCurveThree.pdf",width=8,height=6)
plot(mutaRateV,ResPrb,ylim=c(0,0.006),type="n",
     main="Likelihood curve for 10 singletons and 2 doubletons",
     xlab="Mutation rate",ylab="Probability",cex.lab=1.2,
     col="red",lty=2,lwd=2,cex.lab=1.2)
points(mutaRateV,ResPrb,type="l",col="blue",lwd=3)
#dev.off()

