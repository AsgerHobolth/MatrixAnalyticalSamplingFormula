##------------------------------------------------
## Dependencies
library('PhaseTypeR')
library('partitions')
source("BetaModel.R")
source("EvoModelToPrbModel.R")
source("PrbOfSample.R")
##------------------------------------------------
## Define sample size and expected number of mutations
n <- 4
nMuta <- 10
## Specify alphas in the Beta(alpha,2-alpha)-coalescent
## alpha=1 is the Bolthausen-Sznitman coalescent 
## alpha=2 is the Kingman coalescent
nModels <- 21  # Number of models (should be odd number)
alphaVec <- seq( 1.01,1.99,len=nModels )
##-------------------------------------------------------
## Determine the mutation rate in the Beta-models such that  
## they all have the same number of expected mutations
##------------------------------------------------------
mutaRateVec <- rep( 0,nModels )
for (i in 1:nModels){
  BetaEvoModel <- BetaModel( nSmpl=n,alphaPrm=alphaVec[i],mutaRate=1 )
  ## Tree HEIGHT in the Beta-coalescent as PH object
  BetaHeightPH <- PH( BetaEvoModel$SMat,BetaEvoModel$initVec )
  ## Total tree LENGTH in the Beta-coalescent as PH object
  BetaTotalPH <- 
    reward_phase_type( BetaHeightPH,rowSums(BetaEvoModel$RMat) )
  ## Mean of total tree length
  BetaMutaRate <- nMuta/mean( BetaTotalPH )
  mutaRateVec[i] <- BetaMutaRate
}
## Check that for the Kingman model, the mutation rate is correct
cat("Kingman expected number of mutations:",
    mutaRateVec[nModels]*2*sum(1/1:(n-1)),"\n")
##--------------------------------------------------------------
## We consider all SFS compositions for fixed sample size n 
## and number of mutations nMuta
SFScomps <- as.matrix( t(compositions( nMuta,(n-1) )) )
nComps <- nrow( SFScomps )
cat("Number of compositions:",nComps,"\n")
##---------------------------------------------------------------
## Loop across the models to determine the probabilities for the 
## compositions
PrbMat <- matrix( 0,nrow=nComps,ncol=nModels )  # Store results
for (i in 1:nModels){
  ## Evolutionary Beta-coalescent
  BetaEvoModel <- 
    BetaModel( nSmpl=n,alphaPrm=alphaVec[i],mutaRate=mutaRateVec[i] )
  ## Convert the evolutionary model to probability model
  BetaPrbMdl <- EvoMdlToPrbMdl( BetaEvoModel )
  for (j in 1:nComps){
    PrbMat[j,i] <- SmplPrb( SFScomps[j,],PrbMdl=BetaPrbMdl )
  }
  cat(i,"")
}
## Normalize the probabilities
PrbMat <- PrbMat %*% diag( 1/colSums(PrbMat) )
##--------------------------------------------------------------
## Plot the result in two ways:
## As a function of the parameter alpha (first plot below)
## Cumulative probability (second plot below)
##--------------------------------------------------------------
## First plot
cf <- 0.04 # cut off value
plotIndx <- union( which(PrbMat[,1]>cf),
                   which(PrbMat[,nModels]>cf) ) 
length(plotIndx)
colrs <- rainbow(length(plotIndx))
ParamFun <- function(a){
  return( paste(c( "(",a,")" ) , collapse="" ) )
}
#pdf("BetaCoalescent.pdf",height=6,width=8)
plot(alphaVec,PrbMat[1,],col=colrs[1],
     xlab="Parameter alpha",ylab="Sampling Probability",
     main="Beta-coalescent",
     type="l",ylim=c(0,max(PrbMat)),lwd=2)
for (i in 2:length(plotIndx)){
  points(alphaVec,PrbMat[plotIndx[i],],col=colrs[i],
         type="l",lwd=2)
}
txt <- apply(SFScomps,1,paste,collapse=",")
txt <- sapply(txt,ParamFun,USE.NAMES=FALSE)
legend("topright",txt[plotIndx],lty=1,lwd=2,col=colrs,
       cex=0.6,bty="n")
abline(h=cf,lty="dotted",col="gray")
#dev.off()
## Second plot
#pdf(file="BetaCumPlot.pdf",height=6,width=8)
plot(cumsum( sort(PrbMat[,1],decreasing=TRUE) ),
     ylim=c(0,1),type="l",
     xlab="Compositions sorted after probability (high to low)",
     main="Cumulative probability of compositions",
     ylab="Cumulative probability",lwd=3,col="red",
     cex.lab=1.2,cex.main=1.2,cex.axis=1.1)
points(1:(dim(PrbMat)[1]),
       cumsum( sort(PrbMat[,round(nModels/2)],decreasing=TRUE) ),
       type="l",lwd=3)
points(1:(dim(PrbMat)[1]),
       cumsum( sort(PrbMat[,nModels],decreasing=TRUE) ),
       type="l",lwd=3,col="blue")
abline( h=c(0,0.2,0.5,0.8,1),col="gray",lty="dotted")
abline( v=c( 10*(1:floor(nComps/10)) ),col="gray",lty="dotted" )
legend("bottomright",col=c("red","black","blue"),
       bty="n",lty=1,lwd=3,cex=1.2,
       c("Bolthausen-Sznitman (alpha=1)",
         "alpha=1.5","Kingman (alpha=2)"))
#dev.off()
