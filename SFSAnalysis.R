##---------------------------------------------------
## Analysis of the sampling probabilities of the SFS
## for a fixed number of samples and mutations 
##---------------------------------------------------
library("partitions")
source("EvoModelToPrbModel.R")
source("PrbOfSample.R")
source("SFSModel.R")
##-------------------------------
nSmpls <- 4   ## Sample size
nMuta <- 5    ## Number of mutations
## SFS compositions
SFScomps <- as.matrix( t(compositions( nMuta,(nSmpls-1) )) )
nComps <- nrow( SFScomps )
cat("Number of compositions:",nComps,"\n")
## Average pairwise distance
AvrgPairDist <- rep( 0,nComps )
Wgt <- (1:(nSmpls-1)*((nSmpls-1):1))/choose(nSmpls,2)
for (i in 1:nComps){
  AvrgPairDist[i] <- sum( SFScomps[i,]*Wgt )
}
## Mean total tree length
meanTreeLen <- 2*sum(1/(1:(nSmpls-1)))
## Mutation rate (Watterson's theta)
mutaRate <- nMuta/meanTreeLen
## Evolutionary model
SFSMdl <- SFSModel( nSmpls,mutaRate=mutaRate )
## Convert evolutionary model to probability model 
PrbMdl <- EvoMdlToPrbMdl( SFSMdl )
## Probability for each composition
prbSmpl <- rep(0,nComps)
for (i in 1:nComps){
  prbSmpl[i] <- SmplPrb( SFScomps[i,],PrbMdl=PrbMdl )
  cat(i,"")
}
## Normalize sample probabilities
prbSmpl <- prbSmpl/sum(prbSmpl)
##-------------------------------------------------------------
## Plot average pairwise distance versus probability of sample
##-------------------------------------------------------------
lnPrb <- log(prbSmpl)
plot( AvrgPairDist,lnPrb,
      pch=as.character(SFScomps[,1]),#pch=19,
      cex=0.8,
      xlab="Average pairwise distance",ylab="log-probability",
      main=paste(nMuta,"mutations,",
                 nSmpls,"samples,",
                 nComps,"SFS compositions"))
##-----------------------------------------------------
## Plot the sorted cumulative probability 
##------------------------------------------------------
srtPRB <- sort( prbSmpl,decreasing=TRUE,index.return=TRUE )
plot( cumsum(prbSmpl[srtPRB$ix]) ,
      ylim=c(0,1),
      type="l",xlab="Compositions sorted after probability (high to low)",
      main="Cumulated probability for the compositions",
      ylab="Cumulated probability",lwd=3,col="red")
abline(h=c(0,0.2,0.5,0.8,1),col="gray")
abline( v=c( 10*(1:floor(nComps/10)) ) )
##-----------------------------------------------------------
## Print the sorted compositions 
##-----------------------------------------------------------
hghPrb <- srtPRB$ix
HghPrbComp <- 
  cbind( SFScomps[hghPrb,],
         round(100*prbSmpl[hghPrb],digits=2),
         round(AvrgPairDist[hghPrb],digits=2) )
colnames(HghPrbComp) <- 
  c("singl","doubl","tripl","prb*100","av.pair.dist")
cat("Sorted compositions:","\n")
print(HghPrbComp)



