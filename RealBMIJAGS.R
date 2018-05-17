## R wrapper code to estimate distance-to-food model ##

## TO DO: If there is further push-back need to determine if residual correlation can be more accurately estimated.
##        Is the problem due to the way the model is specified in JAGS?
##        If can't resolve, might be best to drop error variance and serial correlation parameters from the simulation figure and results!

#Set-up: Load in libraries
library(rjags)
coda.options(combine.plots=TRUE,combine.stats=TRUE)

#Specify filename of JAGS code for Bayesian analysis and data in and output directories
rsource="/Volumes/STORAGE/JOMalley/Dartmouth/Biostatistics/DistanceToFood/Code" #If run on a Linux server
setwd(rsource)
datdir="../Data/" #Directory to store data
outdir="../Output/" #Directory to store output

#Specify code and initial values to use
IW=1
subtype=1
if (IW==1) {
 if (subtype==1) {
  JAGScode<-"RealBMICodeIntIW.bug"
  #JAGScode<-"RealBMICodeIntIWlog.bug"
 } else if (subtype==2) {
  JAGScode<-"RealBMICodeIntArea0.bug"
 } else if (subtype==3) {
  JAGScode<-"RealBMICodeIntRSlope0.bug"
 } else {
  JAGScode<-"RealBMICodeIntSerial0.bug"
 } 
 IC<-list(be=c(30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), isigma=0.25, u=0.5) #Initial values
} else {
 JAGScode<-"RealBMICodeIntProdNorm.bug"
 IC<-list(be=c(30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), isigma=0.25, ipntauarea=c(1, 0.1), ipntaurand=c(1, 0.1), psirand=0, psiarea=0, u=0.5) #Initial values	
}

#Read in data
datatype=1
if (datatype==1) {
 #BMIdata<-read.table(paste(datdir,"simdata.txt",sep=""),sep=" ",col.names=c("ID", "wave", "Omni", "Tobs", "TractHome", "TractWork", "BMI", "yob", "smokes", "male", "married", "educ", "tractpov", "etractpov", "unemploy", "DistHome", "DistWork", "DriveDist"), skip=1) #Read in data for analysis
 BMIdata<-read.table(paste(datdir,"datawbugs4town.txt",sep=""),sep=" ",col.names=c("ID", "wave", "Omni", "Tobs", "TractHome", "TractWork", "BMI", "yob", "smokes", "male", "married", "educ", "tractpov", "etractpov", "unemploy", "DistHome", "DistWork", "DriveDist"), skip=1) #Read in data for analysis
} else if (datatype==2) {
 BMIdata<-read.table(paste(datdir,"datawbugs14town.txt",sep=""),sep=" ",col.names=c("ID", "wave", "Omni", "Tobs", "TractHome", "TractWork", "BMI", "yob", "smokes", "male", "married", "educ", "tractpov", "etractpov", "unemploy", "DistHome", "DistWork", "DriveDist"), skip=1) #Read in data for analysis
} else {
 BMIdata<-read.table(paste(datdir,"datawbugs4towncf.txt",sep=""),sep=" ",col.names=c("ID", "wave", "Omni", "Tobs", "TractHome", "TractWork", "BMI", "yob", "smokes", "male", "married", "educ", "tractpov", "etractpov", "unemploy", "DistHome", "DistWork", "DriveDist"), skip=1) #Read in data for analysis
}

#Define parameters and operating characteristics of MCMC procedure
parameters = c("be", "sigma", "rho", "taurand", "tauarea", "corrand", "corarea", "Rsquare")    # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 20000           # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=150000          # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.

#Make an object for fitting model using JAGS
jagsModel = model<-jags.model(JAGScode,data=BMIdata,inits=IC,n.chains=nChains)

# Burn-in phase of model estimation:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )

# Main phase of MCMC model: draws from the posterior distribution are saved in codaSamples
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters ,
                            n.iter=nIter , thin=thinSteps )
# resulting codaSamples object has these indices:
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

mcmcChain = as.matrix( codaSamples )
summary(mcmcChain)


## Prepare for and implement post model-fitting computations to check whether chain converged ##

beSample0 = mcmcChain[,"be[1]"] # Put sampled values in a vector.
beSample1 = mcmcChain[,"be[2]"] # Put sampled values in a vector.
sigmaSample = mcmcChain[,"sigma"] # Put sampled values in a vector.
rhoSample = mcmcChain[,"rho"] # Put sampled values in a vector.
tau1Sample = mcmcChain[,"tauarea[1,1]"] # Put sampled values in a vector.
corrandSample = mcmcChain[,"corrand"] # Put sampled values in a vector.
corareaSample = mcmcChain[,"corarea"] # Put sampled values in a vector.

#Plots of sequences of draws
par(mfrow=c(3,2), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
plot(beSample0[1:nIter],type="l",main="MCMC Chain for Intercept",xlab="iteration",ylab="estimate")
plot(beSample1[1:nIter],type="l",main="MCMC Chain for Wave 2",xlab="iteration",ylab="estimate")
plot(sigmaSample[1:nIter],type="l",main="MCMC Chain for error variance",xlab="iteration",ylab="estimate")
plot(tau1Sample[1:nIter],type="l",main="MCMC Chain for home variance",xlab="iteration",ylab="estimate")
plot(corrandSample[1:nIter],type="l",main="MCMC Chain for individual RE correlation",xlab="iteration",ylab="estimate")
plot(corareaSample[1:nIter],type="l",main="MCMC Chain for area RE correlation",xlab="iteration",ylab="estimate")
dev.copy2eps(file=paste(outdir,"TraceSeqIW4.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Trace plots of long-run averages
par(mfrow=c(1,1), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
trpSample=cbind(corareaSample[1:nIter],corareaSample[(nIter+1):(2*nIter)],corareaSample[(2*nIter+1):(3*nIter)])
for (i in 2:nIter) {
 trpSample[i,1]=((i-1)*trpSample[i-1]+corareaSample[i])/i	
 trpSample[i,2]=((i-1)*trpSample[nIter+i-1]+corareaSample[nIter+i])/i
 trpSample[i,3]=((i-1)*trpSample[2*nIter+i-1]+corareaSample[2*nIter+i])/i
}
mn=mean(corareaSample); sd=sqrt(var(corareaSample))
plot(trpSample[,1],type='l',ylim=c(mn-0.25*sd,mn+0.25*sd),xlim=c(0,nIter),
 main='Trace plot of posterior mean for 3 chains',
 xlab='Draw',
 ylab='Mean')
lines(trpSample[,2],type='l',col='red')
lines(trpSample[,3],type='l',col='green')
dev.copy2eps(file=paste(outdir,"TracePlotIW4.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Autocorrelation plots
par(mfrow=c(1,1), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
autocorr.plot(corareaSample[1:nIter],main="MCMC Chain for p",auto.layout=FALSE)
dev.copy2eps(file=paste(outdir,"AutoCorrPlotIW4.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Gelman-Rubin MCMC convergence statistic
MCMCdiag=gelman.diag(codaSamples,confidence=0.95,autoburnin=FALSE,multivariate=FALSE)$psrf

#Residual analysis
rsquare=mcmcChain[,"Rsquare"]
print(summary(rsquare))

## Generate statistical inferences ##

#Plots of drawn values
par(mfrow=c(3,2), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
plot(density(mcmcChain[,"be[18]"]),xlab="",ylab="Density",
 main="Home distance to FF effect for females")
plot(density(mcmcChain[,"be[21]"]),xlab="",ylab="Density",
 main="Work distance to FF effect for males")
plot(density(mcmcChain[,"tauarea[1,1]"]),xlab="",ylab="Density",
 main="Residential Std. Dev.")
plot(density(mcmcChain[,"tauarea[2,2]"]),xlab="",ylab="Density",
 main="Workplace Std. Dev.")
plot(density(mcmcChain[,"corrand"]),xlab="",ylab="Density",
 main="Correlation of individual random effects")
plot(density(mcmcChain[,"corarea"]),xlab="",ylab="Density",
 main="Correlation of area random effects")
dev.copy2eps(file=paste(outdir,"FixedSerialIW4.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Summary statistics of draws (the types of things you might put in a paper)
mn=apply(mcmcChain,2,mean)
stdev=sqrt(apply(mcmcChain,2,var))
lowl=apply(mcmcChain,2,quantile,0.025)
uppl=apply(mcmcChain,2,quantile,0.975)
pg0=apply(mcmcChain,2,function(x) mean(x>0))
sumdata=cbind(mn,stdev,lowl,uppl,pg0)
write.table(t(sumdata),file=paste(outdir,"BayesSummaryStatisticsIW4.txt",sep=""))

#Generate estimates for joint hypothesis test of female home and male workplace proximity
prHypW=mean(mcmcChain[,"be[18]"]<0)
prHypM=mean(mcmcChain[,"be[21]"]<0)
prJointHyp=mean(mcmcChain[,"be[18]"]<0 & mcmcChain[,"be[21]"]<0)
print(c(prHypW,prHypM,prJointHyp))

#Posterior correlation between error variance and serial correlation
cor(mcmcChain[,"sigma"],mcmcChain[,"rho"])



 