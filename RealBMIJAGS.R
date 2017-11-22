## R wrapper code to estimate distance-to-food model ##

#Set-up: Load in libraries
library(rjags)
coda.options(combine.plots=TRUE,combine.stats=TRUE)

#Specify filename of JAGS code for Bayesian analysis and data in and output directories
rsource="Enter directory where code is stored" #If run on a Linux server
setwd(rsource)
datdir="../Data/" #Sub-directory to store data
outdir="../Output/" #Sub-directory to store output

#Specify code and initial values to use
JAGScode<-"RealBMICodeIntIW.bug"
IC<-list(be=c(30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), isigma=0.25, u=0.5) #Initial values

#Read in data
BMIdata<-read.table(paste(datdir,"simdata.txt",sep=""),sep=" ",col.names=c("ID", "wave", "Omni", "Tobs", "TractHome", "TractWork", "BMI", "yob", "smokes", "male", "married", "educ", "tractpov", "etractpov", "unemploy", "DistHome", "DistWork", "DriveDist"), skip=1) #Read in data for analysis

#Define parameters and operating characteristics of MCMC procedure
parameters = c("be", "sigma", "rho", "taurand", "tauarea", "corrand", "corarea")    # The parameter(s) to be monitored.
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
dev.copy2eps(file=paste(outdir,"TraceSeqIW4sim.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Trace plots of long-run averages
par(mfrow=c(1,1), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
trpSample=cbind(corareaSample[1:nIter],corareaSample[(nIter+1):(2*nIter)],corareaSample[(2*nIter+1):(3*nIter)])
for (i in 2:nIter) {
 trpSample[i,1]=((i-1)*trpSample[i-1]+corareaSample[i])/i	
 trpSample[i,2]=((i-1)*trpSample[nIter+i-1]+corareaSample[nIter+i])/i
 trpSample[i,3]=((i-1)*trpSample[2*nIter+i-1]+corareaSample[2*nIter+i])/i
}
mn=mean(corareaSample); sd=sqrt(var(corareaSample))
plot(trpSample[,1],type='l',ylim=c(mn-sd,mn+sd),xlim=c(0,nIter),
 main='Trace plot of posterior mean for 3 chains',
 xlab='Draw',
 ylab='Mean')
lines(trpSample[,2],type='l',col='red')
lines(trpSample[,3],type='l',col='green')
dev.copy2eps(file=paste(outdir,"TracePlotIW4sim.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Autocorrelation plots
par(mfrow=c(1,1), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
autocorr.plot(corareaSample[1:nIter],main="MCMC Chain for p",auto.layout=FALSE)
dev.copy2eps(file=paste(outdir,"AutoCorrPlotIW4sim.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file


## Generate statistical inferences and outputs ##

#Plots of drawn values
par(mfrow=c(3,2), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
plot(density(mcmcChain[,"be[1]"]),xlab="",ylab="Density",
 main="Intercept")
plot(density(mcmcChain[,"be[2]"]),xlab="",ylab="Density",
 main="Wave 2 effect")
plot(density(mcmcChain[,"be[3]"]),xlab="",ylab="Density",
 main="Wave 3 effect")
plot(density(mcmcChain[,"be[4]"]),xlab="",ylab="Density",
 main="Wave 4 effect")
plot(density(mcmcChain[,"corrand"]),xlab="",ylab="Density",
 main="Correlation for individual random effects")
plot(density(mcmcChain[,"corarea"]),xlab="",ylab="Density",
 main="Correlation for area random effects")
dev.copy2eps(file=paste(outdir,"FixedIW4sim.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file

#Summary statistics of draws (the types of things you might put in a paper)
mn=apply(mcmcChain,2,mean)
stdev=sqrt(apply(mcmcChain,2,var))
lowl=apply(mcmcChain,2,quantile,0.025)
uppl=apply(mcmcChain,2,quantile,0.975)
sumdata=cbind(mn,stdev,lowl,uppl)
write.table(t(sumdata),file=paste(outdir,"BayesSummaryStatisticsIW4sim.txt",sep=""))

 
