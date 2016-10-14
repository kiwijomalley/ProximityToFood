##########################################################################################
# Runs analysis in R using Brugs
# Must install Brugs package first!
##########################################################################################

library(BRugs)

#Set up: Place this script and the associated data in the appropriate directories

#Conditions under which to run model
restrict4=1 #Restrict to 4-town region
tractcf=0 #Tract carried forward from first wave
prodnorm=0 #Use of a product-normal prior
if (prodnorm==1) {
 restrict4=1 
 tractcf=0 
}

modelCheck("RealBMICodeInt.r")   # Check Model written in WinBUGS
modelData("simdata2wbugs.txt")    # Load data
modelData("RealBMIConst4.txt") #Constants
modelCompile(numChains=1)	 # Compile with specified number of chains

if (prodnorm==1) {
 modelInits("RealBMIInitsPN.txt") # Initial values under product normal prior
} else {
 modelInits("RealBMIInits.txt") # Initial values under inverse-Wishart prior
}
#modelInits (c("AnyPsychInParam2.txt", "AnyPsychInParam2.txt"))	# If use two chains

modelGenInits() #Generate initial values for remaining parameters

samplesClear('*')

# Samples
samplesSetBeg(20000)		# Start Storing after burn-in of 20000
samplesSet(c("be","sigma","rho","taurand","corrand","tauarea","corarea","iccind","bothfdist","onefdist"))
samplesSetThin(1)			# Thin by 1
samplesSetFirstChain(1)
samplesSetLastChain(1)
modelUpdate(100000)

# Results for Chain 1
results.chain1<-samplesStats(c("be","sigma","rho","taurand","corrand","tauarea","corarea","iccind","bothfdist","onefdist"), 
	firstChain=samplesGetFirstChain(),
	lastChain=samplesGetFirstChain(), thin=samplesGetThin())

# Results for Chain 2
#results.chain2<-samplesStats(c("be","sigma","rho","taurand","corrand","tauarea","corarea","iccind","bothfdist","onefdist"), 
	#firstChain=samplesGetLastChain(),
	#lastChain=samplesGetLastChain(), thin=samplesGetThin())

# Results for Both Chains Combined 
#results<-samplesStats(c("be","sigma","rho","taurand","corrand","tauarea","corarea","iccind","bothfdist","onefdist"), 
	#firstChain=samplesGetFirstChain(),
	#lastChain=samplesGetLastChain(), thin=samplesGetThin())

#Make my own transformations and posterior inferences
#Need to load in data set etc.
twosamp <- 0
if (twosamp==1) {
 output1 <- data.frame(samplesHistory("*", plot=F,firstChain=samplesGetFirstChain(),
	lastChain=samplesGetFirstChain()))
 output2 <- data.frame(samplesHistory("*", plot=F,firstChain=samplesGetLastChain(),
	lastChain=samplesGetLastChain()))
 output <- rbind(output1,output2)
} else {
 output <- data.frame(samplesHistory("*", plot=F))
}
be <- as.matrix(output[,1:22])
sigma <- as.matrix(output[,29])
rho <- as.matrix(output[,28])
stdint <- sqrt(as.matrix(output[,34]))
stdslope <- sqrt(as.matrix(output[,37]))
corsubj <- as.matrix(output[,25])
stdhome <- sqrt(as.matrix(output[,30]))
stdwork <- sqrt(as.matrix(output[,33]))
corarea <- as.matrix(output[,24])
iccind <- as.matrix(output[,26]) #Individual to observation ICC
bothfooddist <- as.matrix(output[,23]) #Both residential and employment exposure
onefooddist <- as.matrix(output[,27]) #Either residential or employment exposure

par(mfrow=c(3,2), srt=0, mai=c(0.6, 0.6, 0.4, 0.2), mgp=c(2,1,0))
plot(density(be[,17]),xlab="",ylab="Density",
 main="Home distance to FF Effect for Females")
plot(density(be[,20]),xlab="",ylab="Density",
 main="Work distance to FF Effect For Males")
plot(density(stdhome),xlab="",ylab="Density",
 main="Residential Standard Deviation")
plot(density(stdwork),xlab="",ylab="Density",
 main="Workplace Standard Deviation")
plot(density(corsubj),xlab="",ylab="Density",
 main="Intercept and Slope Effect Correlation")
plot(density(corarea),xlab="",ylab="Density",
 main="Home and Work Effect Correlation")

#Output results to a file
if (restrict4==1) {
 dev.copy2eps(file=paste(outdir,"RealAnalPostPlots4.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file
 dev.copy2pdf(file=paste(outdir,"RealAnalPostPlots4.pdf",sep=""), width=6, height=6, out.type="pdf") #to file
 dput(output,file=paste(outdir,"RealBMIResults4.txt",sep=""))
} else {
 dev.copy2eps(file=paste(outdir,"RealAnalPostPlots14.eps",sep=""), width=6, height=6, horizontal=FALSE) #to file
 dev.copy2pdf(file=paste(outdir,"RealAnalPostPlots14.pdf",sep=""), width=6, height=6, out.type="pdf") #to file
 dput(output,file=paste(outdir,"RealBMIResults14.txt",sep=""))
}

#Summary statistics
mn=apply(output,2,mean)
se=apply(output,2,function(x) sqrt(var(x)))
cr=apply(output,2,function(x) quantile(x,c(0.025,0.975)))
estout=data.frame(rbind(mn,se,cr))
names(estout)[1:22]=c("Intercept","Wave 2","Wave 3","Wave 4","Wave 5","Wave 6","Wave 7","Wave 8","YOB",
 "Smokes","Married","Education 1","Education 2","Employ","TractPoverty","EmpTractPoverty",
 "Dist1HomeFem","Dist1HomeMale","Dist1WorkFem","Dist1WorkMale","DriveDistFem","DriveDistMale")
names(estout)[23:37]=c("BothFoodDist","CorHomeWork","CorIndiv","ICCIndiv","OneFoodDist","rho","sigma","StdHome",
 "StdHomeWork","StdHomeWork","StdWork","StdInt","StdIntSlope","StdIntSlope","StdSlope")
