model
{
   BMI[1] ~ dnorm(la[1],isigma.e[1]);
   la[1] <- pred[1];
   pred[1] <- be[1] + waveterms[1] + be[9]*cyob[1] + be[10]*smokes[1] + be[11]*married[1] + be[12]*educ1[1] 
       + be[13]*educ2[1] + be[14]*employ[1] + be[15]*ctractpov[1] + be[16]*cetractpov[1]*employ[1]
       + (be[17]*female[1] + be[18]*male[1])*cDistHome[1] + (be[19]*female[1] + be[20]*male[1])*cDistWork[1] 
       + (be[21]*female[1] + be[22]*male[1])*cDriveDist[1]
       + mu[ID[1], 1] + mu[ID[1], 2]*wave[1] + th[TractHome[1], 1] + th[TractWork[1], 2]*employ[1];
  isigma.e[1] <- isigma;

   waveterms[1] <- be[2]*wave2[1] + be[3]*wave3[1] + be[4]*wave4[1] + be[5]*wave5[1] 
                   + be[6]*wave6[1] + be[7]*wave7[1] + be[8]*wave8[1];
   wave2[1] <- equals(wave[1], 2);
   wave3[1] <- equals(wave[1], 3);
   wave4[1] <- equals(wave[1], 4);
   wave5[1] <- equals(wave[1], 5);
   wave6[1] <- equals(wave[1], 6);
   wave7[1] <- equals(wave[1], 7);
   wave8[1] <- equals(wave[1], 8);
   educ1[1] <- equals(educ[1], 1);
   educ2[1] <- equals(educ[1], 2);
   female[1] <- 1 - male[1];
   employ[1] <- 1 - unemploy[1];
   cyob[1] <- yob[1] - mean(yob[]);
   ctractpov[1] <- tractpov[1] - mean(tractpov[]);
   cetractpov[1] <- etractpov[1] - mean(etractpov[]); #Restrict to employed individuals in specification of model below
   cDistHome[1] <- DistHome[1] - mean(DistHome[]);
   cDistWork[1] <- DistWork[1]*employ[1] - mean(DistWork[])*mean(employ[]); #Only varies for employed individuals
   cDriveDist[1] <- DriveDist[1]*employ[1] - mean(DriveDist[])*mean(employ[]); #Only varies for employed individuals
    
   for (i in 2:Nobs) {
      BMI[i] ~ dnorm(la[i],isigma.e[i]);
      la[i] <- pred[i] + rho*(BMI[i-1] - pred[i-1])*nfst[i];
      pred[i] <- be[1] + waveterms[i] + be[9]*cyob[i] + be[10]*smokes[i] + be[11]*married[i] + be[12]*educ1[i] 
          + be[13]*educ2[i] + be[14]*employ[i] + be[15]*ctractpov[i] + be[16]*cetractpov[i]*employ[i]
          + (be[17]*female[i] + be[18]*male[i])*cDistHome[i] + (be[19]*female[i] + be[20]*male[i])*cDistWork[i] 
          + (be[21]*female[i] + be[22]*male[i])*cDriveDist[i]
          + mu[ID[i], 1] + mu[ID[i], 2]*wave[i] + th[TractHome[i], 1] + th[TractWork[i], 2]*employ[i];
      isigma.e[i] <- (1+ pow(rho, 2)*nfst[i])*isigma;
      nfst[i] <- step(Tobs[i]-1);

      waveterms[i] <- be[2]*wave2[i] + be[3]*wave3[i] + be[4]*wave4[i] + be[5]*wave5[i] 
                      + be[6]*wave6[i] + be[7]*wave7[i] + be[8]*wave8[i];
      wave2[i] <- equals(wave[i], 2);
      wave3[i] <- equals(wave[i], 3);
      wave4[i] <- equals(wave[i], 4);
      wave5[i] <- equals(wave[i], 5);
      wave6[i] <- equals(wave[i], 6);
      wave7[i] <- equals(wave[i], 7);
      wave8[i] <- equals(wave[i], 8);
      educ1[i] <- equals(educ[i], 1);
      educ2[i] <- equals(educ[i], 2);
      female[i] <- 1 - male[i];
      employ[i] <- 1 - unemploy[i];
      cyob[i] <- yob[i] - mean(yob[]);
      ctractpov[i] <- tractpov[i] - mean(tractpov[]);
      cetractpov[i] <- etractpov[i] - mean(etractpov[]);
      cDistHome[i] <- DistHome[i] - mean(DistHome[]);
      cDistWork[i] <- DistWork[i]*employ[i] - mean(DistWork[])*mean(employ[]);
      cDriveDist[i] <- DriveDist[i]*employ[i] - mean(DriveDist[])*mean(employ[]);
   }

   #Distribution for random effects
   for (j in 1:Npat) {
      mu[j,1:2] ~ dmnorm(mn1[], itaurand[,]);
   }

   #Distribution for area effects
   for (k in 1:Narea) {
      th[k,1:2] ~ dmnorm(mn2[],itauarea[,]);
   }

   #Prior for fixed effects
   for (k in 1:P) {
      be[k] ~ dnorm(0,1.0E-6);
   }
   
   #Hyper-priors
   isigma ~ dgamma(1.0E-3,1.0E-3);
   itaurand[1:2,1:2] ~ dwish(Omrand[,],2); 
   itauarea[1:2,1:2] ~ dwish(Omarea[,],2); #Make 2nd param = df
   df <- 4 #2, 4 or 13 for second parameter
   u ~ dbeta(0.5,0.5); #dbeta(1,1) = unif
   rho <- 2*u - 1;

   mn1[1] <- 0; mn1[2] <- 0;
   mn2[1] <- 0; mn2[2] <- 0;
   Omrand[1,1] <- 1; Omrand[1,2] <- 0; Omrand[2,1] <- 0; Omrand[2,2] <- 1;
   Omarea[1,1] <- OmareaMn; Omarea[1,2] <- 0; Omarea[2,1] <- 0; Omarea[2,2] <- OmareaMn;
   OmareaMn <- 1*(df - 3); #To ensure prior mean is 0.1, supply inverse of scaled mean to Wishart

   #Inverse of variances and covariance matrices
   sigma <- 1.0/isigma;
   taurand[1:2,1:2] <- inverse(itaurand[,]);
   tauarea[1:2,1:2] <- inverse(itauarea[,]);
   corrand <- taurand[1,2]/sqrt(taurand[1,1]*taurand[2,2]);
   corarea <- tauarea[1,2]/sqrt(tauarea[1,1]*tauarea[2,2]);

   #Derived quantities of interest
   iccind <- taurand[1,1]/(taurand[1,1]+sigma);
   bothfdist <- (1 - step(be[17]))*(1 - step(be[19]));
   onefdist <- 1 - step(be[17])*step(be[19]);
}
