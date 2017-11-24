# Proximity To Food Synthetic Data Analysis

Data from the Framingham Heart Study is not freely available due to patient confidentiality. The Framingham Heart Study has a formal process for requesting data, and anyone can request data directly from this via: https://www.framinghamheartstudy.org/researchers/application-review.php. In lieu of the actual data, we provide a synthetic data set named simdata.txt. The size and structure of these data are identical to the actual data set, allowing the exact same code to be used to perform the Bayesian analysis on it. The files and data are described in the following and in the file DualNHoodBiosSuppAppend.pdf in PDF format.

Relevant files

 1) R script: RealBMIJAGS.R
 2) JAGS model code: RealBMICodeIntIW.bug
 3) Synthetic data set: simdata.txt
 4) DualNHoodBiosSuppAppend.pdf

Data dictionary

 1) ID: The identification number of the individual. Ranges from 1 to 2,889.
 2) wave: The Framingham Heart Study exam wave; ranges from 1 to 8 for the Offspring cohort and 6 to 8 for the Omni cohort, which only started in wave 6 of the offspring cohort
 3) Tobs: The total number of exams the individual attended
 4) TractHome: The census tract where the individual lived
 5) TractWork: The census tract where an individual was employed (if employed)
 6) BMI: Body Mass Index
 7) yob: Year of birth of individual
 8) smokes: Whether individual is a current smoker
 9) male: Whether individual has male gender
 10) married: Whether individual is currently married
 11) educ: Educational level of individual (0 = High-school or less, 1 = Completed high school; 2 = Other)
 12) tractpov: Percent of households below the poverty line in census tract where individual lived
 13) etractpov: Percent of households below the poverty line in census tract where individual was employed
 14) unemploy: Whether or not individual was unemployed
 15) DistHome: Distance in kilometers from individual's home to the nearest fast food establishment
 16) DistWork: Distance in kilometers from individual's workplace to the nearest fast food establishment (if employed)
 17) DriveDist: Number of fast food restaurants within a 60 meter buffer of the shortest commute between work and home

The first row of the data set contains the above variable names. Workplace distance, driving distance and the poverty of an individual's workplace census tract are only available for employed individuals. For simplicity, real numerical values are given for these variables on such cases. However, the JAGS script is coded so that only observations in which the individual is employed directly contribute to the estimation of model parameters involving these variables.

Instructions

The R script contains the commands to setup the analysis, call JAGS, and perform post-estimation analysis and presentation of the results. The JAGS code performs the model estimation. The entire analysis is performed by running the R code as JAGS is called automatically. The steps are as follows:
 1) Place the R script and the JAGS file in the directory from which you want to perform the analysis
 2) Ensure that the synthetic data set (simdata.txt) is in the directory where the R script calls it
 3) Run the R script by typing source('RealBMIJAGS.R') or using an equivalent method
