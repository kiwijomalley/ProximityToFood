# Proximity To Food Synthetic Data Analysis

Data from the Framingham Heart Study is not freely available due to patient confidentiality. The Framingham Heart Study has a formal process for requesting data, and anyone can request data directly from this via: https://www.framinghamheartstudy.org/researchers/application-review.php.  

In lieu of the actual data, we provide a synthetic data set named simdata.txt. The size and structure of these data are identical to the actual data set, allowing the exact same code to be used to perform the Bayesian analysis on it. The files and data are described in the following.

Relevant files

 R script: RealBMIJAGS.R
 JAGS model code: RealBMICodeIntIW.bug
 Synthetic data set: simdata.txt

Data dictionary

 ID: The identification number of the individual. Ranges from 1 to 2,889.
 wave: The Framingham Heart Study exam wave; ranges from 1 to 8 for the Offspring cohort and 6 to 8 for the Omni cohort, which only started in wave 6 of the offspring cohort
 Tobs: The total number of exams the individual attended
 TractHome: The census tract where the individual lived
 TractWork: The census tract where an individual was employed (if employed)
 BMI: Body Mass Index
 yob: Year of birth of individual
 smokes: Whether individual is a current smoker
 male: Whether individual has male gender
 married: Whether individual is currently married
 educ: Educational level of individual (0 = High-school or less, 1 = Completed high school; 2 = Other)
 tractpov: Percent of households below the poverty line in census tract where individual lived
 etractpov: Percent of households below the poverty line in census tract where individual was employed
 unemploy: Whether or not individual was unemployed
 DistHome: Distance in kilometers from individual's home to the nearest fast food establishment
 DistWork: Distance in kilometers from individual's workplace to the nearest fast food establishment (if employed)
 DriveDist: Number of fast food restaurants within a 60 meter buffer of the shortest commute between work and home

The first row of the data set contains the above variable names. Workplace distance, driving distance and the poverty of an individual's workplace census tract are only available for employed individuals. For simplicity, real numerical values are given for these variables on such cases. However, the JAGS script is coded so that only observations in which the individual is employed directly contribute to the estimation of model parameters involving these variables.

Instructions:

The R script contains the commands to setup the analysis, call JAGS, and perform post-estimation analysis and presentation of the results. The JAGS code performs the model estimation. The entire analysis is performed by running the R code as JAGS is called automatically. The steps are as follows:

 1) Place the R script and the JAGS file in the directory from which you want to perform the analysis
 2) Ensure that the synthetic data set (simdata.txt) is in the directory where the R script calls it
 3) Run the R script by typing source('RealBMIJAGS.R') or using an equivalent method
