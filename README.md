# Proximity To Food Synthetic Data Analysis

Data from the Framingham Heart Study is not freely available due to patient confidentiality. The Framingham Heart Study has a formal process for requesting data, and anyone can request data directly from this via: https://www.framinghamheartstudy.org/researchers/application-review.php.  

In lieu of the actual data, we provide a pseudo data set named simdata4wbugs.txt. The size and structure of these data are identical to the actual data set, allowing the exact same code to be used to perform the Bayesian analysis on it. The data is described in the following.

The following files are relevant:
R script: RealBMIBRugs.r;
WinBUGS model code: RealBMICodeInt.r;
WinBUGS constant data: RealBMIConst4.txt;
WinBUGS Synthetic data set: simdata4wbugs.txt;
WinBUGS initial values: RealBMIInits.txt;

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

The first row of the data set contains the above variable names. Workplace distance, driving distance and the poverty of an individual's workplace census tract are only available for employed individuals. For simplicity, real numerical values are given for these variables on such cases. However, the WinBUGS script is coded so that only observations in which the individual is employed directly contribute to the estimation of model parameters involving these variables.

Instructions for running code:

The R script contains the commands to setup the analysis, call WinBUGS, and perform post-estimation analysis and presentation of the results. The WinBUGS code performs for the model estimation. The entire analysis is performed by running the R code as WinBUGS is called automatically. The steps are as follows:

1) Place the R script and WinBUGS files in the directory from which you want to perform the analysis
2) Ensure that the synthetic data set is in the directory where the R script calls it
3) Run the R script by typing source('RealBMIBRugs.r') or using an equivalent method.
