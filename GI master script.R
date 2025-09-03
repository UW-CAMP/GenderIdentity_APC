## Code for all analysis producing each table and figure 
## of the following manuscript and its supplement:
##
##    Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025.
##      Demographic trends in gender identity among adults in the United States,
##      2014-2021. International Journal of Transgender Health. 
##      Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.
## 
##  Script authors: Barry MP, Godwin J, Goodreau SM
##
##  Notes:
##    - Packages: 00_packages.R must be run before any single other script
##    - BRFSS cleaning: 01-02 must be run before the analysis scripts
##    - Survey-weighted APC: 20 must be run before scripts 21-23
##    - Bayesian smoothing models: 21-23 can be run independently after 20
##    - Sex, sex assigned at birth, sexual orientation by gender identity:
##      25 can be run directly after 02, and does not depend on 20
# Setup ####
rm(list = ls())

# Libraries ####
source("00_packages.R") # call in packages; you may have to open this file and install new packages 

# Load & clean data ####
source("01_BRFSS_index_variable.R") # !! SLOW !! creates manageable versions of BRFSS data (.rds rather than .XPT) - only run once to generate these files
source("02_prepare_BRFSS.R") # run this file to prepare a single file `brfss_final.rds`

# Analyses ####
source("20_survey_weighted.R") # Get survey-weighted GI prevalence estimates
source("21_descriptive_plots.R") # Generate Figs 1, S1, and S2
source("22_one_two_factor_GI.R") # Generate Fig 2, Table 2
source("23_mAPC_GI.R") # Generate Fig 3
source("24_one_two_factor_other.R") # Generate Fig 4
source("25_sex_SAB_SO_GI.R") # Generate Figs 5, 6, Table S2
