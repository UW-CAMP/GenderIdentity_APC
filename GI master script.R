### This code runs the entire YRBS/BRFSS sexual orientation, gender identity, and sex-of-sex-partners analysis (CAMP 60).
# Goodreau SM, Barry MP
# University of Washington

# Abbreviations:
# SO: sexual orientation
# GI: gender identity
# SSP: sex(es) of sex partners

# SAS analyses, address in the corresponding ReadMe file must precede this program steps
rm(list=ls())
set.seed(1)

# Prepare project files
source("00_packages.R") # call in packages; you may have to open this file and install new packages 
source("01_BRFSS_index_variable.R") # !! SLOW !! creates manageable versions of BRFSS data (.rds rather than .XPT) - only run once to generate these files
source("02_prepare_BRFSS.R") # run this file to prepare a single file `brfss_final.rds`

# GI Analyses
source("21_survey_weighted.R") # Get design based prev estimates
source("22_descriptive_plots.R") # Generage Figs 1, S1, and S2
source("23_analysis_1.R") # Fit null, A, P, C, AP, AC, PC models for tw, tm, nbgnc
source("24_analysis_2.R") # Fit Bayesian mAPC models for tw, tm, nbgnc
source("25_analysis_3.R") # Fit null, A, P, C, AP, AC, PC models for dkns, dta
source("26_analysis_4.R") # Sex, Sex assigned at birth, and So by GI
