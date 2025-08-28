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
source("20_survey_weighted.R") # Get design based prev estimates
source("21_descriptive_plots.R") # Generate Figs 1, S1, and S2
source("22_analysis_1.R") # Generate Fig 2, Table 2
source("23_analysis_2.R") # Generate Fig 3
source("24_analysis_3.R") # Generate Fig 4
source("25_analysis_4.R") # Generate Figs 5, 6, Table S2
