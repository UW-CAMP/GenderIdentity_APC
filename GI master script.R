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
source("SOGI_00_packages.R") # call in packages; you may have to open this file and install new packages 
source("SOGI_01_generate_BRFSS_index_variable.R") # !! SLOW !! creates manageable versions of BRFSS data (.rds rather than .XPT) - only run once to generate these files
source("SOGI_02_prepare_BRFSS.R") # run this file to prepare a single file `brfss_final.rds`

# GI Analyses
# source("SOGI_21_GI_prepare.R") # prepare objects for GI analyses
# source("SOGI_22_GI_analyses.R") # GI analyses
source("SOGI_21_GI_design_based.R") # Get design based prev estimates
source("SOGI_22_GI_Q1.R") # Fit null, A, P, C, AP, AC, PC models
source("SOGI_22_GI_Q2_mAPC.R") # Fit Bayesian mAPC models for tw, tm, nbgnc