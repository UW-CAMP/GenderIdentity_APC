##### This script contains all packages used by GI Analyses -----

############################################################### 
# Un-comment and run these lines for any packages not already downloaded
############################################################### 
#install.packages("visdat") #used to visualize missing data patterns in script 01
#install.packages("tidyverse") #used for data management
#install.packages("haven") #used for data management
#install.packages("dplyr") #used for data management
#install.packages("tidyr") #used for data management
#install.packages("tidyselect") #used for data management
#install.packages("DescTools") # for obtaining 95% CI for multinomial proportions ("comparing 18 YO in BRFSS, YRBS")
#install.packages("magrittr") # used for piping
#install.packages("questionr") # used for applying weights to tables in BRFSS and YRBS
#install.packages("colorspace") # used for colors in graphs
#install.packages("dichromat") # used for colors in graphs
#install.packages("labelled") # used for data management
#install.packages("weights") # for using chi squared test with weighted data 
#install.packages("survey") # for weighted prev estimates and 95% CI ("comparing 18 YO")
#install.packages("ggplot2")
#install.packages("ggstream") # for proportions over time plots
#install.packages("trend") # for Mann-Kendall test of monotonicty
#install.packages("multiCA") # multinomial trend test
#install.packages("gridExtra")

############################################################### 
# Load packages
############################################################### 
library(visdat) #used to visualize missing data patterns in script 01
library(tidyverse) #used for data management
library(haven) #used for data management
library(dplyr) #used for data management
library(tidyr) #used for data management
library(tidyselect) #used for data management
library(DescTools) # for obtaining 95% CI for multinomial proportions
library(magrittr) # used for piping
library(questionr) # used for applying weights to tables in BRFSS and YRBS
library(colorspace) # used for colors in graphs
library(dichromat) # used for colors in graphs
library(labelled) # used for data management
library(weights) # for using chi squared test with weighted data 
library(survey) # for weighted prev estimates and 95% CI (comparing 18 YO)
library(ggplot2)
library(ggstream) # for proportions over time plots
library(trend) # for Mann-Kendall test of monotonicty
library(multiCA) # multinomial trend test
library(gridExtra)
