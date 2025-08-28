## Code loading packages necessary for all analyses of the following manuscript:
##
##    Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025.
##      Demographic trends in gender identity among adults in the United States,
##      2014-2021. International Journal of Transgender Health. 
##      Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.
## 
##  Script authors: Barry MP, Godwin J, Goodreau SM
##
##  This script takes large .XPT BRFSS files, assigns index numbers
##  to each BRFSS observation, and generates .rds files for input to 
##  02_prepare_BRFSS.R
##
##  Inputs:
##    - data - raw/LLCP20YY.XPT for YY in 14:21
##
##   Outputs:
##    - data - clean/brfss.rds
##    - data - clean/brfssYY.rds for YY in 14:21

# Setup ####
rm(list = ls())

## Functions ####
tableNA <- function(x, ...){
   table(x, useNA = "ifany", ...)  
}

# Load raw BRFSS ####
brfss <- list()
for (yr in 14:21) {
   cat('Loading year', yr, '\n')
   filename <- paste("data - raw/LLCP20", yr, ".XPT ", sep="")
   brfss[[yr]] <- read_xpt(filename)
   #note: the spaces after ".XPT" should be removed if running in a Windows Environment.
   brfss[[yr]]$index <- as.numeric(paste0(1:nrow(brfss[[yr]]), yr))
}

# Save outputs ####
write_rds(brfss, "data - clean/brfss.rds")
write_rds(brfss[[14]], "data - clean/brfss14.rds")
write_rds(brfss[[15]], "data - clean/brfss15.rds")
write_rds(brfss[[16]], "data - clean/brfss16.rds")
write_rds(brfss[[17]], "data - clean/brfss17.rds")
write_rds(brfss[[18]], "data - clean/brfss18.rds")
write_rds(brfss[[19]], "data - clean/brfss19.rds")
write_rds(brfss[[20]], "data - clean/brfss20.rds")
write_rds(brfss[[21]], "data - clean/brfss21.rds")