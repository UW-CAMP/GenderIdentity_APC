### Prepare for BRFSS GI Analyses

# This script's code uses the combined BRFSS  dataset, previously prepared in script "SOGI_03_Prepare Combined Dataset"
### and prepares objects for the sexual orientation analyses in script "SOGI_07_GI Analyze"

### Prepare workspace ----- 
# clear environment
rm(list = ls())

# packages
source("SOGI_00_packages.R")

# define functions
tableNA <- function(x, ...){
   table(x, useNA = "ifany", ...)  
}

# call in data
brfss <- readRDS("data - clean/brfss_final.rds")
sab <- readRDS("data - clean/brfss_sab.rds")
sab_comp <- readRDS("data - clean/brfss_sab_compare.rds")

### Generate analytic objects -----
# * Multinomial GI arrays --------
# prepare objects for proportion arrays
gi_values <- levels(as.factor(brfss$gender))
sex_values <- levels(as.factor(brfss$sex))

cohort_n <- length(min(brfss$cohort):max(brfss$cohort))
period_n <- length(min(brfss$year):max(brfss$year))
cohort_v <- min(brfss$cohort):max(brfss$cohort) 
period_v <- min(brfss$year):max(brfss$year)

cohort <- list()
for (yr in cohort_v) {
   cohort[[yr]] <- subset(brfss, cohort == yr)  # if no rows, nothing to add
}


sab_values <- levels(as.factor(sab$sab))
sab_bin_values <- levels(as.factor(sab$sab_bin))

sab_period_n <- length(min(sab$year):max(sab$year))
sab_period_v <- min(sab$year):max(sab$year)
sab_cohort_n <- length(unique(sab$cohort))
sab_cohort_v <- sort(unique(sab$cohort))

sab_cohort <- list()
for (yr in sab_cohort_v) {
   sab_cohort[[yr]] <- subset(sab, cohort == yr)  # if no rows, nothing to add
}

sab_comp_cohort <- list()
for (yr in sab_cohort_v) {
   sab_comp_cohort[[yr]] <- subset(sab_comp, cohort == yr)  # if no rows, nothing to add
}

# create proportion array objects (A = all, F = females, M = males)
gi_props_A <- gi_props_F <- gi_props_M <- array(NA, dim=c(cohort_n, period_n, length(unique(brfss$gender))))

# fill the empty arrays with population-level estimates
for(cohort_xi in 1:cohort_n) {
   cohort_vi <- cohort_v[cohort_xi]
   period_x <- which(period_v %in% cohort[[cohort_vi]]$year)
   for (period_xi in period_x) {
      temp1 <- cohort[[cohort_vi]] %>% filter(year == period_v[period_xi]) 
      temp2 <- round(prop.table(questionr::wtd.table(x=temp1$age,
                                                     y=temp1$gender,
                                                     weights=temp1$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      gi_props_A[cohort_xi, period_xi,] <- temp2
      
      temp3 <- temp1 %>% filter(sex == "Female")
      temp4 <- round(prop.table(questionr::wtd.table(x=temp3$age,
                                                     y=factor(temp3$gender, levels = gi_values),
                                                     weights=temp3$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      gi_props_F[cohort_xi, period_xi,] <- temp4
      
      temp5 <- temp1 %>% filter(sex == "Male")
      temp6 <- round(prop.table(questionr::wtd.table(x=temp5$age,
                                                     y=factor(temp5$gender, levels = gi_values),
                                                     weights=temp5$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      gi_props_M[cohort_xi, period_xi,] <- temp6
   }
}

# * Multinomial sex-at-birth arrays ------
# * * Group A: X=gender, Y=sab (categ) -----
sab_props_tw <- sab_props_tm <- sab_props_nb <- sab_props_dkns <- sab_props_ref <- array(NA, dim=c(sab_cohort_n, sab_period_n, length(unique(sab$sab))))

for(sab_cohort_xi in 1:sab_cohort_n) {
   sab_cohort_vi <- sab_cohort_v[sab_cohort_xi]
   sab_period_x <- which(sab_period_v %in% sab_cohort[[sab_cohort_vi]]$year)
   for (sab_period_xi in sab_period_x) {
      temp1 <- sab_cohort[[sab_cohort_vi]] %>% filter(year == sab_period_v[sab_period_xi]) 

   # trans women
      temp3 <- temp1 %>% filter(gender == "1_transwoman")
      temp4 <- round(prop.table(questionr::wtd.table(x=temp3$age,
                                                   y=factor(temp3$sab, levels = sab_values),
                                                   weights=temp3$weight,
                                                   digits=1),
                              margin=1)*100, 3)
      if(nrow(temp4)>0) {
         sab_props_tw[sab_cohort_xi, sab_period_xi,] <- temp4
      } else {
         sab_props_tw[sab_cohort_xi, sab_period_xi,] <- NA
      }
   
   # trans men
      temp5 <- temp1 %>% filter(gender == "2_transman")
      temp6 <- round(prop.table(questionr::wtd.table(x=temp5$age,
                                                     y=factor(temp5$sab, levels = sab_values),
                                                     weights=temp5$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp6)>0) {
         sab_props_tm[sab_cohort_xi, sab_period_xi,] <- temp6
      } else {
         sab_props_tm[sab_cohort_xi, sab_period_xi,] <- NA
      }
   
   # nonbinary people
      temp7 <- temp1 %>% filter(gender == "3_nbgnc")
      temp8 <- round(prop.table(questionr::wtd.table(x=temp7$age,
                                                     y=factor(temp7$sab, levels = sab_values),
                                                     weights=temp7$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp8)>0) {
         sab_props_nb[sab_cohort_xi, sab_period_xi,] <- temp8
      } else {
         sab_props_nb[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
   # those who are unsure of their GI
      temp9 <- temp1 %>% filter(gender == "5_DNKS")
      temp10 <- round(prop.table(questionr::wtd.table(x=temp9$age,
                                                     y=factor(temp9$sab, levels = sab_values),
                                                     weights=temp9$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp10)>0) {
         sab_props_dnks[sab_cohort_xi, sab_period_xi,] <- temp10
      } else {
         sab_props_dkns[sab_cohort_xi, sab_period_xi,] <- NA
      }
   
   # those who did not disclose their GI
      temp11 <- temp1 %>% filter(gender == "6_ref")
      temp12 <- round(prop.table(questionr::wtd.table(x=temp11$age,
                                                      y=factor(temp11$sab, levels = sab_values),
                                                      weights=temp11$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp12)>0) {
         sab_props_ref[sab_cohort_xi, sab_period_xi,] <- temp12
      } else {
         sab_props_ref[sab_cohort_xi, sab_period_xi,] <- NA
      }
   }
}

# * * Group A: X=gender, Y=sab (binary) -----
sab_bin_props_tw <- sab_bin_props_tm <- sab_bin_props_nb <- sab_bin_props_dkns <- sab_bin_props_ref <- array(NA, dim=c(sab_cohort_n, sab_period_n, 2))

for(sab_cohort_xi in 1:sab_cohort_n) {
   sab_cohort_vi <- sab_cohort_v[sab_cohort_xi]
   sab_period_x <- which(sab_period_v %in% sab_cohort[[sab_cohort_vi]]$year)
   for (sab_period_xi in sab_period_x) {
      temp1 <- sab_cohort[[sab_cohort_vi]] %>% filter(year == sab_period_v[sab_period_xi]) 
      
      # trans women
      temp3 <- temp1 %>% filter(gender == "1_transwoman")
      temp4 <- round(prop.table(questionr::wtd.table(x=temp3$age,
                                                     y=factor(temp3$sab_bin, levels = sab_bin_values),
                                                     weights=temp3$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp4)>0) {
         sab_bin_props_tw[sab_cohort_xi, sab_period_xi,] <- temp4
      } else {
         sab_bin_props_tw[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # trans men
      temp5 <- temp1 %>% filter(gender == "2_transman")
      temp6 <- round(prop.table(questionr::wtd.table(x=temp5$age,
                                                     y=factor(temp5$sab_bin, levels = sab_bin_values),
                                                     weights=temp5$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp6)>0) {
         sab_bin_props_tm[sab_cohort_xi, sab_period_xi,] <- temp6
      } else {
         sab_bin_props_tm[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # nonbinary people
      temp7 <- temp1 %>% filter(gender == "3_nbgnc")
      temp8 <- round(prop.table(questionr::wtd.table(x=temp7$age,
                                                     y=factor(temp7$sab_bin, levels = sab_bin_values),
                                                     weights=temp7$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp8)>0) {
         sab_bin_props_nb[sab_cohort_xi, sab_period_xi,] <- temp8
      } else {
         sab_bin_props_nb[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # those who are unsure of their GI
      temp9 <- temp1 %>% filter(gender == "5_DNKS")
      temp10 <- round(prop.table(questionr::wtd.table(x=temp9$age,
                                                      y=factor(temp9$sab_bin, levels = sab_bin_values),
                                                      weights=temp9$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp10)>0) {
         sab_bin_props_dnks[sab_cohort_xi, sab_period_xi,] <- temp10
      } else {
         sab_bin_props_dkns[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # those who did not disclose their GI
      temp11 <- temp1 %>% filter(gender == "6_ref")
      temp12 <- round(prop.table(questionr::wtd.table(x=temp11$age,
                                                      y=factor(temp11$sab_bin, levels = sab_bin_values),
                                                      weights=temp11$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp12)>0) {
         sab_bin_props_ref[sab_cohort_xi, sab_period_xi,] <- temp12
      } else {
         sab_bin_props_ref[sab_cohort_xi, sab_period_xi,] <- NA
      }
   }
}


# * * Group B: X=gender, Y=sex among those who got SAB module -----
sab_module_props_tw <- sab_module_props_tm <- sab_module_props_nb <- sab_module_props_dkns <- sab_module_props_ref <- array(NA, dim=c(sab_cohort_n, sab_period_n, length(unique(sab$sex))))

for(sab_cohort_xi in 1:sab_cohort_n) {
   sab_cohort_vi <- sab_cohort_v[sab_cohort_xi]
   sab_period_x <- which(sab_period_v %in% sab_cohort[[sab_cohort_vi]]$year)
   for (sab_period_xi in sab_period_x) {
      temp1 <- sab_cohort[[sab_cohort_vi]] %>% filter(year == sab_period_v[sab_period_xi]) 
      
      # trans women
      temp3 <- temp1 %>% filter(gender == "1_transwoman")
      temp4 <- round(prop.table(questionr::wtd.table(x=temp3$age,
                                                     y=factor(temp3$sex, levels = sex_values),
                                                     weights=temp3$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp4)>0) {
         sab_module_props_tw[sab_cohort_xi, sab_period_xi,] <- temp4
      } else {
         sab_module_props_tw[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # trans men
      temp5 <- temp1 %>% filter(gender == "2_transman")
      temp6 <- round(prop.table(questionr::wtd.table(x=temp5$age,
                                                     y=factor(temp5$sex, levels = sex_values),
                                                     weights=temp5$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp6)>0) {
         sab_module_props_tm[sab_cohort_xi, sab_period_xi,] <- temp6
      } else {
         sab_module_props_tm[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # nonbinary people
      temp7 <- temp1 %>% filter(gender == "3_nbgnc")
      temp8 <- round(prop.table(questionr::wtd.table(x=temp7$age,
                                                     y=factor(temp7$sex, levels = sex_values),
                                                     weights=temp7$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp8)>0) {
         sab_module_props_nb[sab_cohort_xi, sab_period_xi,] <- temp8
      } else {
         sab_module_props_nb[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # those who are unsure of their GI
      temp9 <- temp1 %>% filter(gender == "5_DNKS")
      temp10 <- round(prop.table(questionr::wtd.table(x=temp9$age,
                                                      y=factor(temp9$sex, levels = sex_values),
                                                      weights=temp9$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp10)>0) {
         sab_module_props_dnks[sab_cohort_xi, sab_period_xi,] <- temp10
      } else {
         sab_module_props_dkns[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # those who did not disclose their GI
      temp11 <- temp1 %>% filter(gender == "6_ref")
      temp12 <- round(prop.table(questionr::wtd.table(x=temp11$age,
                                                      y=factor(temp11$sex, levels = sex_values),
                                                      weights=temp11$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp12)>0) {
         sab_module_props_ref[sab_cohort_xi, sab_period_xi,] <- temp12
      } else {
         sab_module_props_ref[sab_cohort_xi, sab_period_xi,] <- NA
      }
   }
}

# * * Group C: X=gender, Y=sex (not sab, among those who did not receive SAB module) -----
sab_comp_props_tw <- sab_comp_props_tm <- sab_comp_props_nb <- sab_comp_props_dkns <- sab_comp_props_ref <- array(NA, dim=c(sab_cohort_n, sab_period_n, length(unique(sab_comp$sex))))

for(sab_cohort_xi in 1:sab_cohort_n) {
   sab_cohort_vi <- sab_cohort_v[sab_cohort_xi]
   sab_period_x <- which(sab_period_v %in% sab_comp_cohort[[sab_cohort_vi]]$year)
   for (sab_period_xi in sab_period_x) {
      temp1 <- sab_comp_cohort[[sab_cohort_vi]] %>% filter(year == sab_period_v[sab_period_xi]) 
      
      # trans women
      temp3 <- temp1 %>% filter(gender == "1_transwoman")
      temp4 <- round(prop.table(questionr::wtd.table(x=temp3$age,
                                                     y=factor(temp3$sex, levels = sex_values),
                                                     weights=temp3$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp4)>0) {
         sab_comp_props_tw[sab_cohort_xi, sab_period_xi,] <- temp4
      } else {
         sab_comp_props_tw[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # trans men
      temp5 <- temp1 %>% filter(gender == "2_transman")
      temp6 <- round(prop.table(questionr::wtd.table(x=temp5$age,
                                                     y=factor(temp5$sex, levels = sex_values),
                                                     weights=temp5$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp6)>0) {
         sab_comp_props_tm[sab_cohort_xi, sab_period_xi,] <- temp6
      } else {
         sab_comp_props_tm[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # nonbinary people
      temp7 <- temp1 %>% filter(gender == "3_nbgnc")
      temp8 <- round(prop.table(questionr::wtd.table(x=temp7$age,
                                                     y=factor(temp7$sex, levels = sex_values),
                                                     weights=temp7$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      if(nrow(temp8)>0) {
         sab_comp_props_nb[sab_cohort_xi, sab_period_xi,] <- temp8
      } else {
         sab_comp_props_nb[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # those who are unsure of their GI
      temp9 <- temp1 %>% filter(gender == "5_DNKS")
      temp10 <- round(prop.table(questionr::wtd.table(x=temp9$age,
                                                      y=factor(temp9$sex, levels = sex_values),
                                                      weights=temp9$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp10)>0) {
         sab_comp_props_dnks[sab_cohort_xi, sab_period_xi,] <- temp10
      } else {
         sab_comp_props_dkns[sab_cohort_xi, sab_period_xi,] <- NA
      }
      
      # those who did not disclose their GI
      temp11 <- temp1 %>% filter(gender == "6_ref")
      temp12 <- round(prop.table(questionr::wtd.table(x=temp11$age,
                                                      y=factor(temp11$sex, levels = sex_values),
                                                      weights=temp11$weight,
                                                      digits=1),
                                 margin=1)*100, 3)
      if(nrow(temp12)>0) {
         sab_comp_props_ref[sab_cohort_xi, sab_period_xi,] <- temp12
      } else {
         sab_comp_props_ref[sab_cohort_xi, sab_period_xi,] <- NA
      }
   }
}


# * * X=sab, Y=gender ------
sab_props_F <- sab_props_M <- array(NA, dim=c(sab_cohort_n, sab_period_n, length(unique(sab$gender))))

for(sab_cohort_xi in 1:sab_cohort_n) {
   sab_cohort_vi <- sab_cohort_v[sab_cohort_xi]
   sab_period_x <- which(sab_period_v %in% sab_cohort[[sab_cohort_vi]]$year)
   for (sab_period_xi in sab_period_x) {
      temp1 <- sab_cohort[[sab_cohort_vi]] %>% filter(year == sab_period_v[sab_period_xi]) 
      tempX <- temp1 %>% filter (sab == "2_female")
      temp2 <- round(prop.table(questionr::wtd.table(x=tempX$age,
                                                     y=factor(tempX$gender, levels = gi_values),
                                                     weights=tempX$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      sab_props_F[sab_cohort_xi, sab_period_xi,] <- temp2
      
      tempY <- temp1 %>% filter(sab == "1_male") 
      temp4 <- round(prop.table(questionr::wtd.table(x=tempY$age,
                                                     y=factor(tempY$gender, levels = gi_values),
                                                     weights=tempY$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      sab_props_M[sab_cohort_xi, sab_period_xi,] <- temp4
      
   }
}

# * * X=sex (not sab), Y=gender ----
sab_comp_props_F <- sab_comp_props_M <- array(NA, dim=c(sab_cohort_n, sab_period_n, length(unique(sab_comp$gender))))

for(sab_cohort_xi in 1:sab_cohort_n) {
   sab_cohort_vi <- sab_cohort_v[sab_cohort_xi]
   sab_period_x <- which(sab_period_v %in% sab_comp_cohort[[sab_cohort_vi]]$year)
   for (sab_period_xi in sab_period_x) {
      temp1 <- sab_comp_cohort[[sab_cohort_vi]] %>% filter(year == sab_period_v[sab_period_xi]) 
      tempX <- temp1 %>% filter (sex == "Female")
      temp2 <- round(prop.table(questionr::wtd.table(x=tempX$age,
                                                     y=factor(tempX$gender, levels = gi_values),
                                                     weights=tempX$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      sab_comp_props_F[sab_cohort_xi, sab_period_xi,] <- temp2
      
      tempY <- temp1 %>% filter(sex == "Male") 
      temp4 <- round(prop.table(questionr::wtd.table(x=tempY$age,
                                                     y=factor(tempY$gender, levels = gi_values),
                                                     weights=tempY$weight,
                                                     digits=1),
                                margin=1)*100, 3)
      sab_comp_props_M[sab_cohort_xi, sab_period_xi,] <- temp4
      
   }
}

### Save analytic objects for analysis as .rds files --------
# among all respondents (2014:2021), gender identity props by sex (not SAB)
write_rds(gi_props_A, "data - clean/gi_props_A.rds")
write_rds(gi_props_F, "data - clean/gi_props_F.rds")
write_rds(gi_props_M, "data - clean/gi_props_M.rds")

# among Rs who in states who got SAB module, gender identity props by SAB
write_rds(sab_props_M, "data - clean/sab_props_M.rds")
write_rds(sab_props_F, "data - clean/sab_props_F.rds")

# among Rs who in states who did not get SAB module, gender identity props by sex (not SAB)
write_rds(sab_comp_props_M, "data - clean/sab_comp_props_M.rds")
write_rds(sab_comp_props_F, "data - clean/sab_comp_props_F.rds")

# among Rs in states who got SAB module, SAB props (categ), by gender identity
write_rds(sab_props_tw, "data - clean/sab_props_tw.rds")
write_rds(sab_props_tm, "data - clean/sab_props_tm.rds")
write_rds(sab_props_nb, "data - clean/sab_props_nb.rds")
write_rds(sab_props_dkns, "data - clean/sab_props_dkns.rds")
write_rds(sab_props_ref, "data - clean/sab_props_ref.rds")

# among Rs in states who got SAB module, SAB props (binary), by gender identity
write_rds(sab_bin_props_tw, "data - clean/sab_bin_props_tw.rds")
write_rds(sab_bin_props_tm, "data - clean/sab_bin_props_tm.rds")
write_rds(sab_bin_props_nb, "data - clean/sab_bin_props_nb.rds")
write_rds(sab_bin_props_dkns, "data - clean/sab_bin_props_dkns.rds")
write_rds(sab_bin_props_ref, "data - clean/sab_bin_props_ref.rds")

# among Rs in states who got SAB module, sex (*not* SAB) props, by gender identity
write_rds(sab_module_props_tw, "data - clean/sab_module_props_tw.rds")
write_rds(sab_module_props_tm, "data - clean/sab_module_props_tm.rds")
write_rds(sab_module_props_nb, "data - clean/sab_module_props_nb.rds")
write_rds(sab_module_props_dkns, "data - clean/sab_module_props_dkns.rds")
write_rds(sab_module_props_ref, "data - clean/sab_module_props_ref.rds")

# among Rs in states who did not get SAB module, sex props, by gender identity
write_rds(sab_comp_props_tw, "data - clean/sab_comp_props_tw.rds")
write_rds(sab_comp_props_tm, "data - clean/sab_comp_props_tm.rds")
write_rds(sab_comp_props_nb, "data - clean/sab_comp_props_nb.rds")
write_rds(sab_comp_props_dkns, "data - clean/sab_comp_props_dkns.rds")
write_rds(sab_comp_props_ref, "data - clean/sab_comp_props_ref.rds")