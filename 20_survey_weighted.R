## Code for getting Horvitz-Thompson/Hajek estimators for
## prevalence by age x period (survey) x cohort, and age x period (survey)
## of the following manuscript:
##
##    Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025.
##      Demographic trends in gender identity among adults in the United States,
##      2014-2021. International Journal of Transgender Health. 
##      Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.
## 
##  Script authors: Godwin J
##
##  Inputs:
##     -  data - clean/brfss_final.rds
##
##  Outputs:
##     -  data - clean/BRFSS_HT_GI_prevs.csv
##     -  data - clean/BRFSS_HT_AP_GI_prevs.csv
## 


# Setup ####
rm(list = ls())

## Functions ####
tableNA <- function(x, ...){
  table(x, useNA = "ifany", ...)  
}

logit <- function(prob){
  log(prob/(1-prob))
}

expit <- function(logodds){
  exp(logodds)/(1 + exp(logodds))
}

# Load Data ####

brfss <- readRDS(file = "data - clean/brfss_final.rds")

str(brfss)
summary(brfss)
head(brfss)

tableNA(brfss$gender)
tableNA(brfss$sab)
  
# Prep for estimation ####
## Filter to analysis ages ####

data_1x1 <- brfss %>%  
  select(year, age, cohort,
         sex, gender, gender, so, so_new, sab, sex,
         contains("cohort_"), contains("pd_"),
         index, weight, strata) %>% 
  mutate(tw_bin = ifelse(gender == "1_transwoman", 1, 0),
         tm_bin = ifelse(gender == "2_transman", 1, 0),
         nbgnc_bin = ifelse(gender == "3_nbgnc", 1, 0),
         cis_bin = ifelse(gender == "4_cisgender", 1, 0),
         dkns_bin = ifelse(gender == "5_DKNS", 1, 0),
         ref_bin = ifelse(gender == "6_ref", 1, 0)) %>%
  mutate(across(where(is.factor), ~as.character(.x))) %>% 
  rename("period" = "year", "stratum" = "strata")
  
## Check out object
summary(data_1x1)
tableNA(data_1x1$age, data_1x1$cohort)
tableNA(data_1x1$period, data_1x1$cohort)
tableNA(data_1x1$age, data_1x1$period)

## Overall TNB proportion - simple point estimate by year

data_1x1 <- data_1x1 %>% 
  mutate(tnb = ifelse(gender %in% c("1_transwoman",
                                    "2_transman",
                                    "3_nbgnc"), 1, 0)
  )
tableNA(data_1x1$gender, data_1x1$tnb)
periods <- unique(data_1x1$period)
sapply(periods, function(x) 
  sum(data_1x1$weight[data_1x1$period==x & data_1x1$tnb==1]) / 
  sum(data_1x1$weight[data_1x1$period==x])
  )
              
## Specify design object ####

brfss_des <- svydesign(ids = ~1, strata = ~period + stratum,
                       weights = ~weight, data = data_1x1)

brfss_des

# Estimate prevalence ####

prevs <- svyby(~tw_bin + tm_bin + nbgnc_bin + dkns_bin + ref_bin,
               by = ~sex + period + age + cohort,
               design = brfss_des,
               svymean)

prevs_age_per <- svyby(~tw_bin + tm_bin + nbgnc_bin + dkns_bin + ref_bin,
                       by = ~sex + period + age,
                       design = brfss_des,
                       svymean)


# Clean up output & Save ####

## Fix names ####
names(prevs)[names(prevs) %in% paste(c("tw","tm", "nbgnc", "dkns", "ref"),
                                    "bin", sep = "_")] <- 
  paste0("mean.", c("tw","tm", "nbgnc", "dkns", "ref"), "_bin")

names(prevs)

names(prevs_age_per)[names(prevs_age_per) %in% 
                           paste(c("tw","tm", "nbgnc", "dkns", "ref"), "bin",
                                 sep = "_")] <- 
  paste0("mean.", c("tw","tm", "nbgnc", "dkns", "ref"), "_bin")

names(prevs_age_per)

## pivot_longer & save ####
prevs_long <- prevs %>% 
  pivot_longer(cols = contains("_bin"),
               names_to = c(".value", "gender"),
               names_pattern = "(.*)\\.(.*)_bin") %>% 
  ## Calculate logit(phat)
  ## And var(logit(phat)) = var(phat)/[phat^2*(1-phat)^2]
  ## via delta method
  mutate(logit_pijk = logit(mean),
         var = se^2,
         var_logit_pijk = var/(mean^2*(1-mean)^2))

write.csv(prevs_long, row.names = FALSE,
          file = paste0("data - clean/", "BRFSS_HT_GI_prevs.csv"))


prevs_age_per_long <- prevs_age_per %>% 
  pivot_longer(cols = contains("_bin"),
               names_to = c(".value", "gender"),
               names_pattern = "(.*)\\.(.*)_bin") %>% 
  ## Calculate logit(phat)
  ## And var(logit(phat)) = var(phat)/[phat^2*(1-phat)^2]
  ## via delta method
  mutate(logit_pijk = logit(mean),
         var = se^2,
         var_logit_pijk = var/(mean^2*(1-mean)^2))

write.csv(prevs_age_per_long, row.names = FALSE,
          file = paste0("data - clean/", "BRFSS_HT_AP_GI_prevs.csv"))
