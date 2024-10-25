## Fit Multivariate APC 
## Prevalence by age x period (survey) x cohort
## Jessica Godwin

# Setup ####
rm(list = ls())

## Directories ####
## Create tables and plots directory for results if they don't exist
if(!dir.exists("tables/GI_Q2")){
  if(!dir.exists("tables/")){
    dir.create("tables/")
  }
  dir.create("tables/GI_Q2")
}

if(!dir.exists("plots/GI_Q2")){
  if(!dir.exists("plots/")){
    dir.create("plots/")
  }
  dir.create("plots/GI_Q2")
}

## Libraries ####
# install.packages("INLA",
#                  repos=c(getOption("repos"),
#                          INLA="https://inla.r-inla-download.org/R/stable"))
library(INLA) 
library(tidyverse)


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

prev_data <- read.csv("data - clean/BRFSS_HT_GI_prevs.csv")
str(prev_data)
summary(prev_data)
head(prev_data)


# Prep for estimation ####

mod_data <- prev_data %>%  
  arrange(age, period, cohort) %>% 
  mutate(period_fac = factor(period, levels = 2014:2021),
         age_fac = factor(age, levels = 18:79),
         cohort_fac = factor(cohort, levels = 1935:2003),
         sex_fac = factor(sex),
         period_idx = as.numeric(period_fac),
         age_idx = as.numeric(age_fac),
         cohort_idx = as.numeric(cohort_fac),
         sex_idx = as.numeric(sex_fac),
         prec_logit_pijk = 1/var_logit_pijk) %>% 
  rename("gender" = "so") %>% 
  filter(gender %in% c("tw", "tm", "nbgnc")) %>% 
  mutate(gender_idx = case_when(gender == "tw" ~ 1,
                                gender == "tm" ~ 2,
                                TRUE ~ 3)) 

## Check work
coh_summary <- mod_data %>% 
  group_by(cohort) %>% tally()
age_summary <- mod_data %>% 
  group_by(age) %>% tally()
per_summary <- mod_data %>% 
  group_by(period) %>% tally()

# Fit models ####

## Formulas ####
## Create formula objects suitable for INLA for the following model
## logit_pijks =μ_strata +θ_age +φ_period,strata +ψ_cohort,strata.

pc_prec <- list(theta = list(prior = 'pc.prec', param=c(1, .5)))
mod_formula_shareco <- logit_pijk ~ (-1) + gender + 
  f(age_idx, model = "rw2", replicate = gender_idx,
    scale.model = TRUE, constr = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(period_idx, model = "rw2", replicate = gender_idx,
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(cohort_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec)

mod_formula_shareper <- logit_pijk ~ (-1) + gender + 
  f(age_idx, model = "rw2", replicate = gender_idx,
    scale.model = TRUE, constr = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(period_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(cohort_idx, model = "rw2", replicate = gender_idx,
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec)

mod_formula_shareage <- logit_pijk ~ (-1) + gender + 
  f(age_idx, model = "rw2", 
    scale.model = TRUE, constr = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(period_idx, model = "rw2", replicate = gender_idx,
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(cohort_idx, model = "rw2", replicate = gender_idx,
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec)


## Lincombs ####

### Shared Cohort ####

#### Age ####
lincombs_info_age_shareco <-
  data.frame(Index = 1:(max(mod_data$age_idx)*3),
             Gender = rep(c(1,2,3), each = max(mod_data$age_idx)),
             Term = rep(1:max(mod_data$age_idx), 3),
             TermName = "Age",
             SharedName = "Cohort",
             LC_Index = NA)
row_idx <- 0
for(row in 1:(nrow(lincombs_info_age_shareco))){
  if(lincombs_info_age_shareco$Term[row] %in%
     2:(max(lincombs_info_age_shareco$Term) -1)){
    row_idx <- row_idx + 1
    lincombs_info_age_shareco$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    age <- rep(NA, max(mod_data$age_idx)*3)
    age_element <- lincombs_info_age_shareco$Index[row]
    age[age_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(age_idx = age))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareco <- get(object.name)
      names(lincombs_shareco)[row_idx] <- object.name
    }else{
      lincombs_shareco <- c(lincombs_shareco, get(object.name))
      names(lincombs_shareco)[row_idx] <- object.name
    }
    
    ## Keep your environment clean
    rm(object.name)
  }
}

#### Period ####
lincombs_info_per_shareco <-
  data.frame(Index = 1:(max(mod_data$period_idx)*3),
             Gender = rep(c(1,2,3), each = max(mod_data$period_idx)),
             Term = rep(1:max(mod_data$period_idx), 3),
             TermName = "Period",
             SharedName = "Cohort",
             LC_Index = NA)

for(row in 1:(nrow(lincombs_info_per_shareco))){
  if(lincombs_info_per_shareco$Term[row] %in%
     2:(max(lincombs_info_per_shareco$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_per_shareco$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    period <- rep(NA, max(mod_data$period_idx)*3)
    period_element <- lincombs_info_per_shareco$Index[row]
    period[period_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(period_idx = period))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareco <- get(object.name)
      names(lincombs_shareco)[row_idx] <- object.name
    }else{
      lincombs_shareco <- c(lincombs_shareco, get(object.name))
      names(lincombs_shareco)[row_idx] <- object.name
    }
    
    rm(object.name)
  }
}

#### Cohort ####
lincombs_info_coh_shareco <-
  data.frame(Index = 1:max(mod_data$cohort_idx),
             Gender = rep(NA, each = max(mod_data$cohort_idx)),
             Term = rep(1:max(mod_data$cohort_idx), 1),
             TermName = "Cohort",
             SharedName = "Cohort",
             LC_Index = NA)

for(row in 1:(nrow(lincombs_info_coh_shareco))){
  if(lincombs_info_coh_shareco$Term[row] %in% 
     2:(max(lincombs_info_coh_shareco$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_coh_shareco$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    cohort <- rep(NA, max(mod_data$cohort_idx)*1)
    cohort_element <- lincombs_info_coh_shareco$Index[row]
    cohort[cohort_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(cohort_idx = cohort))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareco <- get(object.name)
      names(lincombs_shareco)[row_idx] <- object.name
    }else{
      lincombs_shareco <- c(lincombs_shareco, get(object.name))
      names(lincombs_shareco)[row_idx] <- object.name
    }
    
    rm(object.name)
  }
}

lincombs_info_shareco <- bind_rows(lincombs_info_age_shareco,
                                   lincombs_info_per_shareco,
                                   lincombs_info_coh_shareco)

### Shared Age ####
#### Cohort ####
lincombs_info_coh_shareage <-
  data.frame(Index = 1:(max(mod_data$cohort_idx)*3),
             Gender = rep(c(1,2,3), each = max(mod_data$cohort_idx)),
             Term = rep(1:max(mod_data$cohort_idx), 3),
             TermName = "Cohort",
             SharedName = "Age",
             LC_Index = NA)

row_idx <- 0
for(row in 1:(nrow(lincombs_info_coh_shareage))){
  if(lincombs_info_coh_shareage$Term[row] %in%
     2:(max(lincombs_info_coh_shareage$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_coh_shareage$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    cohort <- rep(NA, max(mod_data$cohort_idx)*3)
    coh_element <- lincombs_info_coh_shareage$Index[row]
    cohort[coh_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(cohort_idx = cohort))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareage <- get(object.name)
      names(lincombs_shareage)[row_idx] <- object.name
    }else{
      lincombs_shareage <- c(lincombs_shareage, get(object.name))
      names(lincombs_shareage)[row_idx] <- object.name
    }
  }
}  
  
#### Period ####
lincombs_info_per_shareage <-
  data.frame(Index = 1:(max(mod_data$period_idx)*3),
             Gender = rep(c(1,2,3), each = max(mod_data$period_idx)),
             Term = rep(1:max(mod_data$period_idx), 3),
             TermName = "Period",
             SharedName = "Age",
             LC_Index = NA)

for(row in 1:(nrow(lincombs_info_per_shareage))){
  if(lincombs_info_per_shareage$Term[row] %in%
     2:(max(lincombs_info_per_shareage$Term) -1)){
    row_idx <- row_idx + 1
    lincombs_info_per_shareage$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    period <- rep(NA, max(mod_data$period_idx)*3)
    per_element <- lincombs_info_per_shareage$Index[row]
    period[per_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(period_idx = period))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareage <- get(object.name)
      names(lincombs_shareage)[row_idx] <- object.name
    }else{
      lincombs_shareage <- c(lincombs_shareage, get(object.name))
      names(lincombs_shareage)[row_idx] <- object.name
    }
  }
}  

#### Age ####

lincombs_info_age_shareage <-
  data.frame(Index = 1:(max(mod_data$age_idx)*1),
             Gender = rep(NA, each = max(mod_data$age_idx)),
             Term = rep(1:max(mod_data$age_idx), 1),
             TermName = "Age",
             SharedName = "Age",
             LC_Index = NA)

for(row in 1:(nrow(lincombs_info_age_shareage))){
  if(lincombs_info_age_shareage$Term[row] %in%
     2:(max(lincombs_info_age_shareage$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_age_shareage$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    age <- rep(NA, max(mod_data$age_idx)*1)
    age_element <- lincombs_info_age_shareage$Index[row]
    age[age_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(age_idx = age))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareage <- get(object.name)
      names(lincombs_shareage)[row_idx] <- object.name
    }else{
      lincombs_shareage <- c(lincombs_shareage, get(object.name))
      names(lincombs_shareage)[row_idx] <- object.name
    }
  }
}  

lincombs_info_shareage <- bind_rows(lincombs_info_age_shareage,
                                    lincombs_info_per_shareage,
                                    lincombs_info_coh_shareage)

### Shared Period ####
#### Cohort ####
lincombs_info_coh_shareper<-
  data.frame(Index = 1:(max(mod_data$cohort_idx)*3),
             Gender = rep(c(1,2,3), each = max(mod_data$cohort_idx)),
             Term = rep(1:max(mod_data$cohort_idx), 3),
             TermName = "Cohort",
             SharedName = "Period",
             LC_Index = NA)

row_idx <- 0
for(row in 1:(nrow(lincombs_info_coh_shareper))){
  if(lincombs_info_coh_shareper$Term[row] %in%
     2:(max(lincombs_info_coh_shareper$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_coh_shareper$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    cohort <- rep(NA, max(mod_data$cohort_idx)*3)
    coh_element <- lincombs_info_coh_shareper$Index[row]
    cohort[coh_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(cohort_idx = cohort))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareper <- get(object.name)
      names(lincombs_shareper)[row_idx] <- object.name
    }else{
      lincombs_shareper <- c(lincombs_shareper, get(object.name))
      names(lincombs_shareper)[row_idx] <- object.name
    }
  }
}  

#### Age ####
lincombs_info_age_shareper <-
  data.frame(Index = 1:(max(mod_data$age_idx)*3),
             Gender = rep(c(1,2,3), each = max(mod_data$age_idx)),
             Term = rep(1:max(mod_data$age_idx), 3),
             TermName = "Age",
             SharedName = "Period",
             LC_Index = NA)

for(row in 1:(nrow(lincombs_info_age_shareper))){
  if(lincombs_info_age_shareper$Term[row] %in%
     2:(max(lincombs_info_age_shareper$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_age_shareper$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    age <- rep(NA, max(mod_data$age_idx)*3)
    age_element <- lincombs_info_age_shareper$Index[row]
    age[age_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(age_idx = age))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareper <- get(object.name)
      names(lincombs_shareper)[row_idx] <- object.name
    }else{
      lincombs_shareper <- c(lincombs_shareper, get(object.name))
      names(lincombs_shareper)[row_idx] <- object.name
    }
  }
}  

#### Period ####
lincombs_info_per_shareper <-
  data.frame(Index = 1:(max(mod_data$period_idx)*1),
             Gender = rep(NA, each = max(mod_data$period_idx)),
             Term = rep(1:max(mod_data$period_idx), 1),
             TermName = "Period",
             SharedName = "Period",
             LC_Index = NA)

for(row in 1:(nrow(lincombs_info_per_shareper))){
  if(lincombs_info_per_shareper$Term[row] %in%
     2:(max(lincombs_info_per_shareper$Term) - 1)){
    row_idx <- row_idx + 1
    lincombs_info_per_shareper$LC_Index[row] <- row_idx
    ## Put linear combination weights 
    ## in for N/As
    period <- rep(NA, max(mod_data$period_idx)*1)
    per_element <- lincombs_info_per_shareper$Index[row]
    period[per_element + c(-1, 0, 1)] <- c(1, -2, 1)
    
    ## Make lincomb in INLA and store
    object.name <- paste("lc", row_idx, sep = "")
    assign(object.name, inla.make.lincomb(period_idx = period))
    
    ## Add to lincombs list
    if(row_idx == 1){
      lincombs_shareper <- get(object.name)
      names(lincombs_shareper)[row_idx] <- object.name
    }else{
      lincombs_shareper <- c(lincombs_shareper, get(object.name))
      names(lincombs_shareper)[row_idx] <- object.name
    }
  }
}  


lincombs_info_shareper <- bind_rows(lincombs_info_age_shareper,
                                    lincombs_info_per_shareper,
                                    lincombs_info_coh_shareper)

## Fit Models ####

### Shared: Age ####
mod_shareage <- inla(mod_formula_shareage,
                     data = mod_data %>% 
                       filter(!is.na(prec_logit_pijk)), 
                     family = "gaussian",
                     control.family = list(hyper = list(
                       prec = list(initial = log(1),
                                   fixed=TRUE))), 
                     control.predictor = list(compute=TRUE),
                     control.compute = list(config = TRUE),
                     lincomb = lincombs_shareage,
                     scale = prec_logit_pijk,
                     verbose = TRUE)

age_info <- unique(mod_data[, c("age_idx", "age")])
per_info <- unique(mod_data[, c("period_idx", "period")])
coh_info <- unique(mod_data[, c("cohort_idx", "cohort")])

#### Plots ####
##### A ####
age_shareage <- mod_shareage$summary.lincomb.derived %>%
  left_join(lincombs_info_shareage,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                         Gender == 2 ~ "Transgender man",
                         Gender == 3 ~ "NB/GNC",
                         TRUE ~ "Shared")) %>%
  filter(TermName == "Age") %>% 
  left_join(age_info,
            by = c("Term" = "age_idx")) %>% 
  mutate(Age = age) %>% 
  ggplot(aes(x = Age, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  ggtitle("Age") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/age_curvatures_shareage.png", age_shareage)

age_shareage_std <- age_shareage + ylim(c(-2, 2))
ggsave(filename = "plots/GI_Q2/age_curvatures_shareage_std.png",
       age_shareage_std)
##### P ####
per_shareage <- mod_shareage$summary.lincomb.derived %>%
  left_join(lincombs_info_shareage,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Period") %>% 
  left_join(per_info,
            by = c("Term" = "period_idx")) %>% 
  mutate(Period = period) %>% 
  ggplot(aes(x = Period, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  ggtitle("Period") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/per_curvatures_shareage.png", per_shareage)

per_shareage_std <- per_shareage + ylim(c(-2, 2))
ggsave(filename = "plots/GI_Q2/per_curvatures_shareage_std.png", 
       per_shareage_std)
##### C ####
coh_shareage <- mod_shareage$summary.lincomb.derived %>%
  left_join(lincombs_info_shareage,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Cohort") %>% 
  left_join(coh_info,
            by = c("Term" = "cohort_idx")) %>% 
  mutate(Cohort = cohort) %>% 
  ggplot(aes(x = Cohort, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  scale_x_reverse() +
  ggtitle("Cohort") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/coh_curvatures_shareage.png", coh_shareage)

coh_shareage_std <- coh_shareage + ylim(c(-2, 2))
ggsave(filename = "plots/GI_Q2/coh_curvatures_shareage_std.png",
       coh_shareage_std)

### Shared: Cohort ####

mod_shareco <- inla(mod_formula_shareco,
                    data = mod_data %>% 
                      filter(!is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.predictor = list(compute=TRUE),
                    control.compute = list(config = TRUE),
                    lincomb = lincombs_shareco,
                    scale = prec_logit_pijk,
                    verbose = TRUE)

#### Plots ####
##### A ####
age_shareco <- mod_shareco$summary.lincomb.derived %>%
  left_join(lincombs_info_shareco,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Age") %>% 
  left_join(age_info,
            by = c("Term" = "age_idx")) %>% 
  mutate(Age = age) %>% 
  ggplot(aes(x = Age, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  ggtitle("Age") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/age_curvatures_shareco.png", age_shareco)

age_shareco_std <- age_shareco + ylim(c(-1,1))
ggsave(filename = "plots/GI_Q2/age_curvatures_shareco_std.png",
       age_shareco_std)
##### P ####
per_shareco <- mod_shareco$summary.lincomb.derived %>%
  left_join(lincombs_info_shareco,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Period") %>% 
  left_join(per_info,
            by = c("Term" = "period_idx")) %>% 
  mutate(Period = period) %>% 
  ggplot(aes(x = Period, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  ggtitle("Period") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/per_curvatures_shareco.png", per_shareco)

per_shareco_std <- per_shareco + ylim(c(-2,2))
ggsave(filename = "plots/GI_Q2/per_curvatures_shareco_std.png",
       per_shareco_std)

##### C ####
coh_shareco <- mod_shareco$summary.lincomb.derived %>%
  left_join(lincombs_info_shareco,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Cohort") %>% 
  left_join(coh_info,
            by = c("Term" = "cohort_idx")) %>% 
  mutate(Cohort = cohort) %>% 
  ggplot(aes(x = Cohort, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  scale_x_reverse() +
  ggtitle("Cohort") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/coh_curvatures_shareco.png", coh_shareco)

coh_shareco_std <- coh_shareco + ylim(c(-1,1))
ggsave(filename = "plots/GI_Q2/coh_curvatures_shareco_std.png",
       coh_shareco_std)

### Shared: Period ####
mod_shareper <- inla(mod_formula_shareper,
                     data = mod_data %>%
                       filter(!is.na(prec_logit_pijk)), 
                     family = "gaussian",
                     control.family = list(hyper = list(
                       prec = list(initial = log(1),
                                   fixed=TRUE))), 
                     control.predictor = list(compute=TRUE),
                     control.compute = list(config = TRUE),
                     lincomb = lincombs_shareper,
                     scale = prec_logit_pijk,
                     verbose = TRUE)

#### Plots ####
##### A ####
age_shareper<- mod_shareper$summary.lincomb.derived %>%
  left_join(lincombs_info_shareper,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Age") %>% 
  left_join(age_info,
            by = c("Term" = "age_idx")) %>% 
  mutate(Age = age) %>% 
  ggplot(aes(x = Age, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  ggtitle("Age") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/age_curvatures_shareper.png", age_shareper)

age_shareper_std <- age_shareper + ylim(c(-1,1))
ggsave(filename = "plots/GI_Q2/age_curvatures_shareper_std.png",
       age_shareper_std)

##### P ####
per_shareper <- mod_shareper$summary.lincomb.derived %>%
  left_join(lincombs_info_shareper,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Period") %>% 
  left_join(per_info,
            by = c("Term" = "period_idx")) %>% 
  mutate(Period = period) %>% 
  ggplot(aes(x = Period, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  ggtitle("Period") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/per_curvatures_shareper.png", per_shareper)

per_shareper_std <- per_shareper + ylim(c(-1,1))
ggsave(filename = "plots/GI_Q2/per_curvatures_shareper_std.png",
       per_shareper_std)

##### C ####
coh_shareper <- mod_shareper$summary.lincomb.derived %>%
  left_join(lincombs_info_shareper,
            by = c("ID" = "LC_Index")) %>% 
  drop_na(Term) %>%  
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Cohort") %>% 
  left_join(coh_info,
            by = c("Term" = "cohort_idx")) %>% 
  mutate(Cohort = cohort) %>% 
  ggplot(aes(x = Cohort, y = `0.5quant`, group = Gender,
             color = Gender, fill = Gender)) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax =  `0.975quant`),
              alpha = 0.5, color = NA) +
  scale_color_manual(values =
                       c("Shared" = "goldenrod",
                         "Transgender woman" = "navy",
                         "Transgender man" = "forestgreen",
                         "NB/GNC" = "firebrick")) +
  scale_fill_manual(values =
                      c("Shared" = "goldenrod",
                        "Transgender woman" = "navy",
                        "Transgender man" = "forestgreen",
                        "NB/GNC" = "firebrick")) +
  theme_classic() +
  facet_wrap(~Gender) +
  scale_x_reverse() +
  ggtitle("Cohort") +
  ylab("Second Differences") 

ggsave(filename = "plots/GI_Q2/coh_curvatures_shareper.png", coh_shareper)

coh_shareper_std <- coh_shareper + ylim(c(-1,1))
ggsave(filename = "plots/GI_Q2/coh_curvatures_shareper_std.png",
       coh_shareper_std)

### Save Model Outputs ####
sink(file = "tables/GI_Q2/mod_shareper.txt")
summary(mod_shareper)
head(mod_shareper$summary.random$age_idx)
head(mod_shareper$summary.linear.predictor)
sink(file = NULL)

sink(file = "tables/GI_Q2/mod_shareco.txt")
summary(mod_shareco)
head(mod_shareco$summary.random$age_idx)
head(mod_shareco$summary.linear.predictor)
sink(file = NULL)

sink(file = "tables/GI_Q2/mod_shareage.txt")
summary(mod_shareage)
head(mod_shareage$summary.random$age_idx)
head(mod_shareage$summary.linear.predictor)
sink(file = NULL)
#
# Manuscript plots ####

## Shared Age ####

### Prep Age ####

age_shareage <- mod_shareage$summary.lincomb.derived %>%
  left_join(lincombs_info_shareage,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Age") %>%
  left_join(age_info,
            by = c("Term" = "age_idx")) %>%
  mutate(Age = age) %>%
  arrange(Term)

### Prep Period ####

per_shareage <- mod_shareage$summary.lincomb.derived %>%
  left_join(lincombs_info_shareage,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Period") %>%
  left_join(per_info,
            by = c("Term" = "period_idx")) %>%
  mutate(Period = period) %>%
  arrange(Term)

### Prep Cohort ####

coh_shareage <- mod_shareage$summary.lincomb.derived %>%
  left_join(lincombs_info_shareage,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Cohort") %>%
  left_join(coh_info,
            by = c("Term" = "cohort_idx")) %>%
  mutate(Cohort = cohort) %>%
  arrange(Term)

## Get limits
age_shareage %>%
  dplyr::summarize(lower = min(`0.025quant`),
            upper = max(`0.975quant`))
per_shareage %>%
  dplyr::summarize(lower = min(`0.025quant`),
            upper = max(`0.975quant`))
coh_shareage %>%
  dplyr::summarize(lower = min(`0.025quant`),
            upper = max(`0.975quant`))

A_xlims <- c(14, 79)
P_xlims <- c(2013.75, 2021.25) #c(2014, 2021)
C_xlims <- c(1935, 2007)

shareage_ylims <- c(-1.2, 1.2)


### Plot ####
png("plots/GI_Q2/SharedAge_panel_curvature.png")
    # width = 300*3, height = 300*3, res=300)
{
  par(mfrow = c(3,3), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))

  ## Transwoman
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)

      mtext("Transgender woman", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)

      ## title
      text("a. Age", x= 15, y= 1.05, cex=.9, adj=0)
      
      ## x-axis
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)

      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)

      ## y = 0
      abline(h = 0)

      ## Ests
      lines(age_shareage$Age[age_shareage$Gender == "Shared"],
            age_shareage$`0.5quant`[age_shareage$Gender == "Shared"],
            lwd = 2, col = "goldenrod")
      ## Intervals
      polygon(x = c(age_shareage$Age[age_shareage$Gender == "Shared"],
                    rev(age_shareage$Age[age_shareage$Gender == "Shared"])),
              y = c(age_shareage$`0.025quant`[age_shareage$Gender == "Shared"],
                    rev(age_shareage$`0.975quant`[age_shareage$Gender ==
                                                    "Shared"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }

    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("b. Period", x= 2014, y= 1.05, cex=.9, adj=0)

      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)

      ## y = 0
      abline(h = 0)

      ## Ests
      lines(per_shareage$Period[per_shareage$Gender == "Transgender woman"],
            per_shareage$`0.5quant`[per_shareage$Gender == "Transgender woman"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareage$Period[per_shareage$Gender == "Transgender woman"],
                    rev(per_shareage$Period[per_shareage$Gender == "Transgender woman"])),
              y = c(per_shareage$`0.025quant`[per_shareage$Gender == "Transgender woman"],
                    rev(per_shareage$`0.975quant`[per_shareage$Gender ==
                                                    "Transgender woman"])),
              col = alpha("navy", 0.35), border = FALSE)
    }

    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("c. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)

      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)

      ## y = 0
      abline(h = 0)

      ## Ests
      lines(coh_shareage$Cohort[coh_shareage$Gender == "Transgender woman"],
            coh_shareage$`0.5quant`[coh_shareage$Gender == "Transgender woman"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareage$Cohort[coh_shareage$Gender == "Transgender woman"],
                    rev(coh_shareage$Cohort[coh_shareage$Gender == "Transgender woman"])),
              y = c(coh_shareage$`0.025quant`[coh_shareage$Gender == "Transgender woman"],
                    rev(coh_shareage$`0.975quant`[coh_shareage$Gender == 
                                                    "Transgender woman"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }

  }

  ## Transman
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)

      ## title
      text("d. Age", x= 15, y= 1.05, cex=.9, adj=0)
      mtext("Transgender man", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)

      ## x-axis 
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)

      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)

      ## y = 0
      abline(h = 0)

      ## Ests
      lines(age_shareage$Age[age_shareage$Gender == "Shared"],
            age_shareage$`0.5quant`[age_shareage$Gender == "Shared"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(age_shareage$Age[age_shareage$Gender == "Shared"],
                    rev(age_shareage$Age[age_shareage$Gender == "Shared"])),
              y = c(age_shareage$`0.025quant`[age_shareage$Gender == "Shared"],
                    rev(age_shareage$`0.975quant`[age_shareage$Gender ==
                                                    "Shared"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }

    {
      ## Period
       plot(NA, xlim = P_xlims, ylim = shareage_ylims,
            xlab = "", ylab = "",
            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
       ## title
       text("e. Period", x= 2014, y= 1.05, cex=.9, adj=0)

       ## Period axis
       mtext("Period", side = 1, line = 2, cex=0.6)
       axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
       mtext(side = 1, 2014:2021, at = 2014:2021,
             cex=0.5, line=0.5)

       ## y = 0
       abline(h = 0)

      ## Ests
      lines(per_shareage$Period[per_shareage$Gender == "Transgender man"],
            per_shareage$`0.5quant`[per_shareage$Gender == "Transgender man"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareage$Period[per_shareage$Gender == "Transgender man"],
                    rev(per_shareage$Period[per_shareage$Gender == "Transgender man"])),
              y = c(per_shareage$`0.025quant`[per_shareage$Gender == "Transgender man"],
                    rev(per_shareage$`0.975quant`[per_shareage$Gender == 
                                                    "Transgender man"])),
              col = alpha("navy", 0.35), border = FALSE)
    }

    {
      ## Cohort
       plot(NA, xlim = C_xlims, ylim = shareage_ylims,
            xlab = "", ylab = "",
            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
       ## title
       text("f. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)

       ## Cohort axis
       mtext("Cohort", side = 1, line = 2, cex=0.6)
       axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
       mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
             cex=0.5, line=0.5)

       ## y = 0
       abline(h = 0)

      ## Ests
      lines(coh_shareage$Cohort[coh_shareage$Gender == "Transgender man"],
            coh_shareage$`0.5quant`[coh_shareage$Gender == "Transgender man"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareage$Cohort[coh_shareage$Gender == "Transgender man"],
                    rev(coh_shareage$Cohort[coh_shareage$Gender == "Transgender man"])),
              y = c(coh_shareage$`0.025quant`[coh_shareage$Gender == "Transgender man"],
                    rev(coh_shareage$`0.975quant`[coh_shareage$Gender == 
                                                    "Transgender man"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }

  }
  
  ## Nonbinary
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("d. Age", x= 15, y= 1.05, cex=.9, adj=0)
      mtext("NB/GNC", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## x-axis 
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareage$Age[age_shareage$Gender == "Shared"],
            age_shareage$`0.5quant`[age_shareage$Gender == "Shared"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(age_shareage$Age[age_shareage$Gender == "Shared"],
                    rev(age_shareage$Age[age_shareage$Gender == "Shared"])),
              y = c(age_shareage$`0.025quant`[age_shareage$Gender == "Shared"],
                    rev(age_shareage$`0.975quant`[age_shareage$Gender ==
                                                    "Shared"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareage$Period[per_shareage$Gender == "NB/GNC"],
            per_shareage$`0.5quant`[per_shareage$Gender == "NB/GNC"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareage$Period[per_shareage$Gender == "NB/GNC"],
                    rev(per_shareage$Period[per_shareage$Gender == "NB/GNC"])),
              y = c(per_shareage$`0.025quant`[per_shareage$Gender == "NB/GNC"],
                    rev(per_shareage$`0.975quant`[per_shareage$Gender == 
                                                    "NB/GNC"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareage$Cohort[coh_shareage$Gender == "NB/GNC"],
            coh_shareage$`0.5quant`[coh_shareage$Gender == "NB/GNC"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareage$Cohort[coh_shareage$Gender == "NB/GNC"],
                    rev(coh_shareage$Cohort[coh_shareage$Gender == "NB/GNC"])),
              y = c(coh_shareage$`0.025quant`[coh_shareage$Gender == "NB/GNC"],
                    rev(coh_shareage$`0.975quant`[coh_shareage$Gender == 
                                                    "NB/GNC"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
}
dev.off()

## Shared Period ####


### Prep Age ####

age_shareper <- mod_shareper$summary.lincomb.derived %>%
  left_join(lincombs_info_shareper,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Age") %>%
  left_join(age_info,
            by = c("Term" = "age_idx")) %>%
  mutate(Age = age) %>%
  arrange(Term)

### Prep Period ####

per_shareper <- mod_shareper$summary.lincomb.derived %>%
  left_join(lincombs_info_shareper,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Period") %>%
  left_join(per_info,
            by = c("Term" = "period_idx")) %>%
  mutate(Period = period) %>%
  arrange(Term)

### Prep Cohort ####

coh_shareper <- mod_shareper$summary.lincomb.derived %>%
  left_join(lincombs_info_shareper,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Cohort") %>%
  left_join(coh_info,
            by = c("Term" = "cohort_idx")) %>%
  mutate(Cohort = cohort) %>%
  arrange(Term)

## Get limits
age_shareper %>%
  dplyr::summarize(lower = min(`0.025quant`),
                   upper = max(`0.975quant`))
per_shareper %>%
  dplyr::summarize(lower = min(`0.025quant`),
                   upper = max(`0.975quant`))
coh_shareper %>%
  dplyr::summarize(lower = min(`0.025quant`),
                   upper = max(`0.975quant`))

### Plot ####
png("plots/GI_Q2/SharedPer_panel_curvature.png")
     # width = 900, height = 900, res=300)
{
  par(mfrow = c(3,3), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))
  
  ## Transwoman
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      mtext("Transgender woman", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## title
      text("a. Age", x= 15, y= 1.05, cex=.9, adj=0)
      
      ## x-axis
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareper$Age[age_shareper$Gender == "Transgender woman"],
            age_shareper$`0.5quant`[age_shareper$Gender == "Transgender woman"],
            lwd = 2, col = "goldenrod")
      ## Intervals
      polygon(x = c(age_shareper$Age[age_shareper$Gender == "Transgender woman"],
                    rev(age_shareper$Age[age_shareper$Gender == "Transgender woman"])),
              y = c(age_shareper$`0.025quant`[age_shareper$Gender == 
                                                "Transgender woman"],
                    rev(age_shareper$`0.975quant`[age_shareper$Gender ==
                                                    "Transgender woman"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("b. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareper$Period[per_shareper$Gender == "Shared"],
            per_shareper$`0.5quant`[per_shareper$Gender == "Shared"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareper$Period[per_shareper$Gender == "Shared"],
                    rev(per_shareper$Period[per_shareper$Gender == "Shared"])),
              y = c(per_shareper$`0.025quant`[per_shareper$Gender == "Shared"],
                    rev(per_shareper$`0.975quant`[per_shareper$Gender ==
                                                    "Shared"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("c. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareper$Cohort[coh_shareper$Gender == "Transgender woman"],
            coh_shareper$`0.5quant`[coh_shareper$Gender == "Transgender woman"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareper$Cohort[coh_shareper$Gender == "Transgender woman"],
                    rev(coh_shareper$Cohort[coh_shareper$Gender == "Transgender woman"])),
              y = c(coh_shareper$`0.025quant`[coh_shareper$Gender == "Transgender woman"],
                    rev(coh_shareper$`0.975quant`[coh_shareper$Gender == 
                                                    "Transgender woman"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
  
  ## Transman
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("d. Age", x= 15, y= 1.05, cex=.9, adj=0)
      mtext("Transgender man", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## x-axis 
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareper$Age[age_shareper$Gender == "Transgender man"],
            age_shareper$`0.5quant`[age_shareper$Gender == "Transgender man"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(age_shareper$Age[age_shareper$Gender == "Transgender man"],
                    rev(age_shareper$Age[age_shareper$Gender == "Transgender man"])),
              y = c(age_shareper$`0.025quant`[age_shareper$Gender == "Transgender man"],
                    rev(age_shareper$`0.975quant`[age_shareper$Gender ==
                                                    "Transgender man"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareper$Period[per_shareper$Gender == "Shared"],
            per_shareper$`0.5quant`[per_shareper$Gender == "Shared"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareper$Period[per_shareper$Gender == "Shared"],
                    rev(per_shareper$Period[per_shareper$Gender == "Shared"])),
              y = c(per_shareper$`0.025quant`[per_shareper$Gender == "Shared"],
                    rev(per_shareper$`0.975quant`[per_shareper$Gender == 
                                                    "Shared"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareper$Cohort[coh_shareper$Gender == "Transgender man"],
            coh_shareper$`0.5quant`[coh_shareper$Gender == "Transgender man"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareper$Cohort[coh_shareper$Gender == "Transgender man"],
                    rev(coh_shareper$Cohort[coh_shareper$Gender == "Transgender man"])),
              y = c(coh_shareper$`0.025quant`[coh_shareper$Gender == "Transgender man"],
                    rev(coh_shareper$`0.975quant`[coh_shareper$Gender == 
                                                    "Transgender man"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
  
  ## Nonbinary
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("d. Age", x= 15, y= 1.05, cex=.9, adj=0)
      mtext("NB/GNC", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## x-axis 
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareper$Age[age_shareper$Gender == "NB/GNC"],
            age_shareper$`0.5quant`[age_shareper$Gender == "NB/GNC"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(age_shareper$Age[age_shareper$Gender == "NB/GNC"],
                    rev(age_shareper$Age[age_shareper$Gender == "NB/GNC"])),
              y = c(age_shareper$`0.025quant`[age_shareper$Gender == "NB/GNC"],
                    rev(age_shareper$`0.975quant`[age_shareper$Gender ==
                                                    "NB/GNC"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareper$Period[per_shareper$Gender == "Shared"],
            per_shareper$`0.5quant`[per_shareper$Gender == "Shared"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareper$Period[per_shareper$Gender == "Shared"],
                    rev(per_shareper$Period[per_shareper$Gender == "Shared"])),
              y = c(per_shareper$`0.025quant`[per_shareper$Gender == "Shared"],
                    rev(per_shareper$`0.975quant`[per_shareper$Gender == 
                                                    "Shared"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareper$Cohort[coh_shareper$Gender == "NB/GNC"],
            coh_shareper$`0.5quant`[coh_shareper$Gender == "NB/GNC"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareper$Cohort[coh_shareper$Gender == "NB/GNC"],
                    rev(coh_shareper$Cohort[coh_shareper$Gender == "NB/GNC"])),
              y = c(coh_shareper$`0.025quant`[coh_shareper$Gender == "NB/GNC"],
                    rev(coh_shareper$`0.975quant`[coh_shareper$Gender == 
                                                    "NB/GNC"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
}
dev.off()

## Shared Cohort####

### Prep Age ####

age_shareco <- mod_shareco$summary.lincomb.derived %>%
  left_join(lincombs_info_shareco,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Age") %>%
  left_join(age_info,
            by = c("Term" = "age_idx")) %>%
  mutate(Age = age) %>%
  arrange(Term)

### Prep Period ####

per_shareco <- mod_shareco$summary.lincomb.derived %>%
  left_join(lincombs_info_shareco,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Period") %>%
  left_join(per_info,
            by = c("Term" = "period_idx")) %>%
  mutate(Period = period) %>%
  arrange(Term)

### Prep Cohort ####

coh_shareco <- mod_shareco$summary.lincomb.derived %>%
  left_join(lincombs_info_shareco,
            by = c("ID" = "LC_Index")) %>%
  drop_na(Term) %>%
  mutate(Gender = case_when(Gender == 1 ~ "Transgender woman",
                            Gender == 2 ~ "Transgender man",
                            Gender == 3 ~ "NB/GNC",
                            TRUE ~ "Shared")) %>%
  filter(TermName == "Cohort") %>%
  left_join(coh_info,
            by = c("Term" = "cohort_idx")) %>%
  mutate(Cohort = cohort) %>%
  arrange(Term)

## Get limits
age_shareco %>%
  dplyr::summarize(lower = min(`0.025quant`),
                   upper = max(`0.975quant`))
per_shareco %>%
  dplyr::summarize(lower = min(`0.025quant`),
                   upper = max(`0.975quant`))
coh_shareco %>%
  dplyr::summarize(lower = min(`0.025quant`),
                   upper = max(`0.975quant`))

### Plot ####
png("plots/GI_Q2/SharedPer_panel_curvature.png")
# width = 900, height = 900, res=300)
{
  par(mfrow = c(3,3), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))
  
  ## Transwoman
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      mtext("Transgender woman", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## title
      text("a. Age", x= 15, y= 1.05, cex=.9, adj=0)
      
      ## x-axis
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareco$Age[age_shareco$Gender == "Transgender woman"],
            age_shareco$`0.5quant`[age_shareco$Gender == "Transgender woman"],
            lwd = 2, col = "goldenrod")
      ## Intervals
      polygon(x = c(age_shareco$Age[age_shareco$Gender == "Transgender woman"],
                    rev(age_shareco$Age[age_shareco$Gender == "Transgender woman"])),
              y = c(age_shareco$`0.025quant`[age_shareco$Gender == 
                                                "Transgender woman"],
                    rev(age_shareco$`0.975quant`[age_shareco$Gender ==
                                                    "Transgender woman"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("b. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareco$Period[per_shareco$Gender == "Transgender woman"],
            per_shareco$`0.5quant`[per_shareco$Gender == "Transgender woman"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareco$Period[per_shareco$Gender == "Transgender woman"],
                    rev(per_shareco$Period[per_shareco$Gender == "Transgender woman"])),
              y = c(per_shareco$`0.025quant`[per_shareco$Gender == "Transgender woman"],
                    rev(per_shareco$`0.975quant`[per_shareco$Gender ==
                                                    "Transgender woman"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("c. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareco$Cohort[coh_shareco$Gender == "Shared"],
            coh_shareco$`0.5quant`[coh_shareco$Gender == "Shared"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareco$Cohort[coh_shareco$Gender == "Shared"],
                    rev(coh_shareco$Cohort[coh_shareco$Gender == "Shared"])),
              y = c(coh_shareco$`0.025quant`[coh_shareco$Gender == "Shared"],
                    rev(coh_shareco$`0.975quant`[coh_shareco$Gender == 
                                                    "Shared"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
  
  ## Transman
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("d. Age", x= 15, y= 1.05, cex=.9, adj=0)
      mtext("Transgender man", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## x-axis 
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareco$Age[age_shareco$Gender == "Transgender man"],
            age_shareco$`0.5quant`[age_shareco$Gender == "Transgender man"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(age_shareco$Age[age_shareco$Gender == "Transgender man"],
                    rev(age_shareco$Age[age_shareco$Gender == "Transgender man"])),
              y = c(age_shareco$`0.025quant`[age_shareco$Gender == "Transgender man"],
                    rev(age_shareco$`0.975quant`[age_shareco$Gender ==
                                                    "Transgender man"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareco$Period[per_shareco$Gender == "Transgender man"],
            per_shareco$`0.5quant`[per_shareco$Gender == "Transgender man"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareco$Period[per_shareco$Gender == "Transgender man"],
                    rev(per_shareco$Period[per_shareco$Gender == "Transgender man"])),
              y = c(per_shareco$`0.025quant`[per_shareco$Gender == "Transgender man"],
                    rev(per_shareco$`0.975quant`[per_shareco$Gender == 
                                                    "Transgender man"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareco$Cohort[coh_shareco$Gender == "Shared"],
            coh_shareco$`0.5quant`[coh_shareco$Gender == "Shared"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareco$Cohort[coh_shareco$Gender == "Shared"],
                    rev(coh_shareco$Cohort[coh_shareco$Gender == "Shared"])),
              y = c(coh_shareco$`0.025quant`[coh_shareco$Gender == "Shared"],
                    rev(coh_shareco$`0.975quant`[coh_shareco$Gender == 
                                                    "Shared"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
  
  ## Nonbinary
  {
    {
      ## Age
      plot(NA, xlim = A_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("d. Age", x= 15, y= 1.05, cex=.9, adj=0)
      mtext("NB/GNC", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## x-axis 
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(15, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(15, 75, 5), at = seq(15, 75, 5),
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(shareage_ylims[1], shareage_ylims[2], .2),
           cex=0.5, labels = FALSE)
      mtext(side = 2, 
            round(seq(shareage_ylims[1], shareage_ylims[2], .2), 1),
            at = seq(shareage_ylims[1], shareage_ylims[2], .2),
            cex=0.5, line=0.5)
      mtext("Second Differences", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(age_shareco$Age[age_shareco$Gender == "NB/GNC"],
            age_shareco$`0.5quant`[age_shareco$Gender == "NB/GNC"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(age_shareco$Age[age_shareco$Gender == "NB/GNC"],
                    rev(age_shareco$Age[age_shareco$Gender == "NB/GNC"])),
              y = c(age_shareco$`0.025quant`[age_shareco$Gender == "NB/GNC"],
                    rev(age_shareco$`0.975quant`[age_shareco$Gender ==
                                                    "NB/GNC"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    {
      ## Period
      plot(NA, xlim = P_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Period", x= 2014, y= 1.05, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, 2014:2021, at = 2014:2021,
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(per_shareco$Period[per_shareco$Gender == "NB/GNC"],
            per_shareco$`0.5quant`[per_shareco$Gender == "NB/GNC"],
            lwd = 2, col = "navy")
      polygon(x = c(per_shareco$Period[per_shareco$Gender == "NB/GNC"],
                    rev(per_shareco$Period[per_shareco$Gender == "NB/GNC"])),
              y = c(per_shareco$`0.025quant`[per_shareco$Gender == "NB/GNC"],
                    rev(per_shareco$`0.975quant`[per_shareco$Gender == 
                                                    "NB/GNC"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    {
      ## Cohort
      plot(NA, xlim = C_xlims, ylim = shareage_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Cohort", x= 1938, y= 1.05, cex=.9, adj=0)
      
      ## Cohort axis
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10),
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(coh_shareco$Cohort[coh_shareco$Gender == "Shared"],
            coh_shareco$`0.5quant`[coh_shareco$Gender == "Shared"],
            lwd = 2, col = "firebrick")
      polygon(x = c(coh_shareco$Cohort[coh_shareco$Gender == "Shared"],
                    rev(coh_shareco$Cohort[coh_shareco$Gender == "Shared"])),
              y = c(coh_shareco$`0.025quant`[coh_shareco$Gender == "Shared"],
                    rev(coh_shareco$`0.975quant`[coh_shareco$Gender == 
                                                    "Shared"])),
              col = alpha("firebrick", 0.35), border = FALSE)
    }
    
  }
}
dev.off()


