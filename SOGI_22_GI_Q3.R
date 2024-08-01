## Carry out chi-squared test for answer to 
## Jessica Godwin

# Setup ####
rm(list = ls())

if(!dir.exists("tables/GI_Q3")){
  if(!dir.exists("tables")){
    dir.create("tables")
  }
  dir.create("tables/GI_Q3")
}

if(!dir.exists("plots/GI_Q3")){
  if(!dir.exists("plots")){
    dir.create("plots")
  }
  dir.create("plots/GI_Q3")
}

## Libraries ####
library(survey)
library(tidyverse)

## Functions ####
tableNA <- function(x, ...){
  table(x, useNA = "ifany", ...)  
}

# Load Data ####

combo <- readRDS(file = "data - clean/combined.rds")

str(combo)
summary(combo)
head(combo)

# Prep for total ests ####

## Filter to analysis ages ####

mod_data <- combo %>% 
  filter(age <= 30 & source == "BRFSS")


# Specify YRBS design ####


brfss_des <- svydesign(ids = ~1, strata = ~year + stratum, weights = ~weight,
                       data = mod_data, nest = TRUE)
brfss_des

# Get totals by age and year ####

## BRFSS ####

brfss_totals <- brfss_E <- list()
for(age_calc in 18:30){
  brfss_totals[[paste0("age_", age_calc)]] <- 
    svyby(~gender, by = ~year, design = subset(brfss_des, age == age_calc),
          svytotal)
  brfss_E[[paste0("age_", age_calc)]] <- 
    svyby(~gender, by = ~year,
          design = subset(brfss_des, age == age_calc & year == 2015),
          svymean)
}

# Get test statistics & p-values ####

## BRFSS ####
brfss_test <- list()
for(age_calc in 18:30){
    observed <- brfss_totals[[paste0("age_", age_calc)]][,-1] %>% 
      select(!contains("se.")) %>% t()
    exp_prop <- brfss_E[[paste0("age_", age_calc)]][,-1] %>% 
      select(!contains("se.")) %>% t()
    expected <- apply(observed, 2, function(x){sum(x)*exp_prop})
    stats <- (observed - expected)^2/expected
    brfss_test[[paste0("age_", age_calc)]] <-
      1 - pchisq(colSums(stats), df = 4)
}


## Save tables ####

### BRFSS ####
brfss_tab_list <- list()
for(age_calc in 18:30){
  brfss_tab_list[[paste0("age_", age_calc)]] <- 
    brfss_test[[paste0("age_", age_calc)]]  %>% 
    t() %>% as.data.frame() %>% 
    select(-`2015`) %>% 
    mutate(age = age_calc) %>% 
    relocate("age", .before = "2014") 
}

brfss_tab <- do.call(rbind.data.frame, brfss_tab_list)
write.csv(brfss_tab, row.names = FALSE,
          file = paste0("tables/GI_Q3/ChiSquared_18_30.csv"))


# Plot changes ####

## Get means by age & year ####
### BRFSS ####

brfss_means <- svyby(~gender, by = ~age + year,
                    design = brfss_des, svymean)
names(brfss_means)[3:8] <- paste0("mean.", names(brfss_means)[3:8])

brfss_prev <- brfss_means %>% 
  pivot_longer(cols = contains("gender"),
               names_to = c("Measure", ".value"),
               names_sep = "\\.") %>% 
  pivot_longer(cols = contains("gender"),
               names_to = "Gender",
               names_prefix = "gender(.)_",
               values_to = "Value") %>% 
  pivot_wider(id_cols = c("age", "year", "Gender"),
              names_from = "Measure",
              values_from = "Value") %>% 
  mutate(Gender = factor(Gender),
         lower = mean - qnorm(.975)*se,
         upper = mean + qnorm(.975)*se)

#### Transwoman ####

##### base R ####

year_cols <- c("navy", "forestgreen", "goldenrod", "firebrick")

plotdat <- brfss_prev %>% 
  filter(Gender == "transwoman" & year %in% seq(2015, 2021, 2)) %>% 
  mutate(age = case_when(year == 2015 ~ age - .2,
                         year == 2017 ~ age - .075,
                         year == 2019 ~ age + 0.075,
                         TRUE ~ age + .2))
y_lims <- c(min(plotdat$lower), max(plotdat$upper) + 0.01)
x_lims <- c(17, 31)

png("plots/GI_Q3/Prev_tw_18to30_baseR.png")
{
  plot(plotdat$age, plotdat$mean,
       xlab = "Age", ylab = "Prevalence",
       ylim = c(0, 0.04),
       type = "n", axes = FALSE,
       main = "")
  axis(1, at = 18:30, labels = as.character(18:30))
  axis(2, at = seq(0, 0.04, 0.005))
  years <- seq(2015, 2021, 2)
  for(yr in years){
    year_idx <- match(yr, years)
    tmpdat <- plotdat %>% 
      filter(year == yr)
    points(tmpdat$age, tmpdat$mean, pch = 19, col = year_cols[year_idx])
    segments(x0 = tmpdat$age, y0 = tmpdat$lower,
             x1 = tmpdat$age, y1 = tmpdat$upper,
             lwd = 2, col = year_cols[year_idx])
  }
  
  legend("topright", pch = 19, col = year_cols,
         legend = years, bty = "n",
         ncol = 1, cex = 1)
  title(main = "Transwoman", adj = 0.5, font.main = 1)
}
dev.off()

#### Transman ####

##### base R ####
plotdat <- brfss_prev %>% 
  filter(Gender == "transman" & year %in% seq(2015, 2021, 2)) %>% 
  mutate(age = case_when(year == 2015 ~ age - .2,
                         year == 2017 ~ age - .075,
                         year == 2019 ~ age + 0.075,
                         TRUE ~ age + .2))
y_lims <- c(min(plotdat$lower), max(plotdat$upper) + 0.01)
x_lims <- c(17, 31)

png("plots/GI_Q3/Prev_tm_18to30_baseR.png")
{
  plot(plotdat$age, plotdat$mean,
       xlab = "Age", ylab = "Prevalence",
       ylim = c(0, 0.04),
       type = "n", axes = FALSE,
       main = "")
  axis(1, at = 18:30, labels = as.character(18:30))
  axis(2, at = seq(0, 0.04, 0.005))
  
  years <- seq(2015, 2021, 2)
  for(yr in years){
    year_idx <- match(yr, years)
    tmpdat <- plotdat %>% 
      filter(year == yr)
    points(tmpdat$age, tmpdat$mean, pch = 19, col = year_cols[year_idx])
    segments(x0 = tmpdat$age, y0 = tmpdat$lower,
             x1 = tmpdat$age, y1 = tmpdat$upper,
             lwd = 2, col = year_cols[year_idx])
  }
  
  legend("topright", pch = 19, col = year_cols,
         legend = years, bty = "n",
         ncol = 1, cex = 1)
  title(main = "Transman", adj = 0.5, font.main = 1)
}
dev.off()

#### Nonbinary ####

##### base R ####

plotdat <- brfss_prev %>% 
  filter(Gender == "nbgnc" & year %in% seq(2015, 2021, 2)) %>% 
  mutate(age = case_when(year == 2015 ~ age - .2,
                         year == 2017 ~ age - .075,
                         year == 2019 ~ age + 0.075,
                         TRUE ~ age + .2))
y_lims <- c(min(plotdat$lower), max(plotdat$upper) + 0.01)
x_lims <- c(17, 31)

png("plots/GI_Q3/Prev_nbgnc_18to30_baseR.png")
{
  plot(plotdat$age, plotdat$mean,
       xlab = "Age", ylab = "Prevalence",
       ylim = c(0, 0.04),
       type = "n", axes = FALSE,
       main = "")
  axis(1, at = 18:30, labels = as.character(18:30))
  axis(2, at = seq(0, 0.04, 0.005))
  
  years <- seq(2015, 2021, 2)
  for(yr in years){
    year_idx <- match(yr, years)
    tmpdat <- plotdat %>% 
      filter(year == yr)
    points(tmpdat$age, tmpdat$mean, pch = 19, col = year_cols[year_idx])
    segments(x0 = tmpdat$age, y0 = tmpdat$lower,
             x1 = tmpdat$age, y1 = tmpdat$upper,
             lwd = 2, col = year_cols[year_idx])
  }
  
  legend("topright", pch = 19, col = year_cols,
         legend = years, bty = "n",
         ncol = 1, cex = 1)
  title(main = "NB/GNC", adj = 0.5, font.main = 1)
}
dev.off()
