## Fit null, A, P, C, AP, AC, PC models
## Prevalence by age x period (survey) x cohort
## Jessica Godwin

# Setup ####
rm(list = ls())

if(!dir.exists("tables/GI_Q1")){
  if(!dir.exists("tables/")){
    dir.create("tables/")
  }
  dir.create("tables/GI_Q1")
}

if(!dir.exists("plots/GI_Q1")){
  if(!dir.exists("plots/")){
    dir.create("plots/")
  }
  dir.create("plots/GI_Q1")
}

# ## Libraries ####
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
  rename("gender" = "so")

coh_summary <- mod_data %>% 
  group_by(cohort) %>% tally()
age_summary <- mod_data %>% 
  group_by(age) %>% tally()
per_summary <- mod_data %>% 
  group_by(period) %>% tally()

# Fit models ####

## Define formulas ####
pc_prec <- list(theta = list(prior = 'pc.prec', param=c(1, .5)))
mod_null <- logit_pijk ~ 1
mod_A <- logit_pijk ~  1 +  
  f(age_idx, model = "rw2", 
    scale.model = TRUE, constr = TRUE, extraconstr = NULL,
    hyper = pc_prec)
mod_P <- logit_pijk ~ 1 + 
  f(period_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec) 

mod_C <- logit_pijk ~ 1 + 
  f(cohort_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec)

mod_AP <- logit_pijk ~ 1 +
  f(age_idx, model = "rw2", 
    scale.model = TRUE, constr = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(period_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec) 

mod_AC <- logit_pijk ~ 1 + 
  f(age_idx, model = "rw2", 
    scale.model = TRUE, constr = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(cohort_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec)

mod_PC <- logit_pijk ~ 1 + 
  f(period_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec) +
  f(cohort_idx, model = "rw2", 
    constr = TRUE, scale.model = TRUE, extraconstr = NULL,
    hyper = pc_prec)



## Transwoman ####

### Models ####

mod_tw_null <- inla(mod_null,
                    data = mod_data %>% 
                      filter(gender == "tw") %>% 
                      filter(!is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE,
                                           dic = TRUE,
                                           waic = TRUE, 
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)
mod_tw_P <- inla(mod_P,
                 data = mod_data %>% 
                   filter(gender == "tw") %>% 
                   filter(!is.na(prec_logit_pijk)), 
                 family = "gaussian",
                 control.family = list(hyper = list(
                   prec = list(initial = log(1),
                               fixed=TRUE))), 
                 control.compute = list(cpo = TRUE,
                                        dic = TRUE,
                                        waic = TRUE, 
                                        config = TRUE),
                 scale = prec_logit_pijk,
                 verbose = TRUE)


mod_tw_C <- inla(mod_C,
                 data = mod_data %>% 
                   filter(gender == "tw") %>% 
                   filter(!is.na(prec_logit_pijk)), 
                 family = "gaussian",
                 control.family = list(hyper = list(
                   prec = list(initial = log(1),
                               fixed=TRUE))), 
                 control.compute = list(cpo = TRUE,
                                        dic = TRUE,
                                        waic = TRUE, 
                                        config = TRUE),
                 scale = prec_logit_pijk,
                 verbose = TRUE)

mod_tw_A <- inla(mod_A,
                 data = mod_data %>% 
                   filter(gender == "tw") %>% 
                   filter(!is.na(prec_logit_pijk)), 
                 family = "gaussian",
                 control.family = list(hyper = list(
                   prec = list(initial = log(1),
                               fixed=TRUE))), 
                 control.compute = list(cpo = TRUE,
                                        dic = TRUE,
                                        waic = TRUE, 
                                        config = TRUE),
                 scale = prec_logit_pijk,
                 verbose = TRUE)          

mod_tw_AP <- inla(mod_AP,
                  data = mod_data %>% 
                    filter(gender == "tw") %>% 
                    filter(!is.na(prec_logit_pijk)), 
                  family = "gaussian",
                  control.family = list(hyper = list(
                    prec = list(initial = log(1),
                                fixed=TRUE))), 
                  control.compute = list(cpo = TRUE,
                                         dic = TRUE,
                                         waic = TRUE, 
                                         config = TRUE),
                  scale = prec_logit_pijk,
                  verbose = TRUE)
mod_tw_AC <- inla(mod_AC,
                  data = mod_data %>% 
                    filter(gender == "tw") %>% 
                    filter(!is.na(prec_logit_pijk)), 
                  family = "gaussian",
                  control.family = list(hyper = list(
                    prec = list(initial = log(1),
                                fixed=TRUE))), 
                  control.compute = list(cpo = TRUE,
                                         dic = TRUE,
                                         waic = TRUE, 
                                         config = TRUE),
                  scale = prec_logit_pijk,
                  verbose = TRUE)
mod_tw_PC <- inla(mod_PC,
                  data = mod_data %>% 
                    filter(gender == "tw") %>% 
                    filter(!is.na(prec_logit_pijk)), 
                  family = "gaussian",
                  control.family = list(hyper = list(
                    prec = list(initial = log(1),
                                fixed=TRUE))), 
                  control.compute = list(cpo = TRUE,
                                         dic = TRUE,
                                         waic = TRUE, 
                                         config = TRUE),
                  scale = prec_logit_pijk,
                  verbose = TRUE)

#### Criteria Tables ####
tw_criteria <- data.frame(Model = c("Null", "A", "P", "C",
                                    "AP", "AC", "PC"),
                          loglik = c(mod_tw_null$mlik[2],
                                     mod_tw_A$mlik[2],
                                     mod_tw_P$mlik[2],
                                     mod_tw_C$mlik[2],
                                     mod_tw_AP$mlik[2],
                                     mod_tw_AC$mlik[2],
                                     mod_tw_PC$mlik[2]),
                          WAIC = c(mod_tw_null$waic$waic,
                                   mod_tw_A$waic$waic,
                                   mod_tw_P$waic$waic,
                                   mod_tw_C$waic$waic,
                                   mod_tw_AP$waic$waic,
                                   mod_tw_AC$waic$waic,
                                   mod_tw_PC$waic$waic),
                          DIC = c(mod_tw_null$dic$dic,
                                  mod_tw_A$dic$dic,
                                  mod_tw_P$dic$dic,
                                  mod_tw_C$dic$dic,
                                  mod_tw_AP$dic$dic,
                                  mod_tw_AC$dic$dic,
                                  mod_tw_PC$dic$dic))

write.csv(tw_criteria, file = paste0("tables/GI_Q1/",
                                     "tw_critera.csv"),
          row.names = FALSE)

#### Model Outputs ####
sink(file = "tables/GI_Q1/mod_tw_summaries.txt")
summary(mod_tw_null)
summary(mod_tw_A)
summary(mod_tw_P)
summary(mod_tw_C)
summary(mod_tw_AP)
summary(mod_tw_AC)
summary(mod_tw_PC)
sink(file = NULL)

### Plots ####

#### AC ####

##### Posterior samples ####
post_tw_AC <- inla.posterior.sample(n = 1000, mod_tw_AC)

age_idx <- mod_data %>% 
  select(age, age_idx) %>% 
  unique()
cohort_idx <- mod_data %>% 
  select(cohort, cohort_idx) %>% 
  unique() %>% 
  arrange(cohort_idx)

## For each draw, filter down to only parameters
##we need to do estimate of prevalence
## & organize it in a usable way
post_tw_AC_df <- lapply(post_tw_AC, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("age_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("age_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

## Turn list of draws into a data.frame
post_tw_AC_df <- do.call(rbind.data.frame, post_tw_AC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

## Combine linear terms in the log odds estimate by draw
## and calculate median & credible interval estimates 
## for log odds and prevalence for each sexual orientation
## along the age and cohort dimensions 
post_tw_AC_plot <- post_tw_AC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tw_AC_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "age", "Age", "Cohort")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####
tw_AC_age_lo <- post_tw_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/TW_AC_age_logodds.png", tw_AC_age_lo)


tw_AC_coh_lo <- post_tw_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_AC_coh_logodds.png", tw_AC_coh_lo)

##### Prevalence ####
tw_AC_age_p <- post_tw_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Prevalence") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_AC_age_prev.png", tw_AC_age_p)


tw_AC_coh_p <- post_tw_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_AC_coh_prev.png", tw_AC_coh_p)

#### PC ####
##### Posteriors ####
post_tw_PC <- inla.posterior.sample(n = 1000, mod_tw_PC)

period_idx <- mod_data %>% 
  select(period, period_idx) %>% 
  unique() %>% 
  arrange(period_idx)

post_tw_PC_df <- lapply(post_tw_PC, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("period_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("period_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tw_PC_df <- do.call(rbind.data.frame, post_tw_PC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tw_PC_plot <- post_tw_PC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tw_PC_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "period", "Period", "Cohort")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Period", period_idx$period[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####

tw_PC_per_lo <- post_tw_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-9, -2, 1), limits = c(-10, -2)) +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_PC_per_logodds.png", tw_PC_per_lo,
       width = 2.66, height = 2.66, units = "in")


tw_PC_coh_lo <- post_tw_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-9, -2, 1), limits = c(-10, -2)) +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_PC_coh_logodds.png", tw_PC_coh_lo,
       width = 2.66, height = 2.66, units = "in")

##### Prevalences ####

tw_PC_per_p <- post_tw_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Prevalence") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_PC_per_prev.png", tw_PC_per_p)


tw_PC_coh_p <- post_tw_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_PC_coh_prev.png", tw_PC_coh_p)

#### AP ####
##### Posteriors ####
post_tw_AP <- inla.posterior.sample(n = 1000, mod_tw_AP)

post_tw_AP_df <- lapply(post_tw_AP, function(draw){
  out <- draw$latent[grepl("period_idx", rownames(draw$latent)) |
                       grepl("age_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("period_idx", rownames(draw$latent)) |
      grepl("age_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tw_AP_df <- do.call(rbind.data.frame, post_tw_AP_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tw_AP_plot <- post_tw_AP_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tw_AP_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "age", "Age", "Period")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                period_idx$period[term_idx]))
##### Log Odds ####
tw_AP_age_lo <- post_tw_AP_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_AP_age_logodds.png", tw_AC_age_lo)


tw_AP_per_lo <- post_tw_AP_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_AP_per_logodds.png", tw_AP_per_lo)



#### C ####

##### Posteriors ####
post_tw_C <- inla.posterior.sample(n = 1000, mod_tw_C)

post_tw_C_df <- lapply(post_tw_C, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tw_C_df <- do.call(rbind.data.frame, post_tw_C_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tw_C_plot <- post_tw_C_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tw_C_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "cohort", "Cohort", term)) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                cohort_idx$cohort[term_idx]))


##### Log Odds ####
tw_C_coh_lo <- post_tw_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_C_coh_logodds.png", tw_C_coh_lo)

##### Prevalences ####
tw_C_coh_p <- post_tw_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Woman") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tw_C_coh_prev.png", tw_C_coh_p)

## Transman ####
### Models ####

mod_tm_null <- inla(mod_null,
                    data = mod_data %>% 
                      filter(gender == "tm") %>% 
                      filter(!is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE,
                                           dic = TRUE,
                                           waic = TRUE, 
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)
mod_tm_P <- inla(mod_P,
                 data = mod_data %>% 
                   filter(gender == "tm") %>% 
                   filter(!is.na(prec_logit_pijk)), 
                 family = "gaussian",
                 control.family = list(hyper = list(
                   prec = list(initial = log(1),
                               fixed=TRUE))), 
                 control.compute = list(cpo = TRUE,
                                        dic = TRUE,
                                        waic = TRUE, 
                                        config = TRUE),
                 scale = prec_logit_pijk,
                 verbose = TRUE)


mod_tm_C <- inla(mod_C,
                 data = mod_data %>% 
                   filter(gender == "tm") %>% 
                   filter(!is.na(prec_logit_pijk)), 
                 family = "gaussian",
                 control.family = list(hyper = list(
                   prec = list(initial = log(1),
                               fixed=TRUE))), 
                 control.compute = list(cpo = TRUE,
                                        dic = TRUE,
                                        waic = TRUE, 
                                        config = TRUE),
                 scale = prec_logit_pijk,
                 verbose = TRUE)


mod_tm_A <- inla(mod_A,
                 data = mod_data %>% 
                   filter(gender == "tm") %>% 
                   filter(!is.na(prec_logit_pijk)), 
                 family = "gaussian",
                 control.family = list(hyper = list(
                   prec = list(initial = log(1),
                               fixed=TRUE))), 
                 control.compute = list(cpo = TRUE,
                                        dic = TRUE,
                                        waic = TRUE, 
                                        config = TRUE),
                 scale = prec_logit_pijk,
                 verbose = TRUE)

mod_tm_AP <- inla(mod_AP,
                  data = mod_data %>% 
                    filter(gender == "tm") %>% 
                    filter(!is.na(prec_logit_pijk)), 
                  family = "gaussian",
                  control.family = list(hyper = list(
                    prec = list(initial = log(1),
                                fixed=TRUE))), 
                  control.compute = list(cpo = TRUE,
                                         dic = TRUE,
                                         waic = TRUE, 
                                         config = TRUE),
                  scale = prec_logit_pijk,
                  verbose = TRUE)

mod_tm_AC <- inla(mod_AC,
                  data = mod_data %>% 
                    filter(gender == "tm") %>% 
                    filter(!is.na(prec_logit_pijk)), 
                  family = "gaussian",
                  control.family = list(hyper = list(
                    prec = list(initial = log(1),
                                fixed=TRUE))), 
                  control.compute = list(cpo = TRUE,
                                         dic = TRUE,
                                         waic = TRUE, 
                                         config = TRUE),
                  scale = prec_logit_pijk,
                  verbose = TRUE)
mod_tm_PC <- inla(mod_PC,
                  data = mod_data %>% 
                    filter(gender == "tm") %>% 
                    filter(!is.na(prec_logit_pijk)), 
                  family = "gaussian",
                  control.family = list(hyper = list(
                    prec = list(initial = log(1),
                                fixed=TRUE))), 
                  control.compute = list(cpo = TRUE,
                                         dic = TRUE,
                                         waic = TRUE, 
                                         config = TRUE),
                  scale = prec_logit_pijk,
                  verbose = TRUE)

#### Criteria Tables ####

tm_criteria <- data.frame(Model = c("Null", "A", "P", "C",
                                      "AP", "AC", "PC"),
                            loglik = c(mod_tm_null$mlik[2],
                                       mod_tm_A$mlik[2],
                                       mod_tm_P$mlik[2],
                                       mod_tm_C$mlik[2],
                                       mod_tm_AP$mlik[2],
                                       mod_tm_AC$mlik[2],
                                       mod_tm_PC$mlik[2]),
                            WAIC = c(mod_tm_null$waic$waic,
                                     mod_tm_A$waic$waic,
                                     mod_tm_P$waic$waic,
                                     mod_tm_C$waic$waic,
                                     mod_tm_AP$waic$waic,
                                     mod_tm_AC$waic$waic,
                                     mod_tm_PC$waic$waic),
                            DIC = c(mod_tm_null$dic$dic,
                                    mod_tm_A$dic$dic,
                                    mod_tm_P$dic$dic,
                                    mod_tm_C$dic$dic,
                                    mod_tm_AP$dic$dic,
                                    mod_tm_AC$dic$dic,
                                    mod_tm_PC$dic$dic))

write.csv(tm_criteria, file = paste0("tables/GI_Q1/",
                                       "tm_critera.csv"),
          row.names = FALSE)

#### Model Outputs ####
sink(file = "tables/GI_Q1/mod_tm_summaries.txt")
summary(mod_tm_null)
summary(mod_tm_A)
summary(mod_tm_P)
summary(mod_tm_C)
summary(mod_tm_AP)
summary(mod_tm_AC)
summary(mod_tm_PC)
sink(file = NULL)

### Plots ####

#### AC ####
##### Posteriors ####
post_tm_AC <- inla.posterior.sample(n = 1000, mod_tm_AC)

post_tm_AC_df <- lapply(post_tm_AC, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("age_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("age_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tm_AC_df <- do.call(rbind.data.frame, post_tm_AC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tm_AC_plot <- post_tm_AC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tm_AC_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "age", "Age", "Cohort")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####

tm_AC_age_lo <- post_tm_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_AC_age_logodds.png", tm_AC_age_lo)


tm_AC_coh_lo <- post_tm_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_AC_coh_logodds.png", tm_AC_coh_lo)

##### Prevalences ####

tm_AC_age_p <- post_tm_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Prevalence") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_AC_age_prev.png", tm_AC_age_p)


tm_AC_coh_p <- post_tm_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_AC_coh_prev.png", tm_AC_coh_p)

#### PC ####

##### Posteriors ####
post_tm_PC <- inla.posterior.sample(n = 1000, mod_tm_PC)

post_tm_PC_df <- lapply(post_tm_PC, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("period_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("period_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tm_PC_df <- do.call(rbind.data.frame, post_tm_PC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tm_PC_plot <- post_tm_PC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tm_PC_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "period", "Period", "Cohort")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Period", period_idx$period[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####

tm_PC_per_lo <- post_tm_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-9, -2, 1), limits = c(-10, -2)) +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_PC_per_logodds.png", tm_PC_per_lo,
       width = 2.66, height = 2.66, units = "in")


tm_PC_coh_lo <- post_tm_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-9, -2, 1), limits = c(-10, -2)) +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_PC_coh_logodds.png", tm_PC_coh_lo,
       width = 2.66, height = 2.66, units = "in")

##### Prevalence ####

tm_PC_per_p <- post_tm_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Prevalence") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_PC_per_prev.png", tm_PC_per_p)


tm_PC_coh_p <- post_tm_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_PC_coh_prev.png", tm_PC_coh_p)


#### AP ####
##### Posteriors ####
post_tm_AP <- inla.posterior.sample(n = 1000, mod_tm_AP)

post_tm_AP_df <- lapply(post_tm_AP, function(draw){
  out <- draw$latent[grepl("period_idx", rownames(draw$latent)) |
                       grepl("age_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("period_idx", rownames(draw$latent)) |
      grepl("age_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tm_AP_df <- do.call(rbind.data.frame, post_tm_AP_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tm_AP_plot <- post_tm_AP_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tm_AP_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "age", "Age", "Period")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                period_idx$period[term_idx]))

##### Log Odds ####
tm_AP_age_lo <- post_tm_AP_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_AP_age_logodds.png", tm_AP_age_lo)


tm_AP_per_lo <- post_tm_AP_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_AP_per_logodds.png", tm_AP_per_lo)


#### C ####
##### Posteriors ####
post_tm_C <- inla.posterior.sample(n = 1000, mod_tm_C)

post_tm_C_df <- lapply(post_tm_C, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_tm_C_df <- do.call(rbind.data.frame, post_tm_C_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_tm_C_plot <- post_tm_C_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tm_C_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "cohort", "Cohort", term)) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####

tm_C_coh_lo <- post_tm_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_C_coh_logodds.png", tm_C_coh_lo)

##### Prevalences ####

tm_C_coh_p <- post_tm_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Transgender Man") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/tm_C_coh_prev.png", tm_C_coh_p)


## Nonbinary/Gender Non-conforming ####

### Models ####
mod_nbgnc_null <- inla(mod_null,
                       data = mod_data %>% 
                         filter(gender == "nbgnc") %>% 
                         filter(!is.na(prec_logit_pijk)), 
                       family = "gaussian",
                       control.family = list(hyper = list(
                         prec = list(initial = log(1),
                                     fixed=TRUE))), 
                       control.compute = list(cpo = TRUE,
                                              dic = TRUE,
                                              waic = TRUE, 
                                              config = TRUE),
                       scale = prec_logit_pijk,
                       verbose = TRUE)
mod_nbgnc_P <- inla(mod_P,
                    data = mod_data %>% 
                      filter(gender == "nbgnc") %>% 
                      filter(!is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE,
                                           dic = TRUE,
                                           waic = TRUE, 
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)


mod_nbgnc_C <- inla(mod_C,
                    data = mod_data %>% 
                      filter(gender == "nbgnc") %>% 
                      filter(!is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE,
                                           dic = TRUE,
                                           waic = TRUE, 
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)


mod_nbgnc_A <- inla(mod_A,
                    data = mod_data %>% 
                      filter(gender == "nbgnc") %>% 
                      filter(!is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE,
                                           dic = TRUE,
                                           waic = TRUE, 
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)

mod_nbgnc_AP <- inla(mod_AP,
                     data = mod_data %>% 
                       filter(gender == "nbgnc") %>% 
                       filter(!is.na(prec_logit_pijk)), 
                     family = "gaussian",
                     control.family = list(hyper = list(
                       prec = list(initial = log(1),
                                   fixed=TRUE))), 
                     control.compute = list(cpo = TRUE,
                                            dic = TRUE,
                                            waic = TRUE, 
                                            config = TRUE),
                     scale = prec_logit_pijk,
                     verbose = TRUE)
mod_nbgnc_AC <- inla(mod_AC,
                     data = mod_data %>% 
                       filter(gender == "nbgnc") %>% 
                       filter(!is.na(prec_logit_pijk)), 
                     family = "gaussian",
                     control.family = list(hyper = list(
                       prec = list(initial = log(1),
                                   fixed=TRUE))), 
                     control.compute = list(cpo = TRUE,
                                            dic = TRUE,
                                            waic = TRUE, 
                                            config = TRUE),
                     scale = prec_logit_pijk,
                     verbose = TRUE)

mod_nbgnc_PC <- inla(mod_PC,
                     data = mod_data %>% 
                       filter(gender == "nbgnc") %>% 
                       filter(!is.na(prec_logit_pijk)), 
                     family = "gaussian",
                     control.family = list(hyper = list(
                       prec = list(initial = log(1),
                                   fixed=TRUE))), 
                     control.compute = list(cpo = TRUE,
                                            dic = TRUE,
                                            waic = TRUE, 
                                            config = TRUE),
                     scale = prec_logit_pijk,
                     verbose = TRUE)


#### Criteria Tables ####
nbgnc_criteria <- data.frame(Model = c("Null", "A", "P", "C",
                                      "AP", "AC", "PC"),
                            loglik = c(mod_nbgnc_null$mlik[2],
                                       mod_nbgnc_A$mlik[2],
                                       mod_nbgnc_P$mlik[2],
                                       mod_nbgnc_C$mlik[2],
                                       mod_nbgnc_AP$mlik[2],
                                       mod_nbgnc_AC$mlik[2],
                                       mod_nbgnc_PC$mlik[2]),
                            WAIC = c(mod_nbgnc_null$waic$waic,
                                     mod_nbgnc_A$waic$waic,
                                     mod_nbgnc_P$waic$waic,
                                     mod_nbgnc_C$waic$waic,
                                     mod_nbgnc_AP$waic$waic,
                                     mod_nbgnc_AC$waic$waic,
                                     mod_nbgnc_PC$waic$waic),
                            DIC = c(mod_nbgnc_null$dic$dic,
                                    mod_nbgnc_A$dic$dic,
                                    mod_nbgnc_P$dic$dic,
                                    mod_nbgnc_C$dic$dic,
                                    mod_nbgnc_AP$dic$dic,
                                    mod_nbgnc_AC$dic$dic,
                                    mod_nbgnc_PC$dic$dic))

write.csv(nbgnc_criteria, file = paste0("tables/GI_Q1/",
                                       "nbgnc_critera.csv"),
          row.names = FALSE)

#### Model Outputs ####
sink(file = "tables/GI_Q1/mod_nbgnc_summaries.txt")
summary(mod_nbgnc_null)
summary(mod_nbgnc_A)
summary(mod_nbgnc_P)
summary(mod_nbgnc_C)
summary(mod_nbgnc_AP)
summary(mod_nbgnc_AC)
summary(mod_nbgnc_PC)
sink(file = NULL)

### Plots ####

#### AC ####
##### Posteriors ####
post_nbgnc_AC <- inla.posterior.sample(n = 1000, mod_nbgnc_AC)

post_nbgnc_AC_df <- lapply(post_nbgnc_AC, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("age_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("age_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_nbgnc_AC_df <- do.call(rbind.data.frame, post_nbgnc_AC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_nbgnc_AC_plot <- post_nbgnc_AC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_nbgnc_AC_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "age", "Age", "Cohort")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####
nbgnc_AC_age_lo <- post_nbgnc_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "firebrick") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo),
              fill = "firebrick",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-9, -2, 1), limits = c(-10, -2)) +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_AC_age_logodds.png", nbgnc_AC_age_lo)


nbgnc_AC_coh_lo <- post_nbgnc_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-9, -2, 1), limits = c(-10, -2)) +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_AC_coh_logodds.png", nbgnc_AC_coh_lo)

##### Prevalences ####
nbgnc_AC_age_p <- post_nbgnc_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Prevalence") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_AC_age_prev.png", nbgnc_AC_age_p)

nbgnc_AC_coh_p <- post_nbgnc_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_AC_coh_prev.png", nbgnc_AC_coh_p)


#### PC ####
##### Posteriors ####
post_nbgnc_PC <- inla.posterior.sample(n = 1000, mod_nbgnc_PC)

post_nbgnc_PC_df <- lapply(post_nbgnc_PC, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("period_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("period_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_nbgnc_PC_df <- do.call(rbind.data.frame, post_nbgnc_PC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_nbgnc_PC_plot <- post_nbgnc_PC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tw_PC_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "period", "Period", "Cohort")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Period", period_idx$period[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####

nbgnc_PC_per_lo <- post_nbgnc_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/\nGender Non-conforming") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-7, -4, 0.5), limits = c(-7.5, -3.75)) +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_PC_per_logodds.png", nbgnc_PC_per_lo,
       width = 2.66, height = 2.66, units = "in")

nbgnc_PC_coh_lo <- post_nbgnc_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/\nGender Non-conforming") +
  ylab("Log Odds") +
  scale_y_continuous(breaks = seq(-7, -4, 0.5), limits = c(-7.5, -3.75)) +
  xlab("Cohort") +
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_PC_coh_logodds.png", nbgnc_PC_coh_lo,
       width = 2.66, height = 2.66, units = "in")

##### Prevalences ####
nbgnc_PC_per_p <- post_nbgnc_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Prevalence") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_PC_per_prev.png", nbgnc_PC_per_p)

nbgnc_PC_coh_p <- post_nbgnc_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_PC_coh_prev.png", nbgnc_PC_coh_p)

#### AP ####
##### Posteriors ####
post_nbgnc_AP <- inla.posterior.sample(n = 1000, mod_nbgnc_AP)

post_nbgnc_AP_df <- lapply(post_nbgnc_AP, function(draw){
  out <- draw$latent[grepl("period_idx", rownames(draw$latent)) |
                       grepl("age_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("period_idx", rownames(draw$latent)) |
      grepl("age_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_nbgnc_AP_df <- do.call(rbind.data.frame, post_nbgnc_AP_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_nbgnc_AP_plot <- post_nbgnc_AP_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_nbgnc_AP_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "age", "Age", "Period")) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                period_idx$period[term_idx]))

##### Log Odds ####

nbgnc_AP_age_lo <- post_nbgnc_AP_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_AP_age_logodds.png", nbgnc_AC_age_lo)

nbgnc_AP_per_lo <- post_nbgnc_AP_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_AP_per_logodds.png", nbgnc_AP_per_lo)

#### C ####
##### Posteriors ####
post_nbgnc_C <- inla.posterior.sample(n = 1000, mod_nbgnc_C)

post_nbgnc_C_df <- lapply(post_nbgnc_C, function(draw){
  out <- draw$latent[grepl("cohort_idx", rownames(draw$latent)) |
                       grepl("(Intercept)", rownames(draw$latent))]
  out_names <- rownames(draw$latent)[
    grepl("cohort_idx", rownames(draw$latent)) |
      grepl("(Intercept)", rownames(draw$latent))]
  out_df <- data.frame(rowname = out_names,
                       draw = out) %>% 
    mutate(term = str_split_i(rowname, "_", 1),
           term_idx = str_split_i(rowname, ":", 2),
           term = ifelse(grepl("Intercept", term), "Intercept", term))
  return(out_df)
})

post_nbgnc_C_df <- do.call(rbind.data.frame, post_nbgnc_C_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_nbgnc_C_plot <- post_nbgnc_C_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_tw_C_df %>% 
              filter(term == "Intercept"),
            by = "draw_no",
            suffix = c("", "_int")) %>% 
  group_by(term, term_idx, draw_no) %>% 
  dplyr::summarize(draw_sum = draw + draw_int,
                   term_idx = as.numeric(term_idx),
                   term = ifelse(term == "cohort", "Cohort", term)) %>% 
  ungroup() %>% 
  group_by(term, term_idx) %>% 
  dplyr::summarize(med_lo = median(draw_sum),
                   lower_lo = quantile(draw_sum, 0.025),
                   upper_lo = quantile(draw_sum, 0.975),
                   med_p = median(expit(draw_sum)),
                   lower_p = quantile(expit(draw_sum), 0.025),
                   upper_p = quantile(expit(draw_sum), 0.975)) %>% 
  mutate(term_idx_name = ifelse(term == "Age", age_idx$age[term_idx],
                                cohort_idx$cohort[term_idx]))

##### Log Odds ####

nbgnc_C_coh_lo <- post_nbgnc_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_C_coh_logodds.png", nbgnc_C_coh_lo)

##### Prevalences ####

nbgnc_C_coh_p <- post_nbgnc_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Nonbinary/Gender Non-conforming") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q1/nbgnc_C_coh_prev.png", nbgnc_C_coh_p)

# Manuscript Plots ####

## Best model ####

## Get min and max endpoints of intervals by sexual orientation
## and age, period, or cohort across all models
best_ranges <- 
  bind_rows(post_tw_PC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Transwomen"),
            post_tm_PC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Transmen"),
            post_nbgnc_AC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Nonbinary/Gender Non-conforming"))

## Across sexual orientation groups 
best_ranges %>% 
  group_by(term) %>% 
  dplyr::summarize(lower = min(lower),
                   upper = max(upper))

## Handset appropriate limits
best_ylims <- c(-10, -2)
A_xlims <- c(18, 79) + c(-1, 1)
C_xlims <- c(1935, 2003) + c(-1,1)
P_xlims <- c(2014, 2021) + c(-1,1)

## AC ####

## Get min and max endpoints of intervals by sexual orientation
## and age, period, or cohort across all models
AC_ranges <- 
  bind_rows(post_tw_AC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Transwomen"),
            post_tm_AC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Transmen"),
            post_nbgnc_AC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Nonbinary/Gender Non-conforming"))
          
## Across sexual orientation groups 
AC_ranges %>% 
  group_by(term) %>% 
  dplyr::summarize(lower = min(lower),
                   upper = max(upper))

## Handset appropriate limits
AC_ylims <- c(-11, -2)
A_xlims <- c(18, 79) + c(-1, 1)
C_xlims <- c(1935, 2003) + c(-1,1)

### png ####
png("plots/GI_Q1/AC_8panel_logodds.png", 
    width = 300*6, height = 300*4, res=300)
{
  par(mfrow = c(2,3), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))
  
  ## Age 
  {
    { 
      ## TW
      plot(NA, xlim = A_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("a. Transwoman", x=20, y=-2.5, cex=.9, adj=0)
      # title("Age", adj = 0, line = 4)
      #mtext("Log Odds", adj = 0.275, side = 2,
      #      line = 2.5, cex = 0.8, outer = TRUE)
      #mtext("Log Odds", adj = 0.85, side = 2,
      #      line = 2.5, cex = 0.8, outer = TRUE)
      mtext("Age", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Age axis
      #axis(1, at = seq(15, 75, 5))
      # mtext(side = 1, text = "Age", line = 2, cex = .9)
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(20, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(20, 75, 5), at = seq(20, 75, 5), 
            cex=0.5, line=0.5)
      axis(2, at = seq(-10, -3, 1), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(-11, -2, 1), at = seq(-11, -2, 1), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tw_AC_plot$term_idx_name[post_tw_AC_plot$term == "Age"],
            post_tw_AC_plot$med_lo[post_tw_AC_plot$term == "Age"],
            lwd = 2, col = "navy")
      polygon(x = c(post_tw_AC_plot$term_idx_name[post_tw_AC_plot$term == "Age"],
                    rev(post_tw_AC_plot$term_idx_name[post_tw_AC_plot$term == 
                                                          "Age"])),
              y = c(post_tw_AC_plot$upper_lo[post_tw_AC_plot$term == "Age"],
                    rev(post_tw_AC_plot$lower_lo[post_tw_AC_plot$term == "Age"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    { 
      ## TM
      plot(NA, xlim = A_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("b. Transman", x=20, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      #axis(1, at = seq(15, 75, 5))
      # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(20, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(20, 75, 5), at = seq(20, 75, 5), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1), )
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tm_AC_plot$term_idx_name[post_tm_AC_plot$term == "Age"],
            post_tm_AC_plot$med_lo[post_tm_AC_plot$term == "Age"],
            lwd = 2, col = "navy")
      polygon(x = c(post_tm_AC_plot$term_idx_name[post_tm_AC_plot$term == "Age"],
                    rev(post_tm_AC_plot$term_idx_name[post_tm_AC_plot$term == 
                                                          "Age"])),
              y = c(post_tm_AC_plot$upper_lo[post_tm_AC_plot$term == "Age"],
                    rev(post_tm_AC_plot$lower_lo[post_tm_AC_plot$term == "Age"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    { 
      ## NB
      plot(NA, xlim = A_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("c. Nonbinary/Gender Non-conforming", x=20, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(20, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(20, 75, 5), at = seq(20, 75, 5), 
            cex=0.5, line=0.5)
      # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1))
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_nbgnc_AC_plot$term_idx_name[post_nbgnc_AC_plot$term == "Age"],
            post_nbgnc_AC_plot$med_lo[post_nbgnc_AC_plot$term == "Age"],
            lwd = 2, col = "navy")
      polygon(x = c(post_nbgnc_AC_plot$term_idx_name[post_nbgnc_AC_plot$term == "Age"],
                    rev(post_nbgnc_AC_plot$term_idx_name[post_nbgnc_AC_plot$term == 
                                                          "Age"])),
              y = c(post_nbgnc_AC_plot$upper_lo[post_nbgnc_AC_plot$term == "Age"],
                    rev(post_nbgnc_AC_plot$lower_lo[post_nbgnc_AC_plot$term == "Age"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
  }
  
  
  ## Cohort 
  {
    { 
      ## TW
      plot(NA, xlim = C_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Transwoman", x=1935, y=-2.5, cex=.9, adj=0)
      mtext("Cohort", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Age axis
      #axis(1, at = seq(1940, 2010, 10))
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
      
      ## Log odds axis
      axis(2, at = seq(-11, -2, 1), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(-11, -2, 1), at = seq(-11, -2, 1), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tw_AC_plot$term_idx_name[post_tw_AC_plot$term ==
                                              "Cohort"],
            post_tw_AC_plot$med_lo[post_tw_AC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_tw_AC_plot$term_idx_name[post_tw_AC_plot$term == 
                                                      "Cohort"],
                    rev(post_tw_AC_plot$term_idx_name[post_tw_AC_plot$term == 
                                                          "Cohort"])),
              y = c(post_tw_AC_plot$upper_lo[post_tw_AC_plot$term == 
                                                 "Cohort"],
                    rev(post_tw_AC_plot$lower_lo[post_tw_AC_plot$term == 
                                                     "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
    
    { 
      ## TM
      plot(NA, xlim = C_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Transman", x=1935, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      #axis(1, at = seq(1940, 2010, 10))
      # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1), )
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tm_AC_plot$term_idx_name[post_tm_AC_plot$term == "Cohort"],
            post_tm_AC_plot$med_lo[post_tm_AC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_tm_AC_plot$term_idx_name[post_tm_AC_plot$term == 
                                                      "Cohort"],
                    rev(post_tm_AC_plot$term_idx_name[post_tm_AC_plot$term == 
                                                          "Cohort"])),
              y = c(post_tm_AC_plot$upper_lo[post_tm_AC_plot$term == 
                                                 "Cohort"],
                    rev(post_tm_AC_plot$lower_lo[post_tm_AC_plot$term ==
                                                     "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
    
    { 
      ## GNC
      plot(NA, xlim = C_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("g. Nonbinary/Gender Non-conforming", x=1935, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      #axis(1, at = seq(1940, 2010, 10))
      # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1))
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_nbgnc_AC_plot$term_idx_name[post_nbgnc_AC_plot$term == "Cohort"],
            post_nbgnc_AC_plot$med_lo[post_nbgnc_AC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_nbgnc_AC_plot$term_idx_name[post_nbgnc_AC_plot$term == 
                                                      "Cohort"],
                    rev(post_nbgnc_AC_plot$term_idx_name[post_nbgnc_AC_plot$term == 
                                                          "Cohort"])),
              y = c(post_nbgnc_AC_plot$upper_lo[post_nbgnc_AC_plot$term == "Cohort"],
                    rev(post_nbgnc_AC_plot$lower_lo[post_nbgnc_AC_plot$term == 
                                                     "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
  }
}
dev.off()

## PC ####
## Get ranges by SO x A/P/C
PC_ranges <- bind_rows(post_tw_PC_plot %>% 
                         group_by(term) %>% 
                         dplyr::summarize(lower = min(lower_lo),
                                          upper = max(upper_lo),
                                          Group = "Bisexual females"),
                       post_tm_PC_plot %>% 
                         group_by(term) %>% 
                         dplyr::summarize(lower = min(lower_lo),
                                          upper = max(upper_lo),
                                          Group = "Bisexual males"),
                       post_nbgnc_PC_plot %>% 
                         group_by(term) %>% 
                         dplyr::summarize(lower = min(lower_lo),
                                          upper = max(upper_lo),
                                          Group = "Lesbians"))
## Across groups 
PC_ranges %>% 
  group_by(term) %>% 
  dplyr::summarize(lower = min(lower),
                   upper = max(upper))

## Hand set limits so they match AC, AP plots
PC_ylims <- c(-11, -2)
P_xlims <- c(2014, 2021) + c(-1, 1)
C_xlims <- c(1935, 2003) + c(-1,1)
P_labs <- c("\'14","\'15","\'16","\'17","\'18","\'19","\'20","\'21")

### png ####

png("plots/GI_Q1/PC_8panel_logodds.png", 
    width = 300*6, height = 300*3, res=300)
{
  par(mfrow = c(2,3), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))
  
  ## Period 
  {
    { 
      ## TW
      
      plot(NA, xlim = P_xlims, ylim = PC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("a. Transwoman", x=2013, y=-2.5, cex=.9, adj=0)
      mtext("Period", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(-11, -2, 1), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(-11, -2, 1), at = seq(-11, -2, 1), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tw_PC_plot$term_idx_name[post_tw_PC_plot$term == "Period"],
            post_tw_PC_plot$med_lo[post_tw_PC_plot$term == "Period"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(post_tw_PC_plot$term_idx_name[post_tw_PC_plot$term == "Period"],
                    rev(post_tw_PC_plot$term_idx_name[post_tw_PC_plot$term == 
                                                          "Period"])),
              y = c(post_tw_PC_plot$upper_lo[post_tw_PC_plot$term ==
                                                 "Period"],
                    rev(post_tw_PC_plot$lower_lo[post_tw_PC_plot$term ==
                                                     "Period"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    { 
      ## TM
      
      plot(NA, xlim = P_xlims, ylim = PC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("b. Transman", x=2013, y=-2.5, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      
      ## Ests
      lines(post_tm_PC_plot$term_idx_name[post_tm_PC_plot$term == "Period"],
            post_tm_PC_plot$med_lo[post_tm_PC_plot$term == "Period"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(post_tm_PC_plot$term_idx_name[post_tm_PC_plot$term == "Period"],
                    rev(post_tm_PC_plot$term_idx_name[post_tm_PC_plot$term == 
                                                          "Period"])),
              y = c(post_tm_PC_plot$upper_lo[post_tm_PC_plot$term == "Period"],
                    rev(post_tm_PC_plot$lower_lo[post_tm_PC_plot$term == "Period"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    { 
      ## Nonbinary/Gender Non-conforming
      plot(NA, xlim = P_xlims, ylim = PC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("c. Nonbinary/Gender Non-conforming", x=2013, y=-2.5, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_nbgnc_PC_plot$term_idx_name[post_nbgnc_PC_plot$term == 
                                              "Period"],
            post_nbgnc_PC_plot$med_lo[post_nbgnc_PC_plot$term == "Period"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(post_nbgnc_PC_plot$term_idx_name[post_nbgnc_PC_plot$term == 
                                                      "Period"],
                    rev(post_nbgnc_PC_plot$term_idx_name[post_nbgnc_PC_plot$term == 
                                                          "Period"])),
              y = c(post_nbgnc_PC_plot$upper_lo[post_nbgnc_PC_plot$term ==
                                                 "Period"],
                    rev(post_nbgnc_PC_plot$lower_lo[post_nbgnc_PC_plot$term == 
                                                     "Period"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
  }
  
  ## Cohort 
  {
    { 
      ## Transwoman
      plot(NA, xlim = C_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("e. Transwoman", x=1935, y=-2.5, cex=.9, adj=0)
      mtext("Cohort", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Cohort axis
      #axis(1, at = seq(1940, 2010, 10))
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
      
      ## Log odds axis
      axis(2, at = seq(-11, -2, 1), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(-11, -2, 1), at = seq(-11, -2, 1), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tw_PC_plot$term_idx_name[post_tw_PC_plot$term ==
                                              "Cohort"],
            post_tw_PC_plot$med_lo[post_tw_PC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_tw_PC_plot$term_idx_name[post_tw_PC_plot$term == 
                                                      "Cohort"],
                    rev(post_tw_PC_plot$term_idx_name[post_tw_PC_plot$term == 
                                                          "Cohort"])),
              y = c(post_tw_PC_plot$upper_lo[post_tw_PC_plot$term == 
                                                 "Cohort"],
                    rev(post_tw_PC_plot$lower_lo[post_tw_PC_plot$term == 
                                                     "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
    
    { 
      ## Transman
      plot(NA, xlim = C_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("f. Transman", x=1935, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      #axis(1, at = seq(1940, 2010, 10))
      # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1), )
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tm_PC_plot$term_idx_name[post_tm_PC_plot$term == "Cohort"],
            post_tm_PC_plot$med_lo[post_tm_PC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_tm_PC_plot$term_idx_name[post_tm_PC_plot$term == 
                                                      "Cohort"],
                    rev(post_tm_PC_plot$term_idx_name[post_tm_PC_plot$term == 
                                                          "Cohort"])),
              y = c(post_tm_PC_plot$upper_lo[post_tm_PC_plot$term == 
                                                 "Cohort"],
                    rev(post_tm_PC_plot$lower_lo[post_tm_PC_plot$term ==
                                                     "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
    
    { 
      ## Nonbinary/ Gender Non-conforming
      plot(NA, xlim = C_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("g. Nonbinary/Gender Non-conforming", x=1935, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      #axis(1, at = seq(1940, 2010, 10))
      # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1))
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_nbgnc_PC_plot$term_idx_name[post_nbgnc_PC_plot$term == "Cohort"],
            post_nbgnc_PC_plot$med_lo[post_nbgnc_PC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_nbgnc_PC_plot$term_idx_name[post_nbgnc_PC_plot$term == 
                                                      "Cohort"],
                    rev(post_nbgnc_PC_plot$term_idx_name[post_nbgnc_PC_plot$term == 
                                                          "Cohort"])),
              y = c(post_nbgnc_PC_plot$upper_lo[post_nbgnc_PC_plot$term == "Cohort"],
                    rev(post_nbgnc_PC_plot$lower_lo[post_nbgnc_PC_plot$term == 
                                                     "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
    
  }
}
dev.off()

## AP ####
## Calculate min and max CI ranges by SO x A/P/C
AP_ranges <- bind_rows(post_tw_AP_plot %>% 
                         group_by(term) %>% 
                         dplyr::summarize(lower = min(lower_lo),
                                          upper = max(upper_lo),
                                          Group = "Transwomen"),
                       post_tm_AP_plot %>% 
                         group_by(term) %>% 
                         dplyr::summarize(lower = min(lower_lo),
                                          upper = max(upper_lo),
                                          Group = "Transmen"),
                       post_nbgnc_AP_plot %>% 
                         group_by(term) %>% 
                         dplyr::summarize(lower = min(lower_lo),
                                          upper = max(upper_lo),
                                          Group = "Nonbinary/Gender Non-conforming"))
## Across groups 
AP_ranges %>% 
  group_by(term) %>% 
  dplyr::summarize(lower = min(lower),
                   upper = max(upper))

## Hand set limits so they match AC, AP plots
AP_ylims <- c(-11, -2)
A_xlims <- c(18, 79) + c(-1, 1)
P_xlims <- c(2014, 2021) + c(-1,1)

### png ####

png("plots/GI_Q1/AP_8panel_logodds.png", 
    width = 300*6, height = 300*4, res=300)
{
  par(mfrow = c(2,3), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))
  
  
  ## Age 
  {
    { 
      ## Transwoman
      plot(NA, xlim = A_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("a. Transwoman", x=20, y=-2.5, cex=.9, adj=0)
      # title("Age", adj = 0, line = 4)
      #mtext("Log Odds", adj = 0.275, side = 2,
      #      line = 2.5, cex = 0.8, outer = TRUE)
      #mtext("Log Odds", adj = 0.85, side = 2,
      #      line = 2.5, cex = 0.8, outer = TRUE)
      mtext("Age", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Age axis
      #axis(1, at = seq(15, 75, 5))
      # mtext(side = 1, text = "Age", line = 2, cex = .9)
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(20, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(20, 75, 5), at = seq(20, 75, 5), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(-11, -2, 1), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(-11, -2, 1), at = seq(-11, -2, 1), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## y = 0
      # abline(h = 0)
      
      ## Ests
      lines(post_tw_AP_plot$term_idx_name[post_tw_AP_plot$term == "Age"],
            post_tw_AP_plot$med_lo[post_tw_AP_plot$term == "Age"],
            lwd = 2, col = "navy")
      polygon(x = c(post_tw_AP_plot$term_idx_name[post_tw_AP_plot$term == "Age"],
                    rev(post_tw_AP_plot$term_idx_name[post_tw_AP_plot$term == 
                                                          "Age"])),
              y = c(post_tw_AP_plot$upper_lo[post_tw_AP_plot$term == "Age"],
                    rev(post_tw_AP_plot$lower_lo[post_tw_AP_plot$term == "Age"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    { 
      ## Transman
      plot(NA, xlim = A_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("b. Transmen", x=20, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      #axis(1, at = seq(15, 75, 5))
      # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(20, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(20, 75, 5), at = seq(20, 75, 5), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1), )
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tm_AP_plot$term_idx_name[post_tm_AP_plot$term == "Age"],
            post_tm_AP_plot$med_lo[post_tm_AP_plot$term == "Age"],
            lwd = 2, col = "navy")
      polygon(x = c(post_tm_AP_plot$term_idx_name[post_tm_AP_plot$term == "Age"],
                    rev(post_tm_AP_plot$term_idx_name[post_tm_AP_plot$term == 
                                                          "Age"])),
              y = c(post_tm_AP_plot$upper_lo[post_tm_AP_plot$term == "Age"],
                    rev(post_tm_AP_plot$lower_lo[post_tm_AP_plot$term == "Age"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    { 
      ##  Nonbinary
      plot(NA, xlim = A_xlims, ylim = AC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("c. Nonbinary/Gender Non-conforming", x=20, y=-2.5, cex=.9, adj=0)
      
      ## Age axis
      mtext("Age", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(20, 75, 5), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(20, 75, 5), at = seq(20, 75, 5), 
            cex=0.5, line=0.5)
      # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
      
      ## Log odds axis
      # axis(2, at = seq(-7, 1, 1))
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_nbgnc_AP_plot$term_idx_name[post_nbgnc_AP_plot$term == "Age"],
            post_nbgnc_AP_plot$med_lo[post_nbgnc_AP_plot$term == "Age"],
            lwd = 2, col = "navy")
      polygon(x = c(post_nbgnc_AP_plot$term_idx_name[post_nbgnc_AP_plot$term == "Age"],
                    rev(post_nbgnc_AP_plot$term_idx_name[post_nbgnc_AP_plot$term == 
                                                          "Age"])),
              y = c(post_nbgnc_AP_plot$upper_lo[post_nbgnc_AP_plot$term == "Age"],
                    rev(post_nbgnc_AP_plot$lower_lo[post_nbgnc_AP_plot$term == "Age"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
  }
  
  
  ## Period 
  {
    { 
      plot(NA, xlim = P_xlims, ylim = PC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("e. Transwoman", x=2013, y=-2.5, cex=.9, adj=0)
      mtext("Period", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(-11, -2, 1), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(-11, -2, 1), at = seq(-11, -2, 1), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tw_AP_plot$term_idx_name[post_tw_AP_plot$term == "Period"],
            post_tw_AP_plot$med_lo[post_tw_AP_plot$term == "Period"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(post_tw_AP_plot$term_idx_name[post_tw_AP_plot$term == "Period"],
                    rev(post_tw_AP_plot$term_idx_name[post_tw_AP_plot$term == 
                                                          "Period"])),
              y = c(post_tw_AP_plot$upper_lo[post_tw_AP_plot$term ==
                                                 "Period"],
                    rev(post_tw_AP_plot$lower_lo[post_tw_AP_plot$term ==
                                                     "Period"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    { 
      ## Transman
      plot(NA, xlim = P_xlims, ylim = PC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("f. Transman", x=2013, y=-2.5, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_tm_AP_plot$term_idx_name[post_tm_AP_plot$term == "Period"],
            post_tm_AP_plot$med_lo[post_tm_AP_plot$term == "Period"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(post_tm_AP_plot$term_idx_name[post_tm_AP_plot$term == "Period"],
                    rev(post_tm_AP_plot$term_idx_name[post_tm_AP_plot$term == 
                                                          "Period"])),
              y = c(post_tm_AP_plot$upper_lo[post_tm_AP_plot$term == "Period"],
                    rev(post_tm_AP_plot$lower_lo[post_tm_AP_plot$term == "Period"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
    
    { 
      ## NHonbinary
      plot(NA, xlim = P_xlims, ylim = PC_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      
      ## title
      text("g. Nonbinary/Gender Non-conforming", x=2013, y=-2.5, cex=.9, adj=0)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## y = 0
      abline(h = 0)
      
      ## Ests
      lines(post_nbgnc_AP_plot$term_idx_name[post_nbgnc_AP_plot$term == 
                                              "Period"],
            post_nbgnc_AP_plot$med_lo[post_nbgnc_AP_plot$term == "Period"],
            lwd = 2, col = "goldenrod")
      polygon(x = c(post_nbgnc_AP_plot$term_idx_name[post_nbgnc_AP_plot$term == 
                                                      "Period"],
                    rev(post_nbgnc_AP_plot$term_idx_name[post_nbgnc_AP_plot$term == 
                                                          "Period"])),
              y = c(post_nbgnc_AP_plot$upper_lo[post_nbgnc_AP_plot$term ==
                                                 "Period"],
                    rev(post_nbgnc_AP_plot$lower_lo[post_nbgnc_AP_plot$term == 
                                                     "Period"])),
              col = alpha("goldenrod", 0.35), border = FALSE)
    }
  }
}
dev.off()

