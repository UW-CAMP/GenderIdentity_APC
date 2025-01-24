## Fit null, A, P, C, AP, AC, PC models
## Prevalence by age x period (survey) x cohort
## Jessica Godwin

# Setup ####
rm(list = ls())

if(!dir.exists("tables/GI_Q5")){
  if(!dir.exists("tables/")){
    dir.create("tables/")
  }
  dir.create("tables/GI_Q5")
}

if(!dir.exists("plots/GI_Q5")){
  if(!dir.exists("plots/")){
    dir.create("plots/")
  }
  dir.create("plots/GI_Q5")
}

# ## Libraries ####
# install.packages("INLA",
#                  repos=c(getOption("repos"),
#                          INLA="https://inla.r-inla-download.org/R/stable"))
library(INLA) 
library(tidyverse)
library(ggpubr)


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
         age_fac = factor(age, levels = 14:79),
         cohort_fac = factor(cohort, levels = 1935:2007),
         sex_fac = factor(sex),
         period_idx = as.numeric(period_fac),
         age_idx = as.numeric(age_fac),
         cohort_idx = as.numeric(cohort_fac),
         sex_idx = as.numeric(sex_fac),
         prec_logit_pijk = 1/var_logit_pijk)

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



## Don't know/Not sure ####

mod_dkns_null <- inla(mod_null,
                      data = mod_data %>% 
                        filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                      family = "gaussian",
                      control.family = list(hyper = list(
                        prec = list(initial = log(1),
                                    fixed=TRUE))), 
                      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                      scale = prec_logit_pijk,
                      verbose = TRUE)
mod_dkns_P <- inla(mod_P,
                   data = mod_data %>%
                     filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                   family = "gaussian",
                   control.family = list(hyper = list(
                     prec = list(initial = log(1),
                                 fixed=TRUE))), 
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                   scale = prec_logit_pijk,
                   verbose = TRUE)


mod_dkns_C <- inla(mod_C,
                   data = mod_data %>%
                     filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                   family = "gaussian",
                   control.family = list(hyper = list(
                     prec = list(initial = log(1),
                                 fixed=TRUE))), 
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                          config = TRUE),
                   scale = prec_logit_pijk,
                   verbose = TRUE)

mod_dkns_A <- inla(mod_A,
                   data = mod_data %>%
                     filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                   family = "gaussian",
                   control.family = list(hyper = list(
                     prec = list(initial = log(1),
                                 fixed=TRUE))), 
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                          config = TRUE),
                   scale = prec_logit_pijk,
                   verbose = TRUE)

mod_dkns_AP <- inla(mod_AP,
                    data = mod_data %>%
                      filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)

mod_dkns_AC <- inla(mod_AC,
                    data = mod_data %>%
                      filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)

mod_dkns_PC <- inla(mod_PC,
                    data = mod_data %>%
                      filter(so == "dkns" & !is.na(prec_logit_pijk)), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)

#### Tables ####
dkns_criteria <- data.frame(Model = c("Null", "A", "P", "C",
                                      "AP", "AC", "PC"),
                            loglik = c(mod_dkns_null$mlik[2],
                                       mod_dkns_A$mlik[2],
                                       mod_dkns_P$mlik[2],
                                       mod_dkns_C$mlik[2],
                                       mod_dkns_AP$mlik[2],
                                       mod_dkns_AC$mlik[2],
                                       mod_dkns_PC$mlik[2]),
                            WAIC = c(mod_dkns_null$waic$waic,
                                     mod_dkns_A$waic$waic,
                                     mod_dkns_P$waic$waic,
                                     mod_dkns_C$waic$waic,
                                     mod_dkns_AP$waic$waic,
                                     mod_dkns_AC$waic$waic,
                                     mod_dkns_PC$waic$waic),
                            DIC = c(mod_dkns_null$dic$dic,
                                    mod_dkns_A$dic$dic,
                                    mod_dkns_P$dic$dic,
                                    mod_dkns_C$dic$dic,
                                    mod_dkns_AP$dic$dic,
                                    mod_dkns_AC$dic$dic,
                                    mod_dkns_PC$dic$dic))

write.csv(dkns_criteria, file = paste0("tables/GI_Q5/",
                                       "dkns_critera.csv"),
          row.names = FALSE)

#### Save Model Outputs ####
sink(file = "tables/GI_Q5/mod_dkns_summaries.txt")
summary(mod_dkns_null)
summary(mod_dkns_A)
summary(mod_dkns_P)
summary(mod_dkns_C)
summary(mod_dkns_AP)
summary(mod_dkns_AC)
summary(mod_dkns_PC)
sink(file = NULL)

### Plots ####

#### AC ####

post_dkns_AC <- inla.posterior.sample(n = 1000, mod_dkns_AC)

age_idx <- mod_data %>% 
  select(age, age_idx) %>% 
  unique()
cohort_idx <- mod_data %>% 
  select(cohort, cohort_idx) %>% 
  unique() %>% 
  arrange(cohort_idx)

post_dkns_AC_df <- lapply(post_dkns_AC, function(draw){
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

post_dkns_AC_df <- do.call(rbind.data.frame, post_dkns_AC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_dkns_AC_plot <- post_dkns_AC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_dkns_AC_df %>% 
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

dkns_AC_age_lo <- post_dkns_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_AC_age_logodds.png", dkns_AC_age_lo)


dkns_AC_coh_lo <- post_dkns_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_AC_coh_logodds.png", dkns_AC_coh_lo)


dkns_AC_age_p <- post_dkns_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Prevalence") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_AC_age_prev.png", dkns_AC_age_p)


dkns_AC_coh_p <- post_dkns_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_AC_coh_prev.png", dkns_AC_coh_p)

#### PC ####

post_dkns_PC <- inla.posterior.sample(n = 1000, mod_dkns_PC)

period_idx <- mod_data %>% 
  select(period, period_idx) %>% 
  unique() %>% 
  arrange(period_idx)

post_dkns_PC_df <- lapply(post_dkns_PC, function(draw){
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

post_dkns_PC_df <- do.call(rbind.data.frame, post_dkns_PC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_dkns_PC_plot <- post_dkns_PC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_dkns_PC_df %>% 
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


dkns_PC_per_lo <- post_dkns_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_PC_per_logodds.png", dkns_PC_per_lo)


dkns_PC_coh_lo <- post_dkns_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure ") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_PC_coh_logodds.png", dkns_PC_coh_lo)

dkns_PC_per_p <- post_dkns_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Prevalence") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_PC_per_prev.png", dkns_PC_per_p)


dkns_PC_coh_p <- post_dkns_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_PC_coh_prev.png", dkns_PC_coh_p)

#### AP ####

post_dkns_AP <- inla.posterior.sample(n = 1000, mod_dkns_AP)

post_dkns_AP_df <- lapply(post_dkns_AP, function(draw){
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

post_dkns_AP_df <- do.call(rbind.data.frame, post_dkns_AP_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_dkns_AP_plot <- post_dkns_AP_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_dkns_AP_df %>% 
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

dkns_AP_age_lo <- post_dkns_AP_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_AP_age_logodds.png", dkns_AC_age_lo)


dkns_AP_per_lo <- post_dkns_AP_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_AP_per_logodds.png", dkns_AP_per_lo)



#### C ####

post_dkns_C <- inla.posterior.sample(n = 1000, mod_dkns_C)

post_dkns_C_df <- lapply(post_dkns_C, function(draw){
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

post_dkns_C_df <- do.call(rbind.data.frame, post_dkns_C_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_dkns_C_plot <- post_dkns_C_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_dkns_C_df %>% 
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

dkns_C_coh_lo <- post_dkns_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_C_coh_logodds.png", dkns_C_coh_lo)

dkns_C_coh_p <- post_dkns_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Don't Know/Not Sure") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/dkns_C_coh_prev.png", dkns_C_coh_p)

## Decline ####

mod_ref_null <- inla(mod_null,
                      data = mod_data %>%
                        filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                        drop_na(), 
                      family = "gaussian",
                      control.family = list(hyper = list(
                        prec = list(initial = log(1),
                                    fixed=TRUE))), 
                      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                                             config = TRUE),
                      scale = prec_logit_pijk,
                      verbose = TRUE)
mod_ref_P <- inla(mod_P,
                   data = mod_data %>%
                     filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                     drop_na(), 
                   family = "gaussian",
                   control.family = list(hyper = list(
                     prec = list(initial = log(1),
                                 fixed=TRUE))), 
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                                          config = TRUE),
                   scale = prec_logit_pijk,
                   verbose = TRUE)

mod_ref_C <- inla(mod_C,
                   data = mod_data %>% 
                     filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                     drop_na(), 
                   family = "gaussian",
                   control.family = list(hyper = list(
                     prec = list(initial = log(1),
                                 fixed=TRUE))), 
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                                          config = TRUE),
                   scale = prec_logit_pijk,
                   verbose = TRUE)

mod_ref_A <- inla(mod_A,
                   data = mod_data %>%
                     filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                     drop_na(), 
                   family = "gaussian",
                   control.family = list(hyper = list(
                     prec = list(initial = log(1),
                                 fixed=TRUE))), 
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                                          config = TRUE),
                   scale = prec_logit_pijk,
                   verbose = TRUE)

mod_ref_AP <- inla(mod_AP,
                    data = mod_data %>%
                      filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                      drop_na(), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)

mod_ref_AC <- inla(mod_AC,
                    data = mod_data %>%
                      filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                      drop_na(), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)

mod_ref_PC <- inla(mod_PC,
                    data = mod_data %>%
                      filter(so == "ref" & !is.na(prec_logit_pijk)) %>% 
                      drop_na(), 
                    family = "gaussian",
                    control.family = list(hyper = list(
                      prec = list(initial = log(1),
                                  fixed=TRUE))), 
                    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                           config = TRUE),
                    scale = prec_logit_pijk,
                    verbose = TRUE)


#### Tables ####
ref_criteria <- data.frame(Model = c("Null", "A", "P", "C",
                                      "AP", "AC", "PC"),
                            loglik = c(mod_ref_null$mlik[2],
                                       mod_ref_A$mlik[2],
                                       mod_ref_P$mlik[2],
                                       mod_ref_C$mlik[2],
                                       mod_ref_AP$mlik[2],
                                       mod_ref_AC$mlik[2],
                                       mod_ref_PC$mlik[2]),
                            WAIC = c(mod_ref_null$waic$waic,
                                     mod_ref_A$waic$waic,
                                     mod_ref_P$waic$waic,
                                     mod_ref_C$waic$waic,
                                     mod_ref_AP$waic$waic,
                                     mod_ref_AC$waic$waic,
                                     mod_ref_PC$waic$waic),
                            DIC = c(mod_ref_null$dic$dic,
                                    mod_ref_A$dic$dic,
                                    mod_ref_P$dic$dic,
                                    mod_ref_C$dic$dic,
                                    mod_ref_AP$dic$dic,
                                    mod_ref_AC$dic$dic,
                                    mod_ref_PC$dic$dic))

write.csv(ref_criteria, file = paste0("tables/GI_Q5/",
                                       "ref_critera.csv"),
          row.names = FALSE)

#### Save Model Outputs ####
sink(file = "tables/GI_Q5/mod_ref_summaries.txt")
summary(mod_ref_null)
summary(mod_ref_A)
summary(mod_ref_P)
summary(mod_ref_C)
summary(mod_ref_AP)
summary(mod_ref_AC)
summary(mod_ref_PC)
sink(file = NULL)

### Plots ####

#### AC ####

post_ref_AC <- inla.posterior.sample(n = 1000, mod_ref_AC)

post_ref_AC_df <- lapply(post_ref_AC, function(draw){
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

post_ref_AC_df <- do.call(rbind.data.frame, post_ref_AC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_ref_AC_plot <- post_ref_AC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_ref_AC_df %>% 
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


ref_AC_age_lo <- post_ref_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_AC_age_logodds.png", ref_AC_age_lo)


ref_AC_coh_lo <- post_ref_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_AC_coh_logodds.png", ref_AC_coh_lo)


ref_AC_age_p <- post_ref_AC_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Prevalence") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_AC_age_prev.png", ref_AC_age_p)


ref_AC_coh_p <- post_ref_AC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_AC_coh_prev.png", ref_AC_coh_p)


#### PC ####

post_ref_PC <- inla.posterior.sample(n = 1000, mod_ref_PC)

post_ref_PC_df <- lapply(post_ref_PC, function(draw){
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

post_ref_PC_df <- do.call(rbind.data.frame, post_ref_PC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_ref_PC_plot <- post_ref_PC_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_dkns_PC_df %>% 
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


ref_PC_per_lo <- post_ref_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_PC_per_logodds.png", ref_PC_per_lo)


ref_PC_coh_lo <- post_ref_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_PC_coh_logodds.png", ref_PC_coh_lo)


ref_PC_per_p <- post_ref_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Prevalence") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_PC_per_prev.png", ref_PC_per_p)

ref_PC_coh_p <- post_ref_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_PC_coh_prev.png", ref_PC_coh_p)

#### AP ####

post_ref_AP <- inla.posterior.sample(n = 1000, mod_ref_AP)

post_ref_AP_df <- lapply(post_ref_AP, function(draw){
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

post_ref_AP_df <- do.call(rbind.data.frame, post_ref_AP_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_ref_AP_plot <- post_ref_AP_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_ref_AP_df %>% 
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

ref_AP_age_lo <- post_ref_AP_plot %>%
  filter(term == "Age") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "navy") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "navy",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Age")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_AP_age_logodds.png", ref_AC_age_lo)


ref_AP_per_lo <- post_ref_AP_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_AP_per_logodds.png", ref_AP_per_lo)



#### C ####

post_ref_C <- inla.posterior.sample(n = 1000, mod_ref_C)

post_ref_C_df <- lapply(post_ref_C, function(draw){
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

post_ref_C_df <- do.call(rbind.data.frame, post_ref_C_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_ref_C_plot <- post_ref_C_df %>% 
  filter(term != "Intercept") %>% 
  left_join(post_dkns_C_df %>% 
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

ref_C_coh_lo <- post_ref_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Lesbian") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_C_coh_logodds.png", ref_C_coh_lo)

ref_C_coh_p <- post_ref_C_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_p)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_p, ymax = upper_p), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Lesbian") +
  ylab("Prevalence") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/GI_Q5/ref_C_coh_prev.png", ref_C_coh_p)

# Manuscript Plots ####

## Best model ####

## Get min and max endpoints of intervals by sexual orientation
## and age, period, or cohort across all models
best_ranges <- 
  bind_rows(post_dkns_PC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Don't Know/Not Sure"),
            post_ref_PC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Declined to Answer"))
## Across gender groups 
best_ranges %>% 
  group_by(term) %>% 
  dplyr::summarize(lower = min(lower),
                   upper = max(upper))

## Handset appropriate limits
best_ylims <- c(-10, -2)
A_xlims <- c(18, 79) + c(-1, 1)
C_xlims <- c(1935, 2003) + c(-1,1)
P_xlims <- c(2014, 2021) + c(-1,1)
P_labs <- c("\'14","\'15","\'16","\'17","\'18","\'19","\'20","\'21")

### png ####
png("plots/GI_Q5/BestModel_4panel_logodds.png", 
    width = 300*5, height = 300*5, res=300)
{
  par(mfrow = c(2,2), lend = 1)
  par(mar=c(5,0,2,0))
  par(oma=c(0,4,2,1))
  
  
  ## Cohort 
  {
    { 
      ## DKNS
      plot(NA, xlim = C_xlims, ylim = best_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("A) PC Model", x=1935, y=-2.5, cex=.9, adj=0)
      mtext("Don't Know/Not Sure", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Cohort axis
      #axis(1, at = seq(1940, 2010, 10))
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(best_ylims[1], best_ylims[2], 2), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(best_ylims[1], best_ylims[2], 2),
            at = seq(best_ylims[1], best_ylims[2], 2), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## Ests
      lines(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term ==
                                            "Cohort"],
            post_dkns_PC_plot$med_lo[post_dkns_PC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
                                                    "Cohort"],
                    rev(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
                                                        "Cohort"])),
              y = c(post_dkns_PC_plot$upper_lo[post_dkns_PC_plot$term == 
                                               "Cohort"],
                    rev(post_dkns_PC_plot$lower_lo[post_dkns_PC_plot$term == 
                                                   "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }
    
    { 
      ## DTA
      plot(NA, xlim = C_xlims, ylim = best_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("B) PC Model", x=1935, y=-2.5, cex=.9, adj=0)
      mtext("Declined to Answer", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Age axis
      #axis(1, at = seq(1940, 2010, 10))
      # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
      mtext("Cohort", side = 1, line = 2, cex=0.6)
      axis(1, at = seq(1940, 2000, 10), cex.axis=0.6, labels=FALSE)
      mtext(side = 1, seq(1940, 2000, 10), at = seq(1940, 2000, 10), 
            cex=0.5, line=0.5)
      
      ## Ests
      lines(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == "Cohort"],
            post_ref_PC_plot$med_lo[post_ref_PC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
                                                    "Cohort"],
                    rev(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
                                                        "Cohort"])),
              y = c(post_ref_PC_plot$upper_lo[post_ref_PC_plot$term == 
                                               "Cohort"],
                    rev(post_ref_PC_plot$lower_lo[post_ref_PC_plot$term ==
                                                   "Cohort"])),
              col = alpha("forestgreen", 0.35), border = FALSE)
    }

  }
  
  ## Period 
  {
    { 
      ## DKNS
      plot(NA, xlim = P_xlims, ylim = best_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      text("C) PC Model", x=2013, y=-2.5, cex=.9, adj=0)
      mtext("Don't Know/Not Sure", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## Log odds axis
      axis(2, at = seq(best_ylims[1], best_ylims[2], 2), cex=0.5, labels = FALSE)
      mtext(side = 2, seq(best_ylims[1], best_ylims[2], 2),
            at = seq(best_ylims[1], best_ylims[2], 2), 
            cex=0.5, line=0.5)
      mtext("Log Odds", side = 2, line = 2, cex = 0.6)
      
      ## Ests
      lines(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term ==
                                            "Period"],
            post_dkns_PC_plot$med_lo[post_dkns_PC_plot$term == "Period"],
            lwd = 2, col = "navy")
      polygon(x = c(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
                                                    "Period"],
                    rev(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
                                                        "Period"])),
              y = c(post_dkns_PC_plot$upper_lo[post_dkns_PC_plot$term == 
                                               "Period"],
                    rev(post_dkns_PC_plot$lower_lo[post_dkns_PC_plot$term == 
                                                   "Period"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
    
    { 
      ## DTA
      plot(NA, xlim = P_xlims, ylim = best_ylims,
           xlab = "", ylab = "",
           main = "", type = "n", frame.plot = TRUE, axes = FALSE)
      ## title
      text("D) PC Model", x=2013, y=-2.5, cex=.9, adj=0)
      mtext("Declined to Answer", adj = 0, side = 3, line = 0.5, outer = FALSE, cex=0.8)
      
      ## Period axis
      mtext("Period", side = 1, line = 2, cex=0.6)
      axis(1, at = 2014:2021, cex.axis=0.6, labels=FALSE)
      mtext(side = 1, P_labs, at = 2014:2021, 
            cex=0.5, line=0.5)
      
      ## Ests
      lines(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == "Period"],
            post_ref_PC_plot$med_lo[post_ref_PC_plot$term == "Period"],
            lwd = 2, col = "navy")
      polygon(x = c(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
                                                    "Period"],
                    rev(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
                                                        "Period"])),
              y = c(post_ref_PC_plot$upper_lo[post_ref_PC_plot$term == 
                                               "Period"],
                    rev(post_ref_PC_plot$lower_lo[post_ref_PC_plot$term ==
                                                   "Period"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
  }
}
dev.off()




# ## AC ####
# AC_ranges <- bind_rows(post_dkns_AC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Bisexual females"),
#                        post_dkns_AC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Bisexual males"),
#                        post_ref_AC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Lesbians"),
#                        post_ref_AC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Gay males"))
# ## Across groups 
# AC_ranges %>% 
#   group_by(term) %>% 
#   dplyr::summarize(lower = min(lower),
#             upper = max(upper))
# 
# AC_ylims <- c(-8.5, 2.5)
# A_xlims <- c(14, 79) + c(-1, 1)
# C_xlims <- c(1935, 2007) + c(-1,1)
# 
# 
# # par(mfrow = c(2,4),lend = 1)
# 
# png("plots/GI_Q5/AC_8panel_logodds.png", 
#     width = 300*8, height = 300*6, res=300)
# {
#   par(mfrow = c(2,4), lend = 1)
#   par(mar=c(5,0,0,0))
#   par(oma=c(5,4,4,1))
#   
#   ## Age 
#   {
#     { 
#       ## Other Response F
#       plot(NA, xlim = A_xlims, ylim = AC_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("a. Other response - female", x= 35, y=2, cex=.9)
#       # title("Age", adj = 0, line = 4)
#       
#       mtext("Log Odds", adj = 0.275, side = 2,
#             line = 2.5, cex = 0.8, outer = TRUE)
#       
#       mtext("Log Odds", adj = 0.85, side = 2,
#             line = 2.5, cex = 0.8, outer = TRUE)
#       
#       mtext("Age", adj = 0, side = 3, line = 0.5, outer = TRUE)
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = .9)
#       
#       ## Log odds axis
#       axis(2, at = seq(-8, 2, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == "Age"],
#             post_dkns_AC_plot$med_lo[post_dkns_AC_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == "Age"],
#                     rev(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == 
#                                                           "Age"])),
#               y = c(post_dkns_AC_plot$upper_lo[post_dkns_AC_plot$term == "Age"],
#                     rev(post_dkns_AC_plot$lower_lo[post_dkns_AC_plot$term == "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Other response M
#       plot(NA, xlim = A_xlims, ylim = AC_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("b. Other response - male", x=33.5, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1), )
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == "Age"],
#             post_dkns_AC_plot$med_lo[post_dkns_AC_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == "Age"],
#                     rev(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == 
#                                                           "Age"])),
#               y = c(post_dkns_AC_plot$upper_lo[post_dkns_AC_plot$term == "Age"],
#                     rev(post_dkns_AC_plot$lower_lo[post_dkns_AC_plot$term == "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer F
#       plot(NA, xlim = A_xlims, ylim = AC_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("c. Declined to answer - female", x = 37.5, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == "Age"],
#             post_ref_AC_plot$med_lo[post_ref_AC_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == "Age"],
#                     rev(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == 
#                                                           "Age"])),
#               y = c(post_ref_AC_plot$upper_lo[post_ref_AC_plot$term == "Age"],
#                     rev(post_ref_AC_plot$lower_lo[post_ref_AC_plot$term == "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer M
#       plot(NA, xlim = A_xlims, ylim = AC_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("d. Declined to answer - male", x=36, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == "Age"],
#             post_ref_AC_plot$med_lo[post_ref_AC_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term ==
#                                                       "Age"],
#                     rev(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == 
#                                                           "Age"])),
#               y = c(post_ref_AC_plot$upper_lo[post_ref_AC_plot$term == "Age"],
#                     rev(post_ref_AC_plot$lower_lo[post_ref_AC_plot$term == 
#                                                      "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#   }
#   
#   
#   ## Cohort 
#   {
#     { 
#       ## Other Response F
#       plot(NA, xlim = C_xlims, ylim = AC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("e. Other response - female", x= 18 + 1940, y=2, cex=.9)
#       mtext("Cohort", adj = 0, side = 3, line = 0.5)
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
#       
#       ## Log odds axis
#       axis(2, at = seq(-8, 2, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term ==
#                                               "Cohort"],
#             post_dkns_AC_plot$med_lo[post_dkns_AC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_dkns_AC_plot$upper_lo[post_dkns_AC_plot$term == 
#                                                  "Cohort"],
#                     rev(post_dkns_AC_plot$lower_lo[post_dkns_AC_plot$term == 
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Other Response M
#       plot(NA, xlim = C_xlims, ylim = AC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("f. Other response - male", x=16.5 + 1940, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1), )
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == "Cohort"],
#             post_dkns_AC_plot$med_lo[post_dkns_AC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_dkns_AC_plot$term_idx_name[post_dkns_AC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_dkns_AC_plot$upper_lo[post_dkns_AC_plot$term == 
#                                                  "Cohort"],
#                     rev(post_dkns_AC_plot$lower_lo[post_dkns_AC_plot$term ==
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer F
#       plot(NA, xlim = C_xlims, ylim = AC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("g. Declined to answer - female", x = 21 + 1940, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == "Cohort"],
#             post_ref_AC_plot$med_lo[post_ref_AC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_ref_AC_plot$upper_lo[post_ref_AC_plot$term == "Cohort"],
#                     rev(post_ref_AC_plot$lower_lo[post_ref_AC_plot$term == 
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer M
#       plot(NA, xlim = C_xlims, ylim = AC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("h. Declined to answer - male", x = 19.5 + 1940, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == "Cohort"],
#             post_ref_AC_plot$med_lo[post_ref_AC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_ref_AC_plot$term_idx_name[post_ref_AC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_ref_AC_plot$upper_lo[post_ref_AC_plot$term == "Cohort"],
#                     rev(post_ref_AC_plot$lower_lo[post_ref_AC_plot$term == 
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#   }
# }
# dev.off()
# 
# ## PC ####
# PC_ranges <- bind_rows(post_dkns_PC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Bisexual females"),
#                        post_dkns_PC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Bisexual males"),
#                        post_ref_PC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Lesbians"),
#                        post_ref_PC_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Gay males"))
# ## Across groups 
# PC_ranges %>% 
#   group_by(term) %>% 
#   dplyr::summarize(lower = min(lower),
#             upper = max(upper))
# 
# PC_ylims <- c(-8.5, 2.5)
# P_xlims <- c(2014, 2021) + c(-1, 1)
# C_xlims <- c(1935, 2007) + c(-1,1)
# 
# 
# # par(mfrow = c(2,4),lend = 1)
# 
# png("plots/GI_Q5/PC_8panel_logodds.png", 
#     width = 300*8, height = 300*6, res=300)
# {
#   par(mfrow = c(2,4), lend = 1)
#   par(mar=c(5,0,0,0))
#   par(oma=c(5,4,4,1))
#   
#   ## Period 
#   {
#     { 
#       ## Other Response F
#       plot(NA, xlim = P_xlims, ylim = PC_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("a. Other response - female", x=2016, y=2, cex=.9)
#       # title("Age", adj = 0, line = 4)
#       
#       mtext("Log Odds", adj = 0.275, side = 2,
#             line = 2.5, cex = 0.8, outer = TRUE)
#       
#       mtext("Log Odds", adj = 0.85, side = 2,
#             line = 2.5, cex = 0.8, outer = TRUE)
#       
#       mtext("Period", adj = 0, side = 3, line = 0.5, outer = TRUE)
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = .9)
#       
#       ## Log odds axis
#       axis(2, at = seq(-8, 2, 1))
#       abline(h = 0)
#       ## Ests
#       lines(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == "Period"],
#             post_dkns_PC_plot$med_lo[post_dkns_PC_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == "Period"],
#                     rev(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
#                                                           "Period"])),
#               y = c(post_dkns_PC_plot$upper_lo[post_dkns_PC_plot$term ==
#                                                  "Period"],
#                     rev(post_dkns_PC_plot$lower_lo[post_dkns_PC_plot$term ==
#                                                      "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Other Response M
#       plot(NA, xlim = P_xlims, ylim = PC_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("b. Other response - male", x=2015.75, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1), )
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == "Period"],
#             post_dkns_PC_plot$med_lo[post_dkns_PC_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == "Period"],
#                     rev(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
#                                                           "Period"])),
#               y = c(post_dkns_PC_plot$upper_lo[post_dkns_PC_plot$term == "Period"],
#                     rev(post_dkns_PC_plot$lower_lo[post_dkns_PC_plot$term == "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer F
#       plot(NA, xlim = P_xlims, ylim = PC_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("c. Declined to answer - female", x = 2016.25, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                               "Period"],
#             post_ref_PC_plot$med_lo[post_ref_PC_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                       "Period"],
#                     rev(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                           "Period"])),
#               y = c(post_ref_PC_plot$upper_lo[post_ref_PC_plot$term ==
#                                                  "Period"],
#                     rev(post_ref_PC_plot$lower_lo[post_ref_PC_plot$term == 
#                                                      "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answre M
#       plot(NA, xlim = P_xlims, ylim =PC_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("d. Declined to answer - male", x=2016.1, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term ==
#                                               "Period"],
#             post_ref_PC_plot$med_lo[post_ref_PC_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term ==
#                                                       "Period"],
#                     rev(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                           "Period"])),
#               y = c(post_ref_PC_plot$upper_lo[post_ref_PC_plot$term == "Period"],
#                     rev(post_ref_PC_plot$lower_lo[post_ref_PC_plot$term == 
#                                                      "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#   }
#   
#   
#   ## Cohort 
#   {
#     { 
#       ## Other Response F
#       plot(NA, xlim = C_xlims, ylim = PC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("e. Other response - female", x= 18 + 1940, y=2, cex=.9)
#       mtext("Cohort", adj = 0, side = 3, line = 0.5)
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
#       
#       ## Log odds axis
#       
#       axis(2, at = seq(-8, 2, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term ==
#                                               "Cohort"],
#             post_dkns_PC_plot$med_lo[post_dkns_PC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_dkns_PC_plot$upper_lo[post_dkns_PC_plot$term == 
#                                                  "Cohort"],
#                     rev(post_dkns_PC_plot$lower_lo[post_dkns_PC_plot$term == 
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Other Response M
#       plot(NA, xlim = C_xlims, ylim = PC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("f. Other response - male", x=16.5 + 1940, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Cohort", line = 2, cex = .9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1), )
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == "Cohort"],
#             post_dkns_PC_plot$med_lo[post_dkns_PC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_dkns_PC_plot$term_idx_name[post_dkns_PC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_dkns_PC_plot$upper_lo[post_dkns_PC_plot$term == 
#                                                  "Cohort"],
#                     rev(post_dkns_PC_plot$lower_lo[post_dkns_PC_plot$term ==
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer F
#       plot(NA, xlim = C_xlims, ylim = PC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("g. Declined to answer - female", x = 21 + 1940, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == "Cohort"],
#             post_ref_PC_plot$med_lo[post_ref_PC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_ref_PC_plot$upper_lo[post_ref_PC_plot$term == "Cohort"],
#                     rev(post_ref_PC_plot$lower_lo[post_ref_PC_plot$term == 
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer M
#       plot(NA, xlim = C_xlims, ylim = PC_ylims,
#            xlab = "Cohort", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("h. Declined to answer - male", x = 19.5 + 1940, y=2, cex=.9)  
#       
#       ## Age axis
#       axis(1, at = seq(1940, 2010, 10))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       # y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == "Cohort"],
#             post_ref_PC_plot$med_lo[post_ref_PC_plot$term == "Cohort"],
#             lwd = 2, col = "forestgreen")
#       polygon(x = c(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                       "Cohort"],
#                     rev(post_ref_PC_plot$term_idx_name[post_ref_PC_plot$term == 
#                                                           "Cohort"])),
#               y = c(post_ref_PC_plot$upper_lo[post_ref_PC_plot$term == "Cohort"],
#                     rev(post_ref_PC_plot$lower_lo[post_ref_PC_plot$term == 
#                                                      "Cohort"])),
#               col = alpha("forestgreen", 0.35), border = FALSE)
#     }
#     
#   }
# }
# dev.off()
# 
# 
# ## AP ####
# AP_ranges <- bind_rows(post_dkns_AP_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Bisexual females"),
#                        post_dkns_AP_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Bisexual males"),
#                        post_ref_AP_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Lesbians"),
#                        post_ref_AP_plot %>% 
#                          group_by(term) %>% 
#                          dplyr::summarize(lower = min(lower_lo),
#                                    upper = max(upper_lo),
#                                    Group = "Gay males"))
# ## Across groups 
# AP_ranges %>% 
#   group_by(term) %>% 
#   dplyr::summarize(lower = min(lower),
#             upper = max(upper))
# 
# AP_ylims <- c(-8.5, 2.5)
# A_xlims <- c(14, 79) + c(-1, 1)
# P_xlims <- c(2014, 2021) + c(-1,1)
# 
# 
# # par(mfrow = c(2,4),lend = 1)
# 
# png("plots/GI_Q5/AP_8panel_logodds.png", 
#     width = 300*8, height = 300*6, res=300)
# {
#   par(mfrow = c(2,4), lend = 1)
#   par(mar=c(5,0,0,0))
#   par(oma=c(5,4,4,1))
#   
#   ## Age 
#   {
#     { 
#       ## Other Response F
#       plot(NA, xlim = A_xlims, ylim = AP_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("a. Other response - female", x= 35, y=2, cex=.9)
#       
#       mtext("Log Odds", adj = 0.275, side = 2,
#             line = 2.5, cex = 0.8, outer = TRUE)
#       
#       mtext("Log Odds", adj = 0.85, side = 2,
#             line = 2.5, cex = 0.8, outer = TRUE)
#       # title("Age", adj = 0, line = 4)
#       mtext("Age", adj = 0, side = 3, line = 0.5, outer = FALSE)
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = .9)
#       
#       ## Log odds axis
#       axis(2, at = seq(-8, 2, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Age"],
#             post_dkns_AP_plot$med_lo[post_dkns_AP_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Age"],
#                     rev(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == 
#                                                           "Age"])),
#               y = c(post_dkns_AP_plot$upper_lo[post_dkns_AP_plot$term == "Age"],
#                     rev(post_dkns_AP_plot$lower_lo[post_dkns_AP_plot$term == "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Other Response M
#       plot(NA, xlim = A_xlims, ylim = AP_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("b. Other response - male", x=33.5, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1), )
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Age"],
#             post_dkns_AP_plot$med_lo[post_dkns_AP_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Age"],
#                     rev(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == 
#                                                           "Age"])),
#               y = c(post_dkns_AP_plot$upper_lo[post_dkns_AP_plot$term == "Age"],
#                     rev(post_dkns_AP_plot$lower_lo[post_dkns_AP_plot$term == "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer F
#       plot(NA, xlim = A_xlims, ylim = AP_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("c. Declined to answer - female", x = 37.5, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == "Age"],
#             post_ref_AP_plot$med_lo[post_ref_AP_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == "Age"],
#                     rev(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == 
#                                                           "Age"])),
#               y = c(post_ref_AP_plot$upper_lo[post_ref_AP_plot$term == "Age"],
#                     rev(post_ref_AP_plot$lower_lo[post_ref_AP_plot$term == "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer M
#       plot(NA, xlim = A_xlims, ylim = AP_ylims,
#            xlab = "Age", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("d. Declined to answer - male", x=36, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = seq(15, 75, 5))
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == "Age"],
#             post_ref_AP_plot$med_lo[post_ref_AP_plot$term == "Age"],
#             lwd = 2, col = "navy")
#       polygon(x = c(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term ==
#                                                       "Age"],
#                     rev(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == 
#                                                           "Age"])),
#               y = c(post_ref_AP_plot$upper_lo[post_ref_AP_plot$term == "Age"],
#                     rev(post_ref_AP_plot$lower_lo[post_ref_AP_plot$term == 
#                                                      "Age"])),
#               col = alpha("navy", 0.35), border = FALSE)
#     }
#     
#   }
#   
#   
#   ## Period 
#   {
#     { 
#       ## Other Response F
#       plot(NA, xlim = P_xlims, ylim = AP_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("a. Other response - female", x=2016, y=2, cex=.9)
#       # mtext("Log Odds", adj = 0.25, side = 2,
#       #       line = 2.5, cex = 0.8, outer = TRUE)
#       # 
#       # mtext("Log Odds", adj = 0.75, side = 2,
#       #       line = 2.5, cex = 0.8, outer = TRUE)
#       
#       # title("Age", adj = 0, line = 4)
#       mtext("Period", adj = 0, side = 3, line = 0.5)
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = .9)
#       
#       ## Log odds axis
#       axis(2, at = seq(-8, 2, 1))
#       abline(h = 0)
#       ## Ests
#       lines(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Period"],
#             post_dkns_AP_plot$med_lo[post_dkns_AP_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Period"],
#                     rev(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == 
#                                                           "Period"])),
#               y = c(post_dkns_AP_plot$upper_lo[post_dkns_AP_plot$term ==
#                                                  "Period"],
#                     rev(post_dkns_AP_plot$lower_lo[post_dkns_AP_plot$term ==
#                                                      "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Other Response M
#       plot(NA, xlim = P_xlims, ylim = AP_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("b. Other response - male", x=2015.75, y=2, cex=.9)
#       
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1), )
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Period"],
#             post_dkns_AP_plot$med_lo[post_dkns_AP_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == "Period"],
#                     rev(post_dkns_AP_plot$term_idx_name[post_dkns_AP_plot$term == 
#                                                           "Period"])),
#               y = c(post_dkns_AP_plot$upper_lo[post_dkns_AP_plot$term == "Period"],
#                     rev(post_dkns_AP_plot$lower_lo[post_dkns_AP_plot$term == "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer F
#       plot(NA, xlim = P_xlims, ylim = AP_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("c. Declined to answer - female", x = 2016.25, y=2, cex=.9)
#       
#       
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == 
#                                               "Period"],
#             post_ref_AP_plot$med_lo[post_ref_AP_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == 
#                                                       "Period"],
#                     rev(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == 
#                                                           "Period"])),
#               y = c(post_ref_AP_plot$upper_lo[post_ref_AP_plot$term ==
#                                                  "Period"],
#                     rev(post_ref_AP_plot$lower_lo[post_ref_AP_plot$term == 
#                                                      "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#     { 
#       ## Declined to Answer M
#       plot(NA, xlim = P_xlims, ylim = AP_ylims,
#            xlab = "Period", ylab = "Log Odds",
#            main = "", type = "n", frame.plot = TRUE, axes = FALSE)
#       ## title
#       text("d. Declined to answer - male", x=2016.1, y=2, cex=.9) 
#       
#       ## Age axis
#       axis(1, at = 2014:2021)
#       # mtext(side = 1, text = "Age", line = 2, cex = 0.9)
#       
#       ## Log odds axis
#       # axis(2, at = seq(-7, 1, 1))
#       
#       ## y = 0
#       abline(h = 0)
#       
#       ## Ests
#       lines(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term ==
#                                               "Period"],
#             post_ref_AP_plot$med_lo[post_ref_AP_plot$term == "Period"],
#             lwd = 2, col = "goldenrod")
#       polygon(x = c(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term ==
#                                                       "Period"],
#                     rev(post_ref_AP_plot$term_idx_name[post_ref_AP_plot$term == 
#                                                           "Period"])),
#               y = c(post_ref_AP_plot$upper_lo[post_ref_AP_plot$term == "Period"],
#                     rev(post_ref_AP_plot$lower_lo[post_ref_AP_plot$term == 
#                                                      "Period"])),
#               col = alpha("goldenrod", 0.35), border = FALSE)
#     }
#     
#   }
# }
# dev.off()
# 
# # Grand proportions for #s in MS body
# 
# yrbs21 <- readRDS("data - raw/yrbs21_tomerge.RDS")
# prop.table(tableNA(yrbs21$GI_21, yrbs21$sex), 2)
# 
