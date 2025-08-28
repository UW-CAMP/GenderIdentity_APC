## Code for Analysis 3, Figure 4
## of the following manuscript:
##
##    Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025.
##      Demographic trends in gender identity among adults in the United States,
##      2014-2021. International Journal of Transgender Health. 
##      Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.
## 
##  Script authors: Godwin J
##
##  This script fits  null (intercept-only), A, P, C, AP, AC, PC
##  Bayesian smoothing models to survey-weighted prevalence estimates
##  for don't know/not sure and those who decline to answer 
##  by age x period (survey) x cohort
##
##  Inputs:
##     -  data - clean/data - clean/BRFSS_HT_GI_prevs.csv
##
##  Outputs:
##     -  tables/analysis 3/mod_[gender]_summaries.txt
##     -  tables/analysis 3/mod_[gender]_summaries.txt
##     -  data - clean/[gender]_PC_posteriors.csv (Posterior samples used in plots)
##     -  plots/analysis 3/dkns_PC_coh_logodds.png (single panel Fig 3)
##     -  plots/analysis 3/dkns_PC_per_logodds.png (single panel Fig 3)
##     -  plots/analysis 3/dta_PC_coh_logodds.png (single panel Fig 3)
##     -  plots/analysis 3/dta_PC_per_logodds.png (single panel Fig 3)
##     -  plots/analysis 3/Fig3_BestModel_4panel_logodds.png (Fig 3 as in manuscript)

# Setup ####
rm(list = ls())

## Folders ####
if(!dir.exists("tables/analysis 3")){
  if(!dir.exists("tables/")){
    dir.create("tables/")
  }
  dir.create("tables/analysis 3")
}

if(!dir.exists("plots/analysis 3")){
  if(!dir.exists("plots/")){
    dir.create("plots/")
  }
  dir.create("plots/analysis 3")
}

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

write.csv(dkns_criteria, file = paste0("tables/analysis 3/",
                                       "dkns_critera.csv"),
          row.names = FALSE)

#### Save Model Outputs ####
sink(file = "tables/analysis 3/mod_dkns_summaries.txt")
summary(mod_dkns_null)
summary(mod_dkns_A)
summary(mod_dkns_P)
summary(mod_dkns_C)
summary(mod_dkns_AP)
summary(mod_dkns_AC)
summary(mod_dkns_PC)
sink(file = NULL)

## Plots ####
age_idx <- mod_data %>% 
  select(age, age_idx) %>% 
  unique()
cohort_idx <- mod_data %>% 
  select(cohort, cohort_idx) %>% 
  unique() %>% 
  arrange(cohort_idx)

period_idx <- mod_data %>% 
  select(period, period_idx) %>% 
  unique() %>% 
  arrange(period_idx)

### PC ####
post_dkns_PC <- inla.posterior.sample(n = 1000, mod_dkns_PC)

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

write.csv(post_dkns_PC_plot, file = "data - clean/dkns_PC_posteriors.csv",
          row.names = FALSE)

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

ggsave(filename = "plots/analysis 3/dkns_PC_per_logodds.png", dkns_PC_per_lo)

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

ggsave(filename = "plots/analysis 3/dkns_PC_coh_logodds.png", dkns_PC_coh_lo)

## Decline ####
mod_dta_null <- inla(mod_null,
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
mod_dta_P <- inla(mod_P,
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

mod_dta_C <- inla(mod_C,
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

mod_dta_A <- inla(mod_A,
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

mod_dta_AP <- inla(mod_AP,
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

mod_dta_AC <- inla(mod_AC,
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

mod_dta_PC <- inla(mod_PC,
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
dta_criteria <- data.frame(Model = c("Null", "A", "P", "C",
                                      "AP", "AC", "PC"),
                            loglik = c(mod_dta_null$mlik[2],
                                       mod_dta_A$mlik[2],
                                       mod_dta_P$mlik[2],
                                       mod_dta_C$mlik[2],
                                       mod_dta_AP$mlik[2],
                                       mod_dta_AC$mlik[2],
                                       mod_dta_PC$mlik[2]),
                            WAIC = c(mod_dta_null$waic$waic,
                                     mod_dta_A$waic$waic,
                                     mod_dta_P$waic$waic,
                                     mod_dta_C$waic$waic,
                                     mod_dta_AP$waic$waic,
                                     mod_dta_AC$waic$waic,
                                     mod_dta_PC$waic$waic),
                            DIC = c(mod_dta_null$dic$dic,
                                    mod_dta_A$dic$dic,
                                    mod_dta_P$dic$dic,
                                    mod_dta_C$dic$dic,
                                    mod_dta_AP$dic$dic,
                                    mod_dta_AC$dic$dic,
                                    mod_dta_PC$dic$dic))

write.csv(dta_criteria, file = paste0("tables/analysis 3/",
                                       "dta_critera.csv"),
          row.names = FALSE)

#### Save Model Outputs ####
sink(file = "tables/analysis 3/mod_dta_summaries.txt")
summary(mod_dta_null)
summary(mod_dta_A)
summary(mod_dta_P)
summary(mod_dta_C)
summary(mod_dta_AP)
summary(mod_dta_AC)
summary(mod_dta_PC)
sink(file = NULL)

## Plots ####
### PC ####

post_dta_PC <- inla.posterior.sample(n = 1000, mod_dta_PC)

post_dta_PC_df <- lapply(post_dta_PC, function(draw){
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

post_dta_PC_df <- do.call(rbind.data.frame, post_dta_PC_df) %>% 
  group_by(term, term_idx) %>% 
  mutate(draw_no = 1:length(term))

post_dta_PC_plot <- post_dta_PC_df %>% 
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

write.csv(post_dta_PC_plot, file = "data - clean/dta_PC_posteriors.csv",
          row.names = FALSE)

dta_PC_per_lo <- post_dta_PC_plot %>%
  filter(term == "Period") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "goldenrod") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "goldenrod",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Period")+
  theme_classic()

ggsave(filename = "plots/analysis 3/dta_PC_per_logodds.png", dta_PC_per_lo)


dta_PC_coh_lo <- post_dta_PC_plot %>%
  filter(term == "Cohort") %>% 
  ggplot(aes(x = term_idx_name, y = med_lo)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = lower_lo, ymax = upper_lo), fill = "forestgreen",
              color = "white", alpha = 0.3) +
  ggtitle("Declined to Answer") +
  ylab("Log Odds") +
  xlab("Cohort")+
  theme_classic()

ggsave(filename = "plots/analysis 3/dta_PC_coh_logodds.png", dta_PC_coh_lo)

# Manuscript Plot: Figure 3 ####
## Get min and max endpoints of intervals by sexual orientation
## and age, period, or cohort across all models
best_ranges <- 
  bind_rows(post_dkns_PC_plot %>% 
              group_by(term) %>% 
              dplyr::summarize(lower = min(lower_lo),
                               upper = max(upper_lo),
                               Group = "Don't Know/Not Sure"),
            post_dta_PC_plot %>% 
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

png("plots/analysis 3/BestModel_4panel_logodds.png", 
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
      lines(post_dta_PC_plot$term_idx_name[post_dta_PC_plot$term == "Cohort"],
            post_dta_PC_plot$med_lo[post_dta_PC_plot$term == "Cohort"],
            lwd = 2, col = "forestgreen")
      polygon(x = c(post_dta_PC_plot$term_idx_name[post_dta_PC_plot$term == 
                                                    "Cohort"],
                    rev(post_dta_PC_plot$term_idx_name[post_dta_PC_plot$term == 
                                                        "Cohort"])),
              y = c(post_dta_PC_plot$upper_lo[post_dta_PC_plot$term == 
                                               "Cohort"],
                    rev(post_dta_PC_plot$lower_lo[post_dta_PC_plot$term ==
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
      lines(post_dta_PC_plot$term_idx_name[post_dta_PC_plot$term == "Period"],
            post_dta_PC_plot$med_lo[post_dta_PC_plot$term == "Period"],
            lwd = 2, col = "navy")
      polygon(x = c(post_dta_PC_plot$term_idx_name[post_dta_PC_plot$term == 
                                                    "Period"],
                    rev(post_dta_PC_plot$term_idx_name[post_dta_PC_plot$term == 
                                                        "Period"])),
              y = c(post_dta_PC_plot$upper_lo[post_dta_PC_plot$term == 
                                               "Period"],
                    rev(post_dta_PC_plot$lower_lo[post_dta_PC_plot$term ==
                                                   "Period"])),
              col = alpha("navy", 0.35), border = FALSE)
    }
  }
}
dev.off()