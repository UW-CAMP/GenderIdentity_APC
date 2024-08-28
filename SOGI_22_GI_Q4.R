## Q4: Pair-wise within-cohort tests for decrease in prevalence
## Jessica Godwin

# Setup ####
rm(list = ls())

if(!dir.exists("tables/GI_Q4")){
  if(!dir.exists("tables")){
    dir.create("tables")
  }
  dir.create("tables/GI_Q4")
}

if(!dir.exists("plots/GI_Q4")){
  if(!dir.exists("plots")){
    dir.create("plots")
  }
  dir.create("plots/GI_Q4")
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
  ## Filter to correct ages & periods
  filter(age <= 30 & source == "BRFSS") %>% 
  filter(year %in% seq(2015, 2021, 2)) %>% 
  ## Create binary indicators of SO
  mutate(tw_bin = ifelse(gender == "1_transwoman", 1, 0),
         tm_bin = ifelse(gender == "2_transman", 1, 0),
         nb_bin = ifelse(gender == "3_nbgnc", 1, 0),
         cohort_2_odd = ifelse(cohort %% 2 == 0, cohort - 1, cohort),
         cohort_2_odd = as.character(cohort_2_odd)) 

# Specify design ####

brfss_des <- svydesign(ids = ~1, strata = ~year + stratum, weights = ~weight,
                       data = mod_data, nest = TRUE)
brfss_des

# Get means and contrasts by cohort and period ####

## BRFSS ####

brfss_means <- list()
cohorts <- unique(brfss_des$variables$cohort_2_odd) %>%  sort()
for(coh in cohorts){
  for(gend in c("1_transwoman", "2_transman", "3_nbgnc")){
    GI_init <- case_when(gend == "1_transwoman" ~ "tw",
                         gend == "2_transman" ~ "tm",
                         gend == "3_nbgnc" ~ "nb",
                         TRUE ~ NA)
    brfss_means[[paste0("coh_", coh, "_", GI_init)]] <- 
      svyby(formula(paste0("~", GI_init, "_bin")), by = ~year, 
            design = subset(brfss_des, cohort_2_odd == coh), svymean) %>% 
      mutate(cohort_2_odd = coh,
             gender = gend,
             source = "BRFSS")
    
    names(brfss_means[[paste0("coh_", coh, "_", GI_init)]])[2] <- "GI_prev"
  }
}

brfss_means_df <- do.call(rbind.data.frame, brfss_means)

# Combine estimates ####

means_df <- brfss_means_df %>% 
  pivot_wider(id_cols = c("cohort_2_odd", "gender", "source"),
              names_from = "year",
              names_glue = "{.value}_{year}",
              values_from = c("GI_prev", "se")) %>% 
  mutate(diff_15_17 = GI_prev_2017 - GI_prev_2015,
         diff_17_19 = GI_prev_2019 - GI_prev_2017,
         diff_19_21 = GI_prev_2021 - GI_prev_2019,
         se_15_17 = sqrt(se_2015^2 + se_2017^2),
         se_17_19 = sqrt(se_2017^2 + se_2019^2),
         se_19_21 = sqrt(se_2019^2 + se_2021^2),
         lower_15_17 = diff_15_17 - qnorm(.975)*se_15_17,
         lower_17_19 = diff_17_19 - qnorm(.975)*se_17_19,
         lower_19_21 = diff_19_21 - qnorm(.975)*se_19_21,
         upper_15_17 = diff_15_17 + qnorm(.975)*se_15_17,
         upper_17_19 = diff_17_19 + qnorm(.975)*se_17_19,
         upper_19_21 = diff_19_21 + qnorm(.975)*se_19_21,
         upper_os_15_17 = diff_15_17 + qnorm(.95)*se_15_17,
         upper_os_17_19 = diff_17_19 + qnorm(.95)*se_17_19,
         upper_os_19_21 = diff_19_21 + qnorm(.95)*se_19_21,
         stat_15_17 = diff_15_17/se_15_17,
         stat_17_19 = diff_17_19/se_17_19,
         stat_19_21 = diff_19_21/se_19_21,
         pval_15_17 = pnorm(stat_15_17),
         pval_17_19 = pnorm(stat_17_19),
         pval_19_21 = pnorm(stat_19_21))

# Plots ####

year_cols <- c("navy", "forestgreen", "goldenrod", "firebrick")

### Transwoman ####
plotdat <- means_df %>%
  filter(gender == "1_transwoman") %>%
  mutate(Cohort_num = as.numeric(cohort_2_odd),
         Cohort = cohort_2_odd)  %>%
  # filter(Cohort < 2001) %>%
  group_by(Cohort_num, Cohort, gender) %>%
  # reframe(diff_15_17 = diff_15_17[!is.na(upper_os_15_17)],
  #         diff_17_19 = diff_17_19[!is.na(upper_os_17_19)],
  #         diff_19_21 = diff_19_21[!is.na(upper_os_19_21)],
  #         upper_os_15_17 = upper_os_15_17[!is.na(upper_os_15_17)],
  #         upper_os_17_19 = upper_os_17_19[!is.na(upper_os_17_19)],
  #         upper_os_19_21 = upper_os_19_21[!is.na(upper_os_19_21)])
  reframe(diff_15_17 = ifelse(sum(is.na(upper_os_15_17)) ==
                                length(upper_os_15_17),
                              NA, diff_15_17[!is.na(upper_os_15_17)]),
          diff_17_19 = ifelse(sum(is.na(upper_os_17_19)) ==
                                length(upper_os_17_19),
                              NA, diff_17_19[!is.na(upper_os_17_19)]),
          diff_19_21 = ifelse(sum(is.na(upper_os_19_21)) ==
                                length(upper_os_19_21),
                              NA, diff_19_21[!is.na(upper_os_19_21)]),
          upper_os_15_17 = ifelse(sum(is.na(upper_os_15_17)) ==
                                    length(upper_os_15_17),
                                  NA, upper_os_15_17[!is.na(upper_os_15_17)]),
          upper_os_17_19 = ifelse(sum(is.na(upper_os_17_19)) ==
                                    length(upper_os_17_19),
                                  NA, upper_os_17_19[!is.na(upper_os_17_19)]),
          upper_os_19_21 = ifelse(sum(is.na(upper_os_19_21)) ==
                                    length(upper_os_19_21),
                                  NA, upper_os_19_21[!is.na(upper_os_19_21)]))

## Works for all four plots
y_lims <- c(-0.04, 0.04)
x_lims <- c(1985, 2007)
n_coh <- 10

png("plots/GI_Q4/Q4.png", width = 960, height = 600)
par(mfrow = c(2,2),lend = 1)
par(mar=c(0,0,0,0))
par(oma=c(7,5,4,1))
{
  #par(mfrow = c(1,1),lend = 1)
  plot(plotdat$Cohort_num, plotdat$upper_os_15_17,
       xlab = "", ylab = "",
       ylim = y_lims, xlim = c(1, n_coh*2) + c(-0.5, 0.5),
       type = "n", xaxt = 'n', yaxt='n',
       main = "")
  
  axis(1, at = seq(1, n_coh*2, 2),
       tick=F, labels=F, cex.axis = 0.9, las = 2)
  
  axis(2, at = seq(-0.04, 0.04, 0.01), cex.axis = 0.9)
  abline(h = 0)
  mtext("Prevalence Change", 2, 3, cex = 0.8)
  
  ## 2015 - 2017
  segments(x0 = seq(1, n_coh*2, 2) - .3, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) - .3, y1 = rev(plotdat$upper_os_15_17),
           lwd = 2, col = alpha(year_cols[1], 0.3))
  points(seq(1, n_coh*2, 2) - .3,
         rev(plotdat$diff_15_17), pch = 19, col = year_cols[1])
  segments(x0 = seq(1, n_coh*2, 2) - .3, y0 = rev(plotdat$diff_15_17),
           x1 = seq(1, n_coh*2, 2) - .3, y1 = rev(plotdat$upper_os_15_17),
           lwd = 2, col = year_cols[1])
  
  ## 2017 - 2019
  segments(x0 = seq(1, n_coh*2, 2) - 0, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) - 0, y1 = rev(plotdat$upper_os_17_19),
           lwd = 2, col = alpha(year_cols[2], 0.3))
  points(seq(1, n_coh*2, 2) - 0,
         rev(plotdat$diff_17_19), pch = 19, col = year_cols[2])
  segments(x0 = seq(1, n_coh*2, 2) - 0, y0 = rev(plotdat$diff_17_19),
           x1 = seq(1, n_coh*2, 2) - 0, y1 = rev(plotdat$upper_os_17_19),
           lwd = 2, col = year_cols[2])
  ## 2019 - 2021
  segments(x0 = seq(1, n_coh*2, 2) + .3, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) + .3, y1 = rev(plotdat$upper_os_19_21),
           lwd = 2, col = alpha(year_cols[3], 0.3))
  points(seq(1, n_coh*2, 2) + .3,
         rev(plotdat$diff_19_21), pch = 19, col = year_cols[3])
  segments(x0 = seq(1, n_coh*2, 2) + .3, y0 = rev(plotdat$diff_19_21),
           x1 = seq(1, n_coh*2, 2) + .3, y1 = rev(plotdat$upper_os_19_21),
           lwd = 2, col = year_cols[3])
  
  
  legend("topright", pch = 19, col = year_cols[1:3],
         # lwd = 2
         legend = paste0(seq(2015, 2019, 2), "-",
                         seq(2017, 2021, 2)),
         bty = "n", title = "Period", adj = 0,
         ncol = 1, cex = 1)
  text("Transwoman", x=2, y=0.035, cex=1.2)
}


### Transman ####
plotdat <- means_df %>%
  filter(gender == "2_transman") %>%
  mutate(Cohort = cohort_2_odd,
         Cohort_num = as.numeric(Cohort)) %>%
  # filter(Cohort < 2001) %>%
  group_by(Cohort_num, Cohort, gender) %>%
  reframe(diff_15_17 = ifelse(sum(is.na(upper_os_15_17)) ==
                                length(upper_os_15_17),
                              NA, diff_15_17[!is.na(upper_os_15_17)]),
          diff_17_19 = ifelse(sum(is.na(upper_os_17_19)) ==
                                length(upper_os_17_19),
                              NA, diff_17_19[!is.na(upper_os_17_19)]),
          diff_19_21 = ifelse(sum(is.na(upper_os_19_21)) ==
                                length(upper_os_19_21),
                              NA, diff_19_21[!is.na(upper_os_19_21)]),
          upper_os_15_17 = ifelse(sum(is.na(upper_os_15_17)) ==
                                    length(upper_os_15_17),
                                  NA, upper_os_15_17[!is.na(upper_os_15_17)]),
          upper_os_17_19 = ifelse(sum(is.na(upper_os_17_19)) ==
                                    length(upper_os_17_19),
                                  NA, upper_os_17_19[!is.na(upper_os_17_19)]),
          upper_os_19_21 = ifelse(sum(is.na(upper_os_19_21)) ==
                                    length(upper_os_19_21),
                                  NA, upper_os_19_21[!is.na(upper_os_19_21)]))

{  
  plot(plotdat$Cohort_num, plotdat$diff_15_17,
       xlab = "", ylab = "",
       ylim = y_lims, xlim = c(1, n_coh*2) + c(-0.5, 0.5),
       type = "n", xaxt = 'n', yaxt='n',
       main = "")
  axis(1, at = seq(1, n_coh*2, 2),
       labels = rev(paste0(seq(x_lims[1], x_lims[2]-4, 2), "-",
                           seq(x_lims[1], x_lims[2]-4, 2) + 1)),
       tick=F, cex.axis = 0.9, las = 2)

  abline(h = 0)
  
  ## 2015 - 2017
  segments(x0 = seq(1, n_coh*2, 2) - .3, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) - .3, y1 = rev(plotdat$upper_os_15_17),
           lwd = 2, col = alpha(year_cols[1], 0.3))
  points(seq(1, n_coh*2, 2) - .3,
         rev(plotdat$diff_15_17), pch = 19, col = year_cols[1])
  segments(x0 = seq(1, n_coh*2, 2) - .3, y0 = rev(plotdat$diff_15_17),
           x1 = seq(1, n_coh*2, 2) - .3, y1 = rev(plotdat$upper_os_15_17),
           lwd = 2, col = year_cols[1])
  
  ## 2017 - 2019
  segments(x0 = seq(1, n_coh*2, 2) - 0, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) - 0, y1 = rev(plotdat$upper_os_17_19),
           lwd = 2, col = alpha(year_cols[2], 0.3))
  points(seq(1, n_coh*2, 2) - 0,
         rev(plotdat$diff_17_19), pch = 19, col = year_cols[2])
  segments(x0 = seq(1, n_coh*2, 2) - 0, y0 = rev(plotdat$diff_17_19),
           x1 = seq(1, n_coh*2, 2) - 0, y1 = rev(plotdat$upper_os_17_19),
           lwd = 2, col = year_cols[2])
  ## 2019 - 2021
  segments(x0 = seq(1, n_coh*2, 2) + .3, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) + .3, y1 = rev(plotdat$upper_os_19_21),
           lwd = 2, col = alpha(year_cols[3], 0.3))
  points(seq(1, n_coh*2, 2) + .3,
         rev(plotdat$diff_19_21), pch = 19, col = year_cols[3])
  segments(x0 = seq(1, n_coh*2, 2) + .3, y0 = rev(plotdat$diff_19_21),
           x1 = seq(1, n_coh*2, 2) + .3, y1 = rev(plotdat$upper_os_19_21),
           lwd = 2, col = year_cols[3])
  
  # legend("topright", pch = 19, col = year_cols[1:3],
  #        # lwd = 2
  #        legend = paste0(seq(2015, 2019, 2), "-",
  #                        seq(2017, 2021, 2)),
  #        bty = "n", title = "Period", adj = 0,
  #        ncol = 1, cex = 1)
  text("Transman", x=2, y=0.035, cex=1.2)
}


### NB/GNC ####
plotdat <- means_df %>%
  filter(gender == "3_nbgnc") %>%
  mutate(Cohort = cohort_2_odd,
         Cohort_num = as.numeric(Cohort)) %>%
  # filter(Cohort < 2001) %>%
  group_by(Cohort_num, Cohort, gender) %>%
  reframe(diff_15_17 = ifelse(sum(is.na(upper_os_15_17)) ==
                                length(upper_os_15_17),
                              NA, diff_15_17[!is.na(upper_os_15_17)]),
          diff_17_19 = ifelse(sum(is.na(upper_os_17_19)) ==
                                length(upper_os_17_19),
                              NA, diff_17_19[!is.na(upper_os_17_19)]),
          diff_19_21 = ifelse(sum(is.na(upper_os_19_21)) ==
                                length(upper_os_19_21),
                              NA, diff_19_21[!is.na(upper_os_19_21)]),
          upper_os_15_17 = ifelse(sum(is.na(upper_os_15_17)) ==
                                    length(upper_os_15_17),
                                  NA, upper_os_15_17[!is.na(upper_os_15_17)]),
          upper_os_17_19 = ifelse(sum(is.na(upper_os_17_19)) ==
                                    length(upper_os_17_19),
                                  NA, upper_os_17_19[!is.na(upper_os_17_19)]),
          upper_os_19_21 = ifelse(sum(is.na(upper_os_19_21)) ==
                                    length(upper_os_19_21),
                                  NA, upper_os_19_21[!is.na(upper_os_19_21)]))

{
  par(lend = 1)
  
  plot(plotdat$Cohort_num, plotdat$diff_15_17,
       xlab = "", ylab = "",
       ylim = y_lims, xlim = c(1, n_coh*2) + c(-0.5, 0.5),
       type = "n", xaxt = 'n', yaxt='n',
       main = "")
  axis(1, at = seq(1, n_coh*2, 2),
       labels = rev(paste0(seq(x_lims[1], x_lims[2]-4, 2), "-",
                           seq(x_lims[1], x_lims[2]-4, 2) + 1)),
       tick=F, cex.axis = 0.9, las = 2)
  
  axis(2, at = seq(-0.04, 0.04, 0.01), cex.axis = 0.9)
  abline(h = 0)
  mtext("Cohort", 1, 5, cex = 0.8)
  mtext("Prevalence Change", 2, 3, cex = 0.8)
  
  ## 2015 - 2017
  segments(x0 = seq(1, n_coh*2, 2) - .3, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) - .3, y1 = rev(plotdat$upper_os_15_17),
           lwd = 2, col = alpha(year_cols[1], 0.3))
  points(seq(1, n_coh*2, 2) - .3,
         rev(plotdat$diff_15_17), pch = 19, col = year_cols[1])
  segments(x0 = seq(1, n_coh*2, 2) - .3, y0 = rev(plotdat$diff_15_17),
           x1 = seq(1, n_coh*2, 2) - .3, y1 = rev(plotdat$upper_os_15_17),
           lwd = 2, col = year_cols[1])
  
  ## 2017 - 2019
  segments(x0 = seq(1, n_coh*2, 2) - 0, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) - 0, y1 = rev(plotdat$upper_os_17_19),
           lwd = 2, col = alpha(year_cols[2], 0.3))
  points(seq(1, n_coh*2, 2) - 0,
         rev(plotdat$diff_17_19), pch = 19, col = year_cols[2])
  segments(x0 = seq(1, n_coh*2, 2) - 0, y0 = rev(plotdat$diff_17_19),
           x1 = seq(1, n_coh*2, 2) - 0, y1 = rev(plotdat$upper_os_17_19),
           lwd = 2, col = year_cols[2])
  ## 2019 - 2021
  segments(x0 = seq(1, n_coh*2, 2) + .3, y0 = y_lims[1],
           x1 = seq(1, n_coh*2, 2) + .3, y1 = rev(plotdat$upper_os_19_21),
           lwd = 2, col = alpha(year_cols[3], 0.3))
  points(seq(1, n_coh*2, 2) + .3,
         rev(plotdat$diff_19_21), pch = 19, col = year_cols[3])
  segments(x0 = seq(1, n_coh*2, 2) + .3, y0 = rev(plotdat$diff_19_21),
           x1 = seq(1, n_coh*2, 2) + .3, y1 = rev(plotdat$upper_os_19_21),
           lwd = 2, col = year_cols[3])
  # 
  # legend("topright", pch = 19, col = year_cols[1:3],
  #        # lwd = 2
  #        legend = paste0(seq(2015, 2019, 2), "-",
  #                        seq(2017, 2021, 2)),
  #        bty = "n", title = "Period", adj = 0,
  #        ncol = 1, cex = 1)
  text("NB/GNC", x=2, y=0.035, cex=1.2)
}
dev.off()
