## Q4: Mann-Kendall within-cohort tests for Ho: decrease in prevalence or non-monotonic
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
library(trend)


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
cohorts <- unique(brfss_des$variables$cohort) %>%  sort()
for(coh in cohorts){
  for(gend in c("1_transwoman", "2_transman", "3_nbgnc")){
    GI_init <- case_when(gend == "1_transwoman" ~ "tw",
                         gend == "2_transman" ~ "tm",
                         gend == "3_nbgnc" ~ "nb",
                         TRUE ~ NA)
    brfss_means[[paste0("coh_", coh, "_", GI_init)]] <- 
      svyby(formula(paste0("~", GI_init, "_bin")), by = ~year, 
            design = subset(brfss_des, cohort == coh), svymean) %>% 
      mutate(cohort = coh,
             gender = gend,
             source = "BRFSS")
    
    names(brfss_means[[paste0("coh_", coh, "_", GI_init)]])[2] <- "GI_prev"
  }
}

brfss_means_df <- do.call(rbind.data.frame, brfss_means) %>% 
  select(gender, cohort, year, GI_prev, se, source) %>% 
  arrange(gender, cohort, year) %>% 
  group_by(gender, cohort) %>% 
  mutate(n_per_coh = n()) %>% 
  ungroup() %>% 
  filter(n_per_coh >= 5) 

# Mann-Kendall Test ####

mk_output <- list()

for(coh in unique(brfss_means_df$cohort)){
  coh_dat <- brfss_means_df %>% 
    filter(cohort == coh) %>% 
    arrange(gender, year)
  
  for(gend in unique(coh_dat$gender)){
    gend_dat <- coh_dat %>% 
      filter(gender == gend) %>% 
      arrange(year)
    
    test_out <- mk.test(gend_dat$GI_prev, alternative = "greater")
    mk_output[[paste0(gend, "__", coh)]] <-
      data.frame(cohort = coh, gender = gend,
                 test_stat = test_out$statistic,
                 pval = test_out$p.value,
                 n_points = nrow(gend_dat))
     
  }
}

mk_output_df <- do.call(rbind.data.frame, mk_output) %>% 
  arrange(gender, cohort) %>% 
  relocate("gender", .before = "cohort")

# Plots ####

gend_cols <- c("forestgreen", "goldenrod", "navy")

mk_plot <- mk_output_df %>% 
  mutate(Cohort = cohort,
         Gender = case_when(grepl("1_", gender) ~ "Transgender woman",
                            grepl("2_", gender) ~ "Transgender man",
                            grepl("3_", gender) ~ "NB/GNC",
                            TRUE ~NA)) %>% 
  ggplot(aes(x = Cohort, y = pval, color = Gender, group = Gender)) +
  geom_hline(yintercept = 0.10, linetype = 2) +
  geom_hline(yintercept = 0.05) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("Transgender woman" = gend_cols[1],
                                "Transgender man" = gend_cols[2],
                                "NB/GNC" = gend_cols[3])) +
  scale_x_continuous(breaks = seq(1988, 1998, 2),
                     labels = paste0("'", seq(88, 98, 2))) +
  facet_wrap(~Gender) +
  ylab("p-value") +
  ggtitle("Mann-Kendall Tests by Gender & Cohort") +
  theme_classic()
    
ggsave(filename = "plots/GI_Q4/MannKendall.png",
       plot = mk_plot,
       height = 4, width = 6, units = "in")
