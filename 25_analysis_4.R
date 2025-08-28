## Code for Analysis 3, Figures 5 and 6 and Table S2
## of the following manuscript:
##
##    Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025.
##      Demographic trends in gender identity among adults in the United States,
##      2014-2021. International Journal of Transgender Health. 
##      Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.
## 
##  Script authors: Godwin J
##
##  This script explores the multivariate relationship between
##  the BRFSS variables: sex, sex assigned at birth, 
##  sexual orientation, and gender identity. 
##  First, survey-weighted estimates of proportions by responses are
##  estimated. Then, estimates are visualized, resulting manuscript
##  Figures 5 & 6.
##
##  Inputs:
##     -  data - clean/brfss_final.rds
##
##  Outputs:
##     -  plots/sex-so-gi/Fig5_sex_stream_sab_by_gi.png
##     -  plots/sex-so-gi/Fig6_so_by_gi.png
##     -  tables/sex-gi/multiCA_sex_[gender].txt
##     -  tables/so-gi/multiCA_so_[gender].txt

# Setup ####
rm(list = ls())

## Folders ####
if(!dir.exists("tables/so-gi/")){
  if(!dir.exists("tables/")){
    dir.create("tables/so-gi/")
  }
  dir.create("tables/so-gi/")
}

if(!dir.exists("tables/sex-gi/")){
  if(!dir.exists("tables/")){
    dir.create("tables/sex-gi/")
  }
  dir.create("tables/sex-gi/")
}

if(!dir.exists("plots/sex-so-gi/")){
  if(!dir.exists("plots/")){
    dir.create("plots/")
  }
  dir.create("plots/sex-so-gi/")
}

## Functions ####
tableNA <- function(x, ...){
   table(x, useNA = "ifany", ...)  
}

# Load Data ####
brfss <- readRDS("data - clean/brfss_final.rds")
head(brfss)
str(brfss)
summary(brfss)

# Prep data ####
data_1x1 <- brfss %>%
  select(year, age, cohort,
         sex, gender, gender, so, so_new, sab, sex,
         contains("cohort_"), contains("pd_"),
         index, weight, strata) %>% 
  mutate(tw_bin = ifelse(gender == "1_Transgender woman", 1, 0),
         tm_bin = ifelse(gender == "2_Transgender man", 1, 0),
         nbgnc_bin = ifelse(gender == "3_nbgnc", 1, 0),
         cis_bin = ifelse(gender == "4_cisgender", 1, 0),
         dkns_bin = ifelse(gender == "5_DKNS", 1, 0),
         ref_bin = ifelse(gender == "6_ref", 1, 0)) %>%
  mutate(across(where(is.factor), ~as.character(.x))) %>% 
  rename("period" = "year")

## Check out object
summary(data_1x1)
tableNA(data_1x1$age, data_1x1$cohort)
tableNA(data_1x1$period, data_1x1$cohort)
tableNA(data_1x1$age, data_1x1$period)

# Get estimates ####
## Specify survey design object 
brfss_des <- svydesign(ids = ~1, strata = ~period + strata,
                       weights = ~weight, data = data_1x1)

brfss_des

## Get weighted estimates

sex_prop <- svyby(~sex, by = ~period + gender, brfss_des, svymean)
so_prop <- svyby(~so_new, by = ~period + gender, brfss_des, svymean)
sab_prop <- svyby(~sab, by = ~period + gender, 
                  subset(brfss_des, period >= 2019 & !is.na(sab)),
                  svymean, na.rm.all = TRUE)


# Plots ####
## Sex (Figure 5, Panel A) ####
names(sex_prop)[3:4] <- paste0("mean.", names(sex_prop)[3:4])
names(sex_prop) <- gsub("sex", "", names(sex_prop))
gend_col <- c("navy", "forestgreen", "firebrick",
              "goldenrod", "grey80")
gend_names <- c("Transgender woman", "Transgender man",
                "NB/GNC", "Don't know/Not sure",
                "Declined to Answer")
names(gend_col) <- gend_names

sex_stream_plot <- sex_prop %>% 
  pivot_longer(cols = contains("."),
               names_to = c(".value", "sex"),
               names_pattern = "(.*)\\.(.*)") %>% 
  mutate(gender = case_when(grepl("1_", gender) ~ "Transgender woman",
                            grepl("2_", gender) ~ "Transgender man",
                            grepl("3_", gender) ~ "NB/GNC",
                            grepl("4_", gender) ~ "Cisgender",
                            grepl("5_", gender) ~ "Don't know/Not sure",
                            grepl("6_", gender) ~ "Declined to Answer")) %>% 
  filter(gender != "Cisgender") %>% 
  mutate(gender = factor(gender, levels = gend_names),
         sex = factor(sex, levels = c("Female", "Male"))) %>% 
  ggplot() +
  xlab("Period") +
  ylab("Proportion") +
  ggtitle("A) Sex") +
  geom_stream(aes(x = period, y = mean, fill = sex), type = "proportion",
              show.legend = FALSE) +
  scale_fill_manual(name = "Sex",
                     values = c("Female" = "navy",
                                "Male" = "forestgreen")) +
  facet_wrap(~gender, ncol = 5) +
  theme_classic()
sex_stream_plot

## Sex & SAB (Figure 5)####
names(sab_prop)[3:6] <- paste0("mean.", c("male", "female", "dkns", "ref"))
names(sab_prop) <- gsub("sab1_", "", names(sab_prop))
names(sab_prop) <- gsub("sab2_", "", names(sab_prop))
names(sab_prop) <- gsub("sab3_", "", names(sab_prop))
names(sab_prop) <- gsub("sab4_", "", names(sab_prop))
names(sab_prop)
sab_cols <- c("navy", "forestgreen", "goldenrod", "grey80")
sab_names <-  c("Female", "Male", "Don't know/Not sure",
                "Declined to Answer")
names(sab_cols) <- sab_names

sab_plot <- 
  sab_prop %>% 
  pivot_longer(cols = contains("."),
               names_to = c(".value", "SAB"),
               names_pattern = "(.*)\\.(.*)") %>% 
  mutate(gender = case_when(grepl("1_", gender) ~ "Transgender woman",
                            grepl("2_", gender) ~ "Transgender man",
                            grepl("3_", gender) ~ "NB/GNC",
                            grepl("4_", gender) ~ "Cisgender",
                            grepl("5_", gender) ~ "Don't know/Not sure",
                            grepl("6_", gender) ~ "Declined to Answer",
                            TRUE ~ NA),
         SAB = case_when(SAB == "male" ~ "Male",
                         SAB == "female" ~ "Female",
                         SAB == "dkns" ~ "Don't know/Not sure",
                         SAB == "ref" ~ "Declined to Answer",
                         TRUE ~ NA),
         mean = ifelse(is.na(mean), 0, mean),
         se = ifelse(is.na(se), 0, se)) %>% 
  filter(gender != "Cisgender" & !is.na(SAB)) %>% 
  mutate(gender = factor(gender, levels = c("Transgender woman", "Transgender man",
                                            "NB/GNC", "Don't know/Not sure",
                                            "Declined to Answer")),
         SAB = factor(SAB, levels = sab_names)) %>% 
  ggplot() +
  geom_bar(aes(x = period, y = mean, fill = SAB), stat = "identity",
           position = "stack") +
  scale_fill_manual(name = "Responses",
                    values = sab_cols) +
  xlab("Period") +
  ylab("Proportion") +
  ggtitle("B) Sex Assigned at Birth") +
  facet_wrap(~gender, ncol = 5) +
  theme_classic() +
  theme(legend.position = "bottom")

sex_sab_plot <- grid.arrange(sex_stream_plot, sab_plot, nrow = 2)

ggsave(filename = "plots/sex-so-gi/Fig5_sex_stream_sab_by_gi.png", plot = sex_sab_plot,
       width = 8, height = 9, units = "in")

## SO (Figure 6) ####
names(so_prop)[3:7] <- paste0("mean.", c("straight", "lesgay", "bi",
                                         "dko", "ref"))
names(so_prop) <- gsub("so_new1_", "", names(so_prop))
names(so_prop) <- gsub("so_new2_", "", names(so_prop))
names(so_prop) <- gsub("so_new3_", "", names(so_prop))
names(so_prop) <- gsub("so_new4_", "", names(so_prop))
names(so_prop) <- gsub("so_new5_", "", names(so_prop))
names(so_prop)
so_names <- c("Straight", "Lesbian/Gay", "Bisexual",
              "Don't know/Not sure",
              "Declined to Answer")
so_cols <- c("navy", "forestgreen", "goldenrod",
             "firebrick", "grey80")
names(so_cols) <- so_names

so_plot <- so_prop %>% 
  pivot_longer(cols = contains("."),
               names_to = c(".value", "SO"),
               names_pattern = "(.*)\\.(.*)") %>% 
  mutate(gender = case_when(grepl("1_", gender) ~ "Transgender woman",
                            grepl("2_", gender) ~ "Transgender man",
                            grepl("3_", gender) ~ "NB/GNC",
                            grepl("4_", gender) ~ "Cisgender",
                            grepl("5_", gender) ~ "Don't know/Not sure",
                            grepl("6_", gender) ~ "Declined to Answer",
                            TRUE ~ NA),
         SO = case_when(SO == "straight" ~ "Straight",
                        SO == "lesgay" ~ "Lesbian/Gay",
                        SO == "bi" ~ "Bisexual",
                        SO == "dko" ~ "Don't know/Not sure",
                        SO == "ref" ~ "Declined to Answer",
                        TRUE ~ NA),
         mean = ifelse(is.na(mean), 0, mean),
         se = ifelse(is.na(se), 0, se)) %>% 
  filter(!is.na(gender) & gender != "Cisgender") %>% 
  mutate(gender = factor(gender, levels = c("Transgender woman", "Transgender man",
                                            "NB/GNC", "Don't know/Not sure",
                                            "Declined to Answer")),
         SO = factor(SO, levels = so_names)) %>% 
  ggplot() +
  geom_stream(aes(x = period, y = mean, fill = SO), type = "proportion") +
  scale_fill_manual(name = "Sexual Orientation",
                    values = so_cols) +
  facet_wrap(~gender, ncol = 5) +
  xlab("Period") +
  ylab("Proportion") +
  ggtitle("Sexual Orientation by Gender Identity") +
  theme_classic()

ggsave(filename = "plots/sex-so-gi/Fig6_so_by_gi.png", plot = so_plot,
       width = 10, height = 8, units = "in")

# Hypothesis tests ####
## Prep data ####
### SO ####
so_to_test<- so_prop %>% 
  pivot_longer(cols = contains("."),
               names_to = c(".value", "SO"),
               names_pattern = "(.*)\\.(.*)") %>% 
  mutate(gender = case_when(grepl("1_", gender) ~ "Transgender woman",
                            grepl("2_", gender) ~ "Transgender man",
                            grepl("3_", gender) ~ "NB/GNC",
                            grepl("4_", gender) ~ "Cisgender",
                            grepl("5_", gender) ~ "Don't know/Not sure",
                            grepl("6_", gender) ~ "Declined to Answer",
                            TRUE ~ NA),
         SO = case_when(SO == "straight" ~ "Straight",
                        SO == "lesgay" ~ "Lesbian/Gay",
                        SO == "bi" ~ "Bisexual",
                        SO == "dko" ~ "Don't know/Not sure",
                        SO == "ref" ~ "Declined to Answer",
                        TRUE ~ NA),
         mean = ifelse(is.na(mean), 0, mean),
         se = ifelse(is.na(se), 0, se)) %>% 
  filter(!is.na(gender) & gender != "Cisgender") %>% 
  mutate(gender = factor(gender, levels = c("Transgender woman", "Transgender man",
                                            "NB/GNC", "Don't know/Not sure",
                                            "Declined to Answer")),
         SO = factor(SO, levels = so_names),
         mean = mean*100) 

### Sex ####
sex_to_test <- sex_prop %>% 
  pivot_longer(cols = contains("."),
               names_to = c(".value", "sex"),
               names_pattern = "(.*)\\.(.*)") %>% 
  mutate(gender = case_when(grepl("1_", gender) ~ "Transgender woman",
                            grepl("2_", gender) ~ "Transgender man",
                            grepl("3_", gender) ~ "NB/GNC",
                            grepl("4_", gender) ~ "Cisgender",
                            grepl("5_", gender) ~ "Don't know/Not sure",
                            grepl("6_", gender) ~ "Declined to Answer")) %>% 
  filter(gender != "Cisgender") %>% 
  mutate(gender = factor(gender, levels = gend_names),
         mean = mean*100)

## multiCA.test ####
for(gend in gend_names){
  sink(file= paste0("tables/so-gi/multiCA_so_",
                    gsub("/", "", gend), ".txt"))
  out <- multiCA.test(SO~period, weights = mean,
               data = so_to_test %>% 
                 filter(gender == gend))
  print(out)
  sink()
  
  sink(file= paste0("tables/sex-gi/multiCA_sex_",
                    gsub("/", "", gend), ".txt"))
  out <- multiCA.test(sex~period, weights = mean,
               data = sex_to_test %>% 
                 filter(gender == gend))
  print(out)
  sink()
}