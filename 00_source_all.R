#01_source_all.R

###################################################################### 

## Code author: Max Fourman

## Description: 
# This code sources the raw csv files from ULTRA for ccp_data, topline, treatment & outcome
# It then joins select columns and calculates the severity index
# The cleaned and prepped file is then saved ready for analysis

###################################################################### 

# Load packages and read csv files

library(remoter)
remoter::client(port=60117)
library("tidyverse")
library("finalfit")

ccp_data <- read.csv("/home/u034/shared/data/wp-4/clinical-data/002_crf_20200904/ccp_data.csv")
topline <- read.csv("/home/u034/shared/data/wp-4/clinical-data/002_crf_20200904/topline.csv")
treatment <- read.csv("/home/u034/shared/data/wp-4/clinical-data/002_crf_20200904/treatment.csv")
outcome <- read.csv("/home/u034/shared/data/wp-4/clinical-data/002_crf_20200904/outcome.csv")

###################################################################### 

#TIER 2 data for cytokines
tier2_topline <- read.csv("/home/u034/shared/data/wp-4/clinical-data/007_crf_20210419/topline.csv")
tier2_topline %>% colnames

tier2_ccp_data <- read.csv("/home/u034/shared/data/wp-4/clinical-data/007_crf_20210419/ccp_data.csv")
tier2_ccp_data %>% colnames


###################################################################### 
### Joins on tier zero data

# join the death variable from topline
ccp_data = ccp_data %>%
  left_join(
    topline %>% select(subjid, death), 
    by = c("subjid"))
tier2_ccp_data = tier2_ccp_data %>%
  left_join(
    topline %>% select(subjid, death), 
    by = c("subjid"))

# join the any_noninvasive, any_oxygen variable from treatment
ccp_data = ccp_data %>% 
  left_join(treatment %>% select(subjid, any_noninvasive, any_oxygen), 
            by = c("subjid"))

###################################################################### 
### ADd severity scale

# Add the severity scale
ccp_data = ccp_data %>% 
  mutate(
    severity2 = case_when(
      death == "Yes" ~ "Death",
      any_invasive == "Yes" ~ "IMV",
      any_noninvasive == "Yes" ~ "NIV/HFNC",
      daily_nasaloxy_cmtrt == "Yes" ~ "NIV/HFNC",
      any_oxygen == "Yes" ~ "Oxygen alone",
      TRUE ~ "Ward"
    ) %>% 
      factor(levels = c("Death", "IMV", "NIV/HFNC", "Oxygen alone", "Ward"))
  )
# Add the severity scale to tier2 data
tier2_ccp_data = tier2_ccp_data %>% 
  mutate(
    severity2 = case_when(
      death == "Yes" ~ "Death",
      any_invasive == "Yes" ~ "IMV",
      any_noninvasive == "Yes" ~ "NIV/HFNC",
      daily_nasaloxy_cmtrt == "Yes" ~ "NIV/HFNC",
      any_oxygen == "Yes" ~ "Oxygen alone",
      TRUE ~ "Ward"
    ) %>% 
      factor(levels = c("Death", "IMV", "NIV/HFNC", "Oxygen alone", "Ward"))
  )

#this is a check step
ccp_data %>% 
  select(subjid, death, any_invasive, any_noninvasive, daily_nasaloxy_cmtrt, any_oxygen, severity2)  %>% 
  head (20)

### fill the missing cestdat values 
ccp_data <- ccp_data %>% 
  group_by(subjid) %>%
  fill(cestdat) %>% 
  ungroup()
### fill the missing cestdat values 
tier2_ccp_data <- tier2_ccp_data %>% 
  group_by(subjid) %>%
  fill(cestdat) %>% 
  ungroup()

###################################################################### 
# Save
save(ccp_data, file = 
       "/home/u034/mfourman/ccp_data_final.rda")
save(tier2_ccp_data, file = 
       "/home/u034/mfourman/tier2_ccp_data_final.rda")

###################################################################### 
