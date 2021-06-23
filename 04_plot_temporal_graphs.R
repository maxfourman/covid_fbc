#004_Plot_temporal_graphs
###################################################################### 

## Code author: Max Fourman

## Description: 
# This code creates the time series graphs with the blood data

###################################################################### 
## Load packages
# May be already done if using previous scripts

install.packages("moonBook", repos='https://cran.ma.imperial.ac.uk/')
library("remoter") 
remoter::client(port=60117)
library(plyr); library(dplyr)
library("tidyverse")
library("finalfit")
library("pROC")
library("ggpubr")
library("ggsci")
###################################################################### 

#load the prepped data from 01_analysis

load("/home/u034/mfourman/ccp_blood_data_long.rda")
###################################################################### 


##to start with work with lymphocyte data
lymphocyte_data <- ccp_blood_data_long %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  distinct()


lymphocyte_data %>% 
  select(!c(dsterm, dsstdtc)) %>% 
  head()

lymphocyte_data %>% 
  distinct(subjid) %>% 
  nrow()

#get the mean of each #requires summarySE function
lymph_plot_data <-
  summarySE(lymphocyte_data, 
            measurevar="value", 
            groupvars=c("outcome", "onset_to_lab"))
lymph_plot_data %>% head(15)



lymph_plot_data %>% 
  filter(onset_to_lab <= 14) %>% 
  ggplot( aes(x=onset_to_lab, y=value, group=outcome, color=outcome)) +
  geom_ribbon(aes(ymin=value-se, ymax=value+se), alpha=0.1, linetype="dotted") +
  geom_line() +
  geom_point() +
  ggtitle("Mean lymphocyte values stratified by outcome") +
  ylab("Mean Lymphocyte count") +
  xlab("Days from symptom onset to lab measurement") +
  scale_color_npg()




#### now some neutrophils
neutrophil_data <- ccp_blood_data_long %>% 
  filter(assay %in% c("daily_neutro_lborres")) %>% 
  drop_na(value) %>%  
  distinct()

neutrophil_data %>% 
  distinct(subjid) %>% 
  nrow()

#get the mean of each #requires summarySE function
neutro_plot_data <-
  summarySE(neutrophil_data, 
            measurevar="value", 
            groupvars=c("outcome", "onset_to_lab"))
neutro_plot_data %>% 
  head()

neutro_plot_data %>% 
  filter(onset_to_lab <= 30) %>% 
  ggplot( aes(x=onset_to_lab, y=value, group=outcome, color=outcome)) +
  geom_ribbon(aes(ymin=value-se, ymax=value+se), alpha=0.1, linetype="dotted") +
  geom_line() +
  geom_point() +
  ggtitle("Mean Neutrophil values stratified by severity2") +
  ylab("Mean Neutrophil count") +
  xlab("Days from symptom onset to lab measurement")

#### now WBC data
wbc_data <- ccp_blood_data_long %>% 
  filter(assay %in% c("daily_wbc_lborres")) %>% 
  distinct()

wbc_data %>% 
  distinct(subjid) %>% 
  nrow()

#get the mean of each #requires summarySE function
wbc_plot_data <-
  summarySE(wbc_data, 
            measurevar="value", 
            groupvars=c("outcome", "onset_to_lab"))
wbc_plot_data %>% head()

wbc_plot_data %>% 
  filter(onset_to_lab <= 14) %>% 
  ggplot( aes(x=onset_to_lab, y=value, group=outcome, color=outcome)) +
  geom_ribbon(aes(ymin=value-se, ymax=value+se), alpha=0.1, linetype="dotted") +
  geom_line() +
  geom_point() +
  ggtitle("Mean WBC values stratified by severity2") +
  ylab("Mean WBC count") +
  xlab("Days from symptom onset to lab measurement")


#### now NLR data
nlr_data <- ccp_blood_data_long %>% 
  filter(assay %in% c("daily_nlr_lborres")) %>% 
  filter(value != 0) %>%                            #there were 22 with value = 0 
  filter(value != Inf) %>% 
  distinct()

nlr_data %>% 
  distinct(subjid) %>% 
  nrow()
nlr_data %>% 
  slice_max(value) %>% 
  select(value)


#get the mean of each #requires summarySE function
nlr_plot_data <-
  summarySE(nlr_data, 
            measurevar="value", 
            groupvars=c("outcome", "onset_to_lab"))
nlr_plot_data %>% head()

nlr_plot_data %>% 
  filter(onset_to_lab <= 14) %>% 
  ggplot( aes(x=onset_to_lab, y=value, group=outcome, color=outcome)) +
  geom_ribbon(aes(ymin=value-se, ymax=value+se), alpha=0.1, linetype="dotted") +
  geom_line() +
  geom_point() +
  ggtitle("Mean NLR values stratified by severity2") +
  ylab("Mean NLR value") +
  xlab("Days from symptom onset to lab measurement")


