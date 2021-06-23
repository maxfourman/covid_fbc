#04_exploratory graphs
###################################################################### 

## Code author: Max Fourman

## Description: 
#

###################################################################### 
## Load packages
# May be already done if using previous scripts

install.packages("moonBook", repos='https://cran.ma.imperial.ac.uk/')
library("remoter") 
remoter::client(port=60117)
library("tidyverse")
library("finalfit")
library("pROC")
library("gridExtra")


###################################################################### 
#load the prepped data from 01_analysis
load("/home/u034/mfourman/ccp_blood_data_long.rda")
load("/home/u034/mfourman/ccp_data_clean.rda")



# PLots 3 x scatter plots to look for outliers
ccp_blood_data_long %>% 
  select(onset_to_lab, assay, value) %>% 
  filter(assay %in% c("daily_wbc_lborres", "daily_neutro_lborres", "daily_lymp_lborres")) %>% 
  ggplot(aes(x=onset_to_lab, y=value)) +
  geom_point()+
  facet_wrap(vars(assay), scales = "free_y") +
  labs(x = "Onset to sample (days)", y = "Value") + 
  theme(legend.position = "top")

# PLot a line plot 
ccp_blood_data_long %>%
  filter(assay %in% c("daily_wbc_lborres", "daily_neutro_lborres", "daily_lymp_lborres")) %>% 
  ggplot(aes(x=onset_to_lab, y=value, group=subjid)) +
  geom_line()+
  facet_wrap(vars(assay), scales = "free_y") +
  labs(x = "Onset to sample (days)", y = "Value") + 
  theme(legend.position = "top")
  
#add a factor for time from symptom onset
ccp_blood_data_long = ccp_blood_data_long %>% 
  mutate(
    onset_to_lab_case = case_when(
      onset_to_lab <= 5 ~ "0-5",
      onset_to_lab <= 10 ~ "5-10",
      onset_to_lab <= 15 ~ "10-15",
      onset_to_lab <= 20 ~ "15-20",
      onset_to_lab <= 25 ~ "20-25",
      TRUE ~ "25-30"
    ) %>% 
      factor(levels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30"))
  )






#Stratify by severity and plot time vs lab measurement

#first clean labels and factors 
ccp_blood_data_long <- ccp_blood_data_long %>% 
  mutate(assay =            
           assay %>%                
           factor() %>% 
           fct_recode(
             "Lymphocytes" = "daily_lymp_lborres",                 
             "Neutrophils"   = "daily_neutro_lborres",
             "NLR" = "daily_nlr_lborres",
             "Platelets" = "daily_plt_lborres",
             "WBC" = "daily_wbc_lborres"))  

ccp_blood_data_long %>% 
  filter(assay %in% c("Lymphocytes")) %>% 
  distinct(subjid, onset_to_lab, value)

# Now plot a faceted grid plot of lymphocyte
ccp_blood_data_long %>% 
  select(subjid, onset_to_lab, value, assay, severity2) %>% 
  filter(assay %in% c("Lymphocytes")) %>% 
  drop_na() %>% 
  ggplot(aes(onset_to_lab, value, colour = severity2)) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(vars(assay, severity2), scales = "fixed", ncol =5) +
  geom_hline(yintercept=1.5, linetype="dashed", color = "gray20") +
  scale_colour_manual(values = c("#3690c0", "#023858", "#c51b7d", "#ff7f00", "#e31a1c")) +
  labs(x = "Onset to sample (days)", 
       y = (expression("Lymphocyte Count (x10"^9*"/L)")), 
       colour = "Severity") + 
  theme(legend.position = "top") 

# Now plot a faceted grid plot of neutrophil
ccp_blood_data_long %>% 
  select(subjid, onset_to_lab, value, assay, severity2) %>% 
  filter(assay %in% c("Neutrophils")) %>% 
  drop_na() %>% 
  ggplot(aes(onset_to_lab, value, colour = severity2)) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(vars(assay, severity2), scales = "fixed", ncol =5) +
  geom_hline(yintercept=7.5, linetype="dashed", color = "gray20") +
  scale_colour_manual(values = c("#3690c0", "#023858", "#c51b7d", "#ff7f00", "#e31a1c")) +
  labs(x = "Onset to sample (days)", 
       y = (expression("Neutrophil Count (x10"^9*"/L)")),
       colour = "Severity") + 
  theme(legend.position = "top") 

# Now plot a faceted grid plot of NLR
ccp_blood_data_long %>% 
  select(subjid, onset_to_lab, value, assay, severity2) %>% 
  filter(assay %in% c("NLR")) %>% 
  drop_na() %>% 
  ggplot(aes(onset_to_lab, value, colour = severity2)) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(vars(assay, severity2), scales = "fixed", ncol =5) +
  scale_colour_manual(values = c("#3690c0", "#023858", "#c51b7d", "#ff7f00", "#e31a1c")) +
  labs(x = "Onset to sample (days)", 
       y = (expression("Neutrophil:Lymphocyte Ratio")),
       colour = "Severity") + 
  theme(legend.position = "top") 











  
  






