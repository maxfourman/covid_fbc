#04_imputation.R

###################################################################### 

## Code author: Max Fourman

## Description: 
# This code does multiple imputation by chained equations for cleaned 
# topline clinical data

###################################################################### 

library(tidyverse)
install.packages("longitudinalData", repos='https://cran.ma.imperial.ac.uk/') 
library(longitudinalData)

###################################################################### 

#load the prepped data from 01_analysis
load("/home/u034/mfourman/ccp_blood_data_long.rda")



##add a factor for time from symptom onset
'''ccp_blood_data_long = ccp_blood_data_long %>% 
  mutate(
    onset_to_lab_case = case_when(
      onset_to_lab <= 3 ~ "0-3",
      onset_to_lab <= 6 ~ "3-6",
      onset_to_lab <= 9 ~ "6-9",
      onset_to_lab <= 12 ~ "9-12",
      onset_to_lab <= 15 ~ "12-15",
      TRUE ~ "15-30"
    ) %>% 
      factor(levels = c("0-3", "3-6", "6-9", "9-12", "12-15", "15-30"))
  )
'''

###################################################################### 
#look at the distribution of the data
table(ccp_blood_data_long$admission_to_lab)
table(ccp_blood_data_long$onset_to_lab)

#check the distribution of the admission to lab data
ccp_blood_data_long %>% 
  filter(onset_to_lab <= 14) %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  ggplot(aes(onset_to_lab)) + 
  geom_histogram()

### check the number of entries per person
ccp_blood_data_long %>% 
  filter(onset_to_lab <= 14) %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  head()

#take 200 rows for a smaller table to test on
test_data <- ccp_blood_data_long %>% 
  filter(onset_to_lab <= 14) %>% 
  head(2000)
test_data


wide_test <- test_data %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  pivot_wider(names_from = c(onset_to_lab), 
              values_from = c(value))
wide_test
colnames(wide_test)
col_order <- c("subjid", "cestdat", "daily_lbdat", "admission_to_lab",
               "severity2", "dsterm", "dsstdtc", "assay",
               "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")
wide_test <- wide_test[, col_order]
wide_test <- as.data.frame(wide_test)
wide_test <- wide_test %>%  select(1, 9:23) 
wide_test %>%  head()

l2 <- longData(wide_test)
long_data_input2           


#make a line plot of each parameter with time from symptom
test_data %>% 
  filter(assay %in% c("daily_wbc_lborres", "daily_neutro_lborres", "daily_lymp_lborres")) %>% 
  ggplot(aes(x=onset_to_lab, y=value, group=subjid)) +
  geom_line()+
  facet_wrap(vars(assay), scales = "free_y") +
  labs(x = "Onset to sample (days)", y = "Value") + 
  theme(legend.position = "top")


ccp_blood_data_long %>% 
  select(subjid, onset_to_lab_case, assay, value) %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  filter(onset_to_lab_case %in% c("0-3", "3-6", "6-9", "9-12")) %>% 
  head()

long_test <- ccp_blood_data_long %>% 
  head(1000)

ccp_blood_data_long %>%
  filter(assay %in% c("daily_wbc_lborres", "daily_neutro_lborres", "daily_lymp_lborres")) %>% 
  ggplot(aes(x=onset_to_lab, y=value, group=subjid)) +
  geom_line(alpha = 0.2)+
  facet_wrap(vars(assay), scales = "free_y") +
  labs(x = "Onset to sample (days)", y = "Value") + 
  theme(legend.position = "top")



wide_test <- test_data %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  pivot_wider(names_from = c(onset_to_lab), 
              values_from = c(value))
              
wide_test


md.pattern(full_wide_test)

full_wide_test %>%  head()

imputed <- mice(full_wide_test,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)
  

sapply(full_wide_test, function(x) sum(is.na(x)))

str(wide_test)


init = mice(wide_test, maxit=0) 
meth = init$method
predM = init$predictorMatrix

predM[, c("subjid", "cestdat", "daily_lbdat", "severity2", "dsterm", "dsstdtc", "days_to_outcome", "assay")]=0


meth[c("7", "9", "12", "15", "0", "5", "11", "16", "19", "22", "24", "2", "4")]="polyreg"

set.seed(103)
imputed = mice(wide_test, method=meth, predictorMatrix=predM, m=5)


###################################################################### 
#Example of using MICE
dat <- read.csv(url("https://goo.gl/4DYzru"), header=TRUE, sep=",")
sapply(dat, function(x) sum(is.na(x)))
original <- dat

dat %>%  head()

set.seed(10)
dat[sample(1:nrow(dat), 20), "Cholesterol"] <- NA
dat[sample(1:nrow(dat), 20), "Smoking"] <- NA
dat[sample(1:nrow(dat), 20), "Education"] <- NA
dat[sample(1:nrow(dat), 5), "Age"] <- NA
dat[sample(1:nrow(dat), 5), "BMI"] <- NA

dat <- dat %>%
  mutate(
    Smoking = as.factor(Smoking),
    Education = as.factor(Education),
    Cholesterol = as.numeric(Cholesterol)
  )
str(dat)


library(mice)
init = mice(dat, maxit=0) 
meth = init$method
predM = init$predictorMatrix

predM[, c("BMI")]=0
meth[c("Age")]=""

meth[c("Cholesterol")]="norm" 
meth[c("Smoking")]="logreg" 
meth[c("Education")]="polyreg"


set.seed(103)
imputed = mice(dat, method=meth, predictorMatrix=predM, m=5)

imputed <- complete(imputed)

sapply(imputed, function(x) sum(is.na(x)))
