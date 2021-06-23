#01b_clean_ccp_TIER2_data

###################################################################### 

## Code author: Max Fourman
## Lines 30-145 adapted from https://github.com/SurgicalInformatics/cocin_ccp/blob/master/03_prep.R

## Description: 
# This code sources the joined ccp file from #00_source_all.R
# It then filters outliers, fills missing data 
# Finally a smaller subset of the data is saved for quicker analysis
# This final table is pivoted 'long' for easier graphing

###################################################################### 
## Load packages


library("remoter") 
remoter::client(port=60117)
library("tidyverse")
library("finalfit")
library("kableExtra")


###################################################################### 
# Load data from 00_source_all.R

load("/home/u034/mfourman/tier2_ccp_data_final.rda")

###################################################################### 
#DATA PREPARATION AND CLEANING
# This section of code used from 
# https://github.com/SurgicalInformatics/cocin_ccp/blob/master/03_prep.R

tier2_ccp_data = tier2_ccp_data %>% 
  mutate(
    ## Units for continuous variables fixed here
    # Potassium has issues, this deals with some. 
    daily_potassium_lborres = case_when(
      daily_potassium_lborres > 100 ~ NA_real_,
      daily_potassium_lborres > 12 ~ daily_potassium_lborres / 10,
      daily_potassium_lborres < 0 ~ abs(daily_potassium_lborres),
      TRUE ~ daily_potassium_lborres),
    # Hb
    # Ignore the units variable as people have unfortunately just got it wrong :(
    daily_hb_lborres = ifelse(daily_hb_lborres < 25, daily_hb_lborres * 10, daily_hb_lborres),
    # WBC
    # Units do not help here either. 
    # We couldn't be sure that those with WBC>100 were definitely factor of 10 wrong or leukaemia. 
    # Spent some time matching up lymph and neut counts, but thought easiest to exclude
    daily_wbc_lborres = ifelse(daily_wbc_lborres > 100, NA_real_, daily_wbc_lborres),
    # Neutrophils
    daily_neutro_lborres = ifelse(daily_neutro_lborres > 100, daily_neutro_lborres / 1000, 
                                  daily_neutro_lborres),
    # Lymphocytes needs looking at: most are x10^9, but some are x10^6
    daily_lymp_lborres = ifelse(daily_lymp_lborres > 100, daily_lymp_lborres / 1000, 
                                daily_lymp_lborres),
    # Platelets
    # Units don't help
    # Few very high have been left in as couldn't be sure not sepsis-related thrombocytosis
    # PT/INR
    # Some PTs are actually INRs
    daily_pt_lborres = ifelse(daily_pt_lborres < 9, daily_pt_lborres * 12, daily_pt_lborres),
    # No good way combining INR so consider using a threshold for abnormal
    # See below for new variable based on INR 1.0 = PT 12.0
    # Bilirubin
    # Don't use daily_bil_lborresu variable as looks incorrect for all daily_bil_lborres values
    # Urea
    # Urea values in different units can't be differentiated by eye. 
    # We changed urea values on basis of units, but wonder if the mg/dL unit has been used incorrectly
    # No pattern across hospitals for units being different. 
    # Decided to leave original values unchanged. 
    # daily_bun_lborres = ifelse(daily_bun_lborresu == "mg/dL", daily_bun_lborres * 2.8, 
    #                            daily_bun_lborres),
    # Creatinine
    # Ignore units variable, mg/dL values are not in the expected range for this unit. 
    # Some high values are left in on the presumption they may be correct. 
    # Glucose
    daily_glucose_lborres = ifelse(daily_glucose_lborres > 100, NA_real_, daily_glucose_lborres),
    # pO2
    daily_pao2_lborres = ifelse(daily_pao2_lborresu == "mmHg", daily_pao2_lborres / 7.5 , 
                                daily_pao2_lborres),
    # Lactate
    # Some very high numbers removed. Patient at 47 died and left in, presumed real.
    daily_lactate_lborres = ifelse(daily_lactate_lborres > 100, NA_real_, daily_lactate_lborres),
    # FiO2 - updated 17/08/2020
    # This may need looked at by hand. L/min have been included by the look of it. 
    daily_fio2_lborres = case_when(
      daily_fio2_lborres <= 1 ~ daily_fio2_lborres, # Presume FiO2
      daily_fio2_lborres <= 2 ~ 0.24,               # Presume these are all L/min
      daily_fio2_lborres <= 3 ~ 0.28,
      daily_fio2_lborres <= 4 ~ 0.32,
      daily_fio2_lborres <= 5 ~ 0.36,
      daily_fio2_lborres <= 6 ~ 0.40,
      daily_fio2_lborres <= 10 ~ 0.50,
      daily_fio2_lborres <= 15 ~ 0.70,
      daily_fio2_lborres > 15 ~ daily_fio2_lborres / 100, # Presume % rather than fraction
      TRUE ~ daily_fio2_lborres),
    # These should all be FiO2%
    daily_fio2b_lborres = case_when(
      daily_fio2b_lborres == 0 ~ 0,
      daily_fio2b_lborres <= 2 ~ 24,               # Presume everything <=15 is actually L/min
      daily_fio2b_lborres <= 3 ~ 28,
      daily_fio2b_lborres <= 4 ~ 32,
      daily_fio2b_lborres <= 5 ~ 36,
      daily_fio2b_lborres <= 6 ~ 40,
      daily_fio2b_lborres <= 10 ~ 50,
      daily_fio2b_lborres <= 15 ~ 70,
      TRUE ~ daily_fio2b_lborres),                   # Presume FiO2%
    # These should all be L/min
    daily_fio2c_lborres_converted = case_when(
      daily_fio2c_lborres == 0 ~ 0,
      daily_fio2c_lborres <= 2 ~ 0.24,               # Presume these are all L/min
      daily_fio2c_lborres <= 3 ~ 0.28,
      daily_fio2c_lborres <= 4 ~ 0.32,
      daily_fio2c_lborres <= 5 ~ 0.36,
      daily_fio2c_lborres <= 6 ~ 0.40,
      daily_fio2c_lborres <= 10 ~ 0.50,
      daily_fio2c_lborres <= 15 ~ 0.70,
      TRUE ~ daily_fio2c_lborres / 100),             # Presume FiO2%
    daily_fio2_combined = case_when(
      !is.na(daily_fio2_lborres) ~ daily_fio2_lborres,
      is.na(daily_fio2_lborres) & !is.na(daily_fio2b_lborres) ~ daily_fio2b_lborres / 100,
      is.na(daily_fio2_lborres) & !is.na(daily_fio2c_lborres) ~ daily_fio2c_lborres_converted),
    # Sa02
    daily_sao2_lborres = case_when(
      daily_sao2_lborres <= 1 ~ daily_sao2_lborres * 100, 
      daily_sao2_lborres <= 5 ~ NA_real_,
      daily_sao2_lborres <= 10 ~ daily_sao2_lborres * 10,
      daily_sao2_lborres <= 50 ~ NA_real_,
      daily_sao2_lborres > 500 ~ daily_sao2_lborres / 10,
      TRUE ~ daily_sao2_lborres), 
    ## Checkbox recodes here ---------------------------------------------
    ethnicity = case_when(
      ethnic___1 == "Checked" ~ "Arab",
      ethnic___2 == "Checked" ~ 	"Black",
      ethnic___3 == "Checked" ~ 	"East Asian",
      ethnic___4 == "Checked" ~ 	"South Asian",
      ethnic___5 == "Checked" ~ 	"West Asian",
      ethnic___6 == "Checked" ~  "Latin American",
      ethnic___7 == "Checked" ~ 	"White",
      ethnic___8 == "Checked" ~ 	"Aboriginal/First Nations",
      ethnic___9 == "Checked" ~ 	"Other"
    ) %>% 
      factor() %>% 
      ff_label("Ethnicity")
  )
###################################################################### 

## Filter outlier lymphocyte counts


# Find the range of values of lymphocytes and neutorphils
tier2_ccp_data %>% 
  drop_na(daily_lymp_lborres) %>% 
  select(daily_lymp_lborres) %>% 
  range()
tier2_ccp_data %>% 
  drop_na(daily_neutro_lborres) %>% 
  select(daily_neutro_lborres) %>% 
  range()

####################################
# Fill missing Age variable based on the topline age 
tier2_ccp_data <- tier2_ccp_data %>% 
  group_by(subjid) %>%
  fill(age) %>% 
  ungroup()
# filter out unreliable ages and children
tier2_ccp_data <- tier2_ccp_data %>% 
  filter(age <=110) %>% 
  filter(age >=18) 


####################################
# check how many outlier lymphocyte counts there are
tier2_ccp_data %>% 
  drop_na(daily_lymp_lborres) %>% 
  mutate(
    neg_daily_lymp_lborres = case_when(
      daily_lymp_lborres < 0 ~ "Negative",
      daily_lymp_lborres > 20 ~ "Outlier",
      TRUE ~ "OK"
    ) %>% 
      factor(levels = c("Negative", "OK", "Outlier"))
  ) %>%
  group_by(neg_daily_lymp_lborres) %>% 
  summarize(count = n())

##Filter out the outliers
tier2_ccp_data <- tier2_ccp_data %>% 
  filter(daily_lymp_lborres >= 0) 

####################################
## Filter outlier neutrophil counts
# check how many outlier neutrophil counts there are
tier2_ccp_data %>% 
  drop_na(daily_neutro_lborres) %>% 
  mutate(
    neg_daily_neutro_lborres = case_when(
      daily_neutro_lborres < 0 ~ "Negative",
      daily_neutro_lborres > 50 ~ "Outlier",
      TRUE ~ "OK"
    ) %>% 
      factor(levels = c("Negative", "OK", "Outlier"))
  ) %>%
  group_by(neg_daily_neutro_lborres) %>% 
  summarize(count = n())

##Filter out the outliers neutrophil
tier2_ccp_data <- tier2_ccp_data %>% 
  filter(daily_neutro_lborres >= 0) 

####################################
# check how many outlier wbc counts there are
tier2_ccp_data %>% 
  drop_na(daily_wbc_lborres) %>% 
  mutate(
    neg_daily_wbc_lborres = case_when(
      daily_wbc_lborres < 0 ~ "Negative",
      daily_wbc_lborres > 70 ~ "Outlier",
      TRUE ~ "OK"
    ) %>% 
      factor(levels = c("Negative", "OK", "Outlier"))
  ) %>%
  group_by(neg_daily_wbc_lborres) %>% 
  summarize(count = n())

##Filter out the outliers wbc count
tier2_ccp_data <- tier2_ccp_data %>% 
  filter(daily_wbc_lborres >= 0) 

# find out how many of the daily_lbdat vs daily_dsstdat dont match
#tier2_ccp_data %>% 
#  filter(daily_lbdat != daily_dsstdat) %>% 
#  select(daily_dsstdat, daily_lbdat, cestdat) %>% 
#  head(100)

# Add in the missing dates from daily_dsstdat if missing from daily_lbdat
tier2_ccp_data <- tier2_ccp_data %>%
  mutate(daily_lbdat = coalesce(daily_lbdat, daily_dsstdat)) 


#drop rows we dont need 
#tier2_ccp_data <- tier2_ccp_data %>%
#  drop_na(daily_lbdat) 



##########################################################################################
#DATA PREPARATION AND CLEANING
## DATES + Drop erroneous dates

#Sort the dates into the correct format
tier2_ccp_data <- tier2_ccp_data %>%
  mutate(cestdat = lubridate::ymd(cestdat),
         daily_lbdat = lubridate::ymd(daily_lbdat),
         dsstdtc = lubridate::ymd(dsstdtc),
         hostdat = lubridate::ymd(hostdat))


# get rid weird crazy dates from 1918
tier2_ccp_data <- tier2_ccp_data %>% 
  filter(daily_lbdat >= "2020-01-01") %>% 
  filter(daily_lbdat <= "2022-01-01") %>% 
  filter(cestdat >= "2020-01-01") %>% 
  filter(cestdat <= "2022-01-01") %>% 
  filter(hostdat >= "2020-01-01") %>% 
  filter(hostdat <= "2022-01-01")

#check that the max and min dates are sensible
tier2_ccp_data %>%
  summarise(
    max_symptom_date = max(cestdat),
    min_symptom_date = min(cestdat),
    max_admission_date = max(hostdat),
    min_admission_date = min(hostdat))

#random graph to plot the situation with the admission dates
tier2_ccp_data %>% 
  filter(hostdat <=  2020-03-01 )
ggplot( aes(hostdat)) + 
  geom_histogram()

# take the difference between the two dates to work out time to symptom onset
tier2_ccp_data <- tier2_ccp_data %>% 
  mutate(onset_to_lab = as.numeric(daily_lbdat - cestdat)) %>% 
  mutate(admission_to_lab = as.numeric(daily_lbdat - hostdat))

#NOT USED calculate the days from lab to outcome (death etc)
#tier2_ccp_data_clean <- tier2_ccp_data_clean %>% 
#  mutate(days_to_outcome = as.numeric(daily_lbdat - dsstdtc)) 


# check for negative dates
tier2_ccp_data %>% 
  drop_na(onset_to_lab) %>% 
  mutate(
    neg_onset_to_lab = case_when(
      onset_to_lab < 0 ~ "Negative",
      onset_to_lab >= 0 ~ "Positive"
    ) %>% 
      factor(levels = c("Negative", "Positive"))
  ) %>%
  group_by(neg_onset_to_lab) %>% 
  summarize(count = n())

# Filter out the negative dates + anyone with hopsital stay over 30 days
tier2_ccp_data <- tier2_ccp_data %>% 
  filter(onset_to_lab >= 0)
#  filter(onset_to_lab <= 30)

### Redundant code
#tier2_ccp_data_clean <- tier2_ccp_data_clean %>% 
#  filter(admission_to_lab >= 0) %>% 
#  filter(admission_to_lab <= 30)

#add NLR column
tier2_ccp_data <- tier2_ccp_data %>% 
  mutate(daily_nlr_lborres = daily_neutro_lborres/daily_lymp_lborres)
#check the NLR column makes sense
tier2_ccp_data %>% 
  select(daily_nlr_lborres, daily_neutro_lborres, daily_lymp_lborres) %>% 
  head()

tier2_ccp_data %>% 
  distinct(subjid) %>% 
  nrow()



#add in the simplifed outcome scale
tier2_ccp_data <- tier2_ccp_data %>%
  mutate(
    outcome = fct_collapse(
      severity2,
      "Death/IMV" = c("Death", "IMV"),
      "Surived" = c("NIV/HFNC", "Oxygen alone", "Ward"))
  )

tier2_ccp_data_clean <- tier2_ccp_data
###################################################################### 
### Save the cleaned tier2_ccp_data file
save(tier2_ccp_data_clean, file = 
       "/home/u034/mfourman/tier2_ccp_data_clean.rda")

