#02_load_and_clean_lothian_data

###################################################################### 

## Code author: Max Fourman

## Description: 
# This code sources CAP data
# It also loads the ICD 10 codes for diagnosis 

###################################################################### 
## Load packages
# May be already done if using previous scripts


library("remoter") 
remoter::client(port=60117)
library("tidyverse")
library("finalfit")



###################################################################### 


# read in the ICD 10 codes 
icd_codes <- read.csv("/home/u034/mfourman/icd10cm_codes_2021.txt",
                        sep = "\t",
                        header = FALSE)

# split the file into two columns
icd_codes <- icd_codes %>% 
  separate(V1, into = c("code", "icd_description"), 
           sep = "\\s", 
           extra = "merge")

#check the above worked
icd_codes %>% head()

###################################################################### 
#load in the trak data
trak_data <- read.csv("/home/u034/mfourman/white_cells_pneumonia.csv")

# NOT USED - if needed to string separate
#trak_data$MAIN_DISCHARGE_DIAGNOSIS <- trak_data$MAIN_DISCHARGE_DIAGNOSIS %>% 
#  str_sub(0,4) %>% 
#  as.factor()

#Join the ICD code descriptions
trak_data <- trak_data %>% 
  left_join(icd_codes, by = c("MAIN_DISCHARGE_DIAGNOSIS" = "code")) 

trak_data %>% head()

#checks which codes do not have a description
trak_data %>% 
  filter(is.na(trak_data$icd_description)) %>% 
  group_by(MAIN_DISCHARGE_DIAGNOSIS) %>% 
  summarise(n= n()) %>% 
  ungroup() %>% 
  arrange(desc(n))




#Use case_when to create a column with a diagnosis based on coalescing similar codes into sensible categories
trak_data <- trak_data %>% 
  mutate(
    diagnosis = case_when(
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J180", "J181", "J189", "J189", "J158", "J156") ~ "unspec_bacterial_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J100", "J101", "J108", "J111") ~ "Influenza",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J128", "J129") ~ "unspec_viral_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J14X") ~ "Hemophilus_influenzae",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J13X") ~ "Streptococcus_pneumoniae",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J157") ~ "Mycoplasma_pneumoniae",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J151") ~ "Pseudomonas",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J152") ~ "staphylococcus_unspecified",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J121") ~ "RSV_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J123") ~ "metapneumovirus_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J122") ~ "Parainfluenza_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J150") ~ "Klebsiella_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J155") ~ "Ecoli_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J154") ~ "other_strep_pneumonia",
      MAIN_DISCHARGE_DIAGNOSIS %in% c("J120") ~ "adenovirus_pneumonia",) %>% 
    factor(levels = c("unspec_bacterial_pneumonia", "Influenza", "unspec_viral_pneumonia", "Hemophilus_influenzae",
                      "Streptococcus_pneumoniae", "Mycoplasma_pneumoniae", "Pseudomonas",
                      "staphylococcus_unspecified", "RSV_pneumonia", "metapneumovirus_pneumonia",
                      "Parainfluenza_pneumonia", "other_viral_pneumonia", "Klebsiella_pneumonia", 
                      "Ecoli_pneumonia", "other_strep_pneumonia", "adenovirus_pneumonia")))

#Check is any diagnoses are missing and then search the code to manually assign them
trak_data %>% 
  filter(is.na(trak_data$diagnosis)) %>% 
  group_by(MAIN_DISCHARGE_DIAGNOSIS) %>% 
  summarise(n= n()) %>% 
  ungroup() %>% 
  arrange(desc(n))

# check the count of the different diagnosis
trak_data %>% 
  group_by(diagnosis) %>% 
  summarise(n= n()) %>% 
  ungroup() %>% 
  arrange(desc(n))

#Group into bacterial, influenza and non_flu_viral

trak_data <- trak_data %>% 
  mutate(
    group_diagnosis = case_when(
      diagnosis == "Influenza" ~ "Influenza",
      diagnosis %in% c("unspec_viral_pneumonia", "RSV_pneumonia", 
                       "metapneumovirus_pneumonia", "adenovirus_pneumonia") ~ "non_flu_viral", 
      is.na(diagnosis) ~ "Other",
      TRUE ~ "bacterial"))

# check the count of the different diagnosis
trak_data %>% 
  group_by(group_diagnosis) %>% 
  summarise(n= n()) %>% 
  ungroup() %>% 
  arrange(desc(n))



trak_data_clean <- trak_data %>% 
  filter(group_diagnosis != "Other") %>% 
  select(-c(MAIN_DIAGNOSIS, MAIN_ADMISSION_DIAGNOSIS, MAIN_DISCHARGE_DIAGNOSIS, icd_description))

trak_data_clean %>% head()


save(trak_data_clean, file = 
       "/home/u034/mfourman/trak_data_clean.rda")