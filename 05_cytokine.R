#05_cytokine analysis


###################################################################### 

## Code author: Max Fourman


## Description: 
# 

###################################################################### 
## Load packages

library("remoter") 
remoter::client(port=60117)
library("tidyverse")
library("ggpubr")
library("finalfit")
library("corrplot")
library("ggsci")
library("RColorBrewer")


###################################################################### 

#if needed extra biochem data
#biochemistry <- read.csv("/home/u034/shared/data/wp-5/data_linkage/release/20210421_v30/wp4_biochemistry_cleaned_20210421_145917.csv")


# Load tier2_topline 
tier2_topline <- read.csv("/home/u034/shared/data/wp-4/clinical-data/007_crf_20210419/topline.csv")

tier2_topline %>% 
  select(subjid, cestdat, hostdat) %>% 
  missing_glimpse()

#add isaric lasso score
tier2_topline <- tier2_topline %>% 
  mutate_at(vars(chrincard:malnutrition_mhyn, dehydration_vsorres, -diabetes_type_mhyn), fct_recode, 
            NULL = "Unknown",
            "No" = "NO",
            "Yes" = "YES") %>%
  mutate(no_comorbid = select(., chrincard:malnutrition_mhyn, -diabetes_type_mhyn) %>% 
           {. == "Yes"} %>% 
           rowSums(na.rm = TRUE)
  ) %>% 
  isaric_lasso(age = age, sex = sex,  comorbid = no_comorbid, rr = rr_vsorres,
               spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
               bun = daily_bun_lborres, crp = daily_crp_lborres,
               output = c("df_vector"), na_to_zeros = FALSE)

#check missing data
tier2_topline %>% 
  select(age, sex, rr_vsorres, oxy_vsorres, 
         daily_gcs_vsorres, daily_neutro_lborres, daily_crp_lborres) %>% 
  missing_glimpse()

###################################################################### 
#load ccp
load("/home/u034/mfourman/tier2_ccp_data_clean.rda")


###################################################################### 
#load cytokines from old unlinked path
#cytokine_assays <- read_csv("/home/u034/shared/data/wp-6/Plasma cytokines/20210305 ISARIC plasma mediators final.csv")
#first find the patient_id we cant match - OLD file
#missing_patient_id_2 <- cytokine_assays %>% 
#  select("patient_id") %>% 
#  anti_join(
#    tier2_topline %>% 
#      select(subjid, cestdat, hostdat), 
#    by = c("patient_id" = "subjid"))
#missing_patient_id_2 %>% distinct(patient_id) %>% nrow()
# drop box patient_id and barcode add in the cestdat and hostdat
#cytokine_assays <- cytokine_assays %>% 
#  select(!c("Box patient_id", "Barcode")) %>% 
    tier2_topline %>% 
##  left_join(
#      select(subjid, cestdat, hostdat), 
#    by = c("patient_id" = "subjid"))

# Load the data linked plasma cytokine data
plasma_cytokines <- read.csv("/home/u034/shared/data/wp-5/data_linkage/release/20210421_v30/wp6_plasma_cytokines_cleaned_20210421_145917.csv")
plasma_cytokines %>% colnames
plasma_cytokines %>% distinct(patient_id) %>% nrow()
plasma_cytokines %>% missing_glimpse()


#first find the patient_id we cant match - NEW file
missing_patient_id <- plasma_cytokines %>% 
  select("patient_id") %>% 
  anti_join(
    tier2_topline %>% 
      select(subjid, cestdat, hostdat), 
    by = c("patient_id" = "subjid"))
missing_patient_id %>% distinct(patient_id) %>% nrow()

#save list of missing in case anyone needs it
write.csv(missing_patient_id, file = 
            "/home/u034/mfourman/missing_patient_ids.csv")

# join in cestdat and hostdat
plasma_cytokines <- plasma_cytokines %>% 
  select("patient_id", "timepoint", "assay", "value") %>% 
  inner_join(
    tier2_topline %>% 
      select(subjid, cestdat, hostdat), 
    by = c("patient_id" = "subjid"))

#check the number of patients remaining
plasma_cytokines %>% distinct(patient_id) %>% nrow
plasma_cytokines %>% head

#filter out only those with a visit day
plasma_cytokines <- plasma_cytokines %>% 
  filter(grepl('Day', timepoint))

# remove the "DAY" to get an integer
plasma_cytokines <- plasma_cytokines %>% 
  mutate(timepoint = str_sub(timepoint, -1))

#cmake factor then get rid of convalescent
plasma_cytokines <- plasma_cytokines %>% 
  mutate(timepoint = as_factor(timepoint)) 

#check it worked
plasma_cytokines$timepoint %>% 
  levels()


# add in the onset_to_admission column
plasma_cytokines <- plasma_cytokines %>% 
  mutate(cestdat = lubridate::ymd(cestdat),
         hostdat = lubridate::ymd(hostdat)) %>% 
  mutate(onset_to_admission = as.numeric(hostdat - cestdat)) %>% 
  rename(visit_day = "timepoint") %>% 
  mutate(visit_day = as.numeric(visit_day)) %>% 
  filter(visit_day %in% c(1:30)) %>% 
  type_convert()


#add in a dailylbdat for the cytokine data for joining
plasma_cytokines <- plasma_cytokines %>% 
  mutate(daily_lbdat = lubridate::ymd(hostdat + (visit_day - 1)))



# add in onset to lab
plasma_cytokines <- plasma_cytokines %>% 
  mutate(onset_to_lab = as.numeric(daily_lbdat - cestdat)) 

#check stepo
plasma_cytokines %>% 
  select(hostdat, daily_lbdat, visit_day, cestdat, onset_to_lab) %>% 
  filter(onset_to_lab == 0) %>% 
  head

#check step
tier2_ccp_data_clean %>% 
  select(hostdat, daily_lbdat, cestdat, onset_to_lab) %>%  
  head

plasma_cytokines %>% 
  missing_glimpse()

plasma_cytokines %>% 
  select(patient_id, onset_to_lab, daily_lbdat, visit_day, assay, value) %>%  
  head(20)


#join ccp with cytokine assay data
cyto_ccp <- plasma_cytokines %>% 
  left_join(
    tier2_ccp_data_clean %>% 
      select(subjid, daily_lbdat, admission_to_lab, 
             daily_wbc_lborres, daily_lymp_lborres, 
             daily_plt_lborres, daily_neutro_lborres, 
             daily_nlr_lborres, daily_crp_lborres, severity2, outcome), 
    by = c("patient_id" = "subjid", "daily_lbdat" = "daily_lbdat"))
###################################################################### 
# check number of patients
cyto_ccp %>% 
  missing_glimpse()
###################################################################### 
cyto_ccp <- cyto_ccp %>% 
  mutate(assay = as_factor(assay)) 


#save and load
save(cyto_ccp, file = 
       "/home/u034/mfourman/cyto_ccp.rda")
###################################################################### 
# load file from above
load("/home/u034/mfourman/cyto_ccp.rda")

cyto_ccp <- cyto_ccp %>% 
  mutate(severity2 = severity2 %>%
           fct_relevel("Ward", "Oxygen alone", "NIV/HFNC", "IMV"))

#Now we have the cyto_ccp dataset joined lets see how big it is etc

cyto_ccp %>%  
  select(patient_id, assay, value, onset_to_lab, daily_wbc_lborres, severity2, outcome) %>% 
  head()

cyto_ccp %>% 
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


cyto_ccp %>% 
  ggplot(aes(onset_to_lab)) +
  geom_histogram()

#censor values for first 50 days
cyto_ccp <- cyto_ccp %>% 
  filter(onset_to_lab >= 0) %>% 
  filter(onset_to_lab <= 50)


#drop days for which we dont have FBC data
cyto_ccp <- cyto_ccp %>% 
  drop_na(daily_neutro_lborres)

#this finds the names of the assays which can be tricky
cyto_ccp$assay %>% 
  levels()

#make a list of cytokines of interest
relevant_cytokines <- c("TNF", "IL-6", "IL-10","IL-17", "GM-CSF...28")

#check missing data now
cyto_ccp %>% 
  select(patient_id, assay, value, onset_to_lab, 
         daily_wbc_lborres, daily_lymp_lborres, 
         daily_plt_lborres, daily_neutro_lborres, 
         daily_nlr_lborres, daily_crp_lborres) %>% 
  missing_glimpse()

cyto_ccp %>%  
  colnames

# create a box plots of day 1 cytokines
cyto_ccp %>% 
  select(patient_id, visit_day, value, assay, severity2) %>% 
  filter(visit_day == 1) %>% 
  filter(assay %in% relevant_cytokines) %>% 
  # filter(visitday_correct.factor == "Day 1-3") %>%
  # group_by(subjid_correct, assay) %>%
  # mutate(value = mean(value, na.rm = TRUE)) %>%
  # distinct(subjid_correct, .keep_all = TRUE) %>% 
  drop_na() %>% 
  ggplot(aes(severity2, value)) +
  geom_jitter(aes(colour = severity2), alpha = 0.3) +
  geom_boxplot(aes(fill = severity2)) +
  facet_wrap(~ assay, scales = "free_y") +
  scale_y_log10() +
  scale_colour_manual(values = c("#3690c0", "#023858", "#c51b7d", "#ff7f00", "#e31a1c")) +
  scale_fill_manual(values = c("#3690c0", "#023858", "#c51b7d", "#ff7f00", "#e31a1c")) + 
  labs(x = "", y = "Value (log scale)", colour = "Severity", fill = "Severity") +
  theme(axis.text.x = element_blank()) + 
  guides(color = FALSE)

# time series mediators
cyto_ccp %>% 
  select(patient_id, onset_to_lab, value, assay, severity2) %>%  
  filter(onset_to_lab >= 0) %>%
  filter(onset_to_lab <= 30) %>%
  filter(assay %in% relevant_cytokines) %>% 
  drop_na() %>% 
  ggplot(aes(onset_to_lab, value, colour = severity2)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  facet_wrap(vars(assay, severity2), scales = "free_y", ncol = 5) +
  scale_y_log10() +
  scale_colour_manual(values = c("#3690c0", "#023858", "#c51b7d", "#ff7f00", "#e31a1c")) +
  labs(x = "Onset to sample (days)", y = "Value (log scale)", colour = "Severity") + 
  theme(legend.position = "top",
        text = element_text(size = 8)) 


cyto_ccp %>%
  filter(assay == "IL-6") %>%
  select(assay, value, daily_crp_lborres ) %>% 
  head()

cyto_ccp %>% 
  select(severity2) %>% 
  head

#Creates a scatterplot of IL6 against CRP 
cyto_ccp %>% 
  drop_na(daily_crp_lborres) %>% 
  select(patient_id, onset_to_lab, assay, value, daily_crp_lborres, severity2) %>% 
  filter(assay == "IL-6") %>% 
  ggplot(aes(x=value, y=daily_crp_lborres)) +
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  stat_cor(method = "pearson") +
  scale_x_log10() +
  ggtitle("CRP vs IL-6") +
  ylab("CRP") +
  xlab("IL-6")

#Creates a scatterplot of IL-6 against Lymphocytes 
cyto_ccp %>% 
  drop_na(daily_crp_lborres) %>% 
  select(patient_id, onset_to_lab, assay, value, daily_lymp_lborres,  outcome) %>% 
  filter(assay == "IL-6") %>% 
  ggplot(aes(x=value, y=daily_lymp_lborres)) +
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  stat_cor(method = "spearman") +
  scale_x_log10() +
  ggtitle("IL-6 vs Lymphocytes") +
  ylab("Lymphocytes") +
  xlab("IL-6 (log scale)")

#Creates a scatterplot of IL-10 against Lymphocytes 
cyto_ccp %>% 
  drop_na(daily_crp_lborres) %>% 
  select(patient_id, onset_to_lab, assay, value, daily_lymp_lborres) %>% 
  filter(assay == "IL-10") %>% 
  ggplot(aes(x=value, y=daily_lymp_lborres)) +
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  stat_cor(method = "spearman") +
  scale_x_log10() +
  ggtitle("IL-10 vs Lymphocytes") +
  ylab("Lymphocytes") +
  xlab("IL-10 (log scale)")

#Creates a scatterplot of  GM-CSF vs Neutrophils
cyto_ccp %>% 
  drop_na(daily_neutro_lborres) %>% 
  select(patient_id, onset_to_lab, assay, value, daily_neutro_lborres, outcome) %>% 
  filter(assay == "GM-CSF...28") %>% 
  ggplot(aes(x=value, y=daily_neutro_lborres)) +
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  stat_cor(method = "pearson") +
  scale_x_log10() +
  ggtitle("GM-CSF vs Neutrophil Count") +
  ylab("Neutrophil Count") +
  xlab("GM-CSF (log scale)")

#Creates a scatterplot of  IL-17 vs Neutrophils
cyto_ccp %>% 
  drop_na(daily_neutro_lborres) %>% 
  select(patient_id, onset_to_lab, assay, value, daily_neutro_lborres, outcome) %>% 
  filter(assay == "IL-17") %>% 
  ggplot(aes(x=value, y=daily_neutro_lborres)) +
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  stat_cor(method = "pearson") +
  scale_x_log10() +
  ggtitle("IL-17 vs Neutrophil Count") +
  ylab("Neutrophil Count") +
  xlab("IL-17 (log scale)")

#Creates a scatterplot of  GM-CSF vs NLR
cyto_ccp %>% 
  drop_na(daily_neutro_lborres) %>% 
  select(patient_id, onset_to_lab, assay, value, daily_nlr_lborres, outcome) %>% 
  filter(assay == "GM-CSF...28") %>% 
  filter(outcome == "Death/IMV") %>% 

  ggplot(aes(x=value, y=daily_nlr_lborres)) +
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  stat_cor(method = "pearson") +
  scale_x_log10() +
  ggtitle("GM-CSF vs NLR") +
  ylab("NLR") +
  xlab("GM-CSF")

# try corrplot
###################################################################### 
##Try time series of cytokines
il6_data <- cyto_ccp %>% 
  filter(assay == "IL-6") %>% 
  distinct()

#get the mean of each #requires summarySE function
il6_plot_data <-
  summarySE(il6_data, 
            measurevar="value", 
            groupvars=c("outcome", "onset_to_lab"))
il6_plot_data %>% head(15)

#create a temporal graph of IL6 data over time
il6_plot_data %>% 
  drop_na(outcome) %>% 
  filter(onset_to_lab <= 14) %>% 
  filter(onset_to_lab >= 0) %>%
  ggplot( aes(x=onset_to_lab, y=value, group=outcome, color=outcome)) +
  geom_ribbon(aes(ymin=value-se, ymax=value+se), alpha=0.1, linetype="dotted") +
  geom_line() +
  geom_point() +
  ggtitle("Mean IL-6 values stratified by outcome") +
  ylab("Mean IL-6 value") +
  xlab("Days from symptom onset to lab measurement") +
  scale_color_npg()


###################################################################### 
 #create a scatter of multiple variables
###################################################################### 

cyto_ccp %>% 
  select(patient_id, onset_to_lab, assay, value, daily_neutro_lborres) %>% 
  filter(assay %in% c("GM-CSF...28", "IL-6")) %>% 
  head() %>% 
  ggplot(aes(value, daily_neutro_lborres)) +
  geom_point() +
  facet_wrap(~ assay, scales = "free_y")


(daily_lbdat, onset_to_lab, admission_to_lab, 
  daily_wbc_lborres, daily_lymp_lborres, 
  daily_plt_lborres, daily_neutro_lborres, 
  daily_nlr_lborres, daily_crp_lborres)
#count the number of values per visit day per patient
cyto_ccp %>% 
  group_by(assay) %>% 
  count()

cyto_ccp %>% 
  colnames

cyto_ccp %>% 
  select(subjid_correct, onset2sample, assay, value, severity2, daily_neutro_lborres ) %>% 
  filter(assay == "IL-6") %>% 
  ggplot(aes(value, daily_crp_lborres, colour = severity2)) +
  geom_point()


###################################################################### 


