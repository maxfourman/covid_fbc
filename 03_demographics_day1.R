#03_demographics_day1

###################################################################### 

## Code author: Max Fourman

## Description: 
# This code loads the topline file made in 01_clean_ccp_data
# Then loads the 

###################################################################### 
## Load packages
# May be already done if using previous scripts

install.packages("moonBook", repos='https://cran.ma.imperial.ac.uk/')
library("remoter") 
remoter::client(port=60117)
library("ggsci")
library("tidyverse")
library("finalfit")
library("pROC")
library("ggpubr")



###################################################################### 

# Load data from topline
# Load the cleaned trak data
load("/home/u034/mfourman/topline.rda")
load("/home/u034/mfourman/trak_data_clean.rda")
load("/home/u034/mfourman/ccp_data_clean.rda")



# add in the columns to the topline file with diagnosis for later comparison with trak data
topline <- topline %>% 
  mutate(diagnosis = "covid19") %>% 
  mutate(group_diagnosis = "covid19")

#Take a smaller subset of topline to merge with the trak data
topline_select_variables <- topline %>% 
  select(subjid, sex, calc_age, daily_wbc_lborres, daily_lymp_lborres, daily_plt_lborres, daily_neutro_lborres, diagnosis, group_diagnosis) %>% 
  mutate(Monocytes = NA)


#get rid of reducnant variables from the trak data
trak_data_clean <- trak_data_clean %>% 
  select(!WARD_CODE) 

#check column names and match the order
trak_data_clean %>%  colnames()
trak_data_clean <- trak_data_clean %>% relocate(Monocytes, .after = last_col())
topline_select_variables %>% colnames

#Make the column names the same
colnames(topline_select_variables) <- colnames(trak_data_clean)

#Bind the two data frames
day_one_combined_data <- rbind(topline_select_variables, trak_data_clean)

day_one_combined_data %>% tail()

day_one_combined_data %>% 
  select(SEX) %>% 
  count()


# first add the NLR column 
day_one_combined_data <- day_one_combined_data %>% 
  mutate(NLR = Neutrophils/Lymphocytes)

#relabel variables
day_one_combined_data <- day_one_combined_data %>% 
  mutate(
    WBC = ff_label(WBC, "WBC"),
    AGE = ff_label(AGE, "Age"),
    SEX = ff_label(SEX, "Sex"),
    Lymphocytes = ff_label(Lymphocytes, "Lymphocytes"),
    Platelets = ff_label(Platelets, "Platelets"),
    Neutrophils = ff_label(Neutrophils, "Neutrophils"))

# FIx the column Names
day_one_combined_data <- day_one_combined_data %>% 
  mutate(group_diagnosis =            
           group_diagnosis %>%    
           factor() %>%        
           fct_recode(           
             "Bacterial" = "bacterial",                 
             "COVID-19"   = "covid19",
             "Other Viral" = "non_flu_viral") %>% 
           ff_label("Group Diagnosis"))
  
# check column names
day_one_combined_data %>% colnames


#this creates a combined table of the variable with data for day 1
day_one_combined_data %>% 
  mutate_if(is.factor, fct_explicit_na) %>% 
  summary_factorlist(
    dependent = c("group_diagnosis"), 
    explanatory = c(
      "AGE", "SEX", "WBC", "Lymphocytes", "Neutrophils", "Monocytes", "Platelets"), 
    p = TRUE, 
    na_include = TRUE, 
    column = TRUE, 
    cont = "median",
    total_col = TRUE, 
    add_col_totals = TRUE, 
    add_row_totals = TRUE, 
    include_row_missing_col = FALSE
  ) -> t1
t1



###################################################################### 
### This creates day 1 violin and box plots



#Create a violin + box plot of the Lymphocytes stratified by organism
day_one_combined_data %>% 
  filter(Lymphocytes <= 30) %>% 
  ggplot(aes(x=group_diagnosis, y=Lymphocytes, colour=group_diagnosis)) + 
  geom_violin() + 
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  geom_boxplot(width=.1, outlier.shape = NA) +
  labs(x = element_blank(), 
       y = (expression("Lymphocyte Count (x10"^9*"/L)")),
       colour = "Aetiology") + 
  geom_hline(yintercept=1.5, linetype="dashed", color = "gray20", alpha = 0.5) +
  scale_y_continuous(breaks=seq(0,10,0.5)) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_colour_npg()

#Create a violin + box plot of the Neutrophils stratified by organism
day_one_combined_data %>% 
  filter(Neutrophils <= 50) %>% 
  ggplot(aes(x=group_diagnosis, y=Neutrophils, color=group_diagnosis)) + 
  geom_violin() + 
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  geom_boxplot(width=.1, outlier.shape = NA) +
  labs(x = element_blank(), 
       y = (expression("Neutrophil Count (x10"^9*"/L)")),
       colour = "Aetiology") + 
  geom_hline(yintercept=7.5, linetype="dashed", color = "gray20", alpha = 0.5) +
  scale_y_continuous(breaks=seq(0,30,2.5)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_color_npg()

#Create a violin + box plot of the NLR stratified by organism
day_one_combined_data %>% 
  filter(NLR <= 100) %>% 
  ggplot(aes(x=group_diagnosis, y=NLR, color=group_diagnosis)) + 
    geom_violin() + 
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    geom_boxplot(width=.1, outlier.shape = NA) +
  labs(x = element_blank(), 
       y = (expression("Neutrophil:Lymphocyte Ratio")),
       colour = "Aetiology") +  
    scale_y_continuous(breaks=seq(0,100,05)) +
    coord_cartesian(ylim = c(0, 50)) +
    scale_color_npg()

###################################################################### 
# Make table of percentage outside of reference ranges

# Creates factor outside reference in each group
day_one_combined_data <- day_one_combined_data %>% 
  mutate(Lymphopenia = case_when(
    Lymphocytes <=1.5 ~ "Lymphopenia",
    Lymphocytes >=4.5 ~ "Lymphocytosis",
    TRUE ~ "Normal"
  ) %>% 
    factor(levels = c("Lymphocytosis", "Normal", "Lymphopenia"))
  ) %>% 
  mutate(Neutropenia = case_when(
    Neutrophils <=2.0 ~ "Neutropenia",
    Neutrophils >=7.5 ~ "Neutrophilia",
    TRUE ~ "Normal"
  ) %>% 
    factor(levels = c("Neutrophilia", "Normal", "Neutropenia"))
  ) 

explanatory_t4 = c("Lymphopenia", "Neutropenia")
#make a table of who is neutropenic
day_one_combined_data %>% 
  summary_factorlist("group_diagnosis", 
                     explanatory_t4, 
                     p = TRUE, 
                     column = TRUE, 
                     cont = "median",
                     add_col_totals = TRUE, 
                     include_row_missing_col = FALSE) -> t4
t4


ccp_data %>% 
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

topline %>% colnames

#create explanatory factor list
explanatory_t2 = c("age.factor", "sex", "ethnicity",
                "chrincard", "renal_mhyn", 
                "malignantneo_mhyn", "modliv", "obesity_mhyn", "chronicpul_mhyn",
                "diabetes_mhyn")
explanatory_t3 = c("daily_wbc_lborres", 
                "daily_neutro_lborres", "daily_lymp_lborres", "daily_plt_lborres",
                "daily_crp_lborres")

topline %>% 
  distinct(subjid, .keep_all = TRUE) %>% 
  summary_factorlist("severity2", 
                     explanatory_t2, 
                     p = FALSE, 
                     column = TRUE, 
                     cont = "median",
                     total_col = TRUE, 
                     add_col_totals = TRUE, 
                     add_row_totals = TRUE, 
                     include_row_missing_col = FALSE) -> t2
t2

topline %>% 
  distinct(subjid, .keep_all = TRUE) %>% 
  summary_factorlist("outcome", 
                     explanatory_t3, 
                     p = TRUE, 
                     column = TRUE, 
                     cont = "median",
                     total_col = TRUE, 
                     add_col_totals = TRUE, 
                     add_row_totals = TRUE, 
                     include_row_missing_col = FALSE) -> t3
t3



save(t1, t2, t3, file = "out.rda")

# check for missing data - also creates a nice table
topline %>% select("daily_wbc_lborres", "daily_lymp_lborres",
                   "daily_neutro_lborres", "daily_haematocrit_lborres",
                   "daily_plt_lborres", "daily_aptt_lborres",
                   "daily_pt_lborres", "daily_esr_lborres",
                   "daily_alt_lborres",
                   "daily_ast_lborres", "daily_glucose_lborres",
                   "daily_lactate_lborres", "daily_ldh_lborres",
                   "daily_crp_lborres") %>% 
  ff_glimpse()

#this creates the biochem table (TABLE 2)
topline %>% 
  mutate_if(is.factor, fct_explicit_na) %>% 
  summary_factorlist(
    explanatory = c(
      "daily_wbc_lborres", "daily_lymp_lborres", "daily_neutro_lborres", "daily_haematocrit_lborres",
      "daily_plt_lborres", "daily_aptt_lborres", "daily_pt_lborres", "daily_esr_lborres", "daily_alt_lborres",
      "daily_ast_lborres", "daily_glucose_lborres", "daily_lactate_lborres", "daily_ldh_lborres",
      "daily_crp_lborres"), 
    na_include = TRUE, 
    add_col_totals = TRUE
  ) -> t2
t2

###################################################################### 
#Attempts to make a ROC 

#if needed to find lothian only data
topline_lothian <- topline %>% 
  filter(place_name.x %in% c("Western General Hospital", 
                                        "Royal Infirmary Of Edinburgh At Little France", 
                                        "St John's Hospital"))

day_one_combined_data %>% 
  select(group_diagnosis) %>% 
  head()

#First make a covid vs non covid 
day_one_combined_data <- day_one_combined_data %>% 
  mutate(covid_status = case_when(
    group_diagnosis == "COVID-19" ~ "COVID-19", 
    TRUE ~ "non_covid_pneumonia") %>% 
      factor(levels = c("COVID-19", "non_covid_pneumonia")))

#check this worked
day_one_combined_data %>% 
  select(covid_status) %>% 
  summary()

nlr_severity_roc <- roc(outcome ~ daily_nlr_lborres, topline)
auc(nlr_severity_roc)

#build individual roc curves
lymph_roc <- roc(covid_status ~ Lymphocytes, day_one_combined_data)
neuts_roc <- roc(covid_status ~ Neutrophils, day_one_combined_data)
nlr_roc <- roc(covid_status ~ NLR, day_one_combined_data)

#calculate AUC and CI
auc(lymph_roc)
auc(neuts_roc)
auc(nlr_roc)
#calculate confidence interval
ci.auc(lymph_roc)
ci.auc(neuts_roc)
ci.auc(nlr_roc)

coords(neuts_roc, x="best", best.method="youden")

power.roc.test(lymph_roc)

#create a roc plot with all 3 on one graph
rocs <- roc(covid_status ~ Lymphocytes + Neutrophils + NLR, 
            data = day_one_combined_data,
            print.auc=TRUE)
rocs

combined_roc_plot <- ggroc(rocs, legacy.axes = TRUE)


#Now plot
combined_roc_plot + theme_minimal() +
  scale_colour_manual(values=c("#2A6EBB", "#E37222", "#333333"),
                      name=NULL,
                      labels=c("Lymphocyte count\n(AUC=0.588; 95% CI:0.583-0.593)\n",
                               "Neutrophil count\n(AUC=0.7132; 95% CI:0.709-0.718)\n",
                               "Neutrophil:Lymphocyte Ratio\n(AUC=0.5963 ; 95% CI:0.591-0.601)")) +
  labs(x="1 - Specificity",
       y="Sensitivity") + 
  theme(text=element_text(family = 'sans'),
        legend.position=c(.65, .25),
        legend.background = element_rect(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        aspect.ratio=1) + 
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour="#E3E4E4", linetype="dashed", size=.25) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), colour="#333333", linetype="solid", size=.4) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), colour="#333333", linetype="solid", size=.4)

###################################################################### 
#Roc plot of lymphocyte monocyte ratio 

#first convert to a factor
trak_data_clean$group_diagnosis <- as.factor(trak_data_clean$group_diagnosis)
trak_data_clean_viral_only %>% group_by(group_diagnosis) %>% summary

#add LMR column
trak_data_clean <- trak_data_clean %>% 
  mutate(LMR = Lymphocytes/Monocytes)

trak_data_clean %>%  select(group_diagnosis) %>%  filter(group_diagnosis != "bacterial") %>%  head


#build individual roc curves
lmr_roc <- roc(controls=trak_data_clean$LMR[trak_data_clean$group_diagnosis=="non_flu_viral"], 
               cases=trak_data_clean$LMR[trak_data_clean$group_diagnosis=="Influenza"])

auc(lmr_roc)

lmr_roc_plot <- ggroc(lmr_roc, legacy.axes = TRUE)
lmr_roc_plot + theme_minimal() +
  labs(x="1 - Specificity",
       y="Sensitivity") + 
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour="#E3E4E4", linetype="dashed", size=.25) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), colour="#333333", linetype="solid", size=.4) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), colour="#333333", linetype="solid", size=.4) +
  ggtitle("ROC plot of LMR to differentiate between \n Flu and non_flu viral pneumonia")

#plot a violin plot
trak_data_clean %>% 
  filter(group_diagnosis != "bacterial") %>% 
  ggplot(aes(x=group_diagnosis, y=LMR, color=group_diagnosis)) + 
  geom_violin() +
  geom_boxplot(width=.1, outlier.shape = NA) +
  theme_minimal() + 
  scale_color_brewer(palette = "Pastel1") +
  ylim(0, 25) +
  ggtitle("Violin plot of LMR  \n Flu and non_flu viral pneumonia")

trak_data_clean %>% colnames

#this creates a combined table of the variable with data for day 1
trak_data_clean %>% 
  filter(group_diagnosis != "bacterial") %>% 
  summary_factorlist(
    dependent = c("group_diagnosis"), 
    explanatory = c(
      "AGE", "SEX", "WBC", "Lymphocytes", "Platelets", "Neutrophils", "Monocytes", "LMR"), 
    na_include = TRUE, 
    add_col_totals = TRUE,
    cont = "median", 
    p = TRUE
  ) -> t3
t3

