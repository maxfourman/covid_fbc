# 02_GMM_attempts
install.packages("lcmm", repos='https://cran.ma.imperial.ac.uk/')
install.packages("LCTMtools", repos='https://cran.ma.imperial.ac.uk/')


library(lcmm)
library(LCTMtools)
library(tidyverse)
library(tidyr)

#load the prepped data from 01_analysis
load("/home/u034/mfourman/ccp_blood_data_long.rda")

##to start with work with lymphocyte data _ limited to 14 days
lymphocyte_data <- ccp_blood_data_long %>% 
  filter(assay %in% c("daily_lymp_lborres")) %>% 
  filter(onset_to_lab <= 14 ) %>% 
  distinct()

lymphocyte_data %>% head()

 
#Find out how many entries per patient
data_counts <- ccp_blood_data %>% 
  group_by(subjid) %>%
  dplyr::summarise(count = n())
#check out these numbers
data_counts %>% 
  summary(count)

# create a histogram of x = observations y = count
data_counts %>% 
  filter(count <= 10) %>% 
  ggplot(aes(count)) +
  geom_histogram(bins = 10) +
  scale_x_continuous(name = "number of observations",
                     breaks=seq(0,10,1))



### First look for missing data
missing_table <- ccp_blood_data_long %>% 
  group_by(onset_to_lab_case) %>% 
  dplyr::summarize(count = n()) 

missing_table%>% 
  ggplot(aes(x=onset_to_lab_case, y=count)) +
  geom_bar(stat='identity')





lymphocyte_data %>% head()

set.seed(2222)

### DO THIS FIRST TO GET DATA INTO CORRECT FORMAT FOR LCMM
# convert subjid to factor, create numeric ID and make tibble into dataframe
lymphocyte_data$subjid <- as.factor(lymphocyte_data$subjid)
lymphocyte_data$numeric_id <- as.numeric(factor(lymphocyte_data$subjid, 
                                           levels=unique(lymphocyte_data$subjid)))
lymphocyte_data <- lymphocyte_data %>% as.data.frame()

lcga1 <- hlme(value ~ onset_to_lab_case, mixture = ~ onset_to_lab_case, subject = "numeric_id", ng = 2, data = lymphocyte_data)


###STEP 1 - LCGA with lcmm
###
lcga1 <- hlme(value ~ onset_to_lab_case, subject = "numeric_id", ng = 1, data = lymphocyte_data)
lcga2 <- gridsearch(rep = 10, maxiter = 10, minit = lcga1,
                    hlme(value ~ onset_to_lab_case, subject = "numeric_id",
                         ng = 2, data = lymphocyte_data, mixture = ~ onset_to_lab_case)) 
lcga3 <- gridsearch(rep = 10, maxiter = 10, minit = lcga1, 
                    hlme(value ~ onset_to_lab_case, subject = "numeric_id",
                         ng = 3, data = lymphocyte_data, mixture = ~ onset_to_lab_case)) 

lcga4 <- gridsearch(rep = 10, maxiter = 10, minit = lcga1, 
                    hlme(value ~ onset_to_lab_case, subject = "numeric_id",
                         ng = 4, data = lymphocyte_data, mixture = ~ onset_to_lab_case)) 

lcga_summary <- summarytable(lcga1, lcga2)
summary(lcga1)
save(lcga_summary, file = 
       "/home/u034/mfourman/lcga_summary.rda")
load("/home/u034/mfourman/lcga_summary.rda")
lcga_summary

###STEP 2 - GMM with only a random intercept
set.seed(2021)
gmm1 <- hlme(value ~ onset_to_lab, 
             subject = "numeric_id", 
             random=~1, ng = 1, 
             data = lymphocyte_data)
gmm2 <- gridsearch(rep = 10, 
                   maxiter = 10, 
                   minit = gmm1, 
                   hlme(value ~ onset_to_lab, 
                        subject = "numeric_id", 
                        random=~1,
                        ng = 2, 
                        data = lymphocyte_data, 
                        mixture = ~ onset_to_lab, 
                        nwg=T))
gmm3 <- gridsearch(rep = 100, 
                   maxiter = 10, 
                   minit = gmm1,
                   hlme(value ~ onset_to_lab, 
                        subject = "numeric_id", 
                        random=~1,
                        ng = 3, 
                        data = lymphocyte_data, 
                        mixture = ~ onset_to_lab, nwg=T))
lymphocyte_data %>% head()
###STEP 2 Refine the preliminary working model from step 1 to determine the optimal number of classes, testing K=1,...7. 
### The number of classes chosen may be chosen based on the lowest Bayesian information criteria (BIC).
gmm_1 <- lcmm::hlme(fixed = value ~ onset_to_lab,
                  random = ~ onset_to_lab,
                  ng = 1,
                  idiag = FALSE, 
                  data = data.frame(lymphocyte_data), subject = "numeric_id")

lin <- c(gmm_1$ng, gmm_1$BIC)

for (i in 3:3) {
  mi <- lcmm::hlme(fixed = value ~ onset_to_lab,
                   mixture = ~ onset_to_lab,
                   random = ~ onset_to_lab,
                   ng = i, nwg = TRUE, 
                   idiag = FALSE, 
                   data = data.frame(lymphocyte_data), subject = "numeric_id")
  
  lin <- rbind(lin, c(i, mi$BIC))
}

modelout <- knitr::kable(lin, col.names = c("k", "BIC"), row.names = FALSE, align = "c")
modelout




#### Maaike way
leuk_hlme_1 <- hlme(value ~ onset_to_lab,
                 subject = "numeric_id", ng = 1, data=lymphocyte_data)
leuk_hlme_2 <- hlme(value ~ onset_to_lab,
                 mixture= ~ onset_to_lab,
                 subject = "numeric_id", ng = 2,
                 data=lymphocyte_data, B=leuk_hlme_1)
leuk_hlme_3 <- hlme(value ~ onset_to_lab,
                    mixture= ~ onset_to_lab,
                 subject = "numeric_id", ng = 3,
                 data=lymphocyte_data, B=leuk_hlme_1)
leuk_hlme_4 <- hlme(value ~ onset_to_lab,
                    mixture= ~ onset_to_lab,
                 subject = "numeric_id", ng = 4,
                 data=lymphocyte_data, B=leuk_hlme_1)
leuk_hlme_5 <- hlme(value ~ onset_to_lab,
                    mixture= ~ onset_to_lab,
                 subject = "numeric_id", ng = 5,
                 data=lymphocyte_data, B=leuk_hlme_1)



############################### NON LINEAR MODELS-quadratic ################################
#gridsearch> run each hlme model 100 tmes (rep=100) using different start values
#to avoid local maxima. Start values are based on the 1 class model (sfr_hlme_1)
# Latent Class Growth Analysis and Growth Mixture Modeling using R: A tutorial for two R-packages and a comparison with Mplus
#
nonlinear_1<-hlme(value ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 1, data=lymphocyte_data)
nonlinear_2<-hlme(value ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  mixture = ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 2, data=lymphocyte_data)

nonlinear_3<-hlme(value ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  mixture = ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 3, data=lymphocyte_data)
nonlinear_4<-hlme(value ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  mixture = ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 4, data=lymphocyte_data)
nonlinear_5<-hlme(value ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  mixture = ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 5, data=lymphocyte_data)
nonlinear_6<-hlme(value ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  mixture = ~ poly(onset_to_lab, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 6, data=sfr_long)


summarytable(nonlinear_1, nonlinear_2)
