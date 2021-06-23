toplinefilems<-read.csv("/home/u034/mcswets/sfr_2020/topline.csv")
ccp_data<-read.csv("/home/u034/mcswets/sfr_2020/ccp_data.csv")
colnames(ccp_data)
sum(unique(toplinefile$subjid))
nrow(toplinefile)
df_1<-toplinefile[,c(1,33,35,36,40,77,78,79,111,112,114,116,118,120,125,126,
                     127:152,154,155,158:183,185,186,191,420,542,557,558,
                     559,560,561,597)]

########################################################################################
#data used 15-12-2020
df_sfr<-read.csv("/home/u034/mcswets/sfr_2020/df_sfr20201207-2.csv")

#####################################################################################

#growth mixed modelling  using lcmm and LCTM tools package
library(lcmm)
library(LCTMtools)
library(tidyr)
library(dplyr)
#hlme function for latent class linear mixed models  (=growth mixed models)
#make sfr dataframe into a long dataframe
df_long<-df_sfr[,c(3:8)]
colnames(df_sfr)
df_long$numeric_id<-seq.int(nrow(df_long))
sfr_long<-gather(df_long, measurement_day, sfr_value, day_0:day_8, factor_key=T)
#hlme function for latent class linear mixed models  (=growth mixed models)
#change measurement day into a number
sfr_long$measurement_day<-gsub("day_","",sfr_long$measurement_day)
sfr_long$measurement_day<-as.numeric(sfr_long$measurement_day)
#set seed and try different number of clusters (=ng variable)
set.seed(2222)
head(sfr_long)
sfr_hlme_1<-hlme(sfr_value ~ measurement_day,
                 subject = "numeric_id", ng = 1, data=sfr_long)
sfr_hlme_2<-hlme(sfr_value ~ measurement_day,
                  mixture= ~ measurement_day,
                  subject = "numeric_id", ng = 2,
                  data=sfr_long, B=sfr_hlme_1)
sfr_hlme_3<-hlme(sfr_value ~ measurement_day,
                 mixture= ~ measurement_day,
                 subject = "numeric_id", ng = 3,
                 data=sfr_long, B=sfr_hlme_1)
sfr_hlme_4<-hlme(sfr_value ~ measurement_day,
                 mixture= ~ measurement_day,
                 subject = "numeric_id", ng = 4,
                 data=sfr_long, B=sfr_hlme_1)
sfr_hlme_5<-hlme(sfr_value ~ measurement_day,
                 mixture= ~ measurement_day,
                 subject = "numeric_id", ng = 5,
                 data=sfr_long, B=sfr_hlme_1)
#pick HLME with lowers BIC value + clinically relevant
summarytable(sfr_hlme_1, sfr_hlme_2, sfr_hlme_3, sfr_hlme_4, sfr_hlme_5)
summarytable(sfr_hlme_4)
#check summary for selected number of clusters
summary(sfr_hlme_4)
#https://rstudio-pubs-static.s3.amazonaws.com/522393_3aa7f65898f8426e9c0a92d7971b619d.html
#check APPA for selected number of clusters
LCTMtoolkit(sfr_hlme_4)
#visualize clusters in plot
newdata<-data.frame(measurement_day=c(0,2,4,6,8))
plot_cluster<-predictY(sfr_hlme_4, newdata, var.time="measurement_day", draws=T,
                       interval= "confidence")
plot_cluster
plot(plot_cluster, legend.loc= "topleft", lty=1, 
     xlab="Day of measurement", ylab ="SFR Value")
rpng.off()
#check size of clusters (+more)
lcmm::postprob(sfr_hlme_4)


#add to dataframe (df_sfr)
library(dplyr)
df_sfr<-df_sfr %>%
  dplyr::rename(numeric_id= num_id)
gmmclust_dataframe<-sfr_hlme_4$pprob[,1:2]
df_sfr<-left_join(df_sfr, gmmclust_dataframe, by ="numeric_id")
df_sfr<-df_sfr %>%
  dplyr::rename(gmmclust= class)

############################### NON LINEAR MODELS-quadratic ################################
#gridsearch> run each hlme model 100 tmes (rep=100) using different start values
#to avoid local maxima. Start values are based on the 1 class model (sfr_hlme_1)
# Latent Class Growth Analysis and Growth Mixture Modeling using R: A tutorial for two R-packages and a comparison with Mplus
#
nonlinear_1<-hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 1, data=sfr_long)
nonlinear_2<-hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                             mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                             subject = "numeric_id", ng = 2, data=sfr_long)
nonlinear_3<-hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                             mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                             subject = "numeric_id", ng = 3, data=sfr_long)
nonlinear_4<-hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                             mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                             subject = "numeric_id", ng = 4, data=sfr_long)
nonlinear_5<-hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                             mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                             subject = "numeric_id", ng = 5, data=sfr_long)
nonlinear_6<-hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                  mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                  subject = "numeric_id", ng = 6, data=sfr_long)

#plot BIC numbers
k<-(c(1,2,3,4,5,6))
BIC_value<-(c(866061,828862,819009,811698,806704,801079))
BIC<-cbind(k,BIC_value)
plot(BIC)
rpng.off()
#pick HLME with lowers BIC value + clinically relevant
summarytable(nonlinear_1, nonlinear_2, nonlinear_3, nonlinear_4, nonlinear_5, nonlinear_6,
             which= c("G", "loglik", "conv", "npm",
                      "AIC", "BIC", "SABIC", "entropy", "%class"))
summarytable(nonlinear_6)
#check summary for selected number of clusters
summary(nonlinear_6)
postprob(nonlinear_4)
#https://rstudio-pubs-static.s3.amazonaws.com/522393_3aa7f65898f8426e9c0a92d7971b619d.html
#check APPA for selected number of clusters
LCTMtoolkit(nonlinear_6)
#increase number of repetitions
nonlinear_4<-gridsearch(rep=5, maxiter=5, minit=nonlinear_1,
                        hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                             mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                             subject = "numeric_id", ng = 4, data=sfr_long))
#visualize clusters in plot
newdata<-data.frame(measurement_day=c(0,2,4,6,8))
plot_cluster<-predictY(nonlinear_6, newdata, var.time="measurement_day", draws=T,
                       interval= "confidence")
plot_cluster
plot(plot_cluster, legend.loc= "topleft", lty=1, 
     xlab="Day of measurement", ylab ="SFR Value")
rpng.off()

summarytable(sfr_hlme_1, sfr_hlme_2, sfr_hlme_3, sfr_hlme_4)
############################### NON LINEAR MODELS- cubic ################################
#gridsearch> run each hlme model 100 tmes (rep=100) using different start values
#to avoid local maxima. Start values are based on the 1 class model (sfr_hlme_1)
# Latent Class Growth Analysis and Growth Mixture Modeling using R: A tutorial for two R-packages and a comparison with Mplus
#
nonlinear_1_c<-hlme(sfr_value ~ poly(measurement_day, degree = 3, raw = TRUE),
                  subject = "numeric_id", ng = 1, data=sfr_long)
nonlinear_2_c<-hlme(sfr_value ~ poly(measurement_day, degree = 3, raw = TRUE),
                  mixture = ~ poly(measurement_day, degree = 3, raw = TRUE),
                  subject = "numeric_id", ng = 2, data=sfr_long)
nonlinear_3_c<-hlme(sfr_value ~ poly(measurement_day, degree = 3, raw = TRUE),
                  mixture = ~ poly(measurement_day, degree = 3, raw = TRUE),
                  subject = "numeric_id", ng = 3, data=sfr_long)
nonlinear_4_c<-hlme(sfr_value ~ poly(measurement_day, degree = 3, raw = TRUE),
                  mixture = ~ poly(measurement_day, degree = 3, raw = TRUE),
                  subject = "numeric_id", ng = 4, data=sfr_long)
nonlinear_5_c<-hlme(sfr_value ~ poly(measurement_day, degree = 3, raw = TRUE),
                  mixture = ~ poly(measurement_day, degree = 3, raw = TRUE),
                  subject = "numeric_id", ng = 5, data=sfr_long)
nonlinear_6_c<-hlme(sfr_value ~ poly(measurement_day, degree = 3, raw = TRUE),
                  mixture = ~ poly(measurement_day, degree = 3, raw = TRUE),
                  subject = "numeric_id", ng = 6, data=sfr_long)

summarytable(nonlinear_1_c, nonlinear_2_c, nonlinear_3_c, nonlinear_4_c,
             which= c("G", "loglik", "conv", "npm",
                      "AIC", "BIC", "SABIC", "entropy", "%class"))
summarytable(nonlinear_6)
#check summary for selected number of clusters
summary(nonlinear_6)
postprob(nonlinear_4)
#https://rstudio-pubs-static.s3.amazonaws.com/522393_3aa7f65898f8426e9c0a92d7971b619d.html
#check APPA for selected number of clusters
LCTMtoolkit(nonlinear_6)
#increase number of repetitions
nonlinear_4<-gridsearch(rep=5, maxiter=5, minit=nonlinear_1,
                        hlme(sfr_value ~ poly(measurement_day, degree = 2, raw = TRUE),
                             mixture = ~ poly(measurement_day, degree = 2, raw = TRUE),
                             subject = "numeric_id", ng = 4, data=sfr_long))
#visualize clusters in plot
newdata<-data.frame(measurement_day=c(0,2,4,6,8))
plot_cluster<-predictY(nonlinear_4_c, newdata, var.time="measurement_day", draws=T,
                       interval= "confidence")
plot_cluster
plot(plot_cluster, legend.loc= "topleft", lty=1, 
     xlab="Day of measurement", ylab ="SFR Value")
rpng.off()
########################################################################################
df_sfr including GMM class membership
write.csv(df_sfr, "df_sfr_17122020.csv")
df_sfr<-read.csv("/home/u034/mcswets/sfr_2020/df_sfr_17122020.csv")
########################################################################################

#plot change in BIC with different clusters
lin<- c(sfr_hlme_1$ng, sfr_hlme_1$BIC)
lin<-rbind(lin, c(sfr_hlme_2$ng, sfr_hlme_2$BIC))
lin<-rbind(lin, c(sfr_hlme_3$ng, sfr_hlme_3$BIC))
lin<-rbind(lin, c(sfr_hlme_4$ng, sfr_hlme_4$BIC))
lin<-rbind(lin, c(sfr_hlme_5$ng, sfr_hlme_5$BIC))
head(lin)
plot(lin, xlab= "Number of clusters",
     ylab= "Bayesian information criteria (BIC)")
rpng.off()

#plot samples from time series for each cluster to see if mean predicted trajectory
# is representative of cluster
cluster1<-subset(df_sfr, df_sfr$gmmclust == 1)
cluster2<-subset(df_sfr, df_sfr$gmmclust == 2)
cluster3<-subset(df_sfr, df_sfr$gmmclust == 3)
cluster4<-subset(df_sfr, df_sfr$gmmclust == 4)
library(ggplot2)
library(tidyr)
cluster1_long<-gather(cluster1, measurement_day, sfr_value, day_0:day_8, factor_key=T)
violin_1<-ggplot(cluster1_long, aes(x=measurement_day, y=sfr_value))
v1<-violin_1 + geom_violin()+ 
  theme_light()+
  ggtitle("Cluster 1")+
  xlab("Measurement day")+
  ylab("SFR value")
cluster2_long<-gather(cluster2, measurement_day, sfr_value, day_0:day_8, factor_key=T)
violin_2<-ggplot(cluster2_long, aes(x=measurement_day, y=sfr_value))
v2<-violin_2 + geom_violin()+ 
  theme_light()+
  ggtitle("Cluster 2")+
  xlab("Measurement day")+
  ylab("SFR value")
cluster3_long<-gather(cluster3, measurement_day, sfr_value, day_0:day_8, factor_key=T)
violin_3<-ggplot(cluster3_long, aes(x=measurement_day, y=sfr_value))
v3<-violin_3 + geom_violin()+ 
  theme_light()+
  ggtitle("Cluster 3")+
  xlab("Measurement day")+
  ylab("SFR value")
cluster4_long<-gather(cluster4, measurement_day, sfr_value, day_0:day_8, factor_key=T)
violin_4<-ggplot(cluster4_long, aes(x=measurement_day, y=sfr_value))
v4<-violin_4 + geom_violin()+ 
  theme_light()+
  ggtitle("Cluster 4")+
  xlab("Measurement day")+
  ylab("SFR value")
#make a grid of the 4 violin plots using cowplot
install.packages("cowplot",repos="https://cran.rstudio.com/")
library(cowplot)
plot_grid(v1,v2,v3,v4)


#transparent graphs for each cluster
p_clus1<- ggplot(data=cluster1_long, aes(x=measurement_day, y=sfr_value))
p1<-p_clus1 + 
  geom_line(aes(group=subjid), alpha=0.1)+
  xlab("Measurement day")+
  ylab("SFR value")+
  ggtitle("Cluster 1")
p_clus2<- ggplot(data=cluster2_long, aes(x=measurement_day, y=sfr_value))
p2<-p_clus2 + 
  geom_line(aes(group=subjid), alpha=0.1)+
  xlab("Measurement day")+
  ylab("SFR value")+
  ggtitle("Cluster 2")
p_clus3<- ggplot(data=cluster3_long, aes(x=measurement_day, y=sfr_value))
p3<-p_clus3 + 
  geom_line(aes(group=subjid), alpha=0.1)+
  xlab("Measurement day")+
  ylab("SFR value")+
  ggtitle("Cluster 3")
p_clus4<- ggplot(data=cluster4_long, aes(x=measurement_day, y=sfr_value))
p4<-p_clus4 + 
  geom_line(aes(group=subjid), alpha=0.1)+quit
xlab("Measurement day")+
  ylab("SFR value")+
  ggtitle("Cluster 4")

library(cowplot)
plot_grid(p1,p2,p3,p4)

#density
d<-density(cluster1_long$sfr_value)
plot(d)
d2<-density(cluster2_long$sfr_value)
plot(d2)
d3<-density(cluster3_long$sfr_value)
plot(d3)
d4<-density(cluster4_long$sfr_value)
plot(d4)
plot_grid(d,d2,d3,d4)
rpng.off()
########################################################################################
#gmmclust characteristics
tapply(df_sfr$age, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$AUC_value, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$worst_sfr, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$onset2admission, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$crp, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$rr_vsorres, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$hodur, df_sfr$gmmclust, quantile, na.rm = TRUE)
tapply(df_sfr$hosp_dur, df_sfr$gmmclust, quantile, na.rm = TRUE)
table(cluster4$resp_support)
table(cluster4$sex)
colnames(cluster1)
table(cluster3$infiltrates_faorres)
table(cluster4$stroke_ceterm)
table(cluster4$clinical_frailty)
table(cluster4$symptom_cluster)
table(cluster4$ethnicity)
sum(is.na(cluster1$ethnicity))
summary(cluster4_long$sfr_value)

stdErr <- function(x) {sd(x)/ sqrt(length(x))}
day0<-tapply(df_sfr$day_0, df_sfr$gmmclus, mean)
day0se<-tapply(df_sfr$day_0, df_sfr$gmmclus, stdErr)
day2<-tapply(df_sfr$day_2, df_sfr$gmmclus, mean)
day4<-tapply(df_sfr$day_4, df_sfr$gmmclus, mean)
day6<-tapply(df_sfr$day_6, df_sfr$gmmclus, mean)
day8<-tapply(df_sfr$day_8, df_sfr$gmmclus, mean)
day2se<-tapply(df_sfr$day_2, df_sfr$gmmclus, stdErr)
day4se<-tapply(df_sfr$day_4, df_sfr$gmmclus, stdErr)
day6se<-tapply(df_sfr$day_6, df_sfr$gmmclus, stdErr)
day8se<-tapply(df_sfr$day_8, df_sfr$gmmclus, stdErr)
mean_sfr_day<-rbind(day0, day2, day4, day6, day8)
se_sfr_day<-rbind(day0se, day2se, day4se, day6se, day8se)
clusters<-c("1", "2", "3", "4")

#switch columns and rows
mean_sfr_day<-t(mean_sfr_day)
se_sfr_day<-t(se_sfr_day)
mean_sfr_day<-data.frame(mean_sfr_day)
se_sfr_day<-data.frame(se_sfr_day)
mean_sfr_day<-cbind(mean_sfr_day, clusters)
se_sfr_day<-cbind(se_sfr_day, clusters)
meansfrday_long<-gather(mean_sfr_day, 
                        measurement_day, sfr_value, day0:day8, factor_key=T)
sesfrday_long<-gather(se_sfr_day, 
                      measurement_day, SE_value, day0se:day8se, factor_key=T)

library(dplyr)
library(ggplot2)
library(tidyr)
sfr_long_df<-cbind(meansfrday_long, sesfrday_long)
sfr_long_df<-subset(sfr_long_df[,c(1:3,6)])
sfr_long_df

ggplot(sfr_long_df, aes(x=measurement_day, 
                        y=sfr_value, colour=clusters, group= clusters)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=sfr_value-SE_value,
                    ymax=sfr_value+SE_value))
#plot mean+ CI non linear function
df_sfr_long<-df_sfr[,c(4:8,150)]
df_sfr_long<-gather(df_sfr_long, measurement_day, sfr_value, day_0:day_8, factor_key=T)
head(df_sfr_long)
ggplot(df_sfr_long, aes(x=measurement_day, 
                        y=sfr_value, colour=clusters, group= clusters))+
  stat_summary(geom="ribbon", fun.data=mean_cl_normal, 
               fun.args=list(conf.int=0.95), fill="lightblue")+
  stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun.y=mean, color="red")


#change 0 to NO and 1 to YES for all comorbidities and complcations
df_sfr[,c(64:68,71:74,77:80,114:142,146)][df_sfr[,c(64:68,71:74,77:80,114:142,146)]== 1]<- "YES"
df_sfr[,c(64:68,71:74,77:80,114:142,146)][df_sfr[,c(64:68,71:74,77:80,114:142,146)]== 0]<- "NO"
df_sfr[,c(64:68,71:74,77:80,114:142,146)][is.na(df_sfr[,c(64:68,71:74,77:80,114:142,146)])]<- "UNKNOWN"
#complications + comorbidities
gmm_comp<-(colSums(cluster1[,c(64:68,71:74,77:80,114:142,146)]== "YES"))
gmm_comp<-(colSums(cluster2[,c(64:68,71:74,77:80,114:142,146)]== "YES"))
gmm_comp<-(colSums(cluster3[,c(64:68,71:74,77:80,114:142,146)]== "YES"))
gmm_comp<-(colSums(cluster4[,c(64:68,71:74,77:80,114:142,146)]== "YES"))
gmm_comp
colnames(cluster1)
#number of comorbidities
gmm_com_num<-(rowSums(cluster4[,c(64:68,71:74,77:80,146)]== "YES"))
table(gmm_com_num)
#number of comorbidities
gmm_compl_num<-(rowSums(cluster4[,c(114:116,118,119,122:125,129:134,136,137,139:142)]== "YES"))
table(gmm_compl_num)

####################################significance###########################################
install.packages("chisq.posthoc.test",repos="https://cran.rstudio.com/")
library(chisq.posthoc.test)
install.packages("BiocManager",repos="https://cran.rstudio.com/")
BiocManager::install("mixOmics")
install.packages("RVAideMemoire",repos="https://cran.rstudio.com/")
library(RVAideMemoire)



cardiac_arrest<-as.table(rbind(c(312, 54, 86, 35),c(4473,3677, 2249,3195)))
cardiac_arrest1<-(c(312, 54, 86, 35))
dimnames(cardiac_arrest)<-list(value=c("cardiac arrest", "no cardiac arrest"),
                               group=c("1", "2", "3", "4"))
cardiac_arrest
chisq.test(cardiac_arrest, correct=FALSE)
chisq.test(cardiac_arrest[,c(2,4)])
chisq.posthoc.test(cardiac_arrest, method= "bonferroni")
chisq.multcomp(cardiac_arrest1, p.method= "bonferroni")

test<-rbind(c(18,90),c(6,149))
chisq.test(test)

####################################viraemia###########################################
#import data 30-11-2020 viraemia set
vir_301120<-read.csv("/home/u034/v1nrodge/viraemia/20201130_vir_metadata.csv")
library(dplyr)
#==============================================================================
head(vir_301120)
#remove some irrelevant columns
vir_301120<-vir_301120[,c(1,5:7)]
#rename subjid and Day column
vir_301120<- vir_301120%>%
  dplyr::rename(subjid= canonical_isaric_id)
vir_301120<- vir_301120%>%
  dplyr::rename(viraemia_day= Day)
#remove rows with missing subjid (=minus 4 rows)
vir_301120<-subset(vir_301120, !is.na(vir_301120$subjid))
#keep highest value for subjects in database twice
vir_301120<- vir_301120 %>%
  group_by (subjid) %>%
  slice(which.max(Vir_Copies))

#join viraemia data with master dataframe
df_sfr<-left_join(df_sfr,vir_301120, by= "subjid")
colnames(df_sfr)

tapply(df_sfr$Vir_Copies, df_sfr$gmmclust, quantile, na.rm = TRUE)
summary(df_sfr$Vir_Copies)
table(df_sfr$Vir_Positive, df_sfr$gmmclust)

#cytokine data from Clark
library(dplyr)
cytokine_data<-read.csv("/Users/Maaike/Rstudio/SFratio/cytokine_data_clark.csv")
c2s(cytokine_data)
#rename 'X' column to subjectid
cytokine_data<-cytokine_data %>%
  dplyr::rename(subjid= X)
#make df with just cluster assignment and subject id
subject_cluster<-df_sfr[,c(3,150)]
#add to cytokine dataframe
cytokine_data<-left_join(cytokine_data,subject_cluster, by= "subjid")
#number of subjects with cytokine and cluster data = 156
cytokine_data<-subset(cytokine_data, !is.na(gmmclust))

#significance 
kruskal.test(cytokine_data$IL.6Ra, cytokine_data$gmmclust)


#co-infection data
library(tidyr)
co_infectiondata<-read.csv("/Users/Maaike/Rstudio/SFratio/resp_infection_data_subjids_first_wave.csv")
c2s(co_infectiondata)
#reshape into long format
co_infection<-reshape(co_infectiondata,
                      direction = "long",
                      varying = list(c("community_acquired_resp_infection", 
                                       "hospital_acquired_resp_infection",
                                       "hospital_acquired_resp_infection_without_community",
                                       "unknown_date_resp_infection")),
                      v.names = "subjid")
#add a column with HAP/CAP for each subject and remove unnecessary columns
co_infection$hapcap<-"hapcap"
co_infection<-co_infection[,c(2,4)]
#only keep unique subjid
co_infection<-aggregate(hapcap ~ subjid, co_infection, function(x) unique(x))
#merge with larger dataframe
co_infection<-left_join(co_infection,subject_cluster, by= "subjid")
#remove subjects without a cluster assignment
co_infection<-subset(co_infection, !is.na(gmmclust))

table(co_infection$gmmclust)


