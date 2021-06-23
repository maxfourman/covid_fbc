## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

isaric_lasso <- function(.data, age, sex, ethnicity, comorbid, rr, spo2, gcs, crp, bun,
                         output = c("vector", "components", "df_vector", "df_components"),
                         na_to_zeros = TRUE, all_na_rm = TRUE){
  .age = enquo(age)
  .sex = enquo(sex)
  .comorbid = enquo(comorbid)
  .rr = enquo(rr)
  .spo2 = enquo(spo2)
  .gcs = enquo(gcs)
  .bun = enquo(bun)
  .crp = enquo(crp)
  out = .data %>%
    dplyr::mutate_at(vars(!! .age, !! .comorbid, !! .rr, !! .spo2, !! .gcs, !! .crp, !! .bun
    ), as.numeric) %>%
    dplyr::mutate_at(vars(!! .sex), ~ as.character(.) %>% tolower() %>% trimws()) %>%
    mutate(
      isaric_lasso_age = case_when(
        !! .age < 50 ~ 0,
        !! .age < 60 ~ 2,
        !! .age < 70 ~ 4,
        !! .age < 80 ~ 6,
        !! .age >= 80 ~ 7,
        TRUE ~ NA_real_),
      isaric_lasso_sex = case_when(
        !! .sex == "male" ~ 1, # Changed from 0
        !! .sex == "female" ~ 0,
        TRUE ~ NA_real_),
      isaric_lasso_comorbid = case_when(
        !! .comorbid == 0 ~ 0,
        !! .comorbid == 1 ~ 1,
        !! .comorbid >= 2 ~ 2,
        TRUE ~ NA_real_),
      isaric_lasso_rr = case_when(
        !! .rr < 20 ~ 0,
        !! .rr < 30 ~ 1,
        !! .rr >= 30 ~ 2,
        TRUE ~ NA_real_),
      isaric_lasso_spo2 = case_when(
        !! .spo2 < 92 ~ 2,
        !! .spo2 >= 92 ~ 0,
        TRUE ~ NA_real_),
      isaric_lasso_gcs = case_when(
        !! .gcs <  15 ~ 2,
        !! .gcs == 15 ~ 0,
        TRUE ~ NA_real_),
      isaric_lasso_bun = case_when(
        !! .bun <= 7 ~ 0,
        !! .bun <= 14 ~ 1,
        !! .bun >  14 ~ 3,
        TRUE ~ NA_real_),
      isaric_lasso_crp = case_when(
        !! .crp < 50 ~ 0,
        !! .crp < 100 ~ 1,
        !! .crp >=100 ~ 2,
        TRUE ~ NA_real_)
    ) %>%
    mutate(
      isaric_lasso = rowSums(dplyr::select(., dplyr::starts_with("isaric_lasso_")),
                             na.rm = na_to_zeros)
    ) %>%
    {if(all_na_rm){
      dplyr::mutate(., isaric_lasso = dplyr::if_else(
        is.na(isaric_lasso_age) &
          is.na(isaric_lasso_sex) &
          is.na(isaric_lasso_comorbid) &
          is.na(isaric_lasso_rr) &
          is.na(isaric_lasso_spo2) &
          is.na(isaric_lasso_gcs) &
          is.na(isaric_lasso_bun) &
          is.na(isaric_lasso_crp), NA_real_, isaric_lasso))
    } else {
      .
    }}
  if(output == "vector"){
    out %>%
      pull(isaric_lasso)
  } else if(output == "components"){
    out %>%
      select(starts_with("isaric_lasso"))
  } else if(output == "df_vector"){
    out %>%
      pull(isaric_lasso) %>%
      bind_cols(.data, isaric_lasso = .)
  } else if(output == "df_components"){
    out
  }
}

