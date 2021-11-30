#### Data Simulation ####
###simulate some exemplary data for testing models
library(MASS)
library(tidyverse)
library(BayesFactor)
library(lme4)


simncalc <- function(n,
                     coeffs_bl, covmatrix, usable, mean_diff_y, sd_diff_y,
                     effect_age_change, sd_age_change, 
                     effect_SBP_baseline, effect_WHR_baseline, 
                     effect_SBP_change, effect_WHR_change, 
                     sd_reff, sd_error){
  
  # This function simulates data for the effects of baseline SBP and WHR and 
  # change in SBP and WHR on WML load.
  # It uses cross-sectional associations from LIFE baseline data and
  # estimates of longitudinal effects from the literature.
  
  #Derive joint distribution of Age, SBP, WHR, ICV with same properties (mean and sd) as baseline data
  datasim <- as.data.frame(mvrnorm(n=n, mu=c(mean(usable$Age_all), mean(usable$ADULT_BP_SBP.1, na.rm=T), 
                                             mean(usable$waist2hip), mean(usable$icv)), Sigma=covmatrix, empirical=TRUE)) # simulate data
  datasim$sex <- rbinom(nrow(datasim),1,table(usable$sex)[1]/table(usable$sex)[2])  
  
  
  #Select sample
  sub <- tibble(
    id = 1:n,
    sub_intercept  = rnorm(n, 0, sd_reff), # random intercept
    age_base = datasim[,"Age_all"],
    sex = datasim[,"sex"],
    SBP_base = datasim[,"ADULT_BP_SBP.1"],
    WHR_base = datasim[, "waist2hip"],
    icv=datasim[,'icv'],
    age_change = rnorm(n, mean_diff_y,sd_diff_y),#based on what is known about followup
    ###########################################################################################
    #from longitudinal studies -> increase in SBP also in older age, supported by Framingham heart (Cheng et al, 2012)
    #and Baltimore Longitudinal Study of Aging: 8.5 mmHg/10 years for men, 4.4. mmHg/decade for women at age 60
    #and Whitehall 2: 1 mmHg/y for older men/women (60 - 70 years)
    #weighted average of SBP change from these three studies: 0.76 mmHg/y
    SBP_change = rnorm(n, mean=0.76*6, sd=4), 
    ##### from interventional studies -> stronger change in SBP
    #Havenon: 15 mmHg difference between standard and intensive BP control arm.
    #Nasrallah: difference, 14.2 mm Hg [95% CI, 13.1 to 15.3 mm Hg])
    WHR_change = rnorm(n, mean=(0.0047/5)*6, sd=0.005),
    #from Baltimore Longitudinal Study of Aging:
    #0.0073 in men, 0.0021 in women over 5 years
  )
  
  tp_n  <- 2
  tp <- tibble(
    tp = c("bl","fu")
  )
  
  df <- crossing(
    id = sub$id, # get subject IDs from the sub data table
    tp = tp$tp, # get stimulus IDs from the stim data table
  ) %>%
    left_join(sub, by = "id") %>% # includes the intercept and conditin for each subject
    left_join(tp, by = "tp")
  
  df <- df %>%
    group_by(id) %>%
    mutate(
      age_change = case_when(
        tp == "bl" ~ 0,
        TRUE ~ age_change
      ),
      SBP_change = case_when(
        tp == "bl" ~ 0,
        TRUE ~ SBP_change
      ),
      WHR_change = case_when(
        tp == "bl" ~ 0,
        TRUE ~ WHR_change
      )
    )
  df <- df %>%
     mutate(intercept = coeffs_bl[1])

  #have to split the simulation in two parts, so that we have dv in correct units
  df <- df %>%
    plyr::mutate(
      res_error = rnorm(nrow(.), 0, sd_err),
      errors=sub_intercept + res_error,
      cross_effects = exp(intercept + 
        scale(age_base, scale=F)* (coeffs_bl[2]) + 
        coeffs_bl[3] * scale(SBP_base, scale=F)  +
        coeffs_bl[4] * scale(WHR_base, scale=F)  + 
        coeffs_bl[5] * sex + coeffs_bl[6]*scale(icv, scale=F)),
      long_effects= (effect_age_change + effect_SBP_baseline * scale(SBP_base, scale=F) + 
                    effect_WHR_baseline*scale(WHR_base, scale=F)) * age_change+
                    effect_SBP_change * SBP_change + effect_WHR_change * WHR_change,
     
      dv = cross_effects + long_effects + errors,
      dv2 = cross_effects + (effect_age_change + effect_WHR_baseline*scale(WHR_base, scale=F)) * age_change,
      #only for quality control
      change_coeff=effect_age_change + effect_SBP_baseline * scale(SBP_base, scale=F) + 
        effect_WHR_baseline*scale(WHR_base, scale=F),
      
      WHR_effects = effect_WHR_baseline*scale(WHR_base, scale=F)* age_change,
      SBP_effects = effect_SBP_baseline * scale(SBP_base, scale=F)* age_change)
  return(df)
}

