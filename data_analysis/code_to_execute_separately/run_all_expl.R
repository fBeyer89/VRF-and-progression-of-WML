library("lme4")
library("lmerTest") 
library("BayesFactor")
library(broom.mixed)
library(mice)
library(miceadds)
library(readxl) #read xlsx, ods etc.
library(tidyverse) #includes ggplot2, dplyr, tidyr, readr, tibble, stringr ...
library(lubridate) #package to work with data and time
#library(eeptools)
#library(data.table)
library(naniar) # helps to deal with NA
#library(lme4)
#library(stringr)

setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/run_exp_LME.R')

miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/imputed_data_08.5.23/imputed_data_08.5.23.Rdata")

# #M1E2a: modifying effect of gender on age-related WML increase
# res_VRF= run_exp_LME(imp, "E2a_sex", n_it=10)
# write.csv(res_VRF[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2a_age_model_res_freq.csv')
# write.csv(res_VRF[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2a_age_model_res_bayes.csv')
# 
# #M2Ea_DBPchange_sex: modifying effect of gender on DBP change on WML progression
# res_DBPchange_sex= run_exp_LME(imp, "E2a_sex_DBP_change", n_it=10)
# write.csv(res_DBPchange_sex[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2a_sex_DBP_change_model_res_freq.csv')
# write.csv(res_DBPchange_sex[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2a_sex_DBP_change_model_res_bayes.csv')
# 
# #M1E2b: modifying effect of gender on DBP-related WML increase
# res_VRF= run_exp_LME(imp, "E2b_DBP", n_it=10)
# write.csv(res_VRF[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2b_DBP_model_res_freq.csv')
# write.csv(res_VRF[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2b_DBP_model_res_bayes.csv')
# 
# #M1E2c: modifying effect of gender on WHR-related WML increase
# res_VRF= run_exp_LME(imp, "E2c_WHR", n_it=10)
# write.csv(res_VRF[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2c_WHR_model_res_freq.csv')
# write.csv(res_VRF[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E2c_WHR_model_res_bayes.csv')
# 
# #M2E3a: modifying effect of gender on WML increase on executive function
# res_exec= run_exp_LME(imp, "E3a_exfunct", n_it=10)
# write.csv(res_exec[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E3a_exfunct_model_res_freq.csv')
# write.csv(res_exec[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E3a_exfunct_model_res_bayes.csv')
# 
# #M3E3b: modifying effect of gender on WML increase on global cognition
# res_globalcog= run_exp_LME(imp, "E3b_globalcog", n_it=10)
# write.csv(res_globalcog[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E3b_globalcog_model_res_freq.csv')
# write.csv(res_globalcog[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/E3b_globalcog_model_res_bayes.csv')

####################################################
# Run sex-stratified analyses (if interaction not significant)
####################################################
source("run_conf_LME_sex_stratified.R")
#Select only females
imp_s0 <- miceadds::subset_datlist(imp, 
                                   select=c(-22),
                                   expr_subset=expression(sex==0),
                                   toclass="mids" ) #toclass command does not work ?!
i=complete(imp_s0, "long",include=TRUE) #workaround to create mids object which can be fed into mice::with
imp_s0_r=as.mids(i)
res_VRF= run_conf_LME_sexstr(imp_s0, sex=0, model="M1_VRF",n_it=10)
write.csv(res_VRF[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M1_VRF_females_model_res_freq.csv')
write.csv(res_VRF[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M1_VRF_females_model_res_bayes.csv')

imp_s1 <- miceadds::subset_datlist(imp, 
                                   select=c(-22),
                                   expr_subset=expression(sex==1),
                                   toclass="mids" ) #toclass command does not work ?!
i=complete(imp_s1, "long",include=TRUE) #workaround to create mids object which can be fed into mice::with
imp_s1_r=as.mids(i)

res_VRF= run_conf_LME_sexstr(imp_s1, sex=1, model="M1_VRF",n_it=10)
write.csv(res_VRF[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M1_VRF_males_model_res_freq.csv')
write.csv(res_VRF[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M1_VRF_males_model_res_bayes.csv')


########################################################
res_M2_exfunct= run_conf_LME_sexstr(imp_s0, sex=0, model="M2_exfunct",n_it=10)
write.csv(res_M2_exfunct[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M2_exfunct_females_model_res_freq.csv')
write.csv(res_M2_exfunct[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M2_exfunct_females_model_res_bayes.csv')

#Select only males
res_M2_exfunct= run_conf_LME_sexstr(imp_s1, sex=1, model="M2_exfunct",n_it=10)
write.csv(res_M2_exfunct[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M2_exfunct_males_model_res_freq.csv')
write.csv(res_M2_exfunct[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M2_exfunct_males_model_res_bayes.csv')


####################################################
res_M3_globalcog= run_conf_LME_sexstr(imp_s0,sex=0, model="M3_globalcog",n_it=10)
write.csv(res_M3_globalcog[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M3_globalcog_females_model_res_freq.csv')
write.csv(res_M3_globalcog[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M3_globalcog_females_model_res_bayes.csv')

#Select only males
res_M3_globalcog= run_conf_LME_sexstr(imp_s1, sex=1, model="M3_globalcog",n_it=10)
write.csv(res_M3_globalcog[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M3_globalcog_males_model_res_freq.csv')
write.csv(res_M3_globalcog[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/M3_globalcog_males_model_res_bayes.csv')
