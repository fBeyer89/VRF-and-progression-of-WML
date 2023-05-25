library("lme4")
library("lmerTest") 
library("BayesFactor")
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer
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
library(mice)
library(performance)
library(influence.ME)
library(car)

setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('test_LME_assumptions.R')

miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/imputed_data_08.5.23/imputed_data_08.5.23.Rdata")

####################################################
# Test assumptions (saves everything automatically)
####################################################
#VRF
print("run VRF")
res_VRF= test_LME_assumptions(imp, "M1_VRF")

#executive function
res_exec= test_LME_assumptions(imp, "M2_exfunct")

#global cognition
res_globalcog= test_LME_assumptions(imp, "M3_globalcog")

####################################################
# Test collinearity once for VRF and cognition models without interactions
# (done only for one imputation as only predictor variables matter)
####################################################
comp_imp=mice::complete(imp, "long")
imp_data=comp_imp[comp_imp$.imp==1,]

test_vif=lm(asinh_wml ~ age_base + age_change + 
              DBP_base + DBP_change + 
              WHR_base + WHR_change + 
              sex + BPmed + TIV, data=imp_data) 
vifres <- as.data.frame(vif(test_vif))
write.csv(vifres, '/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vif/VRF_model.csv')

test_vif_cog=lm(exfunct ~ age_base + age_change + 
                  asinh_wml_base + asinh_wml_change + 
                  education + cesd + sex, data=imp_data) 
test_vif_cog <- as.data.frame(vif(test_vif_cog))
write.csv(test_vif_cog, '/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vif/Cog_model.csv')
