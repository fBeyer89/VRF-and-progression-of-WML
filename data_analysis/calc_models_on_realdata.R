library("lme4")
library("lmerTest") 
library("BayesFactor")
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer
library(broom.mixed)

setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/analysis") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('run_conf_LME.R')
source('run_exp_LME.R')
source('test_LME_assumptions.R')

#miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/imputed_data/imputed_data.Rdata")
#for testing: load smaller dataset:

####################################################
# Run confirmatory models
####################################################
#M1: effect of baseline DBP on WML change (also calculated E1a - E1c: effects of baseline WHR, change DBP and change WHR)
res_VRF= run_conf_LME(imp, "M1_VRF")

#M2: effect of WML change on executive function
res_exec= run_conf_LME(imp, "M2_exfunct")

#M3: effect of WML change on global cognition
res_globalcog= run_conf_LME(imp, "M3_globalcog")

####################################################
# Test assumptions
####################################################
#VRF
res_VRF= test_LME_assumptions(imp, "VRF")

#executive function
res_exec= test_LME_assumptions(imp, "exfunct")

#global cognition
res_globalcog= test_LME_assumptions(imp, "globalcog")

####################################################
# Test collinearity once for VRF and cognition models
# (only for one imputation)
####################################################
comp_imp=mice::complete(imp, "long")
imp_data=comp_imp[comp_imp$.imp==i,]

test_vif=lm(asinh_wml ~ age_base + age_change + 
              SBP_base + SBP_change + 
              WHR_base + WHR_change + 
              sex + BPmed + TIV, data=imp_data) 
vifres <- as.data.frame(vif(test_vif))
write.csv(vifres, '/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vif/VRF_model.csv')

test_vif_cog=lm(exfunct ~ age_base + age_change + 
                  asinh_wml_base + asinh_wml_change + 
                  education + cesd + sex, data=imp_data) 
test_vif_cog <- as.data.frame(vif(test_vif_cog))
write.csv(test_vif_cog, '/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vif/Cog_model.csv')

####################################################
# Run exploratory analyses
####################################################
#M1E2a: modifying effect of gender on age-related WML increase
res_VRF= run_exp_LME(imp, "E2a_age")

#M1E2b: modifying effect of gender on DBP-related WML increase
res_VRF= run_exp_LME(imp, "E2b_DBP")

#M1E2c: modifying effect of gender on WHR-related WML increase
res_VRF= run_exp_LME(imp, "E2c_DBP")

#M2E3a: modifying effect of gender on WML increase on executive function
res_exec= run_exp_LME(imp, "E3a_exfunct")

#M3E3b: modifying effect of gender on WML increase on global cognition
res_globalcog= run_exp_LME(imp, "E3b_globalcog")

####################################################
# Run sex-stratified analyses (if interaction not significant)
####################################################