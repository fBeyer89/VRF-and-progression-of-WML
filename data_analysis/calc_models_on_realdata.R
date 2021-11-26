library("lme4")
library("lmerTest") 
library("BayesFactor")
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer
library(broom.mixed)

setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/analysis") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('run_LME_realdata.R')
source('test_LME_assumptions.R')

miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/imputed_data/imputed_data.Rdata")
####################################################
# Run models
####################################################
#VRF
res_VRF= run_LME(imp, "VRF")

#executive function
res_exec= run_LME(imp, "exfunct")

#global cognition
res_globalcog= run_LME(imp, "globalcog")

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
####################################################
imp_data=mice::complete(imp, "long")

test_vif=lm(asinh_wml ~ age_base + age_change + 
              SBP_base + SBP_change + 
              #WHR_base + WHR_change + 
              sex + BPmed + TIV, data=imp_data) 
vifres <- as.data.frame(vif(test_vif))
write.csv(vifres, '/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vif/VRF_model.csv')

test_vif_cog=lm(exfunct ~ age_base + age_change + 
                  asinh_wml_base + asinh_wml_change + 
                  education + cesd + sex, data=imp_data) 
test_vif_cog <- as.data.frame(vif(test_vif_cog))
write.csv(test_vif_cog, '/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vif/Cog_model.csv')
