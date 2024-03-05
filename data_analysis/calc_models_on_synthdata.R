library("lme4")
library("lmerTest") 
library("BayesFactor")
library(DiagrammeR)
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer
library(broom.mixed)
library(kableExtra)
library(tidyverse)
library(modelsummary)
library(robustlmm)
library(performance)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)
library(mice)
library(ggmice)
library(flextable)
library(ragg)
library(car)
library(influence.ME)


setwd("data_analysis/") #with path to github repository as origin (./)
source("helper/diagnostics_fcns.r")
source("helper/glmm_stability.r")
source('run_conf_LME.R')
source('run_exp_LME.R')
source('test_LME_assumptions.R')

miceadds::load.Rdata(objname = "imp", "synthetic_imputed_data/synthetic_imputed_data.Rdata")

# ####################################################
# # Define directories to store outputs
# ####################################################
outdir_conf='/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/test_synthetic/' #put a file location of your liking here
outdir_expl='/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/test_synthetic/' #put a file location of your liking here
outdir_assumptions="/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/test_synthetic/"
impdata=
# ####################################################
# # Run confirmatory models
# ####################################################
# #M1: effect of baseline DBP:age_change on WML change (also calculated E1a - E1c: effects of baseline WHR:age change, change DBP and change WHR)
res_VRF= run_conf_LME(imp, "M1_VRF",n_it=1, outdir=outdir_conf)
write.csv(res_VRF[[1]], paste0(outdir_conf, 'M1_VRF_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_conf, 'M1_VRF_model_res_bayes.csv'))

# #M2: effect of WML change on executive function -> 
# exec function not imputed in synthetic dataset
#res_exec= run_conf_LME(imp, "M2_exfunct",n_it=10)
#write.csv(res_exec[[1]], paste0(outdir_conf,'/M2_exfunct_model_res_freq.csv'))
#write.csv(res_exec[[2]], paste0(outdir_conf,'/M2_exfunct_model_res_bayes.csv'))

# #M3: effect of WML change on global cognition
res_globalcog= run_conf_LME(imp, "M3_globalcog",n_it=10)
write.csv(res_globalcog[[1]], paste0(outdir_conf,'/M3_globalcog_model_res_freq.csv'))
write.csv(res_globalcog[[2]], paste0(outdir_conf,'/M3_globalcog_model_res_bayes.csv'))

# ####################################################
# # Run exploratory analyses
# ####################################################
# #M1E2a: modifying effect of gender on age-related WML increase
res_VRF= run_exp_LME(imp, "E2a_sex", n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/E2a_age_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/E2a_age_model_res_bayes.csv'))

#M1E2b: modifying effect of gender on DBP-related WML increase
res_VRF= run_exp_LME(imp, "E2b_DBP", n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/E2b_DBP_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/E2b_DBP_model_res_bayes.csv'))

#M1E2c: modifying effect of gender on WHR-related WML increase
res_VRF= run_exp_LME(imp, "E2c_WHR", n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/E2c_DBP_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/E2c_DBP_model_res_bayes.csv'))

#M2E3a: modifying effect of gender on WML increase on executive function ->
# exec func not imputed in synthetic dataset

#res_exec= run_exp_LME(imp, "E3a_exfunct", n_it=10)
#write.csv(res_VRF[[1]], paste0(outdir_expl,'/E3a_exfunct_model_res_freq.csv'))
#write.csv(res_VRF[[2]], paste0(outdir_expl,'/E3a_exfunct_model_res_bayes.csv'))

#M3E3b: modifying effect of gender on WML increase on global cognition
res_globalcog= run_exp_LME(imp, "E3b_globalcog", n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/E3b_globalcog_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/E3b_globalcog_model_res_bayes.csv'))

#M2E2a:E2a_sex_DBP_change
res_VRF= run_exp_LME(imp, "E2a_sex_DBP_change", n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/E2a_sex_DBP_change_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/E2a_sex_DBP_change_res_bayes.csv'))

####################################################
# Run sex-stratified analyses (if interaction not significant)
####################################################
source("run_conf_LME_sex_stratified.R")
#Select only females
imp_s0 <- miceadds::subset_datlist(imp,
                                   expr_subset=expression(sex==0),
                                   toclass="mids" )


res_VRF= run_conf_LME_sexstr(imp_s0, 0, model="M1_VRF",n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/M1_VRF_females_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/M1_VRF_females_model_res_bayes.csv'))

#res_M2= run_conf_LME_sexstr(imp_s0, 0, model="M2_exfunct",n_it=10)
#write.csv(res_M2[[1]], paste0(outdir_expl,'/M2_exfunct_females_model_res_freq.csv'))
#write.csv(res_M2[[2]], paste0(outdir_expl,'/M2_exfunct_females_model_res_bayes.csv'))

M3_globalcog= run_conf_LME_sexstr(imp_s0, 0, model="M3_globalcog",n_it=10)
write.csv(M3_globalcog[[1]], paste0(outdir_expl,'/M3_globalcog_females_model_res_freq.csv'))
write.csv(M3_globalcog[[2]], paste0(outdir_expl,'/M3_globalcog_females_model_res_bayes.csv'))


#Select only males
imp_s1 <- miceadds::subset_datlist(imp,
                                   expr_subset=expression(sex==1),
                                   toclass="mids" )
res_VRF= run_conf_LME_sexstr(imp_s1,1, model="M1_VRF",n_it=10)
write.csv(res_VRF[[1]], paste0(outdir_expl,'/M1_VRF_males_model_res_freq.csv'))
write.csv(res_VRF[[2]], paste0(outdir_expl,'/M1_VRF_males_model_res_bayes.csv'))

#M2_exfunct= run_conf_LME_sexstr(imp_s1, 1, model="M2_exfunct",n_it=10)
#write.csv(M2_exfunct[[1]], paste0(outdir_expl,'/M2_exfunct_males_model_res_freq.csv'))
#write.csv(M2_exfunct[[2]], paste0(outdir_expl,'/M2_exfunct_males_model_res_bayes.csv'))

M3_globalcog= run_conf_LME_sexstr(imp_s1, 1, model="M3_globalcog",n_it=10)
write.csv(M3_globalcog[[1]], paste0(outdir_expl,'/M3_globalcog_males_model_res_freq.csv'))
write.csv(M3_globalcog[[2]], paste0(outdir_expl,'/M3_globalcogt_males_model_res_bayes.csv'))

####################################################
# E4: SBP -> not imputed in synthetic dataset
####################################################
# miceadds::load.Rdata(objname = "imp", "imputed_data_SBP_6.12.23.Rdata")
# 
# 
# res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
#                                  SBP_base + age_change:SBP_base + SBP_change  + 
#                                  WHR_base + WHR_base:age_change + WHR_change + 
#                                  sex + BPmed + TIV + (1|subj)'),
#             REML=F, na.action = na.omit)
# est=summary(mice::pool(res))
# save(res,file=paste0(outdir_expl,'/SBP_model_res.Rdata'))
# write.csv(est,file=paste0(outdir_expl,'/SBP_model_res_freq.csv'))
# 
# #calculate p-values based on Wald nested models comparisons 
# #for DBP_baseline interaction with time
# res_SBP_baseline <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ 
#                                  age_base + age_change + 
#                                  SBP_base + SBP_change +
#                                  WHR_base + WHR_base:age_change + WHR_change + 
#                                  sex + BPmed + TIV + (1|subj)'),
#                          REML=F, na.action = na.omit)
# d1=D1(res,res_SBP_baseline)
# write.csv(d1$result,file=paste0(outdir_expl,'/SBP_model_res_d1.csv'))

####################################################
# Test assumptions & run robust models (saves everything automatically)
####################################################
#VRF
res_VRF= test_LME_assumptions(imp, "M1_VRF", outdir=outdir_assumptions)

#executive function
res_exec= test_LME_assumptions(imp, "M2_exfunct", outdir=outdir_assumptions)

#global cognition
res_globalcog= test_LME_assumptions(imp, "M3_globalcog", outdir=outdir_assumptions)

####################################################
# Test collinearity once for VRF and cognition models
# (done only for one imputation as only predictor variables matter)
####################################################
comp_imp=mice::complete(imp, "long")
imp_data=comp_imp[comp_imp$.imp==1,]

test_vif=lm(asinh_wml ~ age_base + age_change +
              DBP_base + DBP_change +
              WHR_base + WHR_change +
              sex + BPmed + TIV, data=imp_data)
vifres <- as.data.frame(vif(test_vif))
write.csv(vifres, paste0(outdir_assumptions, '/vif/VRF_model.csv'))

test_vif_cog=lm(exfunct ~ age_base + age_change +
                  asinh_wml_base + asinh_wml_change +
                  education + cesd + sex, data=imp_data)
test_vif_cog <- as.data.frame(vif(test_vif_cog))
write.csv(test_vif_cog, paste0(outdir_assumptions, '/vif/Cog_model.csv'))


####################################################
# Explore associations with components
# (done only for one imputation as only predictor variables matter)
####################################################


##### Some explorations during revision ######
# in models without age_change:WHR_base interaction, age_change is still significant.
res0 <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + TIV + (1|subj)'),
             REML=F, na.action = na.omit)
res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
            REML=F, na.action = na.omit)

res_base <- lm(formula = 'DBP_base ~ WHR_base', data=imp$data)
summary(mice::pool(res0))
est=summary(mice::pool(res))
