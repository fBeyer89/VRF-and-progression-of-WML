library("lme4")
library("lmerTest") 
library("BayesFactor")
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer
library(broom.mixed)
library(dplyr)
setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/")

n_it=10

dat=read.csv("../power_simulation/simulated_data.csv")
colnames(dat)=c("subj", "tp", "age_base", "sex", "DBP_base", "WHR_base", "TIV", "age_change", "DBP_change", "WHR_change", "asinh_wml", "education", "cesd", "exfunct", "globalcog")

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
dat= dat%>%
  mutate(across(education:globalcog, scale2))%>%
  mutate_at("TIV", scale2)%>%
  mutate(asinh_wml_base = ifelse(tp == "bl", asinh_wml, lag(asinh_wml)))%>%
  mutate(asinh_wml_change = ifelse(tp== "bl", 0, asinh_wml - asinh_wml_base))

dat$BPmed=rbinom(n=nrow(dat),size=1,prob=0.2)
dat$subj=as.factor(dat$subj)

#introduce random NAs and impute
dat[floor(runif(10,min=1,max=nrow(dat))),"education"]=NA
dat[floor(runif(10,min=1,max=nrow(dat))),"cesd"]=NA

mising_pattern=mice::md.pattern(dat)

m = 5 
iter = 10
imp <- mice::mice(
  dat,
  m = m,
  seed = 8745
)
imp_l <- mice::complete(imp, "long", include = TRUE)
imp.itt <- mice::as.mids(imp_l)

#Test the models for main hypothesis H1 & E1a - E1c
source("run_conf_LME.R")
res=run_conf_LME(imp.itt, "M1_VRF",n_it=10)

res2=run_conf_LME(imp.itt, "M2_exfunct",100)

res3=run_conf_LME(imp.itt, "M3_globalcog",100)

#Test the models for exploratory hypothesis E2a_sex, E2a_DBP, E2a_WHR & E3a - E3c
source("run_exp_LME.R")
res=run_exp_LME(imp.itt, "E2a_sex",10)

res2=run_exp_LME(imp.itt, "E2b_DBP",100)
res3=run_exp_LME(imp.itt, "E2c_WHR",100)

res4=run_exp_LME(imp.itt, "E3a_exfunct",100)
res5=run_exp_LME(imp.itt, "E3b_globalcog",100)


#Test assumptions
source("test_LME_assumptions.R")
res_VRF= test_LME_assumptions(imp=imp.itt, "VRF")
res_VRF= test_LME_assumptions(imp=imp.itt, "exfunct")
res_VRF= test_LME_assumptions(imp=imp.itt, "globalcog")

##code snippets for trying parts of it
comp_imp=mice::complete(imp.itt, "long") 


for (i in c(1:imp.itt$m)){
  
  
  tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                 DBP_base + age_change:DBP_base + DBP_change +
                                                                 WHR_base + WHR_base:age_change + WHR_change +
                                                                 sex + BPmed + TIV + subj"), 
                          data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                          multicore = T,
                          neverExclude = c("age_base", "^age_change$", "^DBP_base$", 
                                           "^WHR_base$", "sex", "BPmed", "TIV", "subj")
  )
  #model 15 is full model; 
  #model 7: age_change:WHR_base dropped
  #model 11: age_change:DBP_base dropped, 
  #model 13: WHR change dropped
  #model 14: DBP change dropped
  
  
  bf_etmp=extractBF(tmp_bf,logbf = F)
  bf_extracted=data.frame(pred=c("age_change:WHR_base", "age_change:DBP_base",
                                 "WHR_change", "DBP_change"), imp_n=rep(i,4))
  bf_extracted[,"bf"] <- c(bf_etmp[15,1] / bf_etmp[7,1], bf_etmp[15,1] / bf_etmp[11,1],
                           bf_etmp[15,1] / bf_etmp[13,1],bf_etmp[15,1] / bf_etmp[14,1])
  
  chains <- posterior(tmp_bf, 15, iterations = n_it, columnFilter="^id$")#The fifteenth model is the full model with all 10 terms.
  #We expect a positive effect of interaction of baseline DBP and WHR with age change
  #and of change in DBP/WHR
  bf_extracted[,"bf_sf"] <- c(mean(chains[,"age_change.&.WHR_base"]>0), mean(chains[,"age_change.&.DBP_base"]>0),
                              mean(chains[,"WHR_change"]>0),mean(chains[,"DBP_change"]>0))
  
  bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5        
  
  
  
  if (i==1){
    bfall=bf_extracted}
  else{
    bfall=rbind(bfall,bf_extracted)
  }}


  
res <- with(imp.itt, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
            REML=F, na.action = na.omit)
  
  
  #res_list=run_LME('asinh_ll', dat, 'simple', n_it)