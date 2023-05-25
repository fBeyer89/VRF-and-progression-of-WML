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
library(mice)

setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/run_conf_LME.R')

miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/imputed_data_08.5.23/imputed_data_08.5.23.Rdata")

res_exec= run_conf_LME(imp, "M2_exfunct",n_it=10)
write.csv(res_exec[[1]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Conf/M2_exfunct_model_res_freq.csv')
write.csv(res_exec[[2]], '/data/pt_life_whm/Results/VRF_cSVD/LME/results/Conf/M2_exfunct_model_res_bayes.csv')
