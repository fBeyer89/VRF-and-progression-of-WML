test_LME_assumptions<- function(imp, model){
  # Checks assumptions for LME 
  library(influence.ME)
  library(robustlmm)
  library(dplyr)
  library(mice)
  library(lme4)

  source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
  #########################
  # Frequentist statistics
  #########################
  # res assumptions for LME
  # normality and homoscedasticity of residuals
  # normality of random effects
  # influential cases
  # stability of LME
  
  # Save separately for each inputed dataset.
  comp_imp=mice::complete(imp, "long") 
  
  if (model == "M1_VRF"){
    print("running model VRF")
    
    for (i in c(1:imp$m)){
    
      tmp=comp_imp[comp_imp$.imp==i,]
      res <- lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                   DBP_base + age_change:DBP_base + DBP_change  + 
                                   WHR_base + WHR_base:age_change + WHR_change + 
                                   sex + BPmed + TIV + (1|subj)',
                  data=tmp, REML=F, na.action = na.omit)
      
      ##check_model from performance package is much nicer!!
      print("run check:model")
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/check_model_",model,"_imp_",i,".png"))
      check_model(res)
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/qqp_vrf_",model,"_imp_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/homoscedasticity_vrf_",model,"_imp_",i,".png"))
      homoscedas=plot(fitted(res), resid(res))
      dev.off()
      
      r=ranef(res)
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/randomeff_vrf_",model,"_imp_",i,".png"))
      hist(r$subj[,1])
      dev.off()
      
      #Influential cases
      print("check for influential cases")
      print(nrow(tmp))
      
      infl=influence.ME::influence(res, group = "subj", count=TRUE, data=tmp)
      cd=as.data.frame(cooks.distance(infl))
      cd$subj=unique(tmp$subj)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"subj"]
      write.csv(infl_cases, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vrf_infl_cases_",model,"_imp_",i,".csv"),
                row.names = F)
      tmp_woi <- tmp %>% filter(!subj %in% infl_cases)
      res <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                   DBP_base + age_change:DBP_base + DBP_change  + 
                                   WHR_base + WHR_base:age_change + WHR_change + 
                                   sex + BPmed + TIV + (1|subj)',
                            data=tmp_woi, REML=F, na.action = na.omit)
      sum=summary(res)
      coeff=coefficients(sum)
      write.csv(coeff, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vrf_coeffs_wo_infl_",model,"imp_",i,".csv"))
      
      # determine stability of LME by excluding levels of random effects, one at a time
      # for further information please see: 
      # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
      
      res_lme4 <- lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                   DBP_base + age_change:DBP_base + DBP_change  + 
                                   WHR_base + WHR_base:age_change + WHR_change + 
                                   sex + BPmed + TIV + (1|subj)',
                            data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vrf_stabmodel_",model,"imp_",i,".csv"))
      
      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('asinh_wml ~ age_base + age_change + 
                                        DBP_base + age_change:DBP_base + DBP_change  + 
                                        WHR_base + WHR_base:age_change + WHR_change + 
                                        sex + BPmed + TIV + (1|subj)', 
                                      data=tmp)
      sum=summary(res_robust)
      coeffs=sum$coefficients
      write.csv(coeffs, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vrf_robustmodel_",model,"imp_",i,".csv"))
    }
    }
  if (model == "M2_exfunct"){
    print("running model exfunct")
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      res <- lme4::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                            data=tmp, REML=F, na.action = na.omit)
      
      ##check_model from performance package is much nicer!!
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/check_model_",model,"_imp_",i,".png"))
      check_model(res)
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/qqp_exfunct",model,"_imp_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/homoscedasticity_",model,"_imp_",i,".png"))
      homoscedas=plot(fitted(res), resid(res))
      dev.off()
      
      r=ranef(res)
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/randomeff_",model,"_imp_",i,".png"))
      hist(r$subj[,1])
      dev.off()
      
      #Influential cases
      infl=influence.ME::influence(res, group = "subj")
      cd=as.data.frame(cooks.distance(infl))
      cd$subj=unique(tmp$subj)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"subj"]
      write.csv(infl_cases, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/infl_cases_",model,"_imp_",i,".csv"),
                row.names = F)
      tmp_woi <- tmp %>% filter(!subj %in% infl_cases)
      res <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)',
                            data=tmp_woi, REML=F, na.action = na.omit)
      sum=summary(res)
      coeff=coefficients(sum)
      write.csv(coeff, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/coeffs_wo_infl_",model,"_imp_",i,".csv"))
      
      # determine stability of LME by excluding levels of random effects, one at a time
      # for further information please see: 
      # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
      
      res_lme4 <- lme4::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                              data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/stabmodel_",model,"_imp_",i,".csv"))
      
      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)', 
                                   data=tmp)
      sum=summary(res_robust)
      coeffs=sum$coefficients
      write.csv(coeffs, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/robustmodel_",model,"_imp_",i,".csv"))
    }
  }
  else{
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      res <- lme4::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                            data=tmp, REML=F, na.action = na.omit)
      
      ##check_model from performance package is much nicer!!
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/check_model_",model,"_imp_",i,".png"))
      check_model(res)
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/qqp_",model,"_imp_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/homoscedasticity_",model,"_imp_",i,".png"))
      homoscedas=plot(fitted(res), resid(res))
      dev.off()
      
      r=ranef(res)
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/randomeff_",model,"_imp_",i,".png"))
      hist(r$subj[,1])
      dev.off()
      
      #Influential cases
      infl=influence.ME::influence(res, group = "subj")
      cd=as.data.frame(cooks.distance(infl))
      cd$subj=unique(tmp$subj)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"subj"]
      write.csv(infl_cases, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/infl_cases_",model,"_imp_",i,".csv"),
                row.names = F)
      tmp_woi <- tmp %>% filter(!subj %in% infl_cases)
      res <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)',
                            data=tmp_woi, REML=F, na.action = na.omit)
      sum=summary(res)
      coeff=coefficients(sum)
      write.csv(coeff, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/coeffs_wo_infl_",model,"_imp_",i,".csv"))
      
      # determine stability of LME by excluding levels of random effects, one at a time
      # for further information please see: 
      # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
      
      res_lme4 <- lme4::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                              data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/stabmodel_",model,"_imp_",i,".csv"))
      
      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)', 
                                   data=tmp)
      sum=summary(res_robust)
      coeffs=sum$coefficients
      write.csv(coeffs, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/robustmodel_",model,"_imp_",i,".csv"))
    }
  }
}