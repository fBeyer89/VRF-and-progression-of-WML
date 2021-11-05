test_LME_assumptions<- function(imp, model){
  # Checks assumptions for LME 
  library(robustlmm)
  library(influence.ME)
  source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
  #########################
  # Frequentist statistics
  #########################
  # Test assumptions for LME
  # normality and homoscedasticity of residuals
  # normality of random effects
  # influential cases
  # stability of LME
  
  # Save separately for each inputed dataset.
  comp_imp=mice::complete(imp, "long") 
  
  if (model == "VRF"){
    for (i in c(1:imp$m)){
    tmp=comp_imp[comp_imp$.imp==i,]
    res <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|sic)',
                data=tmp, REML=F, na.action = na.omit)
    
    jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/qqp_vrf",i,".png"))
    qqp=car::qqPlot(resid(res))
    dev.off()
    
    jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/homoscedasticity_vrf",i,".png"))
    homoscedas=plot(fitted(test), resid(test))
    dev.off()
    
    r=ranef(test)
    jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/randomeff_vrf",i,".png"))
    hist(r$subj[,1])
    dev.off()
    
    #Influential cases
    infl=influence.ME::influence(res, group = "mrt_pseudonym")
    cd=as.data.frame(cooks.distance(infl))
    cd$mrt_pseudonym=row.names(cd)
    infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"mrt_pseudonym"]
    
    tmp <- filter(tmp, (mrt_pseudonym %in% infl_cases))
    res <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|sic)',
                          data=tmp, REML=F, na.action = na.omit)
    coeff=coefficients(res)
    write.csv(coeff, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vrf_coeffs_wo_infl_",i,".csv"))
    
    # determine stability of LME by excluding levels of random effects, one at a time
    # for further information please see: 
    # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
    
    res_lme4 <- lmer4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|sic)',
                          data=tmp, REML=F, na.action = na.omit)
    stab_results <- glmm.model.stab(res_lme4)
    write.csv(stab_results$summary, file = "/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/vrf_model_",i,".csv") 
    
   #Implement robust LMM if necessary
   res_robust= robustlmm::rlmer('asinh_wml ~ age_base + age_change + 
                                      SBP_base + age_change:SBP_base + SBP_change  + 
                                      WHR_base + WHR_base:age_change + WHR_change + 
                                      sex + BPmed + TIV + (1|sic)', 
                                    data=tmp)
    }
    }
  if (model == "exfunct"){
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      res <- lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)',
                            data=tmp, REML=F, na.action = na.omit)
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/qqp_exfunct",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/homoscedasticity_exfunct",i,".png"))
      homoscedas=plot(fitted(test), resid(test))
      dev.off()
      
      r=ranef(test)
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/randomeff_exfunct",i,".png"))
      hist(r$subj[,1])
      dev.off()
      
      #Influential cases
      infl=influence.ME::influence(res, group = "mrt_pseudonym")
      cd=as.data.frame(cooks.distance(infl))
      cd$mrt_pseudonym=row.names(cd)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"mrt_pseudonym"]
      
      tmp <- filter(tmp, (mrt_pseudonym %in% infl_cases))
      res <- lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)',
                            data=tmp, REML=F, na.action = na.omit)
      coeff=coefficients(res)
      write.csv(coeff, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/exfunct_coeffs_wo_infl_",i,".csv"))
      
      # determine stability of LME by excluding levels of random effects, one at a time
      # for further information please see: 
      # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
      
      res_lme4 <- lmer4::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)',
                              data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = "/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/exfunct_model.csv") 
      
      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)', 
                                   data=tmp)
    }
  }
  else{
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      res <- lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)',
                            data=tmp, REML=F, na.action = na.omit)
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/qqp_globalcog_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()
      
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/homoscedasticity_globalcog_",i,".png"))
      homoscedas=plot(fitted(test), resid(test))
      dev.off()
      
      r=ranef(test)
      jpeg(file=paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/randomeff_globalcog_",i,".png"))
      hist(r$subj[,1])
      dev.off()
      
      #Influential cases
      infl=influence.ME::influence(res, group = "mrt_pseudonym")
      cd=as.data.frame(cooks.distance(infl))
      cd$mrt_pseudonym=row.names(cd)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"mrt_pseudonym"]
      
      tmp <- filter(tmp, (mrt_pseudonym %in% infl_cases))
      res <- lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)',
                            data=tmp, REML=F, na.action = na.omit)
      coeff=coefficients(res)
      write.csv(coeff, file = paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/globalcog_coeffs_wo_infl_",i,".csv"))
      
      # determine stability of LME by excluding levels of random effects, one at a time
      # for further information please see: 
      # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
      
      res_lme4 <- lmer4::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)',
                              data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = "/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/globalcog_model.csv") 
      
      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|sic)', 
                                   data=tmp)
    }
  }
}