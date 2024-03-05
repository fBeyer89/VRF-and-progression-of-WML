test_LME_assumptions<- function(imp, model, outdir="/data/pt_life_whm/Results/VRF_cSVD/LME/assumptions/"){
  # Checks assumptions for LME, assumes to be run from data_analysis

  source("wald_confint_rlmm.R")
  source("helper/glmm_stability.r")
  library(influence.ME)
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
      cm=check_model(res)
      #ggsave(cm, filename = '/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures/orig_asinWMH.png', width = 15, height = 20, unit="cm", dpi=300)
      conf=confint(res)
      write.csv(conf, file = paste0(outdir, "/model_M1_VRF/vrf_conf_",model,"_imp_",i,".csv"),
                row.names = F)
      ##check_model from performance package is much nicer!!
      print("run check:model")
      jpeg(file=paste0(outdir, "model_M1_VRF/check_model_",model,"_imp_",i,".png"))
      check_model(res)
      dev.off()
      
      jpeg(file=paste0(outdir, "model_M1_VRF/qqp_vrf_",model,"_imp_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()
      
      jpeg(file=paste0(outdir, "model_M1_VRF/homoscedasticity_vrf_",model,"_imp_",i,".png"))
      homoscedas=plot(fitted(res), resid(res))
      dev.off()
      
      r=ranef(res)
      jpeg(file=paste0(outdir, "model_M1_VRF/randomeff_vrf_",model,"_imp_",i,".png"))
      hist(r$subj[,1])
      dev.off()
      
      # #Influential cases
      print("check for influential cases")
      print(nrow(tmp))

      infl=influence.ME::influence(res, group = "subj", count=TRUE, data=tmp)
      cd=as.data.frame(cooks.distance(infl))
      cd$subj=unique(tmp$subj)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"subj"]
      write.csv(infl_cases, file = paste0(outdir, "model_M1_VRF/vrf_infl_cases_",model,"_imp_",i,".csv"),
                row.names = F)
      tmp_woi <- tmp %>% filter(!subj %in% infl_cases)
      res <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change +
                                   DBP_base + age_change:DBP_base + DBP_change  +
                                   WHR_base + WHR_base:age_change + WHR_change +
                                   sex + BPmed + TIV + (1|subj)',
                            data=tmp_woi, REML=F, na.action = na.omit)
      sum=summary(res)
      coeff=coefficients(sum)
      write.csv(coeff, file = paste0(outdir, "model_M1_VRF/vrf_coeffs_wo_infl_",model,"imp_",i,".csv"))

      # determine stability of LME by excluding levels of random effects, one at a time
      # for further information please see:
      # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r

      res_lme4 <- lme4::lmer(formula = 'asinh_wml ~ age_base + age_change +
                                   DBP_base + age_change:DBP_base + DBP_change  +
                                   WHR_base + WHR_base:age_change + WHR_change +
                                   sex + BPmed + TIV + (1|subj)',
                            data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = paste0(outdir, "model_M1_VRF/vrf_stabmodel_",model,"imp_",i,".csv"))

      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('asinh_wml ~ age_base + age_change + 
                                        DBP_base + age_change:DBP_base + DBP_change  + 
                                        WHR_base + WHR_base:age_change + WHR_change + 
                                        sex + BPmed + TIV + (1|subj)', 
                                      data=tmp)
      sum=summary(res_robust)
      coeffs=sum$coefficients
      conf=confint.rlmerMod(res_robust)
      write.csv(cbind(coeffs,conf), file = paste0(outdir, "model_M1_VRF/vrf_robustmodel_",model,"imp_",i,".csv"))
    }
    }
  if (model == "M2_exfunct"){
    print("running model exfunct")
    print(imp$m)
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      print(tmp)
      res <- lme4::lmer(formula = 'exfunct ~ age_base + age_change +
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                            data=tmp, REML=F, na.action = na.omit)
      check_model(res)
      ggsave(filename = '/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures/orig_exec.png', width = 15, height = 20, unit="cm", dpi=300)
      
      conf=confint(res)
      write.csv(conf, file = paste0(outdir, "model_", model, "/",model,"_imp_",i,".csv"),
                row.names = F)
      ##check_model from performance package is much nicer!!
      png(file=paste0(outdir, "model_M2_exfunct/check_model_",model,"_imp_",i,".png"))
      check_model(res)
      dev.off()

      jpeg(file=paste0(outdir, "model_M2_exfunct/qqp_exfunct",model,"_imp_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()

      jpeg(file=paste0(outdir, "model_M2_exfunct/homoscedasticity_",model,"_imp_",i,".png"))
      homoscedas=plot(fitted(res), resid(res))
      dev.off()

      r=ranef(res)
      jpeg(file=paste0(outdir, "model_M2_exfunct/randomeff_",model,"_imp_",i,".png"))
      hist(r$subj[,1])
      dev.off()
      # 
      # #Influential cases
      infl=influence.ME::influence(res, group = "subj")
      cd=as.data.frame(cooks.distance(infl))
      cd$subj=unique(tmp$subj)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"subj"]
      write.csv(infl_cases, file = paste0(outdir, "model_M2_exfunct/infl_cases_",model,"_imp_",i,".csv"),
                row.names = F)
      tmp_woi <- tmp %>% filter(!subj %in% infl_cases)
      res <- lmerTest::lmer(formula = 'exfunct ~ age_base + age_change +
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                            data=tmp_woi, REML=F, na.action = na.omit)
      sum=summary(res)
      coeff=coefficients(sum)
      write.csv(coeff, file = paste0(outdir, "model_M2_exfunct/coeffs_wo_infl_",model,"_imp_",i,".csv"))

      # # determine stability of LME by excluding levels of random effects, one at a time
      # # for further information please see: 
      # # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r

      res_lme4 <- lme4::lmer(formula = 'exfunct ~ age_base + age_change +
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                              data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = paste0(outdir, "model_M2_exfunct/stabmodel_",model,"_imp_",i,".csv"))

      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('exfunct ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)', 
                                   data=tmp)
      sum=summary(res_robust)
      coeffs=sum$coefficients
      conf=confint.rlmerMod(res_robust)
      write.csv(cbind(coeffs,conf), file = paste0(outdir, "model_M2_exfunct/robustmodel_",model,"_imp_",i,".csv"))
    }
  }
  else{
    print("running model glob")
    print(imp$m)
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      print(nrow(tmp))
      res <- lme4::lmer(formula = 'globalcog ~ age_base + age_change +
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                            data=tmp, REML=F, na.action = na.omit)
      check_model(res)
      ggsave(filename = '/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures/orig_gc.png', width = 15, height = 20, unit="cm", dpi=300)
      
      conf=confint(res)
      write.csv(conf, file = paste0(outdir, "model_", model, "/",model,"_imp_",i,".csv"),
                row.names = F)
      ##check_model from performance package is much nicer!!
      jpeg(file=paste0(outdir, "model_M3_globalcog/check_model_",model,"_imp_",i,".png"))
      check_model(res)
      dev.off()

      jpeg(file=paste0(outdir, "model_M3_globalcog/qqp_",model,"_imp_",i,".png"))
      qqp=car::qqPlot(resid(res))
      dev.off()

      jpeg(file=paste0(outdir, "model_M3_globalcog/homoscedasticity_",model,"_imp_",i,".png"))
      homoscedas=plot(fitted(res), resid(res))
      dev.off()

      r=ranef(res)
      jpeg(file=paste0(outdir, "model_M3_globalcog/randomeff_",model,"_imp_",i,".png"))
      hist(r$subj[,1])
      dev.off()

      # #Influential cases
      infl=influence.ME::influence(res, group = "subj")
      cd=as.data.frame(cooks.distance(infl))
      cd$subj=unique(tmp$subj)
      infl_cases= cd[cd$V1>(3*sd(cd$V1)+mean(cd$V1)),"subj"]
      write.csv(infl_cases, file = paste0(outdir, "model_M3_globalcog/infl_cases_",model,"_imp_",i,".csv"),
                row.names = F)
      tmp_woi <- tmp %>% filter(!subj %in% infl_cases)
      res <- lmerTest::lmer(formula = 'globalcog ~ age_base + age_change +
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                            data=tmp_woi, REML=F, na.action = na.omit)
      sum=summary(res)
      coeff=coefficients(sum)
      write.csv(coeff, file = paste0(outdir, "model_M3_globalcog/coeffs_wo_infl_",model,"_imp_",i,".csv"))

      # # determine stability of LME by excluding levels of random effects, one at a time
      # # for further information please see: 
      # # https://github.com/keyfm/eva/blob/master/trpm8/src/glmm_stability.r
      # 
      res_lme4 <- lme4::lmer(formula = 'globalcog ~ age_base + age_change +
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)',
                              data=tmp, REML=F, na.action = na.omit)
      stab_results <- glmm.model.stab(res_lme4)
      write.csv(stab_results$summary, file = paste0(outdir, "model_M3_globalcog/stabmodel_",model,"_imp_",i,".csv"))

      #Implement robust LMM if necessary
      res_robust= robustlmm::rlmer('globalcog ~ age_base + age_change + 
                                 asinh_wml_base  + asinh_wml_change +
                                 sex + education + cesd + TIV + (1|subj)', 
                                   data=tmp)
      sum=summary(res_robust)
      coeffs=sum$coefficients
      conf=confint.rlmerMod(res_robust)
      write.csv(cbind(coeffs,conf), file = paste0(outdir, "model_M3_globalcog/robustmodel_",model,"_imp_",i,".csv"))
    }
  }
}