run_exp_LME<- function(imp, model, n_it){
  # Calculates exploratory models results for a MIDS (Multiply Imputed Data Set) from mice.
  # Model can take values "E2a_sex", "E2b_DBP", "E2c_WHR", "E3_exfunct" and "E3_globalcog".
  
  #########################
  # Frequentist statistics
  #########################
  # Pool using the mice standard workflow to collect results from 5 imputed datasets
  if (model == "E2a_sex"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + DBP_base:age_change + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    save(res,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + DBP_base:age_change + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                 REML=F, na.action = na.omit)
    save(plot,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_plot', '.RData'))
  }
  if (model == "E2a_sex_DBP_change"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_change:sex +
                                 DBP_base + DBP_base:age_change + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    save(res,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_change:sex +
                                 DBP_base + DBP_base:age_change + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                 REML=F, na.action = na.omit)
    save(plot,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_plot', '.RData'))
  }
  if (model == "E2b_DBP"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + age_change:DBP_base + 
                                 sex:DBP_base + sex:age_change:DBP_base +
                                 DBP_change  + 
                                 WHR_base + age_change:WHR_base + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)',
                REML=F, na.action = na.omit))
    est=summary(mice::pool(res))    
    save(res,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                   age_change:sex +
                                   DBP_base + age_change:DBP_base + 
                                   sex:DBP_base + sex:age_change:DBP_base +
                                   DBP_change  + 
                                   WHR_base + age_change:WHR_base + WHR_change + 
                                   sex + BPmed + TIV + (1|subj)',
                 REML=F, na.action = na.omit))
    save(plot,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_plot', '.RData'))
  }
  if (model == "E2c_WHR"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + age_change:DBP_base + 
                                 DBP_change  + 
                                 WHR_base + age_change:WHR_base + WHR_change + 
                                 sex:WHR_base + sex:age_change:WHR_base +
                                 sex + BPmed + TIV + (1|subj)', 
                REML=F, na.action = na.omit))
    est=summary(mice::pool(res))
    save(res,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + age_change:DBP_base + 
                                 DBP_change  + 
                                 WHR_base + age_change:WHR_base + WHR_change + 
                                 sex:WHR_base + sex:age_change:WHR_base +
                                 sex + BPmed + TIV + (1|subj)',
                                 REML=F, na.action = na.omit))
    save(plot,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_plot', '.RData'))
  }
  if (model == "E3a_exfunct"){
    res <- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_change:sex +
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)',
                REML=F, na.action = na.omit))
    est=summary(mice::pool(res))
    save(res,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_change:sex +
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)',
                                 REML=F, na.action = na.omit))
    save(plot,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_plot', '.RData'))
  }
  if (model == "E3b_globalcog"){
    res <- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_change:sex + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)',
                   REML=F, na.action = na.omit))
  est=summary(mice::pool(res))
  save(res,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_res', '.RData'))
  #recalculate model with lme4 as to use effect for effect prediction
  plot <- with(imp, lme4::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_change:sex + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)',
                               REML=F, na.action = na.omit))
  save(plot,file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_freq_imp_plot', '.RData'))
  }
  
  
  
  #########################
  # Bayesian statistics
  #########################
  #generalTestBF drops factors one by one and returns BF of the reduced model compared to null model
  #(Order of dropped coefficients: age_change:WHR_base, age_change:DBP_base, WHR_change, DBP_change)
  #Thus to obtain the likelihood of the full compared to null model, the respective Bayes Factors have
  #to be inverted. This is done in when calculating the mean BF.
  
  #Pool by calculating mean and standard deviation of BFs across 5 imputations
  comp_imp=mice::complete(imp, "long")  
  if (model == "E2a_sex"){
    for (i in c(1:imp$m)){
      
      
      tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                age_change:sex +
                                                                DBP_base + age_change:DBP_base + DBP_change +
                                                                WHR_base + WHR_base:age_change + WHR_change +
                                                                sex + BPmed + TIV + subj"), 
                              data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                              multicore = T,
                              neverExclude = c("age_base", "^age_change$", "DBP_base", "DBP_change", "WHR_change",
                                               "WHR_base", "BPmed", "TIV", "subj")
      )
      
      full=2 #model 2 is full model; 
      red=1 #age_change:sex dropped

      bf_etmp=extractBF(tmp_bf,logbf = F)
      bf_extracted=data.frame(pred=c("age_change:sex"), imp_n=rep(i,1))
      bf_extracted[,"bf"] <- c(bf_etmp[full,1] / bf_etmp[red,1])
      
      chains <- posterior(tmp_bf, full, iterations = n_it, columnFilter="^subj$")
      #We expect a stronger effect in women, e.g. a negative effect of interaction of sex (males coded as 1) and age change
      bf_extracted[,"siding"] <- c(mean(chains[,"sex:age_change-age_change"]<0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5        
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_imp_',i, '.RData'))
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_osbf=sd(one_sided_bf), mean_bf = mean(bf), mean_siding=mean(siding))
  }
  if (model == "E2a_sex_DBP_change"){
    for (i in c(1:imp$m)){
      
      
      tmp_bf <- generalTestBF(formula = as.formula('asinh_wml ~ age_base + age_change +
                                                                DBP_change:sex +
                                                                DBP_base + age_change:DBP_base + DBP_change +
                                                                WHR_base + WHR_base:age_change + WHR_change +
                                                                sex + BPmed + TIV + subj'), 
                              data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                              multicore = T,
                              neverExclude = c("age_base", "^age_change$", "DBP_base", "WHR_change",
                                               "WHR_base", "BPmed", "TIV", "subj")
      )
      
      full=4 #model 2 is full model; 
      red=3 #DBP_change:sex dropped
      
      bf_etmp=extractBF(tmp_bf,logbf = F)
      bf_extracted=data.frame(pred=c("DBP_change:sex"), imp_n=rep(i,1))
      bf_extracted[,"bf"] <- c(bf_etmp[full,1] / bf_etmp[red,1])
      
      chains <- posterior(tmp_bf, full, iterations = n_it, columnFilter="^subj$")
      #We expect a stronger effect in women, e.g. a negative effect of interaction of sex (males coded as 1) and age change
      bf_extracted[,"siding"] <- c(mean(chains[,"DBP_change:sex-sex"]<0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5        
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_imp_',i, '.RData'))
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_osbf=sd(one_sided_bf), mean_bf = mean(bf), mean_siding=mean(siding))
  }
  if (model == "E2b_DBP"){
    for (i in c(1:imp$m)){
      
      
      tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                age_change:sex +
                                                                DBP_base + age_change:DBP_base + DBP_change +
                                                                sex:DBP_base + sex:age_change:DBP_base + 
                                                                WHR_base + WHR_base:age_change + WHR_change +
                                                                sex + BPmed + TIV + subj"), 
                              data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                              multicore = T, whichModels = "top",
                              neverExclude = c("age_base", "^age_change$", "^DBP_base", 
                                               "WHR_base","WHR_change", "^sex", "DBP_change", "BPmed", "TIV", "subj") 
      )
      
      bf_extracted=extractBF(tmp_bf,logbf = F)
      bf_extracted=bf_extracted[1,] #only keep first row
      bf_extracted$pred=c("age_change:DBP_base:sex")
      
      
      #calculate only the full model for deriving siding factors
      tmp_bf_chains <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                age_change:sex +
                                                                DBP_base + age_change:DBP_base + DBP_change +
                                                                sex:DBP_base + sex:age_change:DBP_base + 
                                                                WHR_base + WHR_base:age_change + WHR_change +
                                                                sex + BPmed + TIV + subj"),
                                     data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  multicore = T,
                                     neverExclude = c("age_base", "age_change", "DBP_base", "DBP_change",
                                                      "WHR_base", "WHR_change", "sex", "BPmed", "TIV", "subj")
      )
      
      #We expect "in women DBP has a stronger effect than in men". 
      #This translated into a negative coefficient for men (women coded as 0/men as 1) 
      #-> attenuated effect of higher DBPxage_change on WML progression
      siding_factor=vector()
      chains <- posterior(tmp_bf_chains,  iterations = n_it, columnFilter="^subj$")
      bf_extracted[1,"siding"]=mean(chains[,"age_change:sex:DBP_base-DBP_base"]<0)
      
      
      #Because the Bayes factor from WhichModels="top"
      #is inverted (i.e. quantifies evidence for reduced model compared to full model),
      #we have to invert it back and multiply by the siding factor to quantifying evidence
      #for the alternative hypotheses/full model
      bf_extracted[, "bf"] <- 1/bf_extracted[, "bf"]
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5 
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_imp_',i, '.RData'))
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
  
  }
  if (model == "E2c_WHR"){
    for (i in c(1:imp$m)){
      
      
      tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                age_change:sex +
                                                                DBP_base + age_change:DBP_base + DBP_change +
                                                                sex:WHR_base + sex:age_change:WHR_base + 
                                                                WHR_base + WHR_base:age_change + WHR_change +
                                                                sex + BPmed + TIV + subj"), 
                              data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                              multicore = T, whichModels = "top",
                              neverExclude = c("age_base", "^age_change$", "DBP_base", 
                                               "^WHR_base","WHR_change", "^sex", "DBP_change", "BPmed", "TIV", "subj") 
      )
      
      bf_extracted=extractBF(tmp_bf,logbf = F)
      bf_extracted=bf_extracted[1,] #only keep first row
      bf_extracted$pred=c("age_change:sex:WHR_base")
      
      
      #calculate only the full model for deriving siding factors
      tmp_bf_chains <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                age_change:sex +
                                                                DBP_base + age_change:DBP_base + DBP_change +
                                                                sex:WHR_base + sex:age_change:WHR_base + 
                                                                WHR_base + WHR_base:age_change + WHR_change +
                                                                sex + BPmed + TIV + subj"),
                                     data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  multicore = T,
                                     neverExclude = c("age_base", "age_change", "DBP_base", "DBP_change",
                                                      "WHR_base", "WHR_change", "sex", "BPmed", "TIV", "subj")
      )
      
      #We expect "in women WHR has a stronger effect than in men". 
      #This translated into a negative coefficient for men (women coded as 0/men as 1) 
      #-> attenuated effect of higher DBPxage_change on WML progression
      chains <- posterior(tmp_bf_chains,  iterations = n_it, columnFilter="^subj$")
      bf_extracted[1,"siding"]=mean(chains[,"age_change:sex:WHR_base-WHR_base"]<0)
      
      
      #Because the Bayes factor from WhichModels="top"
      #is inverted (i.e. quantifies evidence for reduced model compared to full model),
      #we have to invert it back and multiply by the siding factor to quantifying evidence
      #for the alternative hypotheses/full model
      bf_extracted[, "bf"] <- 1/bf_extracted[, "bf"]
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5 
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_imp_',i, '.RData'))
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
     
  }
  if (model == "E3a_exfunct"){
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      modeldat=tmp[!is.na(tmp$exfunct),] #remove timepoints with NA in dependent variable
      tmp_bf <- generalTestBF(formula = as.formula("exfunct ~ age_base + age_change + 
                                                              asinh_wml_change:sex +
                                                              asinh_wml_base + asinh_wml_change + 
                                                              sex + cesd + TIV + education+ subj"), 
                              data=modeldat, whichRandom = "subj", 
                              multicore = T, whichModels = "top",
                              neverExclude = c("age_base", "^age_change$",
                                               "education", "asinh_wml_base", "cesd", "education",
                                               "TIV", "subj")
      )
      
      
      bf_extracted=extractBF(tmp_bf,logbf = F)
      bf_extracted=bf_extracted[1,] #only keep first row
      bf_extracted$pred=c("asinh_wml_change:sex")
      
      
      tmp_bf_chains <- generalTestBF(formula = as.formula("exfunct ~ age_base + age_change + 
                                                              asinh_wml_change:sex +
                                                              asinh_wml_base + asinh_wml_change + 
                                                              sex + cesd + TIV + education+ subj"), 
                              data=modeldat, whichRandom = "subj", 
                              multicore = T, 
                              neverExclude = c("age_base", "age_change", "sex",
                                               "education", "asinh_wml_base", "cesd", "education",
                                               "TIV", "subj")
      )
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      chains <- posterior(tmp_bf_chains, iterations = n_it, columnFilter="^subj$")
      #We expect a stronger effect in men (so additional negative effect in men coded as 1)...
      bf_extracted[,"siding"] <- c(mean(chains[,"asinh_wml_change:sex-sex"]<0))
      
      bf_extracted[, "bf"] <- 1/bf_extracted[, "bf"]
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5 
 
      
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_imp_',i, '.RData'))
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
  }
  if (model == "E3b_globalcog"){
    for (i in c(1:imp$m)){
      tmp=comp_imp[comp_imp$.imp==i,]
      modeldat=tmp[!is.na(tmp$exfunct),] #remove timepoints with NA in dependent variable
      
      tmp_bf <- generalTestBF(formula = as.formula("globalcog ~ age_base + age_change + 
                                                              asinh_wml_change:sex +
                                                              asinh_wml_base + asinh_wml_change + 
                                                              sex + cesd + TIV + education+ subj"), 
                              data=modeldat, whichRandom = "subj", 
                              multicore = T, whichModels = "top",
                              neverExclude = c("age_base", "^age_change$",
                                               "education", "asinh_wml_base", "cesd", "education",
                                               "TIV", "subj")
      )
      
      
      bf_extracted=extractBF(tmp_bf,logbf = F)
      bf_extracted=bf_extracted[1,] #only keep first row
      bf_extracted$pred=c("asinh_wml_change:sex")
      
      
      tmp_bf_chains <- generalTestBF(formula = as.formula("globalcog ~ age_base + age_change + 
                                                              asinh_wml_change:sex +
                                                              asinh_wml_base + asinh_wml_change + 
                                                              sex + cesd + TIV + education+ subj"), 
                                     data=modeldat, whichRandom = "subj", 
                                     multicore = T, 
                                     neverExclude = c("age_base", "age_change", "sex",
                                                      "education", "asinh_wml_base", "cesd", "education",
                                                      "TIV", "subj")
      )
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      chains <- posterior(tmp_bf_chains, iterations = n_it, columnFilter="^subj$")
      #We expect a stronger effect in men (so additional negative effect in men coded as 1)...
      bf_extracted[,"siding"] <- c(mean(chains[,"asinh_wml_change:sex-sex"]<0))
      
      bf_extracted[, "bf"] <- 1/bf_extracted[, "bf"]
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5 
      
      
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0('/data/pt_life_whm/Results/VRF_cSVD/LME/results/Expl/workspace_',model, '_imp_',i, '.RData'))
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
  }
  
  return(list(est, bf_res))
}