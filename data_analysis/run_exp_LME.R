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
    est=summary(mice::pool(res))}
  if (model == "E2b_DBP"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + age_change:DBP_base + 
                                 sex:DBP_base + sex:age_change:DBP_base +
                                 DBP_change  + 
                                 WHR_base + age_change:WHR_base + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    est=summary(mice::pool(res))}
  if (model == "E2c_WHR"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 age_change:sex +
                                 DBP_base + age_change:DBP_base + 
                                 DBP_change  + 
                                 WHR_base + age_change:WHR_base + WHR_change + 
                                 sex:WHR_base + sex:age_change:WHR_base +
                                 sex + BPmed + TIV + (1|subj)'), 
                REML=F, na.action = na.omit)
    est=summary(mice::pool(res))}
  if (model == "E3a_exfunct"){
    res <- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_change:sex +
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
  }
  if (model == "E3b_globalcog"){
    res <- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_change:sex + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                   REML=F, na.action = na.omit)
  est=summary(mice::pool(res))
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
      #We expect a negative effect of interaction of sex (females coded as 0) and age change
      bf_extracted[,"bf_sf"] <- c(mean(chains[,"sex.&.age_change"]<0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5        
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      mutate(mean_one_sided_bf = mean(one_sided_bf), sd_osbf=sd(one_sided_bf))%>%
      mutate(mean_bf = mean(bf), sd_bf=sd(bf))
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
                              neverExclude = c("age_base", "^age_change$", "DBP_base", 
                                               "^WHR_base$","WHR_change", "DBP_change", "BPmed", "TIV", "subj") 
      )
      
      full=21 #model  is full model; 
      red=20 #model: age_change:sex:DBP_base dropped
      
      bf_etmp=extractBF(tmp_bf,logbf = F)
      bf_extracted=data.frame(pred=c("sex:age_change:DBP_base"), imp_n=rep(i,1))
      bf_extracted[,"bf"] <- c(bf_etmp[3,1] / bf_etmp[2,1])
      
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      chains <- posterior(tmp_bf, full, iterations = n_it, columnFilter="^id$")
      #We expect "in women DBP has a stronger effect than in men". 
      #This translated into a negative coefficient for men (women coded as 0/men as 1) 
      #-> attenuated effect of higher DBPxage_change on WML progression
      bf_extracted[,"bf_sf"] <- c(mean(chains[,"sex.&.age_change.&.DBP_base"]<0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5        
      
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf))
  
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
                              multicore = T,
                              neverExclude = c("age_base", "^age_change$", "^DBP_base$", 
                                               "^WHR_base$","WHR_change", "DBP_change", 
                                               "education","BPmed", "TIV", "subj")  
      )
      
      full=21 #model  is full model; 
      red=19 #model: age_change:sex:WHR_base dropped
      
      bf_etmp=extractBF(tmp_bf,logbf = F)
      bf_extracted=data.frame(pred=c("sex:age_change:WHR_base"), imp_n=rep(i,1))
      bf_extracted[,"bf"] <- c(bf_etmp[3,1] / bf_etmp[2,1])
      
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      chains <- posterior(tmp_bf, full, iterations = n_it, columnFilter="^id$")
      #We expect "in women WHR has a stronger effect than in men". 
      #This translated into a negative coefficient for men (women coded as 0/men as 1) 
      #-> attenuated effect of higher WHRxage_change on WML progression
      bf_extracted[,"bf_sf"] <- c(mean(chains[,"sex.&.age_change.&.WHR_base"]<0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5   
      
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf = mean(bf))
  }
  if (model == "E3a_exfunct"){
    for (i in c(1:imp$m)){

      tmp_bf <- generalTestBF(formula = as.formula("exfunct ~ age_base + age_change + 
                                                              asinh_wml_change:sex +
                                                              asinh_wml_base + asinh_wml_change + 
                                                              sex + cesd + TIV + education+ subj"), 
                              data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                              multicore = T,
                              neverExclude = c("age_base", "^age_change$",
                                               "education", "cesd", "education",
                                               "TIV", "subj")
      )
      
      full=9 #model  is full model; 
      red=7 #model: wml_change:sex dropped
      
      bf_etmp=extractBF(tmp_bf,logbf = F)
      bf_extracted=data.frame(pred=c("asinh_wml_change:sex"), imp_n=rep(i,1))
      bf_extracted[,"bf"] <- c(bf_etmp[full,1] / bf_etmp[red,1])
      
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      chains <- posterior(tmp_bf, full, iterations = n_it, columnFilter="^id$")
      #We expect a positive effect (women coded as 0)...
      bf_extracted[,"bf_sf"] <- c(mean(chains[,"asinh_wml_change.&.sex"]>0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5   
      
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf = mean(bf))
  }
  if (model == "E3b_globalcog"){
    for (i in c(1:imp$m)){

      
      tmp_bf <- generalTestBF(formula = as.formula("globalcog ~ age_base + age_change + 
                                                              asinh_wml_change:sex +
                                                              asinh_wml_base + asinh_wml_change + 
                                                              sex + cesd + education + TIV + subj"), 
                              data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                              multicore = T,
                              neverExclude = c("age_base", "^age_change$",
                                               "education", "cesd", "education",
                                               "TIV", "subj")
      )
      
      full=9 #model  is full model; 
      red=7 #model: wml_change:sex dropped
      
      bf_etmp=extractBF(tmp_bf,logbf = F)
      bf_extracted=data.frame(pred=c("asinh_wml_change:sex"), imp_n=rep(i,1))
      bf_extracted[,"bf"] <- c(bf_etmp[full,1] / bf_etmp[red,1])
      
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      chains <- posterior(tmp_bf, full, iterations = n_it, columnFilter="^id$")
      #We expect a positive effect (women coded as 0)...
      bf_extracted[,"bf_sf"] <- c(mean(chains[,"asinh_wml_change.&.sex"]>0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5   
      
      
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
    }
    
    bf_df=as.data.frame(bf)
    bf_res <- bf_df %>% 
      group_by(pred) %>% 
      summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf = mean(bf))
  }
  return(list(est, bf_res))
}