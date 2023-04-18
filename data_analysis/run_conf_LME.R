run_conf_LME<- function(imp, model,n_it=10){
    # Calculates model results for a MIDS (Multiply Imputed Data Set) from mice.
    # Model can take values "M1_VRF", "M2_exfunct" and "M3_globalcog".
  
    #########################
    # Frequentist statistics
    #########################
    # Pool using the mice standard workflow to collect results from 5 imputed datasets
    if (model == "M1_VRF"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    }
    if (model == "M2_exfunct"){
    res <- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                                 REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    }
    if (model == "M3_globalcog"){
    res <- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                     REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    }

    
    #########################
    # Bayesian statistics
    #########################
    
    #Pool by calculating mean and standard deviation of BFs across all imputed datasets 
    comp_imp=mice::complete(imp, "long")  
    if (model == "M1_VRF"){
      for (i in c(1:imp$m)){
        #extract imputed datasets to run functions on individual datasets
        
        tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                 DBP_base + age_change:DBP_base + DBP_change +
                                                                 WHR_base + WHR_base:age_change + WHR_change +
                                                                 sex + BPmed + TIV + subj"), 
                                        data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                                        multicore = T, whichModels="top",
                                        neverExclude = c("age_base", "age_change^", "^DBP_base$", 
                                                         "^WHR_base$", "sex", "BPmed", "TIV", "subj")
        )
        #model 15 is full model; 
        #model 7: age_change:WHR_base dropped
        #model 11: age_change:DBP_base dropped, 
        #model 14: DBP change dropped
        #model 13: WHR change dropped
        
        bf_etmp=extractBF(tmp_bf,logbf = F)
        bf_extracted=data.frame(pred=c("age_change:WHR_base", "age_change:DBP_base",
                                       "WHR_change", "DBP_change"), imp_n=rep(i,4))
        bf_extracted[,"bf"] <- c(bf_etmp[15,1] / bf_etmp[7,1], bf_etmp[15,1] / bf_etmp[11,1],
                                 bf_etmp[15,1] / bf_etmp[13,1],bf_etmp[15,1] / bf_etmp[14,1])
        
        chains <- posterior(tmp_bf, 15, iterations = n_it, columnFilter="^subj$")
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
        
        
        bf_df=as.data.frame(bfall)
        bf_res <- bf_df %>% 
          group_by(pred) %>% 
          summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf))
        }
    
    if (model == "M2_exfunct"){
      for (i in c(1:imp$m)){
      tmp_bf <- generalTestBF(formula = as.formula("exfunct ~ age_base + age_change +
                                                              asinh_wml_base + asinh_wml_change +
                                                              sex + education + cesd + TIV + subj"), 
                                      data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  
                                      multicore = T,
                                      neverExclude = c("age_base", "^age_change$",
                                                       "sex", "education", "cesd", 
                                                       "TIV", "subj")
      )

      #model 3 is full model
      #model 1: asinh_wml_change  dropped
      #model 2: asinh_wml_base dropped
      
      bf_etmp=extractBF(tmp_bf,logbf = F)

      bf_extracted=data.frame(pred=c("asinh_wml_change", "asinh_wml_base"), imp_n=rep(i,2))
      bf_extracted[,"bf"] <- c(bf_etmp[3,1] / bf_etmp[1,1], bf_etmp[3,1] / bf_etmp[2,1])
      
      chains <- posterior(tmp_bf, 3, iterations = n_it, columnFilter="^id$")#The third model is the full model with all 10 terms.
      
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      #We expect a negative effect of baseline and change in WML
      bf_extracted[,"bf_sf"] <- c(mean(chains[,"asinh_wml_change"]<0), mean(chains[,"asinh_wml_base"]<0))
      
      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5        
      
      
      if (i==1){
        bfall=bf_extracted}
      else{
        bfall=rbind(bfall,bf_extracted)
      }}
    
    
      bf_df=as.data.frame(bfall)
      bf_res <- bf_df %>% 
        group_by(pred) %>% 
        summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf))
      }
    
    if(model == "M3_globalcog"){
      for (i in c(1:imp$m)){
        tmp_bf <- generalTestBF(formula = as.formula("globalcog ~ age_base + age_change +
                                                                   asinh_wml_base + asinh_wml_change +
                                                                   sex + education + cesd + TIV + subj"), 
                                        data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  
                                        multicore = T,
                                        neverExclude = c("age_base", "^age_change$",
                                                         "sex", "education", "cesd", 
                                                         "TIV", "subj")
        )
        
        #model 3 is full model
        #model 1: asinh_wml_change  dropped
        #model 2: asinh_wml_base dropped
        
        bf_etmp=extractBF(tmp_bf,logbf = F)
        
        bf_extracted=data.frame(pred=c("asinh_wml_change", "asinh_wml_base"), imp_n=rep(i,2))
        #calculate the bayes factor in favour of the full model (compared to model without term of interest)
        bf_extracted[,"bf"] <- c(bf_etmp[3,1] / bf_etmp[1,1], bf_etmp[3,1] / bf_etmp[2,1])
        
        chains <- posterior(tmp_bf, 3, iterations = n_it, columnFilter="^id$")#The third model is the full model with all 10 terms.
        
        # take one-sided nature of the alternative hypothesis into account 
        # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
        # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
        #We expect a negative effect of baseline and change in WML
        bf_extracted[,"bf_sf"] <- c(mean(chains[,"asinh_wml_change"]<0), mean(chains[,"asinh_wml_base"]<0))
        
        bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "bf_sf"]/0.5        
        
          
        if (i==1){
          bfall=bf_extracted}
        else{
          bfall=rbind(bfall,bf_extracted)
        }}


        bf_df=as.data.frame(bfall)
        bf_res <- bf_df %>% 
          group_by(pred) %>% 
          summarise(mean_one_sided_bf = mean(one_sided_bf), sd_bf=sd(one_sided_bf),mean_bf=mean(bf))
        }
    
    return(list(est, bf_res))
    }
  
  