run_LME<- function(imp, model){
    # Calculates model results for a MIDS (Multiply Imputed Data Set) from mice.
    # Model can take values "VRF", "exfunct" and "globalcog".
  
    #########################
    # Frequentist statistics
    #########################
    # Pool using the mice standard workflow to collect results from 5 imputed datasets
    if (model == "VRF"){
    res <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|sic)'),
                data=data, REML=F, na.action = na.omit)
    est=summary(mice::pool(res))}
    if (model == "exfunct"){
    res <- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|sic)'),
                                 data=data, REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    }
    else{res <- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|sic)'),
                     data=data, REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    }

    
    #########################
    # Bayesian statistics
    #########################
    #generalTestBF drops factors one by one and returns BF of the reduced model compared to null model
    #(Order of dropped coefficients: age_change:WHR_base, age_change:SBP_base, WHR_change, SBP_change)
    #Thus to obtain the likelihood of the full compared to null model, the respective Bayes Factors have
    #to be inverted. This is done in when calculating the mean BF.
    
    #Pool by calculating mean and standard deviation of BFs across 5 imputations
    comp_imp=mice::complete(imp, "long")  
    if (model == "CVR"){
      for (i in c(1:imp$m)){
        
        
        tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                 SBP_base + age_change:SBP_base + SBP_change +
                                                                 WHR_base + WHR_base:age_change + WHR_change +
                                                                 sex + BPmed + TIV + subj"), 
                                        data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                                        whichModels="top", multicore = T,
                                        neverExclude = c("age_base", "^age_change$", "^SBP_base$", 
                                                         "^WHR_base$", "sex", "BPmed", "TIV", "subj")
        )
        
        bf_extracted=extractBF(tmp_bf,logbf = F)
        bf_extracted$imp=i
        bf_extracted$pred=c("age_change:WHR_base", "age_change:SBP_base", "WHR_change", "SBP_change")
        
        # take one-sided nature of the alternative hypothesis into account 
        # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
        # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
        for (j in c(1:4)){#iterate over predictors
        chains <- posterior(tmp_bf, j, iterations = 1000)
        siding_factor <- mean(chains[,1]>0) 
        #We expect a positive effect of SBP/WHR on age change, 
        # as well as for SBP and WHR change. Because the Bayes factor from WhichModels="top"
        # is inverted, the siding factor also has to be calculated inversely
        bf_extracted[j, "one_sided_bf"] <- bf_extracted[j, "bf"] * 0.5/siding_factor
        }
        

        if (i==1){
          bf=bf_extracted}
        else{
          bf=rbind(bf,bf_extracted)
        }
        }
        
        bf_df=as.data.frame(bf)
        bf_res <- bf_df %>% 
          group_by(pred) %>% 
          mutate(mean_one_sided_bf = 1/mean(one_sided_bf), sd_bf=sd(one_sided_bf))%>%
          mutate(mean_bf = 1/mean(bf), sd_bf=sd(bf))
        }
    
    if (model == "exfunct"){
      for (i in c(1:imp$m)){
      tmp_bf <- generalTestBF(formula = as.formula("exfunct ~ age_base + age_change +
                                                              asinh_wml_base + asinh_wml_change +
                                                              sex + education + cesd + TIV + subj"), 
                                      data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  
                                      whichModels="top", multicore = T,
                                      neverExclude = c("age_base", "^age_change$",
                                                       "sex", "education", "cesd", 
                                                       "TIV", "subj")
      )

      bf_extracted=extractBF(tmp_bf,logbf = F)
      bf_extracted$imp=i
      bf_extracted$pred=c("asinh_wml_base", "asinh_wml_change")
      
      # take one-sided nature of the alternative hypothesis into account 
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      for (j in c(1:2)){#iterate over predictors
        chains <- posterior(tmp_bf, j, iterations = 1000)
        siding_factor <- mean(chains[,1]<0) 
        #We expect a negative effect of baseline and change in WML
        # Because the Bayes factor from WhichModels="top"
        # is inverted, the siding factor also has to be calculated inversely
        bf_extracted[j, "one_sided_bf"] <- bf_extracted[j, "bf"] * 0.5/siding_factor
      }
      
      if (i==1){
        bf=bf_extracted}
      else{
        bf=rbind(bf,bf_extracted)
      }
      }
      
      bf_df=as.data.frame(bf)
      bf_res <- bf_df %>% 
        group_by(pred) %>% 
        mutate(mean_one_sided_bf = 1/mean(one_sided_bf), sd_bf=sd(one_sided_bf))%>%
        mutate(mean_bf = 1/mean(bf), sd_bf=sd(bf))
      }
    
    else{
      for (i in c(1:imp$m)){
        tmp_bf <- generalTestBF(formula = as.formula("globalcog ~ age_base + age_change +
                                                                   asinh_wml_base + asinh_wml_change +
                                                                   sex + education + cesd + TIV + subj"), 
                                        data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  
                                        whichModels="top", multicore = T,
                                        neverExclude = c("age_base", "^age_change$",
                                                         "sex", "education", "cesd", 
                                                         "TIV", "subj")
        )
        
        bf_extracted=extractBF(tmp_bf,logbf = F)
        bf_extracted$imp=i
        bf_extracted$pred=c("asinh_wml_base", "asinh_wml_change")

        # take one-sided nature of the alternative hypothesis into account 
        # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
        # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
        for (j in c(1:2)){#iterate over predictors
          chains <- posterior(tmp_bf, j, iterations = 1000)
          siding_factor <- mean(chains[,1]<0) 
          #We expect a negative effect of baseline and change in WML
          # Because the Bayes factor from WhichModels="top"
          # is inverted, the siding factor also has to be calculated inversely
          bf_extracted[j, "one_sided_bf"] <- bf_extracted[j, "bf"] * 0.5/siding_factor
        }
        
        if (i==1){
          bf=bf_extracted}
        else{
          bf=rbind(bf,bf_extracted)
        }
      }
      
      bf_df=as.data.frame(bf)
      bf_res <- bf_df %>% 
        group_by(pred) %>% 
        mutate(mean_one_sided_bf = 1/mean(one_sided_bf), sd_bf=sd(one_sided_bf))%>%
        mutate(mean_bf = 1/mean(bf), sd_bf=sd(bf))
      }
    
    return(list(res, est, bf_res))
    }
  
  