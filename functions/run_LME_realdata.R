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
                data=data, REML=F, na.action = na.omit))
    est=summary(mice::pool(res))}
    if (model == "exfunct"){
    res <- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|sic)'),
                                 data=data, REML=F, na.action = na.omit))
    est=summary(mice::pool(res))
    }
    else{res <- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|sic)'),
                     data=data, REML=F, na.action = na.omit))
    est=summary(mice::pool(res))
    }

    
    #########################
    # Bayesian statistics
    #########################
    #extractBF drops factors one by one and returns BF of the reduced model compared to null model
    #(Order of dropped coefficients: age_change:WHR_base, age_change:SBP_base, WHR_change, SBP_change)
    #Thus to obtain the likelihood of the full compared to null model, the respective Bayes Factors have
    #to be inverted. This is done in when calculating the mean BF.
    
    #Pool by calculating mean and standard deviation of BFs across 5 imputations
    comp_imp=mice::complete(imp, "long")  
    if (model == "CVR"){
      for (i in c(1:imp$m)){
        tmp <- extractBF((generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                 SBP_base + age_change:SBP_base + SBP_change +
                                                                 WHR_base + WHR_base:age_change + WHR_change +
                                                                 sex + BPmed + icv + sic"), 
                                        data=comp_imp[comp_imp$.imp==i,], whichRandom = "sic",  whichModels="top", multicore = T,
                                        neverExclude = c("age_base", "^age_change$", "^SBP_base$", 
                                                         "^WHR_base$", "sex", "BPmed", "icv", "sic")
        )),logbf = F)
        tmp$imp=i
        tmp$pred=c("age_change:WHR_base", "age_change:SBP_base", "WHR_change", "SBP_change")
        if (i==1){
          bf=tmp}
        else{
          bf=rbind(bf,tmp)
        }
        }
        
        bf_df=as.data.frame(bf)
        bf_res <- bf_df %>% 
          group_by(pred) %>% 
          mutate(mean_bf = 1/mean(bf), sd_bf=sd(bf))}
    if (model == "exfunct"){
      for (i in c(1:imp$m)){
      tmp <- extractBF((generalTestBF(formula = as.formula("exfunct ~ age_base + age_change +
                                                                 asinh_wml_base + asinh_wml_change +
                                                                 sex + education + cesd + icv + sic"), 
                                      data=comp_imp[comp_imp$.imp==i,], whichRandom = "sic",  whichModels="top", multicore = T,
                                      neverExclude = c("age_base", "^age_change$",
                                                       "sex", "education", "cesd", 
                                                       "icv", "sic")
      )),logbf = F)
      tmp$imp=i
      tmp$pred=c("asinh_wml_base", "asinh_wml_change")
      if (i==1){
        bf=tmp}
      else{
        bf=rbind(bf,tmp)
      }
      }
      
      bf_df=as.data.frame(bf)
      bf_res <- bf_df %>% 
        group_by(pred) %>% 
        mutate(mean_bf = 1/mean(bf), sd_bf=sd(bf))}
    else{
      for (i in c(1:imp$m)){
        tmp <- extractBF((generalTestBF(formula = as.formula("globalcog ~ age_base + age_change +
                                                                   asinh_wml_base + asinh_wml_change +
                                                                   sex + education + cesd + icv + sic"), 
                                        data=comp_imp[comp_imp$.imp==i,], whichRandom = "sic",  whichModels="top", multicore = T,
                                        neverExclude = c("age_base", "^age_change$",
                                                         "sex", "education", "cesd", 
                                                         "icv", "sic")
        )),logbf = F)
        tmp$imp=i
        tmp$pred=c("asinh_wml_base", "asinh_wml_change")
        if (i==1){
          bf=tmp}
        else{
          bf=rbind(bf,tmp)
        }
      }
      
      bf_df=as.data.frame(bf)
      bf_res <- bf_df %>% 
        group_by(pred) %>% 
        mutate(mean_bf = 1/mean(bf), sd_bf=sd(bf))
      }
    
    return(list(res, est, bf_res))
    }
  
  