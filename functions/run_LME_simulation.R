run_LME<- function(dv, data, model){
  
  # runs the model either with simple or extended set of covariates
  
  if(dv=='log_ll'|dv=='asinh_ll' && model=="covariates"){
    res <- lmerTest::lmer(formula = paste0(dv, '~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + diabetes + smoking + APOE + HT_medication + icv + (1|id)'),
                data=data, REML=F, na.action = na.omit)
    vifres <- as.data.frame(vif(res))
    

    bf <- extractBF((generalTestBF(formula = as.formula(paste0(dv, "~ age_base + age_change +
                                                               SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change +
                                                               sex + diabetes + smoking + APOE + HT_medication + icv + id")), 
                                   data=data, whichRandom = "id",  whichModels="top", multicore = T,
                                   neverExclude = c("age_base", "^age_change$", "^SBP_base$", "^WHR_base$", "sex", "icv", "id")
    )),
    logbf = F)

  }
  if(dv=='log_ll'|dv=='asinh_ll' && model=="simple"){
    res <- lmerTest::lmer(formula = paste0(dv, '~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + icv + (1|id)'),
                data=data, REML=F, na.action = na.omit)
    sres=summary(res)
    coeffs=sres$coefficients
    vifres <- as.data.frame(vif(res))
    
    #extractBF drops factors one by one and returns BF of the reduced model compared to null model
    #(Order of dropped coefficients: age_change:WHR_base, age_change:SBP_base, WHR_change, SBP_change)
    #Thus to obtain the likelihood of the full compared to null model, the respective Bayes Factors have
    #to be inverted. This in done in evaluate_power.R
    bf <- extractBF((generalTestBF(formula = as.formula(paste0(dv, "~ age_base + age_change +
                                                               SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change +
                                                               sex + icv + id")), 
                                   data=data, whichRandom = "id",  whichModels="top", multicore = T,
                                   neverExclude = c("age_base", "^age_change$", "^SBP_base$", "^WHR_base$", "sex", "icv", "id")
                     )),
                    logbf = F)

  }
  return(list(res, sres, coeffs, bf))
  }
  
  