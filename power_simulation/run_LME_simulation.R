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
    
    #generalTestBF drops factors one by one and returns BF of the reduced model compared to full model
    #(Order of dropped coefficients: age_change:WHR_base, age_change:SBP_base, WHR_change, SBP_change)
    #Thus to obtain the likelihood of the full model compared to the model without the specific term,
    #the respective Bayes Factors have to be inverted. This in done in evaluate_power.R
    tmp_bf <- generalTestBF(formula = as.formula(paste0(dv, "~ age_base + age_change +
                                                               SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change +
                                                               sex + icv + id")), 
                            data=data, whichRandom = "id",  whichModels="top", multicore = T,
                            neverExclude = c("age_base", "^age_change$", "^SBP_base$", "^WHR_base$", "sex", "icv", "id")
    )
    bf_extracted=extractBF(tmp_bf,logbf = F)
    for (j in c(1:4)){#iterate over predictors
      chains <- posterior(tmp_bf, j, iterations = 1000)
      siding_factor <- mean(chains[,1]>0) 
      #We expect a positive effect of interaction of baseline SBP and WHR with age change
      #and of change in SBP/WHR
      #Because the Bayes factor from WhichModels="top"
      #is inverted (i.e. quantifies evidence for null),
      # the siding factor also has to be calculated inversely (quantifying evidence
      # for the opposite direction)
      bf_extracted[j, "one_sided_bf"] <- bf_extracted[j, "bf"] * 0.5/siding_factor
    }

  }
  return(list(res, sres, coeffs, bf_extracted))
  }
  
  