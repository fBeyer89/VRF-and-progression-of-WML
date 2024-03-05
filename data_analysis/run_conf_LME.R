run_conf_LME<- function(imp, model,n_it=10, outdir="/data/pt_life_whm/Results/VRF_cSVD/LME/results/Conf/workspace_"){
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
    
    #calculate p-values based on Wald nested models comparisons 
    #for DBP_baseline interaction with time
    res_DBP_baseline <- with(imp, lmerTest::lmer(formula = 'asinh_wml ~ 
                                 age_base + age_change + 
                                 DBP_base + DBP_change +
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                 REML=F, na.action = na.omit)
    d1=D1(res,res_DBP_baseline)
    write.csv(d1$result,paste0(outdir,model, '_freq_imp_res_d1_DBP_baseline.csv'))
    save(res,file=paste0(outdir,model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    save(plot,file=paste0(outdir,model, '_freq_imp_plot', '.RData'))
    }
    if (model == "M2_exfunct"){
    print("running M2")
    res <- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                                 REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    res0<- with(imp, lmerTest::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + 
                                 sex + education + cesd + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    d1=D1(res,res0)
    write.csv(d1$result,paste0(outdir,model, '_freq_imp_res_d1_asinh_wml_change.csv'))
    save(res,file=paste0(outdir,model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'exfunct ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                 REML=F, na.action = na.omit)
    save(plot,file=paste0(outdir,model, '_freq_imp_plot', '.RData'))
    }
    if (model == "M3_globalcog"){
    res <- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                     REML=F, na.action = na.omit)
    est=summary(mice::pool(res))
    res0<- with(imp, lmerTest::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base + 
                                 sex + education + cesd + TIV + (1|subj)'),
                REML=F, na.action = na.omit)
    d1=D1(res,res0)
    write.csv(d1$result,paste0(outdir,model, '_freq_imp_res_d1_asinh_wml_change.csv'))
    
    save(res,file=paste0(outdir,model, '_freq_imp_res', '.RData'))
    #recalculate model with lme4 as to use effect for effect prediction
    plot <- with(imp, lme4::lmer(formula = 'globalcog ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                 sex + education + cesd + TIV + (1|subj)'),
                                  REML=F, na.action = na.omit,
                 REML=F, na.action = na.omit)
    save(plot,file=paste0(outdir,model, '_freq_imp_plot', '.RData'))
    }

    
    #########################
    # Bayesian statistics
    #########################
    
    #Pool by calculating mean and standard deviation of BFs across all imputed datasets
    #Therefore, we extract five datasets in long format:
    comp_imp=mice::complete(imp, "long")  
    
    if (model == "M1_VRF"){
      for (i in c(1:imp$m)){
        print("running model M1")
        if(!file.exists(paste0(outdir,model, '_imp_',i, '.RData'))){
        #extract imputed datasets to run functions on individual datasets
        
        #generalTestBF drops factors one by one and returns BF of the reduced model compared to full model
        #(Order of dropped coefficients: age_change:WHR_base, age_change:DBP_base, WHR_change, DBP_change)
        #Thus to obtain the likelihood of the full model compared to the model without the specific term,
        #the respective Bayes Factors have to be inverted. 
        tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                 DBP_base + age_change:DBP_base + DBP_change +
                                                                 WHR_base + WHR_base:age_change + WHR_change +
                                                                 sex + BPmed + TIV + subj"), 
                                data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj", 
                                multicore = T, whichModels="top",
                                neverExclude = c("age_base", "age_change$", "^DBP_base$", 
                                                 "^WHR_base$", "sex", "BPmed", "TIV", "subj")
        )
        bf_extracted=extractBF(tmp_bf,logbf = F)
        bf_extracted$pred=rownames(bf_extracted)
        
        
        #calculate only the full model for deriving siding factors
        tmp_bf_chains <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change +
                                                                 DBP_base + age_change:DBP_base + DBP_change +
                                                                 WHR_base + WHR_base:age_change + WHR_change +
                                                                 sex + BPmed + TIV + subj"),
                                       data=comp_imp[comp_imp$.imp==i,], whichRandom = "subj",  multicore = T,
                                       neverExclude = c("age_base", "age_change", "DBP_base", "DBP_change",
                                                        "WHR_base", "WHR_change", "sex", "BPmed", "TIV", "subj")
        )
        
        siding_factor=vector()
        chains <- posterior(tmp_bf_chains,  iterations = n_it, columnFilter="^subj$")
        bf_extracted[1,"siding"]=mean(chains[,"age_change.&.WHR_base"]>0)
        bf_extracted[2,"siding"]=mean(chains[,"age_change.&.DBP_base"]>0)
        bf_extracted[3,"siding"]=mean(chains[,"WHR_change"]>0)
        bf_extracted[4,"siding"]=mean(chains[,"DBP_change"]>0)
        
        #We expect a positive effect of interaction of baseline DBP and WHR with age change
        #and of change in DBP/WHR on WML volume change
        #Because the Bayes factor from WhichModels="top"
        #is inverted (i.e. quantifies evidence for reduced model compared to full model),
        #we have to invert it back and multiply by the siding factor to quantifying evidence
        #for the alternative hypotheses/full model
        
        bf_extracted[, "bf"] <- (1/bf_extracted[, "bf"]) 
        bf_extracted[, "one_sided_bf"] <- (bf_extracted[, "bf"]) * bf_extracted[, "siding"]/0.5
        
        
        list2save=list(tmp_bf,chains, bf_extracted)
        save(list2save, 
             file=paste0(outdir,model, '_imp_',i, '.RData'))
        }
        else {
          load(paste0(outdir,model, '_imp_',i, '.RData'))
          bf_extracted=list2save[[3]]
        }
        
        if (i==1){
          bfall=bf_extracted}
        else{
          bfall=rbind(bfall,bf_extracted)}
      }
        
        bf_df=as.data.frame(bfall)
        bf_res <- bf_df %>% 
          group_by(pred) %>% 
          summarise(mean_one_sided_bf = mean(one_sided_bf), sd_osbf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
        
        
        }
    
    if (model == "M2_exfunct"){
      print("running model M2")
      for (i in c(1:imp$m)){
        if(!file.exists(paste0(outdir,model, '_imp_',i, '.RData'))){
      
      tmp=comp_imp[comp_imp$.imp==i,]
      modeldat=tmp[!is.na(tmp$exfunct),] #remove timepoints with NA in dependent variable
      
      tmp_bf=generalTestBF(formula = as.formula("exfunct ~ age_base + age_change +
                                                               asinh_wml_base + asinh_wml_change +
                                                               sex + education + cesd + TIV + subj"),
                                       data=modeldat, whichRandom = "subj",
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

      chains <- posterior(tmp_bf, 3, iterations = n_it, columnFilter="^subj$")#The third model is the full model with all 10 terms.


      # take one-sided nature of the alternative hypothesis into account
      # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
      # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
      #We expect a negative effect of baseline and change in WML
      bf_extracted[,"siding"] <- c(mean(chains[,"asinh_wml_change"]<0), mean(chains[,"asinh_wml_base"]<0))

      bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5
      
      list2save=list(tmp_bf,chains, bf_extracted)
      save(list2save, 
           file=paste0(outdir,model, '_imp_',i, '.RData'))
      }
      else {
        load(paste0(outdir,model, '_imp_',i, '.RData'))
        bf_extracted=list2save[[3]]
        }
      if (i==1){
        bfall=bf_extracted}
      else{
        bfall=rbind(bfall,bf_extracted)
      }  
        }
    
    
      bf_df=as.data.frame(bfall)
      bf_res <- bf_df %>% 
        group_by(pred) %>% 
        summarise(mean_one_sided_bf = mean(one_sided_bf), sd_osbf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
      
      }
    
    if(model == "M3_globalcog"){

      
      for (i in c(1:imp$m)){
        tmp=comp_imp[comp_imp$.imp==i,]
        modeldat=tmp[!is.na(tmp$globalcog),] #remove timepoints with NA in dependent variable
        tmp_bf <- generalTestBF(formula = as.formula("globalcog ~ age_base + age_change +
                                                                   asinh_wml_base + asinh_wml_change +
                                                                   sex + education + cesd + TIV + subj"), 
                                        data=modeldat, whichRandom = "subj",  
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
        
        chains <- posterior(tmp_bf, 3, iterations = n_it, columnFilter="^subj$")#The third model is the full model with all 10 terms.
        
        # take one-sided nature of the alternative hypothesis into account 
        # adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
        # see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
        #We expect a negative effect of baseline and change in WML
        bf_extracted[,"siding"] <- c(mean(chains[,"asinh_wml_change"]<0), mean(chains[,"asinh_wml_base"]<0))
        
        bf_extracted[, "one_sided_bf"] <- bf_extracted[, "bf"] * bf_extracted[, "siding"]/0.5        
        
          
        if (i==1){
          bfall=bf_extracted}
        else{
          bfall=rbind(bfall,bf_extracted)
        }
        list2save=list(tmp_bf,chains, bf_extracted)
        save(list2save, 
             file=paste0(outdir,model, '_imp_',i, '.RData'))
             bf_extracted=list2save[[3]]
        }


        bf_df=as.data.frame(bfall)
        bf_res <- bf_df %>% 
          group_by(pred) %>% 
          summarise(mean_one_sided_bf = mean(one_sided_bf), sd_osbf=sd(one_sided_bf),mean_bf=mean(bf), mean_siding=mean(siding))
        
        }
    
    return(list(est, bf_res))
    }
  
  