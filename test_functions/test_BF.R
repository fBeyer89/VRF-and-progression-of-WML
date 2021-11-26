
comp_imp=mice::complete(imp, "long") 
comp_imp$subj=as.factor(comp_imp$subj)
tmp_bf <- generalTestBF(formula = as.formula("asinh_wml ~ age_base + age_change"), 
                        data=comp_imp[comp_imp$.imp==1,], whichRandom = "subj", 
                        whichModels="top", multicore = T,
                        neverExclude = c("subj"))

bf_extracted=extractBF(tmp_bf,logbf = F)
bf_extracted$imp=i
bf_extracted$pred=c("age_change:WHR_base", "age_change:SBP_base", "WHR_change", "SBP_change")

# take one-sided nature of the alternative hypothesis into account 
# adjust the BF in favour of the alternative by p(effect size>0)/0.5 by sampling from the posterior
# see https://gist.github.com/richarddmorey/7c1bd06a14384412f2145daee315c036 for more detail
for (j in c(1:2)){#iterate over predictors
  chains <- posterior(tmp_bf, j, iterations = 10)
  siding_factor <- mean(chains[,1]>0) #We expect a positive effect of SBP/WHR on age change, 
  # as well as for SBP and WHR change
  bf_extracted[j, "one_sided_bf"] <- bf_extracted[j, "bf"] * 0.5/siding_factor} 