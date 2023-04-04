#test random slopes

dv="asinh_ll"
res <- lmerTest::lmer(formula = paste0(dv, '~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + (1+age_change|id)'),
                      data=dat, na.action = na.omit)


res <- lmerTest::lmer(formula = paste0(dv, '~ age_base + age_change + 
                                 SBP_base + age_change:SBP_base + SBP_change  + WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + (1|id)'),
                      data=dat, REML=F, na.action = na.omit)
