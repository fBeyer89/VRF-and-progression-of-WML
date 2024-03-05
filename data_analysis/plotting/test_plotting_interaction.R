
miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/CVR_WMH_progression/imputed_data_6.12.23/imputed_data_6.12.23.Rdata")
i=1
comp_imp=mice::complete(imp, "long")
tmp=comp_imp[comp_imp$.imp==i,]

#generate quantiles
DBP_cutpoints = quantile(imp$data$DBP_base, probs = c(0,.33,.66,1),  na.rm = TRUE)
WHR_cutpoints = quantile(imp$data$WHR_base, probs = c(0,.33,.66,1),  na.rm = TRUE)

imp$data= imp$data %>%    
  mutate(cut_DBP = cut(DBP_base, breaks = DBP_cutpoints))%>%    
  mutate(cut_WHR = cut(DBP_base, breaks = WHR_cutpoints))


fit <- lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change + 
                                 sex + BPmed + TIV + (1|subj)',REML=F, na.action = na.omit, data=tmp)

# select only levels 30, 50 and 70 from continuous variable Barthel-Index
mydf <- ggpredict(fit, terms = c("age_change",paste0("DBP_base [",((DBP_cutpoints[1]+DBP_cutpoints[2])/2),",",
                                                     (DBP_cutpoints[2]+DBP_cutpoints[3])/2,",", (DBP_cutpoints[3]+DBP_cutpoints[4])/2,
                                                     "]")))
ggplot(mydf, aes(x, predicted, ymin = conf.low, ymax = conf.high, fill=group, colour = group)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5)


plot(mydf)

################
#Prepare imputed data for plotting effects
#Effects to be plotted from model M1   
fitted_lines <-
  tibble(.imp = 1:imp$m) %>% 
  mutate(p = map(.imp, ~ ggpredict(plot$analyses[[.]], terms = "DBP_base")) %>% 
                   data.frame())
  


fitted_lines <-
  tibble(.imp = 1:imp$m) %>% 
  mutate(p = map(.imp, ~ ggpredict(plot$analyses[[.]], terms = c("age_change",paste0("DBP_base [",((DBP_cutpoints[1]+DBP_cutpoints[2])/2),",",
                                                                                     (DBP_cutpoints[2]+DBP_cutpoints[3])/2,",", (DBP_cutpoints[3]+DBP_cutpoints[4])/2,
                                                                                     "]"))) %>% 
                   data.frame())
  )

fitted_lines <- fitted_lines %>% 
  unnest(p) 

fitted_lines=fitted_lines %>% 
  group_by(x,group) %>% #group by predictor values
  summarise(fit_bar = mean(predicted),
            v_w     = mean(std.error^2),
            v_b     = sum((predicted - fit_bar)^2) / (imp$m - 1),
            v_p     = v_w + v_b * (1 + (1 / imp$m)),
            se_p    = sqrt(v_p))%>% 
  # use the _p suffix to indicate these are pooled
  mutate(lwr_p = fit_bar - se_p * 1.96,
         upr_p = fit_bar + se_p * 1.96) 


outcome="asinh_wml"
xmin=min(fitted_lines$x)
xmax=max(fitted_lines$x)
ymax=max(max(imp$data[,paste0(outcome)],na.rm=T))
col2pred="#fdae61"
x_axis_label="Baseline DBP (mmHg)"
y_axis_label="WMH volume"
asterisiks="*"
scatter = 
  fitted_lines %>% 
  ggplot(aes(x = x)) +
  geom_ribbon(aes(ymin = lwr_p, ymax = upr_p, color=group, group=group),
              alpha = 1/2) +
  geom_line(aes(y = fit_bar, color=group, group=group), 
            size = 1/2) +
  # add the observed data for good measure
  geom_point(data = imp$data,
             aes_string(x=predictor, y = outcome, color="cut_DBP"),  alpha=0.5)+
  #add astericks depending on significance
  geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), label = asterisiks, cex = 4, hjust = 0)  +
  theme_bw() +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(x=x_axis_label, y=y_axis_label)+
  xlim(xmin, xmax) 

ggdraw(scatter) #p1 <- 
