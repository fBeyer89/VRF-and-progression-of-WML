#Calculate Power
library(dplyr)

results_dat=read.csv("/data/pt_life_whm/Results/VRF_cSVD/LME/simulations//n_400_600_800_100sim.csv")
n_sim=100
tmp_power<- results_dat %>%
   filter(value == "p_value")%>%
   group_by(nsample)%>%
   summarize(WHR_base = sum(`age_change.WHR_base` < 0.025),
             SBP_base = sum(`age_change.SBP_base` < 0.025),
             WHR_change = sum(`WHR_change` < 0.025),
             SBP_change = sum(`SBP_change` < 0.025))%>%
   mutate(power_WHRb=WHR_base/n_sim,
          power_SBPb=SBP_base/n_sim,
          power_WHRc=WHR_change/n_sim,
          power_SBPc=SBP_change/n_sim)

tmp_power#
write.csv(tmp_power,"/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/power_400_600_800_100sim.csv")

tmp_bf<- results_dat %>%
   filter(value == "BF")%>%
   group_by(nsample)%>%
   summarize(meanBF_WHRb=1/mean(`age_change.WHR_base`),
          meanBF_SBPb=1/mean(`age_change.SBP_base`),
          meanBF_WHRc=1/mean(`WHR_change`),
          meanBF_SBPc=1/mean(`SBP_change`))%>%

tmp_bf_one_sided<- results_dat %>%  
  filter(value == "one_sided_BF")%>%
  group_by(nsample)%>%
  summarize(meanBF_one_sided_WHRb=1/mean(`age_change.WHR_base`),
            meanBF_one_sided_SBPb=1/mean(`age_change.SBP_base`),
            meanBF_one_sided_WHRc=1/mean(`WHR_change`),
            meanBF_one_sided_SBPc=1/mean(`SBP_change`))
#
write.csv(tmp_bf,"/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/bf_400_600_800_100sim.csv")
write.csv(tmp_bf_one_sided,"/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/bf_one_sided_400_600_800_100sim.csv")
