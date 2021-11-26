library("lme4")
library("lmerTest") 
library(doMC)
library(car)
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer

setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/power_simulation/") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('run_LME_simulation.R')
source('simulate_data.R')


### Simulate power curves for SBP
# Simulation parameters
n_sim = 100
n_sample=c(400,600,800)
results_vec=vector()

##### Baseline effect sizes
bl=read.csv("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/PV168_A1_Pilot_subject_list_inclusion_exclusion29.1.19.csv")
usable=subset(bl,(bl$Age_all>50&(!is.na(bl$ADULT_BP_SBP.1))&(!is.na(bl$waist2hip))&(!is.na(bl$icv))))

#Cross section model to derive cross coefficients (#use log-linked GLM from Gamma family )
model_bl=glm(lesionload/1000 ~ scale(Age_all, scale=F) + scale(ADULT_BP_SBP.1, scale=F) + 
               scale(waist2hip, scale=F) + sex + scale(icv, scale=F), 
             data=usable, family = Gamma(link = "log"))
coeffs_bl=coefficients(model_bl)

#Calculate covariance matrix and sex distribution
covmatrix <- cov(usable[,c("Age_all", "ADULT_BP_SBP.1", "waist2hip", "icv")], use = "complete.obs") # determine (co)variance of relevant variables

#Calculate elapsed time between baseline and followup   
tba=read.csv("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/radiological_assessment/radio_assessment_wide.csv")
tba$difftime=difftime(as.Date(tba$Messung.Datum_TI.fu,format="%d.%m.%Y"),as.Date(tba$Messung.Datum_TI.bl,format="%d.%m.%Y"))
mean_diff_y=mean(as.numeric(tba[!is.na(tba$difftime),"difftime"]))/365
sd_diff_y=sd(as.numeric(tba[!is.na(tba$difftime),"difftime"]))/365


###### Effect sizes from literature #####
# Mean rate of change in WMH (cm3/y)
asps=1.38/6
scharf1=0.54
scharf2=1.04
scharf3=1.6
lothanian=4/3
godin=1.07/4
maillard=0.25
peng=(17.82-13.78)/4
nasrallah=1.53/3.98 #standard treatment group
havenon=0.93/(40/12) #glycemic intervention group

#number of participants
asps_n=243
scharf_n1=247
scharf_n2=186
scharf_n3=121
lothanian_n=439
godin_n=1319
nasrallah_n=200
havenon_n=501
maillard_n=1118
peng_n=294

# from publications above: weighted averaged WML change/longitudinal year: 0.64 cm³
# but often, this also contains effects of modifiable risk factors
weighted_av_wmlc=(asps*asps_n+
                    scharf1*scharf_n1+scharf2*scharf_n2+scharf3*scharf_n3+
                    lothanian*lothanian_n+
                    godin*godin_n+
                    nasrallah*nasrallah_n+
                    havenon*havenon_n+
                    maillard_n*maillard+
                    peng_n*peng)/
  (asps_n+scharf_n1+scharf_n2+scharf_n3+lothanian_n+godin_n+nasrallah_n+havenon_n+maillard_n+peng_n)

#standard deviation of rate of change in WMH (cm3/y)
#-> pooled standard deviation
asps_sd=3.76/6
scharf_1=1.27
scharf_2=1.93
scharf_3=2.4
lothanian=4.3/3
godin=2.76/4
nasrallah=0.71
havenon_sd=1.20/(40/12)#glycemic intervention group
maillard_n=1118

pooled_sd_wmlc=((asps_n-1)*asps_sd**2+(scharf_n1-1)*scharf_1**2+(scharf_n2-1)*scharf_2**2+(scharf_n3-1)*scharf_3**2+
                  (lothanian_n-1)*lothanian**2+(godin_n-1)*godin**2+(nasrallah_n-1)*nasrallah**2+(havenon_n-1)*havenon_sd**2)/
  (asps_n+scharf_n1+scharf_n2+scharf_n3+lothanian_n+godin_n+nasrallah_n+havenon_n-8)

#### "SBP
#Effect of baseline SBP on WML progression (in ml per year per mmHg) (interaction with age_change)
#all publications did not consider log-transformed WML data.
baseline_SBP_Godin=0.04/(4*5) #0.04 (0.02) cm³ per 5mmHG SBP difference at baseline (N=1319) over 4 years -Godin
baseline_SBP_dickie=0.0271/3 #associations between VRF at 73 years and change in WMH volume from 73 years to 76 years in cm³
baseline_SBP_gottesmann=1.1/(20*10) #effect of SBP at baseline MRI on WML progression in 10 years in cm³ /20 mmHG - Gottesmann
baseline_SBP_verhaaren=0.08/18 # systolic blood pressure and WML progression:: 0.05 (95% CI, 0.00; 0.09) mL/y per SD increase in systolic blood pressure (~18)
SBP_baseline=(baseline_SBP_Godin+baseline_SBP_dickie+baseline_SBP_gottesmann+baseline_SBP_verhaaren)/4 #calculated per mmHg baseline difference/year

#Effect of change in SBP on WML progression  (in ml per year per mmHg) 
increase_SBP_Godin=0.05/(4*5) #0.05 (0.02) cm³ per 5mmHG SBP increase between timepoints (N=1319) over 4 years -Godin
SBP_change=increase_SBP_Godin

###WHR
#Work around for baseline on progression: use baseline_SBP effect scaled by WHR/mmHg conversion

#Effect of baseline WHR on WML progression in LIFE-Adult usable sample:
#One WHR unit has cross-sectional effect of by exp(2.15) (compared to exp(0.01) for 1 mmHg SBP)
#Thus the effect of 1 mmHg is approximately the effect 1/8.5 WHR unit
#Or 1 WHR unit has the effect of 8.5 mmHg.

#We half the effect size for a conservative estimate.
WHR_baseline=SBP_baseline*8.5/(2)

#Similar for WHR change on progression
WHR_change=SBP_change*8.5/(2)


#Parameter constraints (random intercept and overall model error)
sd_reff=0.5 #for normal distribution with mean 0, SD = 0.5 cm³
sd_err=1 #for normal distribution with mean 0, SD = 1 cm³  


for (n in n_sample){#,800 #number of subjects
    
    for (effect in c(1)){#0.75,1,1.25

      for (i in c(1:n_sim)){#number of simulations per condition
        #set.seed(i)
        #Draw strength of age effect from normal distribution, with mean at half the literature value
        #don't want negative values of age change because regression of WML is rare
        #pooled_sd_wmlc is very large
        effect_age_change=rnorm(n, mean=0.5*weighted_av_wmlc, sd=0.1) # 
        effect_age_change[effect_age_change<0]=0.01
        
        #Draw effects for SBP baseline, WHR baseline and SBP change from normal distributions,
        #considering effect strengths from literature
        effect_SBP_baseline=rnorm(n, mean=SBP_baseline*effect, sd=0.001)
        #Draw effect of SBP change from normal distribution with literature value
        effect_SBP_change=rnorm(n, mean=SBP_change, sd=0.001)
        #Draw effect of WHR baseline from normal distribution with literature value
        effect_WHR_baseline=rnorm(n, mean=WHR_baseline, sd=0.001)
        #Draw effect of WHR change from normal distribution, estimated based on 
        effect_WHR_change=rnorm(n, mean=WHR_change, sd=0.001)

          
        #Simulate
        dat=simncalc(n, 
                     coeffs_bl, covmatrix, usable, mean_diff_y, sd_diff_y,
                     effect_age_change, sd_age_change, effect_SBP_baseline, WHR_baseline, 
                     effect_SBP_change, WHR_change, sd_reff, sd_err)
        dat=as.data.frame(dat)
        scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm))
        dat <- dat %>% mutate_at(c("age_base", "SBP_base", "WHR_base", "icv"), scale2)
        dat$id=as.factor(dat$id)
        dat[dat$dv<=0,"dv"]=0.001 #do not allow negative WML load
        dat$log_ll=log(dat$dv)
        dat$asinh_ll=asinh(dat$dv)
        
        #Calculate model
        res_list=run_LME('asinh_ll', dat, 'simple')
        
        #Extract p-values, effect estimates, two-sided and one-sided BF for
        #predictors age_change:WHR_base, age_change:SBP_base, WHR_change, SBP_change
        results_vec=append(results_vec,res_list[3][[1]][c(11,10,7,5),5])
        results_vec=append(results_vec,res_list[3][[1]][c(11,10,7,5),1])
        results_vec=append(results_vec,res_list[[4]][c(1:4),c(1,5)])
        
        #Plot the simulated data as visual control
        #dat$tertiles=as.factor(ntile(dat$SBP_base,3))
        #ggplot(dat, aes(x = age_base+age_change, y = asinh_ll, group=tertiles, color = tertiles)) +
        #  geom_point(aes(group=id))+
        #  geom_smooth(method = "lm")
        #ggsave(paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/SBP_baseline_sim_effectoverage_sim_",i,"_n_",n,".png"), width = 8, height=5)
        
        #ggplot(dat, aes(x = tp, y = asinh_ll, group=tertiles, color = tertiles)) +
        #  geom_jitter(width = 0.15, alpha=0.5)+
        #  stat_summary(fun.y = mean, geom = "point", size=4)+
        # stat_summary(fun.y = mean, geom = "line") 
        #ggsave(paste0("/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/SBP_baseline_sim_effectperTP_sim_",i,"_n_",n,".png"), width = 8, height=5)
      }
      
    }
  }
      


results_mat <- matrix(results_vec, byrow=T, nrow=n_sim*length(n_sample)*3)
results_dat=as.data.frame(results_mat)
colnames(results_dat)=c("age_change:WHR_base", "age_change:SBP_base", "WHR_change", "SBP_change")
results_dat$value=rep(c("p_value", "effect_size", "BF", "one_sided_BF"),n_sim)
results_dat$nsim=rep(c(1:n_sim),each=3)
results_dat$nsample=rep(n_sample, each=3*n_sim)
write.csv(results_dat,"/data/pt_life_whm/Results/VRF_cSVD/LME/simulations/n_400_600_800_100sim_onesidedBF.csv")


