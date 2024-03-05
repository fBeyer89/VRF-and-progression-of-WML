library("lme4")
library("lmerTest") 
library("BayesFactor")
library(DiagrammeR)
# ensure that lmerTest doesn't mask lmer, which would cause us multiple problems
lmer <- lme4::lmer
library(broom.mixed)
library(kableExtra)
library(tidyverse)
library(modelsummary)
library(robustlmm)
library(performance)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)
library(mice)
library(ggmice)
library(flextable)
library(ragg)
library(car)
library(readxl)


setwd("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/") #set path to github repository here
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/diagnostic_fcns.r")
source("/data/gh_gr_agingandobesity_share/literature/methods/statistics/linear_models_course_rogermundry_2018/functions/glmm_stability.r")
source('run_conf_LME.R')
source('run_exp_LME.R')
source('test_LME_assumptions.R')
source("/data/pt_life_whm/Analysis/multivariate_risk_svd/LIFE/R/bullseye/create_bullseye_plot.R")
source("/data/pt_life_whm/Analysis/multivariate_risk_svd/LIFE/R/bullseye/create_PCA.R")
library(influence.ME)

miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/CVR_WMH_progression/imputed_data_6.12.23/imputed_data_6.12.23.Rdata")



data_wm=read.csv("/data/pt_life_whm/Results/Tables/CVR_WMH_progression/imputed_data_6.12.23/full_dataset_long_excl_a.csv")
data_wm$MR_y_n_fu=1
data_wm= data_wm %>% 
  mutate(subj=as.integer(as.factor(mrt_pseudonym)))
data_wm$order=rownames(data_wm)

aseg <- read.table("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/FreeSurfer/CortParc_SegResults_baseline/FS_results_subcor_LIFE.txt", header = T, sep = "\t")
#aseg$fu <- ifelse(grepl("_fu", aseg$Measure.volume), 1, 0)
aseg$mrt_pseudonym <- substr(aseg$Measure.volume, 1, 10)
aseg_wo_duplicates <-
  aseg %>%
  filter(!duplicated(.[["mrt_pseudonym"]]))


#LI03082832/9675AE8517
aseg_wo_duplicates[nrow(aseg_wo_duplicates) + 1,"mrt_pseudonym"]="LI03082832"
aseg_wo_duplicates[aseg_wo_duplicates$mrt_pseudonym=="LI03082832", "EstimatedTotalIntraCranialVol"]=1350824.471522

#9B19DF4822 LI05095916
aseg_wo_duplicates[nrow(aseg_wo_duplicates) + 1,"mrt_pseudonym"]="LI05095916"
aseg_wo_duplicates[aseg_wo_duplicates$mrt_pseudonym=="LI05095916", "EstimatedTotalIntraCranialVol"]=1491826.541026

#B283FBCEC6 LI01651812
aseg_wo_duplicates[nrow(aseg_wo_duplicates) + 1,"mrt_pseudonym"]="LI01651812"
aseg_wo_duplicates[aseg_wo_duplicates$mrt_pseudonym=="LI01651812", "EstimatedTotalIntraCranialVol"]=1461721.788508

#6A1A2E4DC7 LI01927058
aseg_wo_duplicates[nrow(aseg_wo_duplicates) + 1,"mrt_pseudonym"]="LI01927058"
aseg_wo_duplicates[aseg_wo_duplicates$mrt_pseudonym=="LI01927058", "EstimatedTotalIntraCranialVol"]=1563415.052604

data_wm=merge(data_wm, aseg_wo_duplicates[,c("mrt_pseudonym", "EstimatedTotalIntraCranialVol")], all.x=T, by.x="sic", by.y="mrt_pseudonym")
colnames(data_wm)[3]="pseudonym" #MRT Pseudonym

#This exclusion was used for Bullseye extraction
fs=read.csv("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/FreeSurfer/QA_followup/freesurfer_qa_baseline_and_followup.csv")
data_wm=merge(data_wm, fs, all.x=T)

##Merge with current QA file
qa=read.csv("/data/pt_life_whm/Results/QA/qa_info_all_w_cohorts.csv")
data_wm_f_ex=merge(data_wm, qa, all.x=T, by="pseudonym")
data_wm_f_ex=data_wm_f_ex[data_wm_f_ex$BL_usable!=0,]




res=create_bullseye(sample=data_wm_f_ex,
                    mode="long",
                    location="/data/pt_life_whm/Data/WMparcellations_indiv/",
                    plot=T,
                    filename="/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures/lifec_old_bullseye.jpeg")

# create change plot
res_LIFE=res[[1]]
res_LIFE = res_LIFE %>% mutate(asinhWMH=asinh(WMH_vol_adj))
res_LIFE_bl=res_LIFE[res_LIFE$tp=="bl", ]
res_LIFE_fu=res_LIFE[res_LIFE$tp=="fu", ]

res_LIFE_bl$undefined_diff=res_LIFE_fu$undefined-res_LIFE_bl$undefined
res_LIFE_bl$valdiff=res_LIFE_fu$val-res_LIFE_bl$val
res_LIFE_bl$asinhdiff=res_LIFE_fu$asinhWMH-res_LIFE_bl$asinhWMH

res_sum = res_LIFE_bl %>% group_by(region,depth) %>%
  summarise(mean_val_bl=mean(val),
            med_valdiff=median(valdiff, na.rm=T),
            mean_valdiff=mean(valdiff, na.rm=T),
            med_asinhdiff=median(asinhdiff, na.rm=T))

levels(res_sum$region)=str_replace(levels(res_sum$region), "_", " ") 

p1=ggplot(res_sum, aes(x = region, y =depth)) +
  geom_tile(aes(fill=mean_valdiff),
            color = "gray50",
            lwd = 0.3,
            linetype = 1)+
  coord_curvedpolar(theta= "x")+
  scale_fill_gradient(low = "white", high = "red") +
  labs(fill="mean WMH volume \nchange in mm³")+
  theme_bw() +
  xlab("frontal")+
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank())
ggsave(p1,filename = '/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures/bullseye_WMH_progression.png', width = 20, height = 15, unit="cm", dpi=300)

p2=ggplot(res_sum, aes(x = region, y =depth)) +
  geom_tile(aes(fill=mean_val_bl),
            color = "gray50",
            lwd = 0.3,
            linetype = 1)+
  coord_curvedpolar(theta= "x")+
  scale_fill_gradient(low = "white", high = "red") +
  labs(fill="mean WMH volume \nin mm³")+
  theme_bw() +
  xlab("frontal")+
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank())
ggsave(p2,filename = '/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures/bullseye_WMH_baseline.png', width = 20, height = 15, unit="cm", dpi=300)

#### LOAD DATA ####
res_LIFE=res[[1]]
res_LIFE = res_LIFE %>% mutate(asinhWMH=asinh(WMH_vol_adj))
res_LIFE_bl=res_LIFE[res_LIFE$tp=="bl", c("pseudonym","asinhWMH", "bullseye_c")]%>%
  pivot_wider(id_cols=pseudonym, names_from=bullseye_c, values_from=asinhWMH)
#%>%mutate_at(.vars = c(2:37), ~(scale(.) %>% as.vector))
res_LIFE_bl$cohort="LIFE"

res_LIFE_fu=res_LIFE[res_LIFE$tp=="fu", c("pseudonym","asinhWMH", "bullseye_c")]%>%
  pivot_wider(id_cols=pseudonym, names_from=bullseye_c, values_from=asinhWMH)
#%>%mutate_at(.vars = c(2:37), ~(scale(.) %>% as.vector))
res_LIFE_fu$cohort="LIFE"

pca_LIFE=create_PCA(res_LIFE[res_LIFE$tp=="bl", c("pseudonym", "bullseye_c", "region", "depth", "WMH_vol_adj")],
                    dep="WMH_vol_adj","pca",nfactors=4, "oblimin", ncol=2)
fa.parallel(pca_LIFE[[1]],fa = "pc", n.obs=length(unique(res_LIFE[res_LIFE$tp=="bl",]$pseudonym)))

pca_LIFE[[3]]
ggsave("/data/pt_life_whm/Analysis/VRF-and-progression-of-WML/data_analysis/manuscript/figures//PCA_4factors_oblimin.png", 
       width=51, height=28, units="cm", dpi=400)
#predict scores for baseline and followup
ppca_fu=psych::predict.psych(pca_LIFE[[2]], data=res_LIFE_fu[,c(2:37)],old.data=res_LIFE_bl[,c(2:37)])

#It does not make a difference to use the predicted values without old values as reference.
#ppca_fu2=psych::predict.psych(pca_LIFE[[2]], data=res_LIFE_fu[,c(2:37)])

pscores=data.frame(scores.bl=pca_LIFE[[2]]$scores,
                   scores.fu=ppca_fu,
                   pseudonym=unique(res_LIFE_bl$pseudonym))
write.csv(pscores,
          row.names=F,
          "/data/pt_life_whm/Results/VRF_cSVD/LME/exploratory/proj_LIFE_on 4C_scores_oblimin_LONG.data_LIFE.csv")

## add unclassified
tmp=res_LIFE %>% select(pseudonym:tp) %>% distinct()
colnames(tmp)[3]="time"

pscores_long=pscores%>%
  pivot_longer(cols = !pseudonym, names_to = c("time", "comp"), names_sep=".T")
pscores_long$time=as.factor(pscores_long$time)
levels(pscores_long$time)=c("bl","fu")

## Merge with originally imputed data
add_pheno=merge(pscores_long,tmp, by=c("time", "pseudonym"), all.y=T)

data_wm_w_p=merge(data_wm, add_pheno, by=c("time", "pseudonym"), all.x=T)



data_ana = data_wm_w_p %>% 
  pivot_wider(id_cols=c("pseudonym","subj", "order", "sex","education","EstimatedTotalIntraCranialVol","time","wml","undefined", "DBP","WHR",
                          "learning","recall","recognition","phon_f","sem_f","TMT","proc","BPmed","cesd",
                          "age","exfunct","memo","globalcog"),
                        names_from = "comp",
                        values_from="value")


mean <- data_ana%>% group_by(time)%>%summarise(mean_C1=mean(C1, na.rm=T),
                                               mean_C2=mean(C2, na.rm=T),
                                                  mean_C3=mean(C3, na.rm=T),
                                                  mean_C4=mean(C4, na.rm=T))
                                               #mean_C5=mean(C5, na.rm=T))

ggplot(data = data_ana, aes(x = age, y=undefined, group=pseudonym, color=time)) + geom_point(aes(shape=time), size=2) + geom_line() +
  xlab("Age") + ylab("undefined") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))+
  labs(shape="Time point", color="Age")

ggplot(data = data_ana, aes(x = time, y=C3, group=pseudonym, color=time)) + 
  geom_point(aes(shape=time), size=2) + 
  geom_line() +
  xlab("Age") + ylab("C3") + 
  geom_point(data= mean, mapping=aes(x=time, group=time, y = mean_C3), 
             color=c(5,5), size=4)+
  geom_line(data=mean, mapping=aes(x=time, group=1, y=mean_C3), col=5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))+
  labs(shape="Time point", color="Age")

ggplot(data = data_ana, aes(x = time, y=C4, group=pseudonym, color=time)) + 
  geom_point(aes(shape=time), size=2) + 
  geom_line() +
  xlab("Age") + ylab("C4") + 
  geom_point(data= mean, mapping=aes(x=time, group=time, y = mean_C4), 
             color=c(5,5), size=4)+
  geom_line(data=mean, mapping=aes(x=time, group=1, y=mean_C4), col=5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))+
  labs(shape="Time point", color="Age")

ggplot(data = data_ana, aes(x = time, y=C1, group=pseudonym, color=time)) + 
  geom_point(aes(shape=time), size=2) + 
  geom_line() +
  xlab("Age") + ylab("C1") + 
  geom_point(data= mean, mapping=aes(x=time, group=time, y = mean_C1), 
             color=c(5,5), size=4)+
  geom_line(data=mean, mapping=aes(x=time, group=1, y=mean_C1), col=5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))+
  labs(shape="Time point", color="Age")
ggplot(data = data_ana, aes(x = age_change, y=C1, group=pseudonym, color=time)) + 
  geom_point(aes(shape=time), size=2) + 
  geom_line() +
  xlab("Age") + ylab("C1") + 
  geom_point(data= mean, mapping=aes(x=time, group=time, y = mean_C1), 
             color=c(5,5), size=4)+
  geom_line(data=mean, mapping=aes(x=time, group=1, y=mean_C1), col=5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))+
  labs(shape="Time point", color="Age")

ggplot(data = data_ana %>% pivot_longer(cols=(C2:C4), names_to = "comp", values_to = "value" ), aes(x = age, y=value, group=pseudonym, color=time)) + geom_point(aes(shape=time), size=2) + geom_line() +
  xlab("Age") + ylab("Value") + facet_wrap(~comp, ncol=2)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))+
  labs(shape="Time point", color="Age")


#Merge with imputed data:
comp_imp=mice::complete(imp, "long")

data_short = data_ana %>%  
  arrange(as.numeric(order))%>% 
  select(c("subj","time","undefined","order",
           "C1", "C2", "C3","C4", "age"))%>%
  mutate(      c1_base = ifelse(time == "bl", C1, lag(C1)),
               c1_change = ifelse(time == "bl", 0, C1- c1_base),
               c3_base = ifelse(time == "bl", C3, lag(C3)),
               c3_change = ifelse(time == "bl", 0, C3- c3_base),
               #c5_base = ifelse(time == "bl", C5, lag(C5)),
               #c5_change = ifelse(time == "bl", 0, C5- c5_base),
               c2_base = ifelse(time == "bl", C2, lag(C2)),
               c2_change = ifelse(time == "bl", 0, C2- c2_base),
               c4_base = ifelse(time == "bl", C4, lag(C4)),
               c4_change = ifelse(time == "bl", 0, C4- c4_base),
               und_base = ifelse(time == "bl", asinh(undefined), lag(asinh(undefined))),
               und_change = ifelse(time == "bl", 0, asinh(undefined)- und_base))

#continue without imputation
imp_bound=mice::cbind(imp, data_short)

miceadds::write.mice.imputation(imp_bound , "imputed_data_withWMHcomponents", mids2spss=FALSE)


#### differences over time
diffC3 <- with(imp_bound, lmerTest::lmer(formula = 'C3 ~ age_base + age_change + TIV + (1|subj)'),
            REML=F, na.action = na.omit)
diffC1 <- with(imp_bound, lmerTest::lmer(formula = 'C1 ~ age_change + TIV + (1|subj)'),
               REML=F, na.action = na.omit)
diffC2 <- with(imp_bound, lmerTest::lmer(formula = 'C2 ~ age_change + TIV + (1|subj)'),
               REML=F, na.action = na.omit)
diffC4 <- with(imp_bound, lmerTest::lmer(formula = 'C4 ~ age_change + TIV + (1|subj)'),
               REML=F, na.action = na.omit)
diffC5 <- with(imp_bound, lmerTest::lmer(formula = 'C5 ~ age_change + TIV + (1|subj)'),
               REML=F, na.action = na.omit)
#WMH volume
res <- with(imp_bound, lmerTest::lmer(formula = 'asinh_wml ~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change +
                                    BPmed + TIV + (1|subj)'),
            REML=F, na.action = na.omit)
est=summary(mice::pool(res))


## Calc M1 with C3/C5
resC=list()
i=1
dbpchange_est=c()
dbpbase_est=c()
for (comp in c(paste0("C", seq(1:5)))){
  
  resC[[i]] <-with(imp_bound, lmerTest::lmer(formula=paste0(comp,'~ age_base + age_change + 
                                 DBP_base + age_change:DBP_base + DBP_change  + 
                                 WHR_base + WHR_base:age_change + WHR_change +
                                    BPmed + TIV + (1|subj)'),
                                        REML=F, na.action = na.omit))
 
  
  est=summary(mice::pool(resC[[i]]))
  i=i+1
}

summary(mice::pool(resC[[1]]))
summary(mice::pool(resC[[2]]))
summary(mice::pool(resC[[3]]))
summary(mice::pool(resC[[4]]))
summary(mice::pool(resC[[5]]))


## C3
resex=list()
resgc=list()
i=1

for (comp in c(paste0("C", seq(1:5)))){
  
  resex[[i]] <-with(imp_bound, lmerTest::lmer(formula = paste0('exfunct ~ age_base + age_change + 
                                 c',i,'_base + c',i,'_change + 
                                  education + cesd + TIV +
                             (1|subj)'),REML=F, na.action = na.omit))
  resgc[[i]] <-with(imp_bound, lmerTest::lmer(formula = paste0('globalcog ~ age_base + age_change + 
                                 c',i,'_base + c',i,'_change + 
                                  education + cesd + TIV +
                             (1|subj)'),REML=F, na.action = na.omit))
  i=i+1
  }

summary(mice::pool(resex[[1]]))
summary(mice::pool(resex[[2]]))
summary(mice::pool(resex[[3]]))
summary(mice::pool(resex[[4]]))
summary(mice::pool(resex[[5]]))

summary(mice::pool(resgc[[1]]))
summary(mice::pool(resgc[[2]]))
summary(mice::pool(resgc[[3]]))
summary(mice::pool(resgc[[4]]))
summary(mice::pool(resgc[[5]]))

resgc <-with(imp_bound, lmerTest::lmer(formula = paste0('globalcog ~ age_base + age_change + 
                                 asinh_wml_base + asinh_wml_change + 
                                  education + cesd + TIV +
                             (1|subj)'),REML=F, na.action = na.omit))

resund <-with(imp_bound, lmerTest::lmer(formula = paste0('globalcog ~ age_base + age_change + 
                                 und_base + und_change + 
                                  education + cesd + TIV +
                             (1|subj)'),REML=F, na.action = na.omit))
