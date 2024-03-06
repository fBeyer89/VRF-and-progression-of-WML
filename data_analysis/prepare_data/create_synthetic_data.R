# Create synthetic data for publication
library(mice)
library(dplyr)
library(miceadds)


#Starting with the original dataset
wml_long_excl_a=read.csv("/data/pt_life_whm/Results/VRF_cSVD/imputed_data/imputed_data_6.12.23/full_dataset_long_excl_a.csv")

#Prepare data frame for imputation
test = wml_long_excl_a %>% 
  mutate(subj=as.integer(as.factor(mrt_pseudonym)),
         time=as.factor(time),
         education=education)%>% 
  select(c("subj","time", "age", "sex","education","TIV", "wml","DBP","WHR","globalcog", "BPmed", "cesd"))

# Set default predictor matrix
ini <- mice::mice(test, maxit=0)
ini$loggedEvents #logged events point out issues we know

# create imputation model with predictor matrix (following 'type')
# indicators:
#-2: cluster variable, 1: fixef with randint, 0: not used
# 2: random effect of predictor which we do not have with two timepoints

#Imputes all variables to preserve correlation structure
colnames(test)

pred <- ini$predictorMatrix
pred["age",] <-        c(-2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
pred["sex",] <-        c(-2, 1, 1, 0, 0, 0 ,0, 1, 1, 0, 0, 0)
pred["DBP",] <-        c(-2, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0)
pred["wml",] <-        c(-2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
pred["WHR",] <-        c(-2, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0)
pred["TIV",] <-        c(-2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
pred["globalcog",] <-  c(-2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
pred["education",] <-  c(-2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
pred["BPmed",] <-  c(-2, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0)
pred["cesd",] <-  c(-2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

# set imputation method 
impmethod <- c("", "",  #do not impute subj and timepoint
               "2l.lmer", #use lmer mixed model for age, DBP, WHR, WMH, global cognition (vary over time)
               "2lonly.pmm", "2lonly.pmm","2lonly.pmm", #second level imputation  of sex, education and TIV only
               "2l.lmer", "2l.lmer", "2l.lmer", "2l.lmer",#use lmer mixed model for  DBP, WHR, WMH, global cognition
               "2l.lmer", "2l.lmer") #second level imputation of BP & cesd
names(impmethod) <- colnames(test)

# where to impute -> everywhere except subj and time
tmp=matrix(TRUE, nrow(test), ncol(test))
tmp[,c(1:2)]=FALSE

# using mice to impute the full dataset -> https://www.mdpi.com/2624-8611/3/4/45
m = 5
iter = 10
imp <- mice::mice(
  test,
  method = impmethod,
  pred = pred,
  maxit = iter,
  where=tmp,
  m = m,
  seed = 8745
)
# where=tmp,

#See results
plot(imp)

#make original data all NA so that individuals cannot be identified
tmp_impdata=imp$data
tmp_impdata[,c(3:ncol(tmp_impdata))]=NA
imp$data=tmp_impdata

#Check the imputed data
imp_l <- mice::complete(imp, action = 'long', include = TRUE) #include = TRUE option to be able to save to mids later

############################################################################
# calculate between and within measures of Age, SBP, WHR, asinh(WML)
############################################################################

imp_l <- imp_l %>%
  group_by(subj, .imp) %>%
  mutate(
    subj=as.factor(subj),
    age_base = ifelse(time == "bl", age, lag(age)),
    DBP_base = ifelse(time == "bl", DBP, lag(DBP)),
    WHR_base = ifelse(time == "bl", WHR, lag(WHR)),
    age_change = ifelse(time == "bl", 0, age - age_base),
    DBP_change = ifelse(time == "bl", 0, DBP - DBP_base),
    WHR_change = ifelse(time == "bl", 0, WHR - WHR_base),
    asinh_wml = asinh(wml),
    asinh_wml_base = ifelse(time == "bl", asinh_wml, lag(asinh_wml)),
    asinh_wml_change = ifelse(time == "bl", 0, asinh_wml - asinh_wml_base))



imp.itt <- mice::as.mids(imp_l)


#write full dataset
setwd("/data/pt_life_whm/Results/VRF_cSVD/imputed_data/")
miceadds::write.mice.imputation(imp.itt , "synthetic_imputed_data", mids2spss=FALSE)


imp_l_NA <- mice::complete(imp.itt, action = 'long', include = FALSE)
