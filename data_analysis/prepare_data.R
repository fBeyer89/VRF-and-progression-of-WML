#library(kit)
library(readxl) #read xlsx, ods etc.
library(tidyverse) #includes ggplot2, dplyr, tidyr, readr, tibble, stringr ...
library(lubridate) #package to work with data and time
#library(eeptools)
#library(data.table)
library(naniar) # helps to deal with NA
#library(lme4)
#library(stringr)

#Prepare data for CVR-cSVD RR
#Collect all variables in wide format and then transform to long

############################################################################
#  Load WML data (as only participants with both time point WML data will be included)
############################################################################
wml=read.csv("/data/pt_life_whm/Results/Tables/longvols_w_pseudonym_qa.csv")
colnames(wml)[1]="mrt_pseudonym"
############################################################################
#  Add PV-specific pseudonyms to the MRI pseudonyms
############################################################################
pv_ids=read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/PV609_mrt_pseudonyme_2023-02-03.xlsx")
wml=merge(wml[,c("mrt_pseudonym", "qa_check", "qa_comment", "wml_vol_bl", "wml_vol_fu")], pv_ids, by="mrt_pseudonym", all.x=T)

############################################################################
#  Add SIC pseudonyms to the MRI pseudonyms
############################################################################
sic_pseudo=read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/pseudo_mrt_20201214.xlsx")
wml=merge(wml, sic_pseudo, all.x=T, by.x="mrt_pseudonym", by.y="pseudonym")

############################################################################
#  Load baseline demographics
############################################################################
# load data on participants' sex and birth date
basic <- readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_R00001.xlsx")

# rename variables to allow merging and merge
basic <- basic %>%
  rename(
    pv_pseudonym = TEILNEHMER_SIC
  )
# define sex: female is 0, male is 1
basic$sex <- ifelse(basic$TEILNEHMER_GESCHLECHT == 2, 0, basic$TEILNEHMER_GESCHLECHT)
# birthmonth and year are available -- turn them into an usable form
basic$birth <- as.Date(parse_date_time(basic$TEILNEHMER_GEB_JJJJMM, "ym"))

wml=merge(wml, basic[,c("pv_pseudonym", "sex", "birth")], by="pv_pseudonym", all.x=T)

############################################################################
#  Load education (contains duplicates)
############################################################################
demo <- read_excel("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_D00140_NODUP.xlsx")
#remove duplicates (assuming that education was the same at both timepoints)
demo = demo[!duplicated(demo$SIC),]
# categorise participants with less than a tertiary degree (score below 3.6) as 1 in a variable called education
demo$education <- ifelse(demo$SES2_SESBLDG < 3.6, 1, 0)
colnames(demo)[1]="pv_pseudonym"

wml=merge(wml, demo[,c("pv_pseudonym", "education")], by = "pv_pseudonym", all.x=T)

##########################################################################
# Load MRI dates (for LST, the earlier scan ("bl") was always selected)
MRI_bl=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00197_NODUP.xlsx")
MRI_fu=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1//PV0609_T01158_NODUP.xlsx")
MRI_bl$MRT_DATUM_bl <- as.Date(MRI_bl$MRT_DATUM)

# For controlling that earlier time point was selected.
MRI_bl_duplicates <-
  MRI_bl  %>% 
  arrange(MRT_DATUM_bl) %>%
  filter(duplicated(.[["MRT_SIC"]]))

MRI_bl_wo_duplicates <-
  MRI_bl  %>% 
  arrange(MRT_DATUM_bl) %>%
  filter(!duplicated(.[["MRT_SIC"]]))

MRI_fu$MRT_DATUM_fu <- as.Date(MRI_fu$MRT_ZW_U_DATUM)

wml=merge(wml, MRI_bl_wo_duplicates[,c("MRT_SIC","MRT_DATUM_bl")], by.x="pv_pseudonym", by.y="MRT_SIC", all.x=T)
wml=merge(wml, MRI_fu[,c("MRT_ZW_U_SIC","MRT_DATUM_fu")], by.x="pv_pseudonym", by.y="MRT_ZW_U_SIC", all.x=T, by="pseudonym")

############################################################################
#  Load medication, medical anamnese and radiological assessment data
#  for Exclusion
############################################################################
medanam_bl=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00173_NODUP.xlsx")
medanam_fu=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T00173_NODUP.xlsx")
#("T01228.xlsx")
kardanam_fu <- read_excel("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_kardanam_20210511.xlsx")

# rename medical anamnesis labels to make them understandable
medanam_bl <- medanam_bl %>%
    rename(
    epilepsy_bl = MEDANAM_F0167, 
    parkinson_bl = MEDANAM_F0171,
    MS_bl = MEDANAM_F0179,
    stroke_bl = MEDANAM_F0024,
    hyp_treat_bl = MEDANAM_F0039
  )
#medanam_bl_duplicates <-
#  medanam_bl  %>% 
#  arrange(MEDANAM_DATUM) %>%
#  filter(duplicated(.[["SIC"]]))
medanam_bl_wo_duplicates <-
  medanam_bl  %>% 
  arrange(MEDANAM_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

medanam_fu <- medanam_fu %>%
  rename(
    other_diseases = mediz_an.f31a, 
    parkinson_fu = mediz_an.f27,
    MS_fu = mediz_an.f26,
    dementia_fu = mediz_an.f30,
    dat_medanam_fu = edat)
medanam_fu <- medanam_fu %>%
  mutate(epilepsy_fu = case_when(
    other_diseases == "EPILEPSIE" ~ 1,
    TRUE ~ 0
  ))

#Duplicates in medanam_fu (use only those from FU_A1, not A1_FU1_ANSCH which has
#been usually performed before the 2nd
#MRI, maybe this refers to postal questionnaires)
medanam_fu_wo_duplicates <-
  medanam_fu  %>% 
  filter(sgroup!="A1_FU1_ANSCH") %>% 
  filter(sgroup!="A1_F1_MRT_VERGLEICH1")

kardanam_fu <- kardanam_fu%>%
  rename(stroke_fu = kard_an.f7,
         hyp_treat_fu = kard_an.f10,
         dat_kardanam_fu = edat)

kardanam_fu_wo_duplicates <-
  kardanam_fu  %>% 
  filter(sgroup!="A1_FU1_ANSCH") %>% 
  filter(sgroup!="A1_F1_MRT_VERGLEICH1")

wml=merge(wml, medanam_bl_wo_duplicates[,c("SIC", "epilepsy_bl", "stroke_bl", "MS_bl", "parkinson_bl", "hyp_treat_bl")], by.x="pv_pseudonym", by.y="SIC", all.x=T)
wml=merge(wml, medanam_fu_wo_duplicates[,c("pseudonym","dat_medanam_fu", "epilepsy_fu", "MS_fu", "parkinson_fu", "dementia_fu")], by.x="pv_pseudonym", by.y="pseudonym",all.x=T)
wml=merge(wml, kardanam_fu_wo_duplicates[,c("pseudonym", "dat_kardanam_fu", "stroke_fu", "hyp_treat_fu")], by.x="pv_pseudonym", by.y="pseudonym",all.x=T)

## MMSE score
sidam_bl=read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00043_NODUP.xlsx")
colnames(sidam_bl)[1]="SIC"
a <- c("SIDAM1_")
b <- c("1", "2", "3", "4", "5", 
       "7", "8", "9", "10", "11", 
       "12A","12B", "12C", "12D", "12E", 
       "14", "15", "6A", "6B", "6C",
       "16A", "16B", "16C", "35A", "35B",
       "36", "37", "40A", "40B", "40C")
sidam_cols <- paste0(rep(a, each = length(b)), b)
sidam_bl$MMSE = rowSums(sidam_bl[,sidam_cols])
sidam_bl = sidam_bl %>%
  mutate(MMSE = case_when(
    MMSE >30 ~ NA_real_,
    TRUE ~ MMSE)) %>%
  mutate(MMSE_exclude_bl = case_when(
    MMSE < 24 ~ 1,
    TRUE  ~ 0))

sidam_bl_wo_duplicates <-
  sidam_bl  %>% 
  arrange(SIDAM1_STARTZEIT) %>%
  filter(!duplicated(.[["SIC"]]))

#Check that derived version and self-calculated version agree:
# sidam_bl_sum=read_xlsx("/data/pt_life/ResearchProjects/LLammer/Data/data/Baseline/PV0573_D00060_NODUP.xlsx")
# pseudo <- read_excel("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/pseudo_mrt_20201214.xlsx")
# colnames(pseudo)[1]="SIC"
# colnames(pseudo)[2]="mrt_pseudonym"
# sidam_bl_pseudo=merge(sidam_bl, pseudo, all.x=T)
# 
# colnames(sidam_bl_sum)[1]="pv_pseudonym"
# pseudo <- read_excel("/data/pt_life/ResearchProjects/LLammer/Data/PV573_PV-MRT-Pseudonymliste.xlsx")
# sidam_bl_sum_pseudo=merge(sidam_bl_sum, pseudo)
# test=merge(sidam_bl_pseudo[,c("mrt_pseudonym", "MMSE")],sidam_bl_sum_pseudo[,c("mrt_pseudonym", "MMST_MMST")])
# plot(test$MMSE, test$MMST_MMST)

###
sidam_fu=read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T00043_NODUP.xlsx")
a <- c("SIDAM1_")
b <- c("1", "2", "3", "4", "5", 
       "7", "8", "9", "10", "11", 
       "12A","12B", "12C", "12D", "12E", 
       "14", "15", "6A", "6B", "6C",
       "16A", "16B", "16C", "35A", "35B",
       "36", "37", "40A", "40B", "40C")
sidam_cols <- paste0(rep(a, each = length(b)), b)
sidam_fu$MMSE = rowSums(sidam_fu[,sidam_cols])
sidam_fu = sidam_fu %>%
  mutate(MMSE = case_when(
    MMSE >30 ~ NA_real_,
    TRUE ~ MMSE)) %>%
  mutate(MMSE_exclude_fu = case_when(
    MMSE < 24 ~ 1,
    TRUE  ~ 0))

sidam_fu_wo_duplicates <-
  sidam_fu  %>% 
  arrange(SIDAM1_STARTZEIT) %>%
  filter(!duplicated(.[["SIDAM1_SIC"]]))

wml=merge(wml, sidam_bl_wo_duplicates[,c("SIC", "MMSE_exclude_bl")], by.x="sic", by.y="SIC", all.x=T)
wml=merge(wml, sidam_fu_wo_duplicates[,c("SIDAM1_SIC", "MMSE_exclude_fu")], by.x="sic", by.y="SIDAM1_SIC", all.x=T)

############################################################################
# load data from radiological assessments (is already in long format, merge later)
############################################################################
radio <- read.csv("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/radiological_assessment/radio_assessment_with_followup.csv", na.strings = c("NA",""))

# exclude participants woith tumors and non-congenital lesions
radio <- radio %>%
  mutate(medBefund_date = as.Date(med_Befund.Datum_Bef, format="%d.%m.%Y")) %>%
  mutate(excludetumor = case_when(
    med_Befund..Tumours == "keine Kontrolle notwendig" ~ 0,
    is.na(med_Befund..Tumours) ~ 0,
    TRUE ~ 1
  )) %>%
  mutate(excludelesion = case_when(
    med_Befund..Lesions == "kongenital" ~ 0,
    is.na(med_Befund..Lesions) ~ 0,
    TRUE ~ 1
  ))%>%
  mutate(exclude_usable = case_when(
    med_Befund.Bewertung_Bef == "verwendungsf\xe4hig nein" ~ 1,
    TRUE ~ 0
  ))

radio_fu <-
  radio %>%
  arrange(medBefund_date) %>%
  filter(duplicated(.[["pseudonym"]])) %>%
  filter(medBefund_date>"2015-01-01")
colnames(radio_fu)=paste0(colnames(radio_fu),"_fu")

radio_bl <-
  radio %>%
  arrange(medBefund_date) %>%
  filter(!duplicated(.[["pseudonym"]])) %>%
  filter(medBefund_date<"2018-01-01")

colnames(radio_bl)=paste0(colnames(radio_bl),"_bl")

wml=merge(wml, radio_bl[,c("pseudonym_bl", "exclude_usable_bl", 
                           "excludelesion_bl", "excludetumor_bl")], 
          by.x="mrt_pseudonym", by.y="pseudonym_bl", all.x=T)
wml=merge(wml, radio_fu[,c("pseudonym_fu", "exclude_usable_fu", 
                           "excludelesion_fu", "excludetumor_fu")], 
          by.x="mrt_pseudonym", by.y="pseudonym_fu", all.x=T)
############################################################################
#Medication Anamnese
############################################################################
# Exclude centrally active drugs begin with the following ATC-codes: 
# M03B - centrally active muscle relexants
# N06D - antidementive drugs but we will not exclude particpants taking Ginko (N06DP)
# N05A - antipsychotic drugs
# N05B - anxiolytic drugs but we will not exclude particpants taking lavendar oil (N05BP03)
# N06B - psychostimulative drugs, ADHD drugs & nootropic drugs
# N05C - sedative and hypnotic drugs but we will not exclude participants taking drugs based on lemon balm, passiflora incarnata, valerian and humulus
# N06A - antidepressants but we will not exclude participants taking homeopathic antidepressants
# N02A - opioids
# N03 - antiepileptic drugs
# N07A - parasympathomimetics
# N02C - migraine remedies but we will not exclude participants taking homeopathic or plant-based migraine remedies 
# N07C - vertigo medications but we will not exclude participants taking homeopathic or anthroposophic antivertiginosa (N07CH)
# N04 - parkinson's disease medication
# R05DA - opium alkaloids and derivatives

medianam_bl<- read_excel("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_D00038_NODUP.xlsx")

# make sure no hashtags are around ACT-codes
medianam_bl$ADULT_MEDA_H_ATC <- str_replace_all(medianam_bl$ADULT_MEDA_H_ATC, "#", "") 

to_match <- c("^M03B", "M03B", 
              "^N02A", "N02A", "^N03", "N03", "^N07A", "N07A", "^N07CA", 
              "^N02CA", "N02CA", "^N02CB", "N02CB", "^N02CC", "N02CC", "^N02CD", "N02CD", "^N02CX", "N02CX", 
              "^N04", "N04",
              "^N06DA", "N06DA", "^N06DX01", "N06DX01", "^N06DX30", "N06DX30", 
              "^N06B", "N06B", 
              "^N06AA", "N06AA", "^N06AB", "N06AB", "^N06AF", "N06AF", "^N06AG", "N06AG", "^N06AP", "N06AP", "^N06AX", "N06AX", 
              "^N05A",  "N05A", "^N05CA", "N05CA", "^N05CB", 
              "^N05BA", "N05BA", "^N05BB", "N05BB", "^N05BC", "N05BC", "^N05BD", "N05BD",
              "^N05BE", "N05BE", "^N05BX", "N05BX", "^N05BP02", "N05BP02",
              "N05CB", "^N05CC", "N05CC", "^N05CD", "N05CD", "^N05CE", 
              "N05CE", "^N05CF", "N05CF", "^N05CH", "N05CH", "^N05CM", 
              "N05CM", "^N05CX", "N05CX", "^N05CP02", "N05CP02", "^N05CP03", ", N05CP03", 
              "N07CA",  
              "^R05DA", "R05DA")
medianam_bl <- medianam_bl %>%
  mutate(centr_act_med_bl = case_when(
    grepl(paste0(to_match, collapse = "|"), ADULT_MEDA_H_ATC) ~ 1,
    TRUE ~ 0
  ))


# in the ATC nomenclature all antihypertensive drugs begin with "C02", all diuretics begin with "C03", 
# all beta-blockers with "C07", all Ca-chanel blockers with "C08" and all RAAS affecting drugs with "C09"  
# all antidiabetic drugs begin with "A10" in the ATC nomenclature
medianam_bl <- medianam_bl %>%
  mutate(BPmed_bl = case_when(
    grepl("^C02", ADULT_MEDA_H_ATC) | grepl(", C02", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C03", ADULT_MEDA_H_ATC) | grepl(", C03", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C07", ADULT_MEDA_H_ATC) | grepl(", C07", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C08", ADULT_MEDA_H_ATC) | grepl(", C08", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C09", ADULT_MEDA_H_ATC) | grepl(", C09", ADULT_MEDA_H_ATC) ~ 1,
    is.na(ADULT_MEDA_H_ATC) ~ NA_real_,
    TRUE ~ 0
  ))

#Remove duplicates
medianam_bl_wo_duplicates <-
  medianam_bl %>%
  arrange(ADULT_MEDA_H_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

# Identify cases to exclude at followup 
medianam_fu <- read_excel("/data/pt_life/ResearchProjects/LLammer/Data/data/Followup_update/PV0573_IDOM_Medikamente_horizontal.xlsx")

medianam_fu$ADULT_MEDA_H_ATC <- str_replace_all(medianam_fu$ATC, ",", ", ")
medianam_fu$ADULT_MEDA_H_ATC <- str_replace_all(medianam_fu$ADULT_MEDA_H_ATC, "#", "") 
to_match <- c("^M03B", "M03B", 
              "^N02A", "N02A", "^N03", "N03", "^N07A", "N07A", "^N07CA", 
              "^N02CA", "N02CA", "^N02CB", "N02CB", "^N02CC", "N02CC", "^N02CD", "N02CD", "^N02CX", "N02CX", 
              "^N04", "N04",
              "^N06DA", "N06DA", "^N06DX01", "N06DX01", "^N06DX30", "N06DX30", 
              "^N06B", "N06B", 
              "^N06AA", "N06AA", "^N06AB", "N06AB", "^N06AF", "N06AF", "^N06AG", "N06AG", "^N06AP", "N06AP", "^N06AX", "N06AX", 
              "^N05A",  "N05A", "^N05CA", "N05CA", "^N05CB", 
              "^N05BA", "N05BA", "^N05BB", "N05BB", "^N05BC", "N05BC", "^N05BD", "N05BD",
              "^N05BE", "N05BE", "^N05BX", "N05BX", "^N05BP02", "N05BP02",
              "N05CB", "^N05CC", "N05CC", "^N05CD", "N05CD", "^N05CE", 
              "N05CE", "^N05CF", "N05CF", "^N05CH", "N05CH", "^N05CM", 
              "N05CM", "^N05CX", "N05CX", "^N05CP02", "N05CP02", "^N05CP03", ", N05CP03", 
              "N07CA",  
              "^R05DA", "R05DA")
medianam_fu<- medianam_fu %>%
  mutate(centr_act_med_fu = case_when(
    grepl(paste0(to_match, collapse = "|"), ADULT_MEDA_H_ATC) ~ 1,
    TRUE ~ 0
  ))

medianam_fu <- medianam_fu %>%
  mutate(BPmed_fu = case_when(
    grepl("^C02", ADULT_MEDA_H_ATC) | grepl(", C02", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C03", ADULT_MEDA_H_ATC) | grepl(", C03", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C07", ADULT_MEDA_H_ATC) | grepl(", C07", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C08", ADULT_MEDA_H_ATC) | grepl(", C08", ADULT_MEDA_H_ATC) ~ 1,
    grepl("^C09", ADULT_MEDA_H_ATC) | grepl(", C09", ADULT_MEDA_H_ATC) ~ 1,
    is.na(ADULT_MEDA_H_ATC) ~ NA_real_,
    TRUE ~ 0
  ))

#Remove duplicates (which are truely duplicates, not repeated measures)
medianam_fu_wo_duplicates <-
  medianam_fu %>%
  arrange(DATUM) %>%
  filter(!duplicated(.[["PSEUDONYM"]]))

wml=merge(wml, medianam_bl_wo_duplicates[,c("SIC", "centr_act_med_bl", "BPmed_bl")], by.x="pv_pseudonym", by.y="SIC",all.x=T)
wml=merge(wml, medianam_fu_wo_duplicates[,c("PSEUDONYM", "centr_act_med_fu", "BPmed_fu")], by.x="pv_pseudonym",by.y="PSEUDONYM", all.x=T)


############################################################################
#Predictors of interest
#Rule for missing/implausible/wrong values: impute in multi-level imputation
#Exclusion of participant if both time points are missing.
############################################################################
#Blood pressure (after revision DBP; from derivative not raw data)
bp_bl=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_D00073_NODUP.xlsx")

# Replace all missing or wrong values with NA
bp_bl<- replace_with_na_at(bp_bl, c("ADULT_BP_DBP", "ADULT_BP_SBP"), ~.x %in% c(998, 996))

#bp_bl$SBP_bl <- (bp_bl$BLUTDRUCKMESS_F0024 + bp_bl$BLUTDRUCKMESS_F0019 + bp_bl$BLUTDRUCKMESS_F0014)/3
#bp_bl$DBP_bl <- (bp_bl$BLUTDRUCKMESS_F0025 + bp_bl$BLUTDRUCKMESS_F0020 + bp_bl$BLUTDRUCKMESS_F0015)/3 #VERIFY

bp_bl= bp_bl %>%
  mutate(ADULT_BP_DBP = case_when(
    ADULT_BP_DBP > 140  ~ NA_real_,
    ADULT_BP_SBP < ADULT_BP_DBP ~ NA_real_,
    TRUE ~ ADULT_BP_DBP))

#Remove duplicates
bp_bl_wo_duplicates <-
  bp_bl %>%
  arrange(ADULT_BP_S010063_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

# SBP followup
bp_fu=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T01170_NODUP.xlsx")

# Replace all missing or wrong values with NA
bp_fu<- replace_with_na_at(bp_fu, c("BP_F0024", "BP_F0019", "BP_F0014",
                                    "BP_F0015", "BP_F0020", "BP_F0025"), ~.x %in% c(998, 996))

bp_fu$SBP_fu <- (bp_fu$BP_F0024 + bp_fu$BP_F0019 + bp_fu$BP_F0014)/3
bp_fu$DBP_fu <- (bp_fu$BP_F0025 + bp_fu$BP_F0020 + bp_fu$BP_F0015)/3 #VERIFY

bp_fu= bp_fu %>%
  mutate(DBP_fu = case_when(
    DBP_fu > 140  ~ NA_real_,
    SBP_fu < DBP_fu ~ NA_real_,
    TRUE ~ DBP_fu))

#Remove duplicates ##NEEDS TO BE CHECKED FOR 2nd BP measurement
bp_fu_wo_duplicates <-
  bp_fu %>%
  arrange(BP_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

#WHR
#Using derivative
whr_bl=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_D00074_NODUP.xlsx")

#Using Raw data
#whr_bl=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2019_PV465_Evelyn/data/PV0465_T00047_NODUP.xlsx")

# Replace all missing or wrong values in waist and hip with NA
#whr_bl<- replace_with_na_at(whr_bl, c("ANTHRO_F0012",  
#                                      "ANTHRO_F0015"), ~.x %in% c(999, 998, 996))

#whr_bl$WHR_bl <- (whr_bl$ANTHRO_F0012)/
#                 (whr_bl$ANTHRO_F0015)

# Identify biologically implausible values to be imputed
whr_bl= whr_bl %>%
  mutate(whr_bl = case_when(
    BMI_WAIST_HIP_RATIO < 0.5  ~ NA_real_,
    BMI_WAIST_HIP_RATIO > 1.5 ~ NA_real_,
    TRUE ~ BMI_WAIST_HIP_RATIO))

whr_bl_wo_duplicates <-
  whr_bl %>%
  arrange(BMI_S010061_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

####
whr_fu=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T01169_NODUP.xlsx")

# Replace all missing or wrong values in waist and hip with NA
whr_fu<- replace_with_na_at(whr_fu, c("ANTHRO_F0027", "ANTHRO_F0028"), ~.x %in% c(999, 996))

whr_fu$WHR_fu <- (whr_fu$ANTHRO_F0027/whr_fu$ANTHRO_F0028)

# Identify biologically implausible values to be imputed
whr_fu= whr_fu %>%
  mutate(WHR_fu = case_when(
    WHR_fu < 0.5  ~ NA_real_,
    WHR_fu > 1.5 ~ NA_real_,
    TRUE ~ WHR_fu))

wml=merge(wml, bp_bl_wo_duplicates[,c("SIC", "SBP_bl")], by.x="pv_pseudonym", by.y="SIC", all.x=T)
wml=merge(wml, bp_fu_wo_duplicates[,c("SIC", "SBP_fu")], by.x="pv_pseudonym", by.y="SIC", all.x=T)
wml=merge(wml, whr_bl_wo_duplicates[,c("SIC", "WHR_bl")], by.x="pv_pseudonym", by.y="SIC", all.x=T)
wml=merge(wml, whr_fu[,c("SIC", "WHR_fu")], by.x="pv_pseudonym", by.y="SIC", all.x=T)

#Cognitive Scores
cerad_bl=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00044_NODUP.xlsx") 
cerad_bl_anim=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00042_NODUP.xlsx") 
cerad_bl_tmt=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00041_NODUP.xlsx")
#cerad_bl_young=readxl::read_xlsx("/data/pt_life/ResearchProjects/LLammer/Data/data/Baseline/PV0573_T00195_NODUP.xlsx")

#remove duplicates
cerad_bl_wo_duplicates <-
  cerad_bl %>%
  arrange(CERAD_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))
cerad_bl_anim_wo_duplicates <-
  cerad_bl_anim %>%
  arrange(CERAD_VF_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))
cerad_bl_tmt_wo_duplicates <-
  cerad_bl_tmt %>%
  arrange(TMT_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))
#cerad_bl_young_wo_duplicates <-
#  cerad_bl_young %>%
#  arrange(WORTLISTE_DATUM) %>%
#  filter(!duplicated(.[["SIC"]]))
#colnames(cerad_bl_young_wo_duplicates)[2]="CERAD_DATUM"

#Only duplicated in cerad_fu (where the earlier one seems to be the correct one)
cerad_fu=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T00044_NODUP.xlsx") 
cerad_fu_anim=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T00042_NODUP.xlsx") 
cerad_fu_tmt=readxl::read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T00041_NODUP.xlsx") 

cerad_fu_wo_duplicates <-
  cerad_fu %>%
  arrange(CERAD_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

############################################################
# Functions to calculate cognition from questionnaires
############################################################

calc_cog_outcomes<-function(cerad, cerad_young, cerad_anim, tmt){
########
# Memory
########
# determine groups of columns belonging to the different cognitive tests 
a <- c("CERAD_WL1_", "CERAD_WL2_", "CERAD_WL3_")
b <- c("BUTTER", "ARM", "STRAND", "BRIEF", "KONIGI", "HUTTE", "STANGE", "KARTE", "GRAS", "MOTOR")
c <- c("KIRCHE", "KAFFEE", "BUTTER", "DOLLAR", "ARM", "STRAND", "FUNF", "BRIEF", "HOTEL", "BERG", "KONIGI", "HUTTE", 
       "PANTOF", "STANGE", "DORF", "BAND", "KARTE", "HEER", "GRAS", "MOTOR")
learning_cols <- paste0(rep(a, each = length(b)), b)
recall_cols <- paste0("CERAD_WL4_", b)
recognition_cols <- paste0("CERAD_WLW_", c)

cerad <- replace_with_na_at(cerad, .vars = c(recognition_cols), 
                            condition =  ~.x %in% c(97, 98, 99))
cerad$learning <- rowSums(cerad[,learning_cols]) 
cerad$recall <- rowSums(cerad[, recall_cols])
cerad$recognition <- rowSums(cerad[, recognition_cols])

# some learning and recall scores have to be turned into NAs, 
# because the test was marked as "does not apply" (97), "unrateable" (98) 
# or "refused to respond" (99) in a separate column
d <- c("97", "98", "99")
learning_na_cols <- paste0(rep(c(a), each = length(d)), d)
recall_na_cols <- paste0(rep(c("CERAD_WL4_"), each = length(d)), d)
cerad$learning_na <- rowSums(cerad[,learning_na_cols], na.rm = T)
cerad$recall_na <- rowSums(cerad[,recall_na_cols], na.rm = T)
cerad$learning <- ifelse(cerad$learning_na > 0, NA, cerad$learning)
cerad$recall <- ifelse(cerad$recall_na > 0, NA, cerad$recall)

if (ncol(cerad_young)>0){
# CERAD scores are named differently for participants aged <= 60
a2 <- c("WORTLISTE_WL1_", "WORTLISTE_WL2_", "WORTLISTE_WL3_")
learning_cols2 <- paste0(rep(a2, each = length(b)), b)
recall_cols2 <- paste0("WORTLISTE_WL4_", b)
recognition_cols2 <- paste0("WORTLISTE_WLW_", c)

cerad_young <- replace_with_na_at(cerad_young, .vars = c(recognition_cols2), 
                            condition =  ~.x %in% c(97, 98, 99))
cerad_young$learning <- rowSums(cerad_young[,learning_cols2]) 
cerad_young$recall <- rowSums(cerad_young[, recall_cols2])
cerad_young$recognition <- rowSums(cerad_young[, recognition_cols2])

# some learning and recall scores have to be turned into NAs, 
# because the test was marked as "does not apply" (97), "unrateable" (98) 
# or "refused to respond" (99) in a separate column
d <- c("97", "98", "99")
learning_na_cols <- paste0(rep(c(a2), each = length(d)), d)
recall_na_cols <- paste0(rep(c("WORTLISTE_WL4_"), each = length(d)), d)
cerad_young$learning_na <- rowSums(cerad_young[,learning_na_cols], na.rm = T)
cerad_young$recall_na <- rowSums(cerad_young[,recall_na_cols], na.rm = T)
cerad_young$learning <- ifelse(cerad_young$learning_na > 0, NA, cerad_young$learning)
cerad_young$recall <- ifelse(cerad_young$recall_na > 0, NA, cerad_young$recall)

cog_complete=merge(cerad[,c("SIC", "CERAD_DATUM", "learning", "learning_na", "recall", "recall_na", "recognition")],
                   cerad_young[,c("SIC", "learning", "CERAD_DATUM", "learning_na", "recall", "recall_na", "recognition")], all=T)
}
else {
  colnames(cerad)[1]="SIC"
  cog_complete=cerad[,c("SIC", "CERAD_DATUM", "learning", "learning_na", "recall", "recall_na", "recognition")]
}
########
#Fluency
########
phon_flu_cols <- paste0("CERAD_S_", c(15, 30, 45, 60), "S_CORR")
sem_flu_cols <- paste0("CERAD_VF_", c(15, 30, 45, 60), "S_CORR")

# replace error codes with NA
cerad <- replace_with_na_at(cerad, .vars = c(phon_flu_cols), condition =  ~.x %in% c(97, 98, 99))
cerad_anim <- replace_with_na_at(cerad_anim, .vars = c(sem_flu_cols), condition = ~.x %in% c(997, 998, 999))

cerad$phon_flu <- rowSums(cerad[,phon_flu_cols])
cerad_anim$sem_flu <- rowSums(cerad_anim[, sem_flu_cols])
colnames(cerad_anim)[1]="SIC"

########
# Executive function
########
# replace error codes with NA
tmt <- replace_with_na_at(tmt, .vars = c("TMT_TIMEA", "TMT_TIMEB"), 
                         condition = ~.x %in% c(997, 998, 999))

# exclude participants with long reaction times
tmt <- replace_with_na_at(tmt, .vars = c("TMT_TIMEA", "TMT_TIMEB"), 
                         condition = ~.x > 300)
                         
# calculate TMT composite score and negate it to facilitate later calculations
tmt$TMT <- -(tmt$TMT_TIMEB/tmt$TMT_TIMEA)
tmt$proc <- -tmt$TMT_TIMEA
colnames(tmt)[1]="SIC"

#Filter out two cases who have participated in pilot and adult
cog_complete <-
  cog_complete %>%
  arrange(CERAD_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

cog_complete=merge(cog_complete,cerad[,c("SIC", "phon_flu")], all=T )
cog_complete=merge(cog_complete,cerad_anim[,c("SIC", "sem_flu")], all=T )
cog_complete=merge(cog_complete,tmt[,c('SIC', 'TMT', 'proc')], all=T )
return(cog_complete)
}

cog_bl=calc_cog_outcomes(cerad_bl_wo_duplicates, cerad_bl_young_wo_duplicates, 
                         cerad_bl_anim_wo_duplicates, cerad_bl_tmt_wo_duplicates)
colnames(cog_bl)[2:length(colnames(cog_bl))]=paste0(colnames(cog_bl)[2:length(colnames(cog_bl))],"_bl")

cog_fu=calc_cog_outcomes(cerad=cerad_fu_wo_duplicates, cerad_young=data.frame(), 
                         cerad_anim=cerad_fu_anim, 
                         tmt=cerad_fu_tmt)

colnames(cog_fu)[2:length(colnames(cog_fu))]=paste0(colnames(cog_fu)[2:length(colnames(cog_fu))],"_fu")

wml=merge(wml, cog_bl, by.x="pv_pseudonym", by.y="SIC", all.x=T)
wml=merge(wml, cog_fu, by.x="pv_pseudonym", by.y="SIC", all.x=T)
############################################################################
# Control variables
############################################################################
####################
#Hypertensive mediation
#combine both definitions
wml = wml %>% mutate(BPmed_bl_n = case_when(
  hyp_treat_bl == TRUE ~ 1,
  TRUE ~ BPmed_bl
))

wml = wml %>% mutate(BPmed_fu_n = case_when(
  hyp_treat_fu == TRUE ~ 1,
  TRUE ~ BPmed_fu
))
####################

####################
#CESD
####################
cesd_bl=read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult Basis/PV0609_T00013_NODUP.xlsx")
CESD_cols <- sprintf("CES_D_%s",seq(1:20))
cesd_bl$CESD_sum_bl <- rowSums(cesd_bl[,CESD_cols])

cesd_bl_wo_duplicates <-
  cesd_bl %>%
  arrange(CES_D_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

#Composite data
#cesd_bl = cesd_bl %>%
#  mutate(cesd_bl =CES_D_SCORE_SUM_CES_D)

cesd_fu=read_xlsx("/data/gh_gr_agingandobesity_share/life_shared/Data/Raw/Questionnaires_raw/2022_PV609_Beyer/Adult FU1/PV0609_T00013_NODUP.xlsx")
CESD_cols <- sprintf("CES_D_%s",seq(1:20))
cesd_fu$CESD_sum_fu <- rowSums(cesd_fu[,CESD_cols])

cesd_fu_wo_duplicates <-
  cesd_fu %>%
  arrange(CES_D_DATUM) %>%
  filter(!duplicated(.[["SIC"]]))

wml=merge(wml, cesd_bl_wo_duplicates[,c("SIC", "CESD_sum_bl")], by.x="pv_pseudonym", by.y="SIC",all.x=T)
wml=merge(wml, cesd_fu_wo_duplicates[,c("SIC", "CESD_sum_fu")], by.x="pv_pseudonym", by.y="SIC",all.x=T)

#temporary fix
wml$CESD_sum_fu=wml$CESD_sum_bl


############################################################################
#TIV
############################################################################
#load aseg data
aseg <- read.table("/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/FreeSurfer/CortParc_SegResults_baseline/FS_results_subcor_LIFE.txt", header = T, sep = "\t")
#aseg$fu <- ifelse(grepl("_fu", aseg$Measure.volume), 1, 0)
aseg$mrt_pseudonym <- substr(aseg$Measure.volume, 1, 10)
aseg$TIV <- scale(aseg$EstimatedTotalIntraCranialVol)[,1]

aseg_wo_duplicates <-
  aseg %>%
  filter(!duplicated(.[["mrt_pseudonym"]]))

wml=merge(wml, aseg_wo_duplicates[,c("mrt_pseudonym", "TIV")], all.x=T)


############################################################################
# Reshaping the data into long format.
############################################################################
varying=list(c(6,7), #wml
             c(11, 12), #mrt datum
             c(13,19), #epilepsy
             c(14, 24), #stroke
             c(15, 20), #MS
             c(16, 21), #parkinson
             c(26,27), #MMSE
             c(28, 31), #radio assessment usable
             c(29,32), #exclude tumour
             c(30, 33), #exclude lesion
             c(34, 36), #central medication
             c(38, 39), #SBP
             c(40, 41), # WHR
             c(43, 53), #learning
             c(44, 54), #learning_na
             c(45, 55), #recall
             c(46, 56), #recall_na
             c(47, 57), #recognition
             c(48, 58), #phonemic fluency
             c(49, 59), #semantic fluency
             c(50, 60), #TMT
             c(51, 61), #processing
             c(62,63), #BPmed
             c(64, 65) #CESD
                         )
varying_names=c("wml", "mrt_datum",
                "epilepsy", "stroke", "MS", "parkison",
                "MMSE", "radio_usable",
                "radio_lesion", "radio_tumour","med_central",
                "SBP", "WHR", 
                "learning", "learning_na", "recall", "recall_na",
                "recognition", "phon_f", "sem_f", "TMT", "proc",
                "BPmed", "cesd")
wml_long=reshape(wml, varying=varying, 
                  v.names =varying_names, 
                  times=(c("bl","fu")),idvar='mrt_pseudonym', direction='long')

# Calculate age in years
wml_long = wml_long %>% 
  mutate(age = as.integer(mrt_datum - birth)/365)

############################################################################
# Apply exclusion criteria
############################################################################
#We will exclude participants with neurological or psychiatric disease at baseline 
#or follow-up (i.e. radiological finding of ischemic, traumatic or hemorrhagic lesion
#in MRI, incidental finding leading to non-usability of participant, multiple sclerosis, 
#Parkinson’s disease, epilepsy, previous stroke, self-reported dementia, intake of 
#centrally active medication or a score of < 24 in the Mini mental state examination, 
#see Supplementary Table 1). If participants lack information on these variables for 
#one or both timepoints, we will not exclude the participant.
#Further, participants for whom the Lesion Segmentation Toolbox did not run 
#correctly or who were labeled 1 or 3 during quality control will be excluded 
#from all analyses (H1 – H3). 

exclude_ids = wml_long %>%
  group_by(mrt_pseudonym) %>%
  filter(epilepsy == 1| stroke ==1 | MS == 1 | parkison == 1 |
           radio_usable == 1 | radio_lesion == 1 | radio_tumour == 1 |
           med_central == 1 | MMSE == 1 | qa_check ==1 | qa_check == 3)
wml_long <- filter(wml_long, !(mrt_pseudonym %in% exclude_ids$mrt_pseudonym))

# Exclude participants with missing/implausible values of SBP/WHR at baseline AND
# followup
#wml_long = wml_long %>% 
#  group_by(mrt_pseudonym) %>%
#  mutate(sum_na_SBP = sum(is.na(SBP)),
#         sum_na_WHR = sum(is.na(WHR))) %>%
#  filter(sum_na_SBP < 2 & sum_na_WHR < 2)

# Exclude participants without cognitive data (e.g. na values for learning,
# recall, fluencies and TMT)
#wml_long = wml_long %>% 
#  group_by(mrt_pseudonym) %>%
#  mutate(sum_na_cognition = sum(is.na(learning)+is.na(recall)+is.na(phon_f)+is.na(sem_f)+is.na(TMT)+is.na(proc))) %>%
#  filter(sum_na_cognition < 12)

# Exclude participants out of the specified age range (45 -85y at baseline)
ids = wml_long %>% 
  filter(time == "bl" & !is.na(age) & age>45 & age < 85)

wml_long <- filter(wml_long, (mrt_pseudonym %in% ids$mrt_pseudonym))


############################################################################
# Standardize variables
############################################################################
wml_long <- wml_long %>%
  mutate_at(.vars = c("learning", "recall", "recognition", "sem_f", "phon_f", "TMT", 
                      "proc"), ~(scale(.) %>% as.vector)) 

# calculate composite scores for cognitive functions depending on number of available tests
wml_long$exfunct <- apply(wml_long[, c("phon_f", "sem_f", "TMT")], 1, mean, na.rm=TRUE)
wml_long$memo <- apply(wml_long[, c("learning", "recognition", "recall")], 1, mean, na.rm=TRUE)
wml_long$globalcog <- apply(wml_long[, c("exfunct", "memo", "proc")], 1, mean, na.rm=TRUE)

# scale outcome variables to improve interpretability
wml_long <- wml_long %>%
  mutate_at(.vars = c("exfunct", "memo", "globalcog"), ~(scale(.) %>% as.vector))

############################################################################
# Missing values
############################################################################
#Prepare data frame for imputation
test = wml_long %>% 
           mutate(subj=as.integer(as.factor(mrt_pseudonym)),
                  time=as.factor(time),
                  sic=as.factor(sic),
                  BPmed = as.factor(BPmed),
                  education=as.factor(education))%>%
          select(c(1,6,8,18,19,20,21,31,32,33,35,37,38,39,40,41,42,43,44,45,46,47,48))%>%
          relocate(mrt_pseudonym,wml, mrt_datum,.after = globalcog)%>%
          relocate(subj, time, age, .before = sex)%>%
          relocate(SBP, WHR, cesd, BPmed, .after = TIV)

#inspect pattern of missingness
mising_pattern=mice::md.pattern(test)

#right now, there are missing values in Age which need to be removed
#(missing MRI data)              
test=test[!is.na(test$age),]

#Impute all NAs in this dataframe
# set default predictor matrix
ini <- mice::mice(test, maxit=0)
ini$loggedEvents #logged events point out issues we know

# create imputation model with predictor matrix (following 'type')
# indicators:
#-2: cluster variable, 1: fixef with randint, 0: not used
# 2: random effect of predictor which we do not have with two timepoints

#Need to add WHR here once it is varying ;)

colnames(test)

pred <- ini$predictorMatrix
pred["education",] <- c(-2, 0, 1, 1, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0)
pred["TIV",] <- c(-2, 0, 1, 1, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
pred["SBP",] <- c(-2, 1, 1, 1, 1, 0, 0, 0 , 1, 1,0,0,0,0,0,0,0,0,0,0,0,0,0)  
#pred["WHR",] <- c(-2, 1, 1, 1, 1, 0, 1, 1, 0) 
pred["cesd",] <- c(-2, 1, 1, 1, 1, 0, 1, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
pred["BPmed",] <- c(-2, 1, 1, 1, 0, 0, 1, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0)  


# set imputation method 

impmethod <- c("", "", "", "", 
               "2l.bin", "2lonly.pmm", #second level imputation  of education and TIV only
               "2l.pan", "", "2l.pan", "2l.bin", # first level imputation of SBP, WHR, cesd, BPmed
               "", "", "", "", "", "", "", "", "", "", "", "", "")
names(impmethod) <- colnames(test)


m = 5 
iter = 10
imp <- mice::mice(
  test,
  method = impmethod,
  pred = pred,
  maxit = iter,
  m = m,
  seed = 8745
)

#See results
plot(imp)

############################################################################
# calculate between and within measures of Age, SBP, WHR, asinh(WML)
############################################################################
imp_l <- mice::complete(imp, "long", include = TRUE)

imp_l <- imp_l %>%
  group_by(mrt_pseudonym, .imp) %>%
  mutate(
    subj=as.factor(subj),
    cesd = asinh(cesd),
    age_base = ifelse(time == "bl", age, lag(age)),
    SBP_base = ifelse(time == "bl", SBP, lag(SBP)),
    #WHR_base = ifelse(time == "bl", WHR, lag(WHR)),
    age_change = ifelse(time == "bl", 0, age - age_base),
    SBP_change = ifelse(time == "bl", 0, SBP - SBP_base),
    #WHR_change = ifelse(time == "bl", WHR - WHR_base),
    asinh_wml = asinh(wml),
    asinh_wml_base = ifelse(time == "bl", asinh_wml, lag(asinh_wml)),
    asinh_wml_change = ifelse(time == "bl", 0, asinh_wml - asinh_wml_base))

imp.itt <- mice::as.mids(imp_l)


#write full dataset
setwd("/data/pt_life_whm/Results/Tables/")
miceadds::write.mice.imputation(imp.itt , "imputed_data" )





