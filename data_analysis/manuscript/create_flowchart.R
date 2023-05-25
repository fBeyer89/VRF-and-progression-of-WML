library(DiagrammeR)

all_wMRI_fu=nrow(LST_overview[!is.na(LST_overview$LPA_long_complete),])
all_complete=length(LST_overview[!is.na(LST_overview$LPA_long_complete)&LST_overview$LPA_long_complete==1,"LPA_long_complete"])
missing_baseline_MRI=8
missing_fu=2
procissue=1
grViz("digraph flowchart {
node [shape = rectangle, width = 8, fillcolor = Antiquewhite, style = filled, fontsize = 30]
a [label = '@@1', width = 8, height = 1.25]
b [label = '@@2']
c [label = '@@3']
d [label = '@@4']

''

a -> b -> c -> d 


}
[1]: paste(paste0('All participants with MRI at follow-up (n = ', all_wMRI_fu, ')'), sep = '\\n')
[2]: paste(paste0('Participants with WML segmentation at baseline and followup (n = ', all_complete, ')'), sep = '\\n')
[3]: paste(paste0('Participants without neurological disease (n = ', after_exclusion, ')'), sep = '\\n')
[4]: paste(paste0('Diagnosed with dementia or MMSE < 24 (n = ', all_wMRI_fu, ')'), sep = '\\n')
[5]: paste(paste0('With neurological disease (MS, epilepsy, Parkinson's disease,stroke) (n = ', excl_path, ')'), sep = '\\n')
[6]: paste(paste0('With dementia (n = ', excl_dement, ')'), sep = '\\n')
[7]: paste(paste0('With radiological findings (n = ', excl_radio, ')'), sep = '\\n')
[8]: paste(paste0('With radiological findings (n = ', excl_radio, ')'), sep = '\\n')
")

#paste(paste0('missing baseline MRI = ', missing_baseline_MRI, ', missing followup MRI = ', missing_fu, ', processing issue = ', procissue, ''), sep = '\\n')

h [label = '@@3']
i [label = '@@4']
j [label = '@@5']
k [label = '@@6']
l [label = '@@7']
m [label = '@@8']
n [label = '@@9']
o [label = '@@10']
p [label = '@@11']
q [label = '@@12']
r [label = '@@13']
s [label = '@@14']
t [label = '@@15']


[4]: paste(paste0('usable MRI (n = ', length(which(df$FS_usable == 1)) , ')'), paste0('BL = ', length(which(df$FS_usable == 1 & df$ fu == 0)), '   |  ', 'FU = ', length(which(df$FS_usable == 1 & df$fu == 1))), sep = '\\n')
[5]: paste(paste0('no HCV outlier (n = ', length(which(df$FS_usable == 1 & df$outlier_HCV != 1)), ')'), paste0('BL = ', length(which(df$FS_usable == 1 & df$outlier_HCV != 1 & df$fu == 0)), '  |  ', 'FU = ', length(which(df$FS_usable == 1 & df$outlier_HCV != 1 & df$fu == 1))), sep = '\\n')
[6]: paste(paste0('all control variables for model 2 (n = ', length(which(df$FS_usable == 1 & df$outlier_HCV != 1 & !is.na(df$CESD_sum))), ')'), paste0('BL = ', length(which(df$FS_usable == 1 & df$outlier_HCV != 1 & !is.na(df$CESD_sum) & df$fu == 0)), '  |  ', 'FU = ', length(which(df$FS_usable == 1 & df$outlier_HCV != 1 & !is.na(df$CESD_sum) & df$fu == 1))), sep = '\\n')
[7]: paste(paste0('executive functions available (n = ', length(which(!is.na(df$exfunct))), ')'), paste0('BL = ', length(which(!is.na(df$exfunct) & df$ fu == 0)), '  |  ', 'FU = ', length(which(!is.na(df$exfunct) & df$fu == 1))), sep = '\\n')
[8]: paste(paste0('no executive functions outlier (n = ', length(which(df$outlier_exfunct != 1)), ')'), paste0('BL = ', length(which(df$outlier_exfunct != 1 & df$fu == 0)), '  |  ', 'FU = ', length(which(df$outlier_exfunct != 1 & df$fu == 1))), sep = '\\n')
[9]: paste(paste0('all control variables for model 2 (n = ', length(which(df$outlier_exfunct != 1 & !is.na(df$CESD_sum))), ')'), paste0('BL = ', length(which(df$outlier_exfunct != 1 & !is.na(df$CESD_sum) & df$ fu == 0)), '  |  ', 'FU = ', length(which(df$outlier_exfunct != 1 & !is.na(df$CESD_sum) & df$fu == 1))), sep = '\\n')
[10]: paste(paste0('memory available (n = ', length(which(!is.na(df$memo))), ')'), paste0('BL = ', length(which(!is.na(df$memo) & df$fu == 0)), '  |  ', 'FU = ', length(which(!is.na(df$memo) & df$fu == 1))), sep = '\\n')
[11]: paste(paste0('no memory outlier (n = ', length(which(df$outlier_memo != 1)), ')'), paste0('BL = ', length(which(df$outlier_memo != 1 & df$ fu == 0)), '  |  ', 'FU = ', length(which(df$outlier_memo != 1 & df$ fu == 1))), sep = '\\n')
[12]: paste(paste0('all control variables for model 2 (n = ', length(which(df$outlier_memo != 1 & !is.na(df$CESD_sum))), ')'), paste0('BL = ', length(which(df$outlier_memo != 1 & !is.na(df$CESD_sum) & df$fu == 0)), '  |  ', 'FU = ', length(which(df$outlier_memo != 1 & !is.na(df$CESD_sum) & df$fu == 1))), sep = '\\n')
[13]: paste(paste0('processing speed available (n = ', length(which(!is.na(df$procspeed))), ')'), paste0('BL = ', length(which(!is.na(df$procspeed) & df$fu == 0)), '  |  ', 'FU = ', length(which(!is.na(df$procspeed) & df$fu == 1))), sep = '\\n')
[14]: paste(paste0('no processing speed outlier (n = ', length(which(df$outlier_procspeed != 1)), ')'), paste0('BL = ', length(which(df$outlier_procspeed != 1 & df$ fu == 0)), '  |  ', 'FU = ', length(which(df$outlier_procspeed != 1 & df$ fu == 1))), sep = '\\n')
[15]: paste(paste0('all control variables for model 2 (n = ', length(which(df$outlier_procspeed != 1 & !is.na(df$CESD_sum))), ')'), paste0('BL = ', length(which(df$outlier_procspeed != 1 & !is.na(df$CESD_sum) & df$fu == 0)), '  |  ', 'FU = ', length(which(df$outlier_procspeed != 1 & !is.na(df$CESD_sum) & df$fu == 1))), sep = '\\n')
