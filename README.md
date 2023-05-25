# VRF-and-progression-of-WML
This repository contains analysis and simulation scripts for the pre-registered project on the association of vascular risk factors and progression of white matter lesions.

## prepare data
- LST used all earlier scans (described as "bl" in bids format) for all participants (see [/data/gh_gr_agingandobesity_share/life_shared/Analysis/MRI/lifebids/README.md]

## data_analysis
- `calc_models_on_realdata.R`: function to run scripts on actual data
- `run_conf_LME.R`: run confirmatory analyses for H1 - H3
- `run_exp_LME.R: run exploratory analyses for E1a-E3b
- `test_LME_assumptions.R`: test assumptions of models except for VIF
    


## power_simulation
- `test_all_functions_with_simulated_data.R`: function which loads simulated data and can be used to test all models implemented in data_analysis folder
