# Simulation scripts

This folder contains necessary scripts for the power analysis using simulated data

1. Simulation & Model (M1) estimation in simulated data
Both steps are performed in `calc_models_on_simdata.R`.  

We simulated the data based on assumptions regarding the direction and size of the association of VR factors and WMH progression.  
This was based on published effect sizes.  
We used three sample sizes (N= 400, 600, 800) and three WHR scalings (0.5, 1, 1.5), which sample the uncertainty with regard to WHR effect size.  
Simulation happens in `/simulate_data.R` called from `calc_models_on_simdata.R`.
Simulated data are saved in `simulated_data.csv`.

Then, `run_LME_simulation.R` fits the model M1 to the simulated data and and results are saved for frequentist and Bayesian inference. 
`evaluate_power.R` calculates the power from the results saved in the previous step.

2. Testing of functions with simulated data
`test_all_functions_with_simulated_data.R` uses the simulated data to check that all functions used for the main analysis work. 


