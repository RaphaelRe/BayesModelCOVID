A Bayesian hierarchical approach to account for reporting uncertainty, variants of concern and vaccination coverage when estimating the effects of non-pharmaceutical interventions on the spread of infectious disease
-------------------------------------------------------------------------------

**Online material**

This repo contains 3 directories:
* data 
* Application
* Simulation study


## Data
Folder contains the used data. This includes all time-shifting distributions, $\pi^{D}$, prevalence of variants of concern and data for application.

## Application
Files used in the application:
* `fit_model_europe.py`: Main file to run the model on 20 European countries
* `prepare_data.R`: Script to build the main data frame for the application. Downloads data, defines NPIs and brings everything together
* `estimate_new_mutants.R`: File which was used to estimate the prevalence of the variants of concern
* `calculate_weighted_ifr.R`: File which was used to calculate $\pi^{D}$

## Simulation Study
Contains all files for the simulation study:
* `fit_model_simulation_study.py`: Main fail to fit the model to the simulation study
* `simulate_data_dynamic.R`: Script to simulate the data for the simulation study
