#!/bin/sh

source virtualenv/bin/activate

# fit models on covid data
echo "Running models on Covid data..."
cd application
python fit_model_europe.py > model_europe.out
python sensitivity_analysis.py > sensitivity.out
wait
# plot results
echo "Plot results..."
mkdir -p plots/main_results
mkdir -p plots/sensitivity_analysis
R < plot_results.R --no-save > plots_europe.out 2>&1
R < plot_sensitivity.R --no-save > plots_sensitivity.out 2>&1
wait



cd ../simulation_study/
# simulated data
echo "Simulate datasets..."
R < simulate_data_dynamic.R --no-save > simulation_standard.out 2>&1
R < simulate_data_dynamic_diffusion.R --no-save > simulation_diffusion.out 2>&1

wait

echo "Running models on simulated data..."
python fit_model_simulation_study.py > simulation.out
python fit_model_simulation_study_diffusion.py > diffusion.out
wait

# plot results
echo "Plot results..."
mkdir -p plots/standard 
mkdir -p plots/diffusion

R < plot_results_standard.R --no-save > plots_standard.out 2>&1
R < plot_results_diffusion.R --no-save > plots_diffusion.out 2>&1
wait

echo "==== Finished ===="
