
run_bin_crct.R and run_bin_adapt_crct.R are the two main files that run the simulations
run_bin_crct.R is for the non-adaptive trial, and run_bin_adapt_crct.R is for the adaptive trial.
save_adapt_sims.R and save_nonadapt_sims.R collate the results.

The following are files that make functions for the trial:
assignCluster.R
makeDecision.R
nonAdaptTrial.R
runSimTrial.R
testInterim.R
testFull.R

runSimTrialLate.R is the 'runSimTrial.R' function for the simulated trial that has 2 interims and has the first interim 50% through the trial.
assignClusterLate.R is the 'assignCluster.R' function for the simulated trial that has 2 interims and has the first interim 50% through the trial.

adapt_sim_convergence.R checks how the convergence of the models in the simulation.
adapt_sim_performance.R calculates the operating characteristics of interest in the trial and graphs them over the trial properties.
adapt_sim_performance_newscen.R is for the sensitivity analysis that has a scenario with a less clear optimal arm.

adapt_arm.stan is the Bayesian hierarchical model used in all simulations.

single_example.R and single_example_null.R are the scripts to run the single adaptive trial examples in additional file 3.
