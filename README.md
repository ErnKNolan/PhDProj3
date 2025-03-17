
run_bin_crct.R and run_bin_adapt_crct.R are the two main files that run and collate the simulation results.
run_bin_crct.R is for the non-adaptive trial, and run_bin_adapt_crct.R is for the adaptive trial.

The following are files that make functions for the trial:
assignCluster.R
makeDecision.R
nonAdaptTrial.R
runSimTrial.R
testInterim.R
testFull.R

adapt_sim_convergence.R checks how the convergence of the models in the simulation.
adapt_sim_performance.R calculates the operating characteristics of interest in the trial
and graphs them over the trial properties.

correct_differences.R is used in adapt_sim_performance.R to compare the operating characteristics
between trial properties.

adapt_arm.stan is the Bayesian hierarchical model used in all simulations.
