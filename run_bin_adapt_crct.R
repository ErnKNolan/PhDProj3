#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 ADAPTIVE SIMULATION
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES

pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan,tidyr)

source("make_clusters.R")
source("testFull.R")
source("testInterim.R")
source("makeDecision.R")
source("runSimTrial.R")
source("assignCluster.R")

#for late interim additional analyses
source("runSimTrialLate.R")
source("assignClusterLate.R")
#PROJECT 3
#The different trial properties
set.seed(656038738)
#sensitivity seed
#set.seed(294761845)
#NEW SCOPE
#properties <- expand.grid(nblock=c(2,3),trt_eff_scen = c(3,6), ctrl_prop = c(0.1), icc = c(0.2,0.05), k = c(15),n_per_k = c(10,25,50))
#new scope 25clust
#properties <- expand.grid(nblock=c(2,3),trt_eff_scen = c(3,6), ctrl_prop = c(0.1), icc = c(0.2,0.05), k = c(25),n_per_k = c(10,25,50))
#new scope late start
#properties <- expand.grid(nblock=c(3),trt_eff_scen = c(3,6), ctrl_prop = c(0.1), icc = c(0.2,0.05), k = c(15,25),n_per_k = c(10,25,50))
#extra scenario
properties <- expand.grid(nblock=c(2,3),trt_eff_scen = c(5), ctrl_prop = c(0.1), icc = c(0.2,0.05), k = c(15,25),n_per_k = c(10,25,50))


#bind to properties
properties <- rbind(properties) %>%
  mutate(t4 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0.4,
                        trt_eff_scen == 5 ~ ctrl_prop+0.1,
                        trt_eff_scen == 6 ~ ctrl_prop+0.3),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0.1,
                        trt_eff_scen == 5 ~ ctrl_prop+0.05,
                        trt_eff_scen == 6 ~ ctrl_prop+0.2),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0,
                        trt_eff_scen == 5 ~ ctrl_prop+0,
                        trt_eff_scen == 6 ~ ctrl_prop+0.1),
         t1 = ctrl_prop,
         interim = floor(k/(nblock)),
         draws = ifelse(n_per_k == 75,1250,750))
drawsdat <- properties$draws
properties <- properties %>% dplyr::select(-draws)
properties2 <- properties %>% mutate(row = row_number()) 

#Put in the paths and options for the trial
plan(multisession,workers=32)
baepath <- "adapt_arm2.stan"
#baepath <- "prior_sens.stan"
set_cmdstan_path(path="/root/.cmdstan/cmdstan-2.33.1")
outdir <- "SimTrash"
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
adaption <- "both" #this can be early_stopping, arm_dropping, or both
drop_cut <- 0.05
stop_cut <- 0.15
t <- 4

#for sensitivity 1, 8, 18, 39, 43
#Run the trial
test <- list()
for(j in 1:24){
  test[[length(test)+1]] <- future_replicate(2500,future.seed=42L,runSimTrial(properties,mod,outdir,j,adaption,drop_cut,stop_cut,t=t,nblock=properties$nblock,draws=drawsdat[j]))
                            #future_replicate(2500,future.seed=42L,runSimTrialLate(properties,mod,outdir,j,adaption,drop_cut,stop_cut,t=t,nblock=properties$nblock,draws=drawsdat[j]))
  saveRDS(test,"adapt_newscen.RDS")
}

