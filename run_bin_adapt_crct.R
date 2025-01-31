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

#PROJECT 3
#The different trial properties
set.seed(656038738)
#The first interim is after either 3 or 5 clusters (out of 5 and 10)
properties <- expand.grid(trt_eff_scen = c(2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(25,50,75), k = c(15,25),nblock=c(2,3)) 
properties2 <- properties %>% mutate(row = row_number()) 
#bind to properties
properties <- rbind(properties) %>%
  mutate(t4 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t1 = ctrl_prop,
         interim = floor(k/(nblock)),
         draws = ifelse(n_per_k == 75,1250,750))

drawsdat <- properties$draws
properties <- properties %>% dplyr::select(-draws)
#Put in the paths and options for the trial
plan(multisession,workers=20)
baepath <- "adapt_arm2.stan"
set_cmdstan_path(path="/root/.cmdstan/cmdstan-2.33.1")
outdir <- "SimTrash"
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
adaption <- "both" #this can be early_stopping, arm_dropping, or both
drop_cut <- 0.05
stop_cut <- 0.15
t <- 4

#Run the trial
test <- list()
for(j in c(1:48)){
  test[[length(test)+1]] <- future_replicate(2500,future.seed=42L,runSimTrial(properties,mod,outdir,j,adaption,drop_cut,stop_cut,t=t,nblock=properties$nblock,draws=drawsdat[j]))
  saveRDS(test,"adaptprob2.RDS")
}
#saveRDS(test,here("Data","pilot_sim.RDS"))
#Take out the trial properties
trial_props <- list()
for(j in 1:48){
  for(i in seq(3,12500,5)){
    trial_props[[length(trial_props)+1]] <- test[[j]][[i]]
    trial_props[[length(trial_props)]]$sim <- (i+2)/5
    trial_props[[length(trial_props)]]$property <- j
  }
}
trial_props <- bind_rows(trial_props)
#saveRDS(trial_props,here("Data","adaptprob_trial_props.RDS"))

#Take out the interim analyses
interim <- list()
for(j in 1:48){
  for(i in seq(1,12500,5)){
    interim[[length(interim)+1]] <- test[[j]][[i]]
    interim[[length(interim)]]$sim <- (i+4)/5
    interim[[length(interim)]]$property <- j
  }
}
interim <- bind_rows(interim)
#saveRDS(interim,here("Data","adapt_interim.RDS"))

#Take out the full analyses
tempd <- list()
for(j in c(1:23)){
  for(i in seq(2,75,5)){
    tempd[[length(tempd)+1]] <- test[[j]][[i]]
    tempd[[length(tempd)]]$sim <- (i+3)/5
    tempd[[length(tempd)]]$property <- j
  }
}
outsim <- bind_rows(tempd)

#pull out the cluster data
clusts <- list()
for(j in c(1:11)){
  for(i in seq(5,12500,5)){
    clusts[[length(clusts)+1]] <- test[[j]][[i]]
  }
}

#clean the cluster data
clusters <- plyr::ldply(clusts, rbind) %>%
  mutate(arm = rep(c("arm2","arm3","arm4"),times=length(test)*2500),`2`=NA,
         property = rep(c(1:length(test)),each=3*2500), #specify the property
         sim = rep(c(1:2500),each=3,times=length(test))) %>% #specify the sim
  rename(interim1 = `1`,
         interim2 = `2`) %>%
  pivot_wider(id_cols=c(property,sim),names_from = c(arm), values_from=c(interim1,interim2)) %>%
  merge(properties2,by.x="property",by.y="row") %>%
  group_by(property,sim) %>%
  mutate(arm2 = sum(c(interim1_arm2,interim2_arm2,interim),na.rm=TRUE),
         arm3 = sum(c(interim1_arm3,interim2_arm3,interim),na.rm=TRUE),
         arm4 = sum(c(interim1_arm4,interim2_arm4,interim),na.rm=TRUE),
         stop_int1 = ifelse(interim1_arm2 == 0 & interim1_arm3 == 0 & interim1_arm4 == 0,1,0),
         stop_int2 = ifelse(interim2_arm2 == 0 & interim2_arm3 == 0 & interim2_arm4 == 0,1,0))

outsim2 <- merge(outsim,clusters,by.x=c("property","sim"),by.y = c("property","sim"))
saveRDS(clusters,here("Data","clusters.RDS"))
saveRDS(outsim2,here("Data","prob_outsim.RDS"))

#sensitivity analysis 1st interim of 2 interim trials
#Take out the interim analyses
interim_s <- list()
for(j in 1:24){
  for(i in seq(1,12500,5)){
    interim_s[[length(interim_s)+1]] <- test[[j]][[i]]
    interim_s[[length(interim_s)]]$sim <- (i+4)/5
    interim_s[[length(interim_s)]]$property <- j
  }
}
interim_s <- bind_rows(interim_s)
#saveRDS(interim_s,here("Data","adapt_interim_s.RDS"))

#nonadaptive
test <- readRDS(here("Data","nonadapt.RDS"))
nonadapt <- list()
for(j in 1:24){
  for(i in 1:2500){
    nonadapt[[length(nonadapt)+1]] <- test[[j]][[i]]
    nonadapt[[length(nonadapt)]]$sim <- i
    nonadapt[[length(nonadapt)]]$property <- j
  }
}
nonadapt <- bind_rows(nonadapt)
properties2 <- properties %>% mutate(row = row_number()) 
nonadapt_out <- merge(nonadapt,properties2,by.y=c("row"),by.x="property")
saveRDS(nonadapt_out,here("Data","nonadapt_out.RDS"))
