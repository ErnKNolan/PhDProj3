#AUTHOR: Erin Nolan
#TITLE: NON-ADAPTIVE SIMULATIONS FOR OPTIMISING IMPLEMENTATION STRATEGIES
#PURPOSE: RUNS THE SIMULATIONS FOR MULTIARM cRCT UNDER VARIOUS PROPERTIES

#Run the functions defined in make_clusters and fit_bae
pacman::p_load(here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan)

source("make_clusters.R")
source("runNonAdaptTrial.R")
source("testInterim.R")

#The different trial properties
set.seed(656038738)
#trt_eff_scen = the treatment scenario
#ctrl_prop = the baseline proportion of the event of interest
#icc = intra-class correlation
#n_per_k = number of participants per cluster
#k = number of clusters
properties <- expand.grid(trt_eff_scen = c(2,3), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(25,50,75), k = c(15,25),nblock=1)

#bind to properties
properties <- rbind(properties) %>%
  mutate(t4 = case_when(trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t3 = case_when(trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t2 = case_when(trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0),
         t1 = ctrl_prop,
         draws = ifelse(n_per_k == 75,1250,750))

drawsdat <- properties$draws
properties <- properties %>% dplyr::select(-draws)


#set how many workers you want to use
plan(multisession,workers=20) 

#Put in the paths for the simulation
baepath <- "adapt_arm2.stan" #the file path for the Stan model code
set_cmdstan_path(path="/root/.cmdstan/cmdstan-2.33.1") #where cmdstan is located
outdir <- "SimTrash" #where the stan files will be output
mod <- cmdstan_model(baepath, pedantic = F, compile=T)

#run the trial
#test <- list() #test is the list of all the output from the simulated trials
for(j in 1:24){ #loops through properties 1 to 24
  test[[length(test)+1]] <- future_replicate(5,runNonAdaptTrial(properties,mod,outdir,j,t=t,draws=drawsdat[j]),
                                             future.seed = 42L)
  saveRDS(test,"nonadapt_test.RDS")
}

nonadapt <- list()
for(j in 1:24){
  for(i in seq(1,12500,5)){
    nonadapt[[length(nonadapt)+1]] <- test[[j]][[i]]
    nonadapt[[length(nonadapt)]]$sim <- (i+4)/5
    nonadapt[[length(nonadapt)]]$property <- j
  }
}
nonadapt <- bind_rows(nonadapt)
properties2 <- properties %>% mutate(row = row_number()) 
nonadapt_out <- merge(nonadapt,properties2,by.y=c("row"),by.x="property")
saveRDS(nonadapt_out,here("Data","nonadapt_out.RDS"))
