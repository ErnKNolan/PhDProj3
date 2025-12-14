#AUTHOR: Erin Nolan
#TITLE: CASE STUDY
#PURPOSE: RUN THE CASE STUDY FOR PROJECT 3

pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan,tidyr)
set.seed(656038738)
plan(multisession,workers=32)
#trt1 = ctrl = 0.172
#trt2 = HD = 0.619
#trt2_2 = HD2 = 0.649
#trt4 = LD = 0.692
#for simplicity, assume all same amount as ctrl (rounded) - 31 clusters, 6 teachers per cluster at time 2
#icc = 0.07

#--------------------------------------------------------------------------------------------------------
#DESIGN 1: Consequtive trials
source("make_clusters.R")
source("testModel.R")

#SETUP
#assign properties to trials
properties <- expand.grid(trt_eff_scen=c("null","alt"),icc = c(0.07), k = c(31),n_per_k = c(6)) %>%
  mutate(t2 = ifelse(trt_eff_scen=="null",0.172,0.619),
         t1 = 0.172,
         draws = 1250)
properties2 <- properties %>% mutate(row = row_number()) 
baepath <- "conseq_trial.stan"
set_cmdstan_path(path="/root/.cmdstan/cmdstan-2.33.1")
mod <- cmdstan_model(baepath, pedantic = F, compile=T)

#RUN THE TRIALS
expdat <- makeClusters(t=2,nid=6,t1=31,t2=31)
conseqTrialAlt <- function(){
  #trial 1
  trial1 <- testModel(expdat=expdat,t=2,mod=mod,outdir=outdir,rho=0.07,
                      t1=0.172,t2=0.619,draws=1250)
  trial <- 1
  #trial 2, if trial 1 successful
  if(trial1[[1]][3,2] > 0.85){
    trial2 <- testModel(expdat=expdat,t=2,mod=mod,outdir=outdir,rho=0.07,
                        t1=0.619,t2=0.649,draws=1250)
    success <- as.numeric(trial2[[1]][3,2] > 0.85)
    trial <- 2
  }else{
    success <- 0
  }
  #trial 3, if trial 2 successful
  if(success == 1){
    trial3 <- testModel(expdat=expdat,t=2,mod=mod,outdir=outdir,rho=0.07,
                        t1=0.649,t2=0.692,draws=1250)
    success <- as.numeric(trial3[[1]][3,2] > 0.85)
    trial <- 3
  }else{
    success <- 0
  }
  res <- cbind(success,trial)
  return(res)
}

conseqTrialNull <- function(){
  #trial 1
  trial1 <- testModel(expdat=expdat,t=2,mod=mod,outdir=outdir,rho=0.07,
                      t1=0.172,t2=0.172,draws=1250)
  trial <- 1
  #trial 2, if trial 1 successful
  if(trial1[[1]][3,2] > 0.85){
    trial2 <- testModel(expdat=expdat,t=2,mod=mod,outdir=outdir,rho=0.07,
                        t1=0.172,t2=0.172,draws=1250)
    success <- as.numeric(trial2[[1]][3,2] > 0.85 | trial2[[1]][3,2] < 0.15)
    trial <- 2
  }else{
    success <- 0
  }
  #trial 3, if trial 2 successful
  if(success == 1){
    trial3 <- testModel(expdat=expdat,t=2,mod=mod,outdir=outdir,rho=0.07,
                        t1=0.172,t2=0.172,draws=1250)
    success <- as.numeric(trial3[[1]][3,2] > 0.85 | trial3[[1]][3,2] < 0.15)
    trial <- 3
  }else{
    success <- 0
  }
  res <- cbind(success,trial)
  return(res)
}

conseq_success_alt <- future_replicate(2500,future.seed=42L,conseqTrialAlt())
saveRDS(conseq_success_alt,"conseq_success_alt.RDS") #0.1228
#X = 1, X = 2, X = 3
#(((0*1*2)+(1703*2*2)+(707*3*2))*31*6)/2500 = 882
conseq_success_null <- future_replicate(2500,future.seed=42L,conseqTrialNull())
saveRDS(conseq_success_null,"conseq_success_null.RDS") #0.016
#2125 = 1, 253 = 2, 122 = 3
#(((2125*1*2)+(253*2*2)+(122*3*2))*31*6)/2500 = 446
#--------------------------------------------------------------------------------------------------------
#DESIGN 2: Multiarm fixed trial

#SETUP
source("make_clusters.R")
source("testInterim.R")
source("nonAdaptTrial.R")
j <- 1
t <- 4
baepath <- "adapt_arm2.stan"
mod <- cmdstan_model(baepath, pedantic = F, compile=T)

#expdat <- makeClusters(t=4,nid=6,t1=31,t2=31,t3=31,t4=31)
# 
# caseNonAdaptAlt <- function(){
#   full <- nonAdaptTrial(properties,mod=mod,outdir,j=j,t=t,draws=1250)
#   
#   #full <- nonAdaptTrial(expdat=expdat,t=t,mod=mod,outdir=outdir,j=j,
#   #                    rho=0.07,t1=0.172,t2=0.619,t3=0.649,t4=0.692,draws=1250)
#   #res <- as.numeric(full[[1]][2,2] > 0.85)
#   return(full)
# }

# caseNonAdaptNull <- function(){
#   
#   #full <- CaseInterim(expdat=expdat,t=3,mod=mod,outdir=outdir,
#   #                    rho=0.07,t1=0.172,t2=0.172,t3=0.172,t4=0.172,draws=1250)
#   #res <- as.numeric(full[[1]][2,2] > 0.85 | full[[1]][1,2] > 0.85)
#   #return(res)
#   full <- nonAdaptTrial(properties,mod=mod,outdir,j=j,t=t,draws=1250)
#   return(full)
# }
#ALT SCENARIO
properties <- expand.grid(icc = 0.07, k = 31, n_per_k = 6, nblock=1) %>%
  mutate(t4=0.692,t3=0.649,t2=0.619,t1=0.172)
nonadapt_alt <- future_replicate(2500,future.seed=42L,nonAdaptTrial(properties,mod=mod,outdir,j=j,t=t,draws=1250))
estimate <- bind_rows(nonadapt_alt) %>% filter(grepl("pp_trt4",variable))
mean(estimate$mean > 0.85) 
saveRDS(nonadapt_alt,"case_nonadapt_alt.RDS") #0.282

#get the output to a readable state
nonadapt_alt2 <- list(nonadapt_alt)
nonadapt <- list()
for(j in 1){
  for(i in 1:2500){
    nonadapt[[length(nonadapt)+1]] <- nonadapt_alt2[[j]][[i]]
    nonadapt[[length(nonadapt)]]$sim <- i
    nonadapt[[length(nonadapt)]]$property <- j
  }
}
nonadapt <- bind_rows(nonadapt)
nonadapt_out <- merge(nonadapt,properties)
nonadapt_out$trt_eff_scen <- "alt"

#NULL SCENARIO
properties <- expand.grid(icc = 0.07, k = 31, n_per_k = 6, nblock=1) %>%
  mutate(t4=0.172,t3=0.172,t2=0.172,t1=0.172)
nonadapt_null <- future_replicate(2500,future.seed=42L,nonAdaptTrial(properties,mod=mod,outdir,j=j,t=t,draws=1250))
estimate <- bind_rows(nonadapt_null) %>% filter(grepl("pp_trt",variable)) %>%
  mutate(group = (row_number() - 1) %/% 3) %>%
  group_by(group) %>%
  mutate(success = max(ifelse(mean > 0.85,1,0))) %>%
  filter(row_number() == 1)
mean(estimate$success)
saveRDS(nonadapt_null,"case_nonadapt_null.RDS") #0.0564

#get the output to a readable state
nonadapt_null2 <- list(nonadapt_null)
nonadaptnull <- list()
for(j in 1){
  for(i in 1:2500){
    nonadaptnull[[length(nonadaptnull)+1]] <- nonadapt_null2[[j]][[i]]
    nonadaptnull[[length(nonadaptnull)]]$sim <- i
    nonadaptnull[[length(nonadaptnull)]]$property <- j
  }
}
nonadaptnull <- bind_rows(nonadaptnull)
nonadaptnull_out <- merge(nonadaptnull,properties)
nonadaptnull_out$trt_eff_scen <- "null"
nonadaptnull_out$interim <- 0
sim_pow <- rbind(nonadapt_out,nonadaptnull_out)

#--------------------------------------------------------------------------------------------------------
#DESIGN 3: Multiarm adaptive trial - 1 interim
#SETUP
source("make_clusters.R")
source("testFull.R")
source("testInterim.R")
source("makeDecision.R")
source("runSimTrial.R")
#source("CaseSimTrialNull.R")
source("assignCluster.R")
baepath <- "adapt_arm2.stan"
mod <- cmdstan_model(baepath, pedantic = F, compile=T)
drop_cut <- 0.05
stop_cut <- 0.15
t <- 4
j <- 1

#ALT SCENARIO
properties <- expand.grid(nblock=c(3), trt_eff_scen=c("alt"),icc = c(0.07), k = c(31),n_per_k = c(6)) %>%
  mutate(t4 = 0.692,t3=0.649,t2 = 0.619,t1 = 0.172,interim = floor(k/(nblock)))

adapt_alt <- future_replicate(2500,future.seed=42L,runSimTrial(properties,mod,outdir,j=j,adaption="both",drop_cut=drop_cut,stop_cut=stop_cut,t=4,nblock=3,draws=1250))
adapt_alt2 <- list(adapt_alt)

#Take out the full analyses
tempd <- list()
for(j in c(1:1)){
  for(i in seq(2,12500,5)){
    tempd[[length(tempd)+1]] <- adapt_alt2[[j]][[i]]
    tempd[[length(tempd)]]$sim <- (i+3)/5
    tempd[[length(tempd)]]$property <- j
  }
}
outsim <- bind_rows(tempd)

#pull out the cluster data
clusts <- list()
for(j in c(1:1)){
  for(i in seq(5,12500,5)){
    clusts[[length(clusts)+1]] <- adapt_alt2[[j]][[i]]
  }
}

#clean the cluster data
clusters <- plyr::ldply(clusts, rbind) %>%
  mutate(arm = rep(c("arm2","arm3","arm4"),times=length(adapt_alt2)*2500),
         property = rep(c(1:length(test)),each=3*2500), #specify the property
         sim = rep(c(1:2500),each=3,times=length(test))) %>% #specify the sim
  rename(interim1 = `1`,
         interim2 = `2`) %>%
  pivot_wider(id_cols=c(property,sim),names_from = c(arm), values_from=c(interim1,interim2)) %>%
  merge(properties) %>%
  group_by(property,sim) %>%
  mutate(arm2 = sum(c(interim1_arm2,interim2_arm2,interim),na.rm=TRUE),
         arm3 = sum(c(interim1_arm3,interim2_arm3,interim),na.rm=TRUE),
         arm4 = sum(c(interim1_arm4,interim2_arm4,interim),na.rm=TRUE),
         stop_int1 = ifelse(interim1_arm2 == 0 & interim1_arm3 == 0 & interim1_arm4 == 0,1,0),
         stop_int2 = ifelse(interim2_arm2 == 0 & interim2_arm3 == 0 & interim2_arm4 == 0,1,0))

outsim2 <- merge(outsim,clusters,by.x=c("property","sim"),by.y = c("property","sim"))

trial_success <- outsim2 %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(sim,interim1_arm2,interim1_arm3,interim1_arm4,
                        interim2_arm2,interim2_arm3,interim2_arm4,
                        stop_int1,stop_int2),names_from=variable,values_from=mean) %>%
  mutate(power = case_when(!is.na(interim2_arm2) & (pp_trt4 >= 0.85 & interim2_arm4 > 0) ~ 1,
                           .default = 0)) 
adapt_alt_power <- mean(trial_success$power) #0.2412
saveRDS(adapt_alt,"case_adapt_alt.RDS")

#NULL SCENARIO
properties <- expand.grid(nblock=c(3), trt_eff_scen=c("null"),icc = c(0.07), k = c(31),n_per_k = c(6)) %>%
  mutate(t4=0.172,t3 = 0.172,t2 = 0.172,t1 = 0.172,interim = floor(k/(nblock)))
adapt_null <- future_replicate(2500,future.seed=42L,runSimTrial(properties,mod,outdir,j=j,adaption="both",drop_cut=drop_cut,stop_cut=stop_cut,t=4,nblock=3,draws=1250))
adapt_null2 <- list(adapt_null)

#Take out the full analyses
tempd <- list()
for(j in c(1:1)){
  for(i in seq(2,12500,5)){
    tempd[[length(tempd)+1]] <- adapt_null2[[j]][[i]]
    tempd[[length(tempd)]]$sim <- (i+3)/5
    tempd[[length(tempd)]]$property <- j
  }
}
outsim <- bind_rows(tempd)

#pull out the cluster data
clusts <- list()
for(j in c(1:1)){
  for(i in seq(5,12500,5)){
    clusts[[length(clusts)+1]] <- adapt_null2[[j]][[i]]
  }
}

#clean the cluster data
clusters <- plyr::ldply(clusts, rbind) %>%
  mutate(arm = rep(c("arm2","arm3","arm4"),times=length(adapt_null2)*2500),
         property = rep(c(1:length(test)),each=3*2500), #specify the property
         sim = rep(c(1:2500),each=3,times=length(test))) %>% #specify the sim
  rename(interim1 = `1`,
         interim2 = `2`) %>%
  pivot_wider(id_cols=c(property,sim),names_from = c(arm), values_from=c(interim1,interim2)) %>%
  merge(properties) %>%
  group_by(property,sim) %>%
  mutate(arm2 = sum(c(interim1_arm2,interim2_arm2,interim),na.rm=TRUE),
         arm3 = sum(c(interim1_arm3,interim2_arm3,interim),na.rm=TRUE),
         arm4 = sum(c(interim1_arm4,interim2_arm4,interim),na.rm=TRUE),
         stop_int1 = ifelse(interim1_arm2 == 0 & interim1_arm3 == 0 & interim1_arm4 == 0,1,0),
         stop_int2 = ifelse(interim2_arm2 == 0 & interim2_arm3 == 0 & interim2_arm4 == 0,1,0))

outsim2_null <- merge(outsim,clusters,by.x=c("property","sim"),by.y = c("property","sim"))

trial_success <- outsim2 %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(sim,interim1_arm2,interim1_arm3,interim1_arm4,
                        interim2_arm2,interim2_arm3,interim2_arm4,
                        stop_int1,stop_int2),names_from=variable,values_from=mean) %>%
  mutate(type1 = case_when(!is.na(interim1_arm2) & !is.na(interim2_arm2) &
                           ((pp_trt2 >= 0.85 & interim2_arm2 > 0) | (pp_trt3 >= 0.85 & interim2_arm3 > 0) | (pp_trt4 >= 0.85 & interim2_arm4 > 0)) ~ 1, #2 interim trials
                           .default = 0)) 
adapt_null_type1 <- mean(trial_success$type1) #0.0268
saveRDS(adapt_null,"case_adapt_null.RDS")

#getting the scaled power for the fixed and adaptive scenarios
outsim2$trt_eff_scen <- "alt"
outsim2_null$trt_eff_scen <- "null"
sim_pow <- rbind(outsim2,outsim2_null)
#add in properties and nblocks
outsim2_null$property <- 2
nonadapt_out$property <- 3
nonadaptnull_out$property <- 4
nonadapt_out$nblock <- 1
nonadaptnull_out$nblock <- 1
#combine the datasets
sim_pow2 <- bind_rows(outsim2,outsim2_null,nonadapt_out,nonadaptnull_out)
#to get power, type 1 error and scaled power
trial_success <- sim_pow2 %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen,k,icc,n_per_k,nblock,
                        interim1_arm2,interim1_arm3,interim1_arm4,
                        interim2_arm2,interim2_arm3,interim2_arm4,
                        stop_int1,stop_int2),names_from=variable,values_from=mean) %>%
  mutate(type1 = case_when(trt_eff_scen == "null" & is.na(interim1_arm2) & (pp_trt2 >= 0.85 | pp_trt3 >= 0.85 | pp_trt4 >= 0.85) ~ 1, #fixed trials
                           trt_eff_scen == "null" & !is.na(interim1_arm2) & is.na(interim2_arm2) & 
                             ((pp_trt2 >= 0.85 & interim1_arm2 > 0) | (pp_trt3 >= 0.85 & interim1_arm3 > 0) | (pp_trt4 >= 0.85 & interim1_arm4 > 0)) ~ 1, #1 interim trials 
                           trt_eff_scen == "null" & !is.na(interim1_arm2) & !is.na(interim2_arm2) &
                             ((pp_trt2 >= 0.85 & interim2_arm2 > 0) | (pp_trt3 >= 0.85 & interim2_arm3 > 0) | (pp_trt4 >= 0.85 & interim2_arm4 > 0)) ~ 1, #2 interim trials
                           trt_eff_scen =="alt" ~ NA,#effect scenario
                           .default = 0), 
         power = case_when(trt_eff_scen == "null" ~ NA,
                           trt_eff_scen %in% c("alt") & is.na(interim1_arm2) & pp_trt4 >= 0.85 ~ 1, #fixed trials
                           trt_eff_scen %in% c("alt") & !is.na(interim1_arm2) & (pp_trt4 >= 0.85 & interim1_arm4 > 0) ~ 1, #1 interim trials 
                           trt_eff_scen %in% c("alt") & !is.na(interim2_arm2) & (pp_trt4 >= 0.85 & interim2_arm4 > 0) ~ 1, #2 interim trials
                           .default = 0),
         correct = ifelse(type1 %in% c(0) | power %in% c(1),1,0)) 

#select variables for scaling the power
scale_type1 <- trial_success %>% filter(trt_eff_scen =="null") %>% group_by(k,icc,n_per_k,nblock) %>% filter(row_number()==1) %>%
  dplyr::select(k,icc,n_per_k,nblock)
#get the pp_trts we need to scale
prop1 <- trial_success %>% filter(trt_eff_scen =="null") %>% group_by(k,icc,n_per_k,nblock)%>% 
  mutate(pp_trt2 = case_when(interim1_arm2 %in% c(0) ~ NA, #dont want to include cutpoints for dropped arms
                             interim2_arm2 %in% c(0) ~ NA,
                             .default=pp_trt2),
         pp_trt3 = case_when(interim1_arm3 %in% c(0) ~ NA,
                             interim2_arm3 %in% c(0) ~ NA,
                             .default=pp_trt3),
         pp_trt4 = case_when(interim1_arm4 %in% c(0) ~ NA,
                             interim2_arm4 %in% c(0) ~ NA,
                             .default=pp_trt4)) %>%
  group_split()
#loop function to find the threshold related to a type 1 error of 0.05
scale <- list()
for(i in 1:2){ #1:the number of rows in prop1
  pp_trt2 <- prop1[[i]]$pp_trt2
  pp_trt3 <- prop1[[i]]$pp_trt3
  pp_trt4 <- prop1[[i]]$pp_trt4
  x <- matrix(nrow=1001,ncol=2)
  ind <- 1
  for(y in seq(0,1,by=0.001)){
    x[ind,1] <- y
    x[ind,2] <- sum(pp_trt2 >= y | pp_trt3 >= y | pp_trt4 >= y,na.rm=T)/2500#2500 = number of sims here
    ind <- ind+1
  }
  #finding the relevant quantile that has the type 1 error closest to 0.05
  quantile <- data.frame(prop=x[which.min(abs(x[,2] - 0.05)),1],
                         icc=prop1[[i]][1,]$icc,
                         k=prop1[[i]][1,]$k,
                         n_per_k=prop1[[i]][1,]$n_per_k,
                         nblock=prop1[[i]][1,]$nblock)
  scale[[i]] <- quantile
}
#add in all the thresholds
scale_type1 <- bind_rows(scale)
#determine the scaled power from the thresholds
trial_success2 <- merge(trial_success,scale_type1,by=c("k","icc","n_per_k","nblock")) %>%
  mutate(scale_power = case_when(trt_eff_scen == "null" ~ NA,
                                 trt_eff_scen %in% c("alt") & is.na(interim1_arm2) & pp_trt4 >= prop ~ 1, #fixed trials
                                 trt_eff_scen %in% c("alt") & !is.na(interim1_arm2) & (pp_trt4 >= prop & interim1_arm4 > 0) ~ 1, #1 interim trials 
                                 trt_eff_scen %in% c("alt") & !is.na(interim2_arm2) & (pp_trt4 >= prop & interim2_arm4 > 0) ~ 1, #2 interim trials
                                 .default = 0),
         scale_type1 = case_when(trt_eff_scen == "null" & is.na(interim1_arm2) & (pp_trt2 >= prop | pp_trt3 >= prop | pp_trt4 >= prop) ~ 1, #fixed trials
                                 trt_eff_scen == "null" & !is.na(interim1_arm2) & is.na(interim2_arm2) & 
                                   ((pp_trt2 >= prop & interim1_arm2 > 0) | (pp_trt3 >= prop & interim1_arm3 > 0) | (pp_trt4 >= prop & interim1_arm4 > 0)) ~ 1, #1 interim trials 
                                 trt_eff_scen == "null" & !is.na(interim1_arm2) & !is.na(interim2_arm2) &
                                   ((pp_trt2 >= prop & interim2_arm2 > 0) | (pp_trt3 >= prop & interim2_arm3 > 0) | (pp_trt4 >= prop & interim2_arm4 > 0)) ~ 1, #2 interim trials
                                 trt_eff_scen %in% c("alt") ~ NA,#effect scenario
                                 .default = 0))

correct <- trial_success2 %>% 
  group_by(k,n_per_k,icc,nblock) %>%
  summarise(correct = sum(correct)/n(),
            type1 = sum(type1,na.rm=TRUE)/(n()/2),
            power = sum(power,na.rm=TRUE)/(n()/2),
            scale_power = sum(scale_power,na.rm=T)/(n()/2),
            scale_type1 = sum(scale_type1,na.rm=T)/(n()/2)) %>%
  mutate(#k2 = ifelse(k == 15,"1_15","0_25"),
    k3 = factor(k, levels=c("25","15")),
    nint = nblock-1,
    n_per_k = factor(n_per_k)) %>%
  ungroup()
  
