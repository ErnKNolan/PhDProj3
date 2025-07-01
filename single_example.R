#Single adaptive trial to show example
set.seed(656038738)
pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan,tidyr)
source(here("Programs","make_clusters.R"))
source(here("Programs","testFull.R"))
source(here("Programs","testInterim.R"))
source(here("Programs","makeDecision.R"))
source(here("Programs","runSimTrial.R"))
source(here("Programs","assignCluster.R"))

plan(multisession,workers=20)
baepath <- "C:/Users/ENolan/OneDrive - HMRI/Documents/PhDProject2/Programs/PhDProject3/Programs/adapt_arm2.stan"
#set_cmdstan_path(path="C:/Users/nolan/Documents/.cmdstan/cmdstan-2.33.1")
#baepath <- "D:/Programs/PhDProject3/Programs/adapt_arm.stan"
set_cmdstan_path(path="C:/Users/ENolan/OneDrive - HMRI/Documents/.cmdstan/cmdstan-2.33.1")
#outdir <- "J:/Sims"
outdir <- "C:/Users/ENolan/Downloads/Sims"
mod <- cmdstan_model(baepath, pedantic = F, compile=T)

properties <- expand.grid(trt_eff_scen = c(2), ctrl_prop = c(0.1), icc = c(0.05), n_per_k = c(25), k = c(15),nblock=c(2)) 
properties2 <- properties %>% mutate(row = row_number()) 
#bind to properties
properties <- rbind(properties) %>%
  mutate(t4 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0.4,
                        trt_eff_scen == 5 ~ ctrl_prop+0.4),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0.1,
                        trt_eff_scen == 5 ~ ctrl_prop+0.35),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0,
                        trt_eff_scen == 5 ~ ctrl_prop+0.25),
         t1 = ctrl_prop,
         interim = floor(k/(nblock)),
         draws = ifelse(n_per_k == 75,1250,750))

drawsdat <- properties$draws
properties <- properties %>% dplyr::select(-draws)

adaption <- "both" #this can be early_stopping, arm_dropping, or both
drop_cut <- 0.05
stop_cut <- 0.15
t <- 4

#Running the trial
j <- 1
#The interim clusters are full clusters divided by block
properties$kt2 <- floor(properties$k[j]/(properties$nblock[j]))
properties$kt3 <- floor(properties$k[j]/(properties$nblock[j]))
properties$kt4 <- floor(properties$k[j]/(properties$nblock[j]))


#number of interims
nint <- properties$nblock[j] - 1
#making the cluster and drop matrices
mat <- matrix(0,nrow = t-1,ncol=properties$nblock[j]-1)
drops <- matrix(1,nrow=t-1,ncol=properties$nblock[j]-1)
#interim values
properties_int <- properties[j,]
ints <- 0

i <- 1
#cluster properties
clusters <- data.frame(interim = properties_int$interim, kt2 = properties_int$kt2, kt3 = properties_int$kt3, kt4 = properties_int$kt4)
#make the cluster dataset
intclusters <- makeClusters(t=4,nid=properties$n_per_k[j],
                            t1=clusters$interim,t2=clusters$kt2,
                            t3=clusters$kt3,t4=clusters$kt4)
#test first interim
prop <- list(t1=properties$t1[j],t2=properties$t1[j],t3=properties$t3[j],t4=properties$t4[j])
comp <- (1:t)
o <- list()
beta <- vector()
for(i in 1:length(comp)){
  o[[i]] <- prop[[comp[i]]]/(1-prop[[comp[i]]])
}
beta <- append(beta,log(o[[1]]))
for(i in 2:length(comp)){
  beta <- append(beta,log(o[[i]]/o[[1]]))
}
#intercept 01, beta = log(01), log(02/01), log(03/01),log(04/01)
#when modelling, factor variable of trt
sigma2 <- (pi ^ 2) / 3
theta <- sqrt((properties$icc[j]*sigma2)/(1-properties$icc[j]))
names(theta)<-c("site.(Intercept)")
#fitting model
results <- vector()
resp <- suppressMessages(simulate.formula( ~ factor(trt) + (1|site), nsim = 1, family = binomial, 
                                           newdata = intclusters,newparams = list(beta=beta, theta=theta)))

#make the data stan ready
resp <- as.vector(resp[,1])
N_obs <- dim(intclusters)[1]

intclusters$resp <- resp
siteunique <- paste0(intclusters$trt,intclusters$site)
ascendsite <- as.integer(factor(siteunique,levels=unique(siteunique)))

N_site <- length(unique(ascendsite))
N_trt_groups <- length(unique(intclusters$trt))
data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
             site = ascendsite, trt = as.numeric(intclusters$trt), resp = resp)

#run the model
res <- mod$sample(
  data = data, 
  init = 0,
  iter_warmup = 750,
  iter_sampling = 750,
  chains = 4, 
  parallel_chains = 1,
  adapt_delta = 0.8,
  refresh = 0, 
  max_treedepth=10,
  output_dir=outdir
)
#model results
res$summary(variables=c("beta_trt"))
res$summary(variables=c("beta_trt"), ~quantile(.x, probs = c(0.025, 0.975)))

res$summary(variables=c("pp_trt2","pp_trt3","pp_trt4"))
#descriptive results
desc_dat <- intclusters
desc_dat$resp <- resp
desc_dat %>% group_by(trt) %>% summarise(prop = mean(resp==1))

#DROP trt 2
i <- 1
drops[1] <- 0
mat[,i] <- assignCluster(mat=mat, k=properties$k[j]-floor(properties$k[j]/(properties$nblock[j])), nblock=properties$nblock[j]-1, drops=drops, i=i)
#make dataset for finish analysis
clusters <- data.frame(k = properties_int$k, kt2 = 7, kt3 = 19, kt4 = 19)
fullclusters <- makeClusters(t=4,nid=properties_int$n_per_k,
                             t1=clusters$k,
                             t2=clusters$kt2,
                             t3=clusters$kt3,
                             t4=clusters$kt4)

prop <- list(t1=properties$t1[j],t2=properties$t1[j],t3=properties$t3[j],t4=properties$t4[j])
comp <- (1:t)
o <- list()
beta <- vector()
for(i in 1:length(comp)){
  o[[i]] <- prop[[comp[i]]]/(1-prop[[comp[i]]])
}
beta <- append(beta,log(o[[1]]))
for(i in 2:length(comp)){
  beta <- append(beta,log(o[[i]]/o[[1]]))
}
#intercept 01, beta = log(01), log(02/01), log(03/01),log(04/01)
#when modelling, factor variable of trt
sigma2 <- (pi ^ 2) / 3
theta <- sqrt((properties$icc[j]*sigma2)/(1-properties$icc[j]))
names(theta)<-c("site.(Intercept)")
#fitting model
results <- vector()

resp <- suppressMessages(simulate.formula( ~ factor(trt) + (1|site), nsim = 1, family = binomial, 
                                           newdata = fullclusters,newparams = list(beta=beta, theta=theta)))
resp_dat <- cbind(fullclusters,resp) 
#write over simulated data with the previous interims data
#only writes over data that would have been made in previous interim
names(resp_dat) <- c("iid","site","trt","resp")
resp_dat <- merge(resp_dat,intclusters,by=c("iid","site","trt"),all.x=TRUE)%>%
  within(., resp <- ifelse(!is.na(resp.y), resp.y, resp.x)) %>%
  dplyr::select(-resp.x,-resp.y) %>% 
  arrange(trt,site) %>%
  mutate(site_unique = paste0(trt,site))

resp <- as.vector(resp_dat[,4])
N_obs <- dim(fullclusters)[1]

siteunique <- paste0(fullclusters$trt,fullclusters$site)
ascendsite <- as.integer(factor(siteunique,levels=unique(siteunique)))

N_site <- length(unique(ascendsite))
N_trt_groups <- length(unique(fullclusters$trt))
data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
             site = ascendsite, trt = as.numeric(fullclusters$trt), resp = resp)

res <- mod$sample(
  data = data, 
  init = 0,
  iter_warmup = 750,
  iter_sampling = 750,
  chains = 4, 
  parallel_chains = 1,
  adapt_delta = 0.9,
  refresh = 0, 
  max_treedepth=10,
  output_dir=outdir
)
#model results
res$summary(variables=c("beta_trt"))
res$summary(variables=c("beta_trt"), ~quantile(.x, probs = c(0.025, 0.975)))

res$summary(variables=c("pp_trt2","pp_trt3","pp_trt4"))
#descriptive results
desc_dat <- resp_dat
desc_dat %>% group_by(trt) %>% summarise(prop = mean(resp==1))

#conclude treatment 4 is optimal
