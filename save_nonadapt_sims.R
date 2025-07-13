#Save the simulations into a readable form for testing power etc
#For the nonadaptive simulations
pacman::p_load(future,here,future.apply,tictoc,car,ggforce,rsimsum,dplyr,cmdstanr,rstan,tidyr)
test <- readRDS(here("Data","sims_nonadapt_newscopek25.RDS"))

#new scope
#properties <- expand.grid(trt_eff_scen = c(3,6), ctrl_prop = c(0.1), icc = c(0.2,0.05), k = c(15),n_per_k = c(10,25,50),nblock=1)

#new scope 25 clusters
properties <- expand.grid(trt_eff_scen = c(3,6), ctrl_prop = c(0.1), icc = c(0.2,0.05), k = c(25),n_per_k = c(10,25,50),nblock=1)


#bind to properties
properties <- rbind(properties) %>%
  mutate(t4 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.5,
                        trt_eff_scen == 2 ~ ctrl_prop+0.4,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0.4,
                        trt_eff_scen == 5 ~ ctrl_prop+0.4,
                        trt_eff_scen == 6 ~ ctrl_prop+0.3),
         t3 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.3,
                        trt_eff_scen == 2 ~ ctrl_prop+0.3,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0.1,
                        trt_eff_scen == 5 ~ ctrl_prop+0.35,
                        trt_eff_scen == 6 ~ ctrl_prop+0.2),
         t2 = case_when(trt_eff_scen == 1 ~ ctrl_prop+0.1,
                        trt_eff_scen == 2 ~ ctrl_prop+0.2,
                        trt_eff_scen == 3 ~ ctrl_prop+0,
                        trt_eff_scen == 4 ~ ctrl_prop+0,
                        trt_eff_scen == 5 ~ ctrl_prop+0.3,
                        trt_eff_scen == 6 ~ ctrl_prop+0.1),
         t1 = ctrl_prop,
         interim = floor(k/(nblock)))
properties2 <- properties %>% mutate(row = row_number()) 

nonadapt <- list()
for(j in 1:12){
  for(i in 1:2500){
    nonadapt[[length(nonadapt)+1]] <- test[[j]][[i]]
    nonadapt[[length(nonadapt)]]$sim <- i
    nonadapt[[length(nonadapt)]]$property <- j
  }
}
nonadapt <- bind_rows(nonadapt)
nonadapt_out <- merge(nonadapt,properties2,by.y=c("row"),by.x="property")
saveRDS(nonadapt_out,here("Data","nonadapt_out_newscopek25.RDS"))
