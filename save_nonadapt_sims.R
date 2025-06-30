#Save the simulations into a readable form for testing power etc
#For the nonadaptive simulations
test <- readRDS(here("Data","simulations.RDS"))

properties <- expand.grid(trt_eff_scen = c(5), ctrl_prop = c(0.1), icc = c(0.05,0.2), n_per_k = c(25,50,75), k = c(15,25),nblock=1)

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
         draws = ifelse(n_per_k == 75,1250,750))

nonadapt <- list()
for(j in 1:21){
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
