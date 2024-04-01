erin.long<-readRDS("C:/Users/JDizon/Downloads/test_dat.RDS")

resp <- as.vector(erin.long[,4])
N_obs <- dim(erin.long)[1]
N_site <- length(unique(erin.long$site))
N_trt_groups <- length(unique(erin.long$trt))
data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
             site = erin.long$site, trt = as.numeric(erin.long$trt), resp = resp, repeats=rep(1:N_obs))


erin.reduced <- erin.long %>%
  group_by(site,trt,sim_1)%>%
  summarise(repeats=n())%>%
  data.frame()

resp.reduced <- as.vector(erin.reduced[,3])
N_obs.reduced <- dim(erin.reduced)[1]
N_site.reduced <- length(unique(erin.reduced$site))
N_trt_groups.reduced <- length(unique(erin.reduced$trt))
data.reduced <- list(N_obs = N_obs.reduced, N_site = N_site.reduced, N_trt_groups = N_trt_groups.reduced, 
             site = unlist(erin.reduced$site), trt = as.numeric(erin.reduced$trt), resp = resp.reduced, repeats=erin.reduced$repeats)


mod.erin<-cmdstan_model("C:/Users/Jdizon/MARS/adapt_arm 1.stan", pedantic = F, compile=T);
fit.erin <<- mod.erin$sample(
  data = data, 
  init = .75,
  iter_warmup = 500,
  iter_sampling = 250,
  chains = 4, 
  parallel_chains = 4,
  adapt_delta = 0.8,
  refresh = 0, 
  max_treedepth=10,
)

mod.erin.reduced<-cmdstan_model("C:/Users/Jdizon/MARS/adapt_arm 1 - reduced.stan", pedantic = F, compile=T);
fit.erin.reduced <<- mod.erin.reduced$sample(
  data = data.reduced, 
  init = .75,
  iter_warmup = 500,
  iter_sampling = 250,
  chains = 4, 
  parallel_chains = 4,
  adapt_delta = 0.8,
  refresh = 0, 
  max_treedepth=10
)

mean(unlist(fit.erin$summary(variables="ypred")[c(2)]))
mean(unlist(fit.erin.reduced$summary(variables="ypred")[c(2)]))
