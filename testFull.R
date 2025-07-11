#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 FIT BAYESIAN MODELS ADAPTIVE
#PURPOSE: FUNCTION TO MAKE THE GENERATED DATASET AND FIT MODEL

#t = number of treatment groups
#expdat = the generated dataset (without response) to read in
#rho = intra-cluster correlation
#mod = model set by cmdstan_model
#... = the proportion for each treatment group i.e t1=0.2, t2=0.5 etc
testFull <- function(t,expdat,rho,mod,outdir,int_dat,draws,...){
  tic()
  prop <- list(...)
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
  theta <- sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  #fitting model
  results <- vector()
  
  resp <- suppressMessages(simulate.formula( ~ factor(trt) + (1|site), nsim = 1, family = binomial, 
                                             newdata = expdat,newparams = list(beta=beta, theta=theta)))
  resp_dat <- cbind(expdat,resp) 
  names(resp_dat) <- c("iid","site","trt","resp")
  resp_dat <- merge(resp_dat,int_dat,by=c("iid","site","trt"),all.x=TRUE) %>%
    within(., resp <- ifelse(!is.na(resp.y), resp.y, resp.x)) %>%
    dplyr::select(-resp.x,-resp.y) %>% 
    arrange(trt,site) %>%
    mutate(site_unique = paste0(trt,site))
  
expdat$siteunique <- paste0(expdat$trt,expdat$site)
expdat$ascendsite <- as.integer(factor(expdat$siteunique,levels=unique(expdat$siteunique)))
  resp <- as.vector(resp_dat[,4])
  N_obs <- dim(expdat)[1]
  N_site <- length(unique(expdat$siteunique))
  N_trt_groups <- length(unique(expdat$trt))
  data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
               site = expdat$ascendsite, trt = as.numeric(expdat$trt), resp = resp)
  
  res <- mod$sample(
    data = data, 
    init = 0,
    iter_warmup = draws,
    iter_sampling = draws,
    chains = 4, 
    parallel_chains = 1,
    adapt_delta = 0.9,
    refresh = 0, 
    max_treedepth=10,
    output_dir=outdir
    
  )
  print(j)
  time <- toc()
  time <- time$toc - time$tic
  resp_dat %>% group_by(trt) %>% summarise(n = n_distinct(site)) %>% print()
  results <- list(data.frame(res$summary(variables=c("pred_prob_trt","pp_trt2","pp_trt3","pp_trt4","ov_fut","beta_trt","probd_trt2","probd_trt3","probd_trt4","sigma_alpha")),time=time),resp=list(resp))
  
  return(results)
  
}

