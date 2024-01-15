#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 FIT BAYESIAN MODELS ADAPTIVE
#PURPOSE: FUNCTION TO MAKE THE GENERATED DATASET AND FIT MODEL

#t = number of treatment groups
#expdat = the generated dataset (without response) to read in
#rho = intra-cluster correlation
#mod = model set by cmdstan_model
#... = the proportion for each treatment group i.e t1=0.2, t2=0.5 etc
testFull <- function(t,expdat,rho,mod,outdir,int_dat,properties_int,int_prior_mu,int_prior_sd,trt_prior_mu,trt_prior_sd,...){
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
  return <- tryCatch({
    
    resp <- suppressMessages(simulate.formula( ~ factor(trt) + (1|site), nsim = 1, family = binomial, 
                                               newdata = expdat,newparams = list(beta=beta, theta=theta)))
    resp_dat <- cbind(expdat,resp) 
    names(resp_dat) <- c("iid","site","trt","resp")
    resp_dat <- merge(resp_dat,int_dat,by=c("iid","site","trt"),all.x=TRUE) %>%
      within(., resp <- ifelse(!is.na(resp.y), resp.y, resp.x)) %>%
      dplyr::select(-resp.x,-resp.y) %>% 
      arrange(site)

    #SECOND HALF OF DATA ONLY
    #WIP
    #adding in the interim properties to determine drops
    testdat <- resp_dat %>% cbind(properties_int)
      
    #if stop, then we're just analysing the interim data
    #TODO - just read out the interim data instead of rerunning
    if(properties_int$stop==1) {
      testdat <- testdat %>%
      mutate(dropt2 = 0,dropt3 = 0,dropt4 = 0)
    #making the rules to drop pre-interim data
    #if an arm if being dropped then we keep the pre-interim data to analyse
    } else {
      testdat <- testdat %>%
      mutate(dropt2 = case_when(drop=="trt2" ~ 0,drop != "trt2" ~ interim),
      dropt3 = case_when(drop=="trt3" ~ 0,drop != "trt3" ~ interim),
      dropt4 = case_when(drop=="trt4" ~ 0,drop != "trt4" ~ interim))
      }
    
    #get rid of the pre interim data and just use the posteriors as priors
    testdat <- testdat %>% 
      group_by(trt) %>% 
      mutate(keep = case_when(trt == 1 & row_number() > interim*n_per_k ~ 1,
                              trt == 2 & row_number() > dropt2*n_per_k ~ 1,
                              trt == 3 & row_number() > dropt3*n_per_k ~ 1,
                              trt == 4 & row_number() > dropt4*n_per_k ~ 1)) %>%
      filter(!is.na(keep)) 
    #%>%
    #  dplyr::select(iid,site,trt,resp)
    
    #if an arm is dropped, dont use its prior as we will re-analyse that data at final
    if(properties_int$drop == "trt2"){
      trt_prior_mu[1] <- 0
      trt_prior_sd[1] <- 2
    } else if(properties_int$drop == "trt3"){
      trt_prior_mu[2] <- 0
      trt_prior_sd[2] <- 2
    } else if(properties_int$drop == "trt4"){
      trt_prior_mu[3] <- 0
      trt_prior_sd[3] <- 2
    } else {
      x <- 1 #if I feel later to add in something for no drops
    }
    
    #END WIP
    
    
    resp <- unlist(as.vector(testdat[,4]))
    N_obs <- dim(testdat)[1]
    N_site <- length(unique(testdat$site))
    N_trt_groups <- length(unique(testdat$trt))
    data <- list(N_obs = N_obs, N_site = N_site, N_trt_groups = N_trt_groups, 
                 site = testdat$site, trt = as.numeric(testdat$trt), resp = resp,
                 int_prior_mu = int_prior_mu, int_prior_sd = int_prior_sd,
                 trt_prior_mu = trt_prior_mu, trt_prior_sd = trt_prior_sd)

      res <- mod_full$sample(
      data = data, 
      init = 0,
      iter_warmup = 250,
      iter_sampling = 250,
      chains = 4, 
      parallel_chains = 1,
      adapt_delta = 0.99,
      refresh = 0, 
      max_treedepth=12,
      output_dir=outdir)
    print(j)
    time <- toc()
    time <- time$toc - time$tic
    results <- list(data.frame(res$summary(variables=c("pred_prob_trt","beta_trt","pp_trt2","pp_trt3","pp_trt4","thetarep")),time=time))
  },
  
  error=function(e) { message(conditionMessage(e)) 
    res <- list(data.frame(variable=NA,mean=NA,median=NA,sd=NA,mad=NA,q5=NA,q95=NA,rhat=NA,ess_bulk=NA,ess_tail=NA,time=NA))
  })
  
  return(results)
  
}

