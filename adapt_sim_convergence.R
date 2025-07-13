#AUTHOR: Erin Nolan
#TITLE: PROJECT 2 SIMULATION CONVERGENCE
#PURPOSE: PLOTTING THE CONVERGENCE OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,gtsummary)

outsim2 <- readRDS(here("Data","prob_outsim_newscope.RDS"))
outsim2 <- readRDS(here("Data","prob_outsim_newscope_k25.RDS"))
outsim2 <- readRDS(here("Data","nonadapt_out_newscope.RDS"))
outsim2 <- readRDS(here("Data","nonadapt_out_newscopek25.RDS"))

#Convergence diagnostics-----------------------------
# ESS and rhat
ESS <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>%
  group_by(property,sim) %>%
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), #currently treating missing ess as 'bad' but that isnt necessarily true
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         bess_bulk300 = ifelse(ess_bulk < 300 | is.na(ess_bulk),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         sess_bulk300 = ifelse(sum(bess_bulk300) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  group_by(property) %>%
  filter(variable == "beta_trt[4]") %>%
  summarise(m_ess_bulk = mean(sess_bulk == 1),
            m_ess_tail = mean(sess_tail == 1),
            m_ess_bulk300 = mean(sess_bulk300 == 1),
            mrhat = mean(srhat == 1))

#by treatment group
ESS_trt <- outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>%
  group_by(property,sim,variable) %>%
  mutate(bess_bulk = ifelse(ess_bulk < 400 | is.na(ess_bulk),1,0), #currently treating missing ess as 'bad' but that isnt necessarily true
         bess_tail = ifelse(ess_tail < 400 | is.na(ess_tail),1,0),
         bess_bulk300 = ifelse(ess_bulk < 300 | is.na(ess_bulk),1,0),
         brhat = ifelse(rhat > 1.05 | is.na(rhat),1,0),
         sess_bulk = ifelse(sum(bess_bulk) >= 1,1,0),
         sess_tail = ifelse(sum(bess_tail) >= 1,1,0),
         sess_bulk300 = ifelse(sum(bess_bulk300) >= 1,1,0),
         srhat = ifelse(sum(brhat) >= 1,1,0)) %>%
  group_by(property,variable) %>%
  summarise(m_ess_bulk = mean(sess_bulk == 1),
            m_ess_tail = mean(sess_tail == 1),
            m_ess_bulk300 = mean(sess_bulk300 == 1),
            mrhat = mean(srhat == 1))

#Plot of ess
outsim2 %>% 
  filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]")) %>% 
  ggplot(aes(x=ess_bulk,group=n_per_k)) + 
  geom_density(aes(fill=n_per_k),alpha=0.3,position="identity") +
  facet_wrap(~k)
