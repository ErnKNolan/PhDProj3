#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 SIMULATION PERFORMANCE
#PURPOSE: PLOTTING THE POWER AND ERROR OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,cowplot)

#read in the data
outsim2 <- readRDS(here("Data","prob_outsim.RDS"))
#properties2 <- readRDS(here("Data","properties2.RDS"))
#interim 1 and interim 2 of the designs
interim <- readRDS(here("Data","adapt_interim.RDS"))
interim_s <- readRDS(here("Data","adapt_interim_s.RDS"))

#keep the betas
interim2 <- interim %>% group_by(sim,property) %>% filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]"))
interim_s2 <- interim_s %>% group_by(sim,property) %>% filter(variable %in% c("beta_trt[2]","beta_trt[3]","beta_trt[4]"))
interim2 %>% group_by(variable) %>% summarise(essb_low = mean(ess_bulk < 400),
                                              esst_low = mean(ess_tail < 400),
                                              rhat_hi = mean(rhat > 1.05))
interim_s2 %>% group_by(variable) %>% summarise(essb_low = mean(ess_bulk < 400),
                                              esst_low = mean(ess_tail < 400),
                                              rhat_hi = mean(rhat > 1.05))


#The convergence for the interims
ESS <- interim %>% 
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
  dplyr::select(sim,property,sess_bulk,sess_tail,srhat)

#the convergence for the 2nd interim sets
ESS_s <- interim_s %>% 
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
  dplyr::select(sim,property,sess_bulk,sess_tail,srhat)

#put the two interims together
ESS2 <- rbind(ESS,ESS_s) %>%
  group_by(property,sim) %>%
  mutate(essb = ifelse(sum(sess_bulk) >= 1,1,0),
         esst = ifelse(sum(sess_tail) >= 1,1,0),
         rhatt = ifelse(sum(srhat) >= 1,1,0)) %>%
  filter(row_number() == 1) %>%
  dplyr::select(sim,property,essb,esst,rhatt)

#proportion of any poor convergence
ESS2 %>% ungroup() %>% summarise(mean = mean(essb == 1 | esst == 1 | rhatt == 1),
                                 essb = mean(essb == 1),
                                 esst = mean(esst == 1),
                                 rhatt = mean(rhatt == 1))

#RUN THIS ONE IF YOU WANT PLOT FOR ONLY POOR CONVERGE
int_sens <- outsim2 %>% 
  merge(ESS2) %>%
  group_by(property,sim) %>%
  filter(essb == 1 | esst == 1 | rhatt == 1)

#RUN THIS ONE IF YOU WANT PLOT FOR ONLY GOOD CONVERGE
int_sens <- outsim2 %>% 
  merge(ESS2) %>%
  group_by(property,sim) %>%
  filter(essb == 0 , esst == 0 , rhatt == 0)

#proportion of trials (total) that had ess bulk < 400, ess tail < 400, rhat > 1.05
tot <- outsim2 %>%
  merge(ESS2) %>%
  filter(variable == "pred_prob_trt[1]") %>%
  summarise(mean_bulk = mean(essb == 1),
            mean_tail = mean(esst == 1),
            mean_rhat = mean(rhatt == 1))


#trial error types
trial_success <- int_sens %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen,k,icc,n_per_k,nblock),names_from=variable,values_from=mean) %>%
  mutate(correct = ifelse(trt_eff_scen == 3 & pp_trt2 < 0.95 & pp_trt3 < 0.95 & pp_trt4 < 0.95,1,
                          ifelse(trt_eff_scen == 2 & pp_trt4 > 0.95,1,0)),
         type1 = ifelse(trt_eff_scen == 3 & (pp_trt2 >= 0.95 | pp_trt3 >= 0.95 | pp_trt4 >= 0.95),1,
                        ifelse(trt_eff_scen == 2, NA, 0)),
         power = ifelse(trt_eff_scen == 2 & pp_trt4 > 0.95,1,
                        ifelse(trt_eff_scen == 3, NA, 0))) 
correct <- trial_success %>% group_by(k,n_per_k,icc,nblock) %>%
  summarise(correct = sum(correct)/n(),
            alts = sum(trt_eff_scen == 2),
            nulls = sum(trt_eff_scen == 3),
            type1 = sum(type1,na.rm=TRUE)/nulls,
            power = sum(power,na.rm=TRUE)/alts) %>%
  mutate(#k2 = ifelse(k == 15,"1_15","0_25"),
    k3 = factor(k, levels=c("25","15")),
    nint = nblock-1) %>%
  ungroup() %>%
  dplyr::select(-k,-nblock,-alts,-nulls) %>%
  rename(`Total correct` = correct,
         `Type 1 error` = type1,
         `1 - Type 2 error` = power,
         `N clusters per arm` = k3,
         ICC = icc,
         `N interims` = nint)

#trial error types by properties
png(filename=here("Output","correct_loop_sens.png"),width=6,height=4,res=300,units="in")
nested_loop_plot(resdf = correct, 
                 x = "n_per_k", steps = c("N clusters per arm","ICC","N interims"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#003331","#9E702B","#B69DE9")) +
  labs(color="Inference",shape="Inference",linetype="Inference",size="Inference")
dev.off()


#stops in a trial by properties
stops <- int_sens %>% group_by(property,nblock,icc,k,n_per_k,trt_eff_scen) %>% 
  summarise(stop1 = mean(stop_int1 > 0, na.rm=TRUE),
            stop2 = mean(stop_int2 > 0 & stop_int1 == 0, na.rm=TRUE)) %>%
  mutate(stop2 = ifelse(nblock < 3, NA, stop2),
         stop = ifelse(!is.na(stop2),stop1+stop2,stop1),
         Scenario = case_when(trt_eff_scen == 2 ~ "Effect",
                              trt_eff_scen == 3 ~ "Null"),
         nint = nblock-1,
         k2 = factor(k, levels=c("25","15"))) %>% 
  ungroup() %>%
  dplyr::select(-property,-trt_eff_scen,-nblock,-k) %>%
  rename(`N clusters per arm` = k2,
         ICC = icc,
         `N interims` = nint,
         `Stopped interim 1` = stop1,
         `Stopped interim 2` = stop2,
         `Stopped total` = stop)

#stops loop plot
png(filename=here("Output","stops_loop_sens.png"),width=7,height=4,res=300,units="in")
nested_loop_plot(resdf = stops, 
                 x = "n_per_k", steps = c("N clusters per arm","ICC","Scenario","N interims"),
                 steps_y_base = -0.04, steps_y_height = 0.03, steps_y_shift = 0.06,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.35,0.45),
                 y_breaks = c(0,0.1,0.2,0.3,0.4),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Inference",shape="Inference",linetype="Inference",size="Inference")
dev.off()

#number of clusters PER ARM
narm2 <- int_sens %>% 
  group_by(icc,k,n_per_k,trt_eff_scen,nblock) %>% 
  summarise(arm2 = mean(arm2),arm3 = mean(arm3), arm4 = mean(arm4)) %>%
  mutate(arm2 = (arm2 / (k*3)),
         arm3 = (arm3 / (k*3)),
         arm4 = (arm4 / (k*3))) %>%
  mutate(Scenario = case_when(trt_eff_scen == 2 ~ "Effect",
                              trt_eff_scen == 3 ~ "Null"),
         nint = nblock-1,
         k2 = factor(k, levels=c("25","15"))) %>%
  ungroup() %>%
  dplyr::select(-trt_eff_scen,-nblock,-k) %>%
  rename(`N clusters per arm` = k2,
         ICC = icc,
         `N interims` = nint,
         `Arm 2` = arm2,
         `Arm 3` = arm3,
         `Arm 4` = arm4)

#number of clusters loop plot
png(filename=here("Output","narm_loop_sens.png"),width=6,height=4,res=300,units="in")
nested_loop_plot(resdf = narm2, 
                 x = "n_per_k", steps = c("ICC","N interims","Scenario","N clusters per arm"),
                 steps_y_base = -0.03, steps_y_height = 0.03, steps_y_shift = 0.05,
                 x_name = "Sample size per cluster", y_name = "Proportion of clusters in each arm",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 4, 
                 hline_intercept = 0, 
                 ylim = c(-0.30,0.65),
                 y_breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 10)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")
dev.off()


#Arms dropped by interim
narmdrop <- int_sens %>% group_by(icc,k,n_per_k,trt_eff_scen,nblock) %>%
  summarise(int1_arm2 = mean(interim1_arm2 == 0),
            int1_arm3 = mean(interim1_arm3 == 0), 
            int1_arm4 = mean(interim1_arm4 == 0),
            int2_arm2 = mean(interim2_arm2 == 0),
            int2_arm3 = mean(interim2_arm3 == 0), 
            int2_arm4 = mean(interim2_arm4 == 0)) %>%
  mutate(nint = nblock-1,
         k2 = factor(k, levels=c("25","15"))) %>%
  ungroup() %>%
  dplyr::select(-nblock,-k) %>%
  pivot_longer(cols = c("int1_arm2","int1_arm3","int1_arm4","int2_arm2","int2_arm3","int2_arm4"),
               names_to = c("int","arm"),names_sep = "_", values_to = c("prop")) %>%
  pivot_wider(id_cols = c("icc","n_per_k","k2","trt_eff_scen","nint","int"),names_from = "arm",values_from = "prop") %>%
  mutate(trt_eff_scen = ifelse(trt_eff_scen == "2","Effect","Null"),
         int = ifelse(int == "int1","1","2")) %>%
  rename(`Arm 2` = arm2,
         `Arm 3` = arm3,
         `Arm 4` = arm4,
         ICC = icc,
         `Scenario` = trt_eff_scen,
         `n clusters per arm` = k2,
         `Interim` = int,
         `Number of interims` = nint)

#which arms are dropped at which interims
narmdrop_i1 <- narmdrop %>% 
  filter(`Number of interims` == 1,`Interim` == "1") %>% 
  dplyr::select(-`Number of interims`,-`Interim`)
narmdrop_i2 <- narmdrop %>% filter(`Number of interims` == 2) %>%
  dplyr::select(-`Number of interims`)

p1 <- nested_loop_plot(resdf = narmdrop_i1, 
                       x = "n_per_k", steps = c("n clusters per arm","ICC","Scenario"),
                       steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                       x_name = " ", y_name = " ",
                       spu_x_shift = 30,
                       line_alpha = 0.6,
                       point_alpha = 0.8,
                       ylim = c(-0.75,1),
                       y_breaks = c(0,0.25,0.5,0.75,1),
                       steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                       hline_intercept = c(0,1),
                       post_processing = list(
                         add_custom_theme = list(
                           axis.text.x = element_text(angle = -90, 
                                                      vjust = 0.5, 
                                                      size = 5)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

p2 <- nested_loop_plot(resdf = narmdrop_i2, 
                       x = "n_per_k", steps = c("n clusters per arm","ICC","Scenario","Interim"),
                       steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                       x_name = "Sample size per cluster", y_name = " ",
                       spu_x_shift = 30,
                       line_alpha = 0.6,
                       point_alpha = 0.8,
                       ylim = c(-0.75,1),
                       y_breaks = c(0,0.25,0.5,0.75,1),
                       steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                       hline_intercept = c(0,1),
                       post_processing = list(
                         add_custom_theme = list(
                           axis.text.x = element_text(angle = -90, 
                                                      vjust = 0.5, 
                                                      size = 5)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")
p2 <- p2 + theme(legend.position="none")

left <- plot_grid(p1,p2,nrow=2)
plots1 <- plot_grid(left,legend,rel_widths = c(1,0.2))

png(filename=here("Output","pdrop_loop_sens.png"),width=10,height=6,res=300,units="in")
plots1 +
  draw_label("Proportion of trials that drop:", x=0, y=0.5, vjust= 1.5, angle=90)
dev.off()

