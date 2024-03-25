#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 SIMULATION PERFORMANCE
#PURPOSE: PLOTTING THE POWER AND ERROR OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,cowplot)

#read in the data
#outsim2 <- readRDS(here("Data","outsim_adapt_opt.RDS"))
#outsim2 <- readRDS(here("Data","outsim_adapt_prob.RDS"))
#outsim2 <- readRDS(here("Data","outsim_nonadapt.RDS"))
#properties2 <- readRDS(here("Data","properties2.RDS"))

#Power for these sims-----------------------------------------------
#Power defined as number of sims that trt4 is max eff more than 95% of time
trial_success <- outsim2 %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen,stop_int1,stop_int2),names_from=variable,values_from=mean) %>%
  mutate(bayesr = ifelse(trt_eff_scen == 3 & (pp_trt2 >= 0.95 | pp_trt3 >= 0.95 | pp_trt4 >= 0.95),1,
                         ifelse(trt_eff_scen %in% c(1,2) & pp_trt4 >= 0.95,1,0))) 
power <- trial_success %>% group_by(property) %>%
  summarise(bayesr = sum(bayesr)/n())


#TEST
trial_success <- outsim2 %>% 
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
            type1 = sum(type1,na.rm=TRUE)/(n()/2),
            power = sum(power,na.rm=TRUE)/(n()/2)) %>%
  mutate(k2 = ifelse(k == 15,"1_15","0_25"),
         nint = nblock-1) %>%
  ungroup() %>%
  dplyr::select(-k,-nblock) %>%
  rename(`Total correct` = correct,
         `Type 1 error` = type1,
         `1 - Type 2 error` = power,
         `Number of clusters` = k2,
         ICC = icc,
         `N interims` = nint)

png(filename=here("Output","correct_loop.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = correct, 
                 x = "n_per_k", steps = c("Number of clusters","ICC","N interims"),
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
#END TEST

#TEST STOPS
stops <- clusters %>% group_by(property,nblock,icc,k,n_per_k,trt_eff_scen) %>% 
  summarise(stop1 = mean(stop_int1 > 0, na.rm=TRUE),
            stop2 = mean(stop_int2 > 0 & stop_int1 == 0, na.rm=TRUE)) %>%
  mutate(stop2 = ifelse(nblock < 3, NA, stop2),
         stop = ifelse(!is.na(stop2),stop1+stop2,stop1)) %>% 
  ungroup() %>%
  dplyr::select(-property)

png(filename=here("Output","stops_loop.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = stops, 
                 x = "n_per_k", steps = c("k","icc","trt_eff_scen","nblock"),
                 steps_y_base = -0.02, steps_y_height = 0.01, steps_y_shift = 0.02,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.125,0.20),
                 y_breaks = c(0,0.05,0.1,0.15,0.2),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Inference",shape="Inference",linetype="Inference",size="Inference")
dev.off()
#END TEST

#TEST N K PER ARM
narm <- clusters %>% group_by(icc,k,n_per_k,trt_eff_scen,nblock) %>% summarise(arm2 = mean(arm2),arm3 = mean(arm3), arm4 = mean(arm4))

png(filename=here("Output","narm_loop.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = narm, 
                 x = "n_per_k", steps = c("icc","nblock","trt_eff_scen","k"),
                 steps_y_base = -2, steps_y_height = 2, steps_y_shift = 2,
                 x_name = "Sample size per cluster", y_name = "Mean number of clusters",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0), 
                 y_breaks = c(0,10,20,30),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")
dev.off()

narmdrop <- clusters %>% group_by(icc,k,n_per_k,trt_eff_scen,nblock) %>%
  summarise(int1_arm2 = mean(interim1_arm2 == 0),
            int1_arm3 = mean(interim1_arm3 == 0), 
            int1_arm4 = mean(interim1_arm4 == 0),
            int2_arm2 = mean(interim2_arm2 == 0),
            int2_arm3 = mean(interim2_arm3 == 0), 
            int2_arm4 = mean(interim2_arm4 == 0)) %>%
  mutate(nint = nblock-1) %>%
  dplyr::select(-nblock) %>%
  pivot_longer(cols = c("int1_arm2","int1_arm3","int1_arm4","int2_arm2","int2_arm3","int2_arm4"),
               names_to = c("int","arm"),names_sep = "_", values_to = c("prop")) %>%
  pivot_wider(id_cols = c("icc","n_per_k","k","trt_eff_scen","nint","int"),names_from = "arm",values_from = "prop") %>%
  rename(`Arm 2` = arm2,
         `Arm 3` = arm3,
         `Arm 4` = arm4,
         ICC = icc,
         `Scenario` = trt_eff_scen,
         `n clusters per arm` = k,
         `Interim` = int,
         `Number of interims` = nint)


narmdrop_i1 <- narmdrop %>% 
  filter(`Number of interims` == 1,`Interim` == "int1") %>% 
  dplyr::select(-`Number of interims`,-`Interim`)
narmdrop_i2 <- narmdrop %>% filter(`Number of interims` == 2) %>%
  dplyr::select(-`Number of interims`)

p1 <- nested_loop_plot(resdf = narmdrop_i1, 
                 x = "n_per_k", steps = c("n clusters per arm","ICC","Scenario"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials that drop:",
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
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

p2 <- nested_loop_plot(resdf = narmdrop_i2, 
                       x = "n_per_k", steps = c("n clusters per arm","ICC","Scenario","Interim"),
                       steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                       x_name = "Sample size per cluster", y_name = "Proportion of trials that drop:",
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
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")

left <- plot_grid(p1,p2,nrow=2)
plots1 <- plot_grid(left,legend,rel_widths = c(1,0.3))

#png(filename=here("Output","pdrop_loop.png"),width=10,height=6,res=300,units="in")
plots1
#dev.off()
#END TEST

#MC Error for power (ref Morris 2019 table 6)
MCSE_power <- data.frame(MCSE = sqrt((power$bayesr*(1-power$bayesr))/ 2500), property=power$property)

#BACK INTO POWER
#merge outsim2 back into outsim
outsim3 <- merge(outsim2,power,by="property") %>%
  group_by(property) %>%
  filter(sim == 1, variable == "pred_prob_trt[4]")

#plot the power over different features
#setting up the data
power_plot_opt <- outsim3 %>%
  mutate(join = paste0("n per k=",n_per_k,", k=",k),
         sample_n = n_per_k*k, #total sample size
         deff = 1+icc*((n_per_k)-1), #design effect
         eff_n = sample_n / deff, #effective sample size
         iccf = paste0("ICC = ",factor(icc)),
         ctrlpf = factor(ctrl_prop),
         test = paste0(icc,ctrl_prop),
         kf = factor(k))

#save the power plot
#saveRDS(power_plot_opt,here("Data","nonadapt_power.RDS"))

#graph comparing adaptive and non-adaptive
#cleaning the two datasets
#power_plot <- readRDS(here("Data","power_plot.RDS"))
power_plot_opt <- readRDS(here("Data","power_plot_opt.RDS"))
power_plot_prob <- readRDS(here("Data","power_plot_prob.RDS"))
nonadapt_power <- readRDS(here("Data","nonadapt_power.RDS"))
nonadapt_plot <- nonadapt_power %>% ungroup() %>% filter(k %in% c(5,10), ctrl_prop == 0.1) %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr) 
nest_plot <- power_plot_opt %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)
nest_plot_prob <- power_plot_prob %>% ungroup() %>% dplyr::select(trt_eff_scen,icc,n_per_k,k,bayesr)

#combine the two datasets
loop_plot <- inner_join(nest_plot,nonadapt_plot,by=c("icc","trt_eff_scen","n_per_k","k")) %>%
  rename(nonadapt = bayesr.y,
         both = bayesr.x) %>%
  inner_join(nest_plot_prob) %>%
  rename(prob = bayesr) %>%
  mutate(k = fct_rev(factor(k))) %>%
  rename(`Optimisation ties` = both,
         `Non-adpative` = nonadapt,
         `Probability ties` = prob,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
  mutate(Scenario = case_when(Scenario == 1 ~ "Strong effect",
                            Scenario == 2 ~ "Mid/Weak effect",
                            Scenario == 3 ~ "Null effect"))
loop_plot$Scenario <- factor(loop_plot$Scenario, levels = c("Strong effect","Mid/Weak effect","Null effect"))
loop_plot_effs <- loop_plot %>% filter(Scenario != "Null effect")
loop_plot_null <- loop_plot %>% filter(Scenario =="Null effect") %>% dplyr::select(-Scenario)

#The graph
png(filename=here("Output","nest_loop.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_effs, 
                     x = "n_per_k", steps = c("ICC","k", "Scenario"),
                     steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                     x_name = "Sample size per cluster", y_name = "Power",
                     spu_x_shift = 25,
                     line_alpha = 0.6,
                     point_alpha = 0.6,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                     hline_intercept = c(0,0.8,0.9,1), 
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90, 
                                                    vjust = 0.5, 
                                                    size = 5))))+
  scale_colour_manual(values=c("#8D4585","#003151","orange"))+
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#Graph for the null scenario
png(filename=here("Output","nest_loop_null.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = loop_plot_null, 
                 x = "n_per_k", steps = c("ICC","k"),
                 steps_y_base = -0.025, steps_y_height = 0.025, steps_y_shift = 0.025,
                 x_name = "Sample size per cluster", y_name = "Type 1 error",
                 spu_x_shift = 50,
                 line_alpha = 0.6,
                 point_alpha = 0.6,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0), 
                 ylim = c(-0.1,0.25),
                 y_breaks = c(0,0.05,0.1,0.15,0.2,0.25),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#8D4585","#003151","orange"))+
  labs(color="Adaptive design",shape="Adaptive design",linetype="Adaptive design",size="Adaptive design")
dev.off()

#Trial properties
#trial_props <- readRDS(here("Data","adaptopt_trial_props.RDS"))
#trial_props <- readRDS(here("Data","adaptprob_trial_props.RDS"))

stops <- trial_props %>% group_by(property) %>% summarise(stops = (sum(stop)/n()))
trt2drps <- trial_props %>% group_by(property) %>% summarise(trt2drps = (sum(drop=="trt2")/n()))
trt3drps <- trial_props %>% group_by(property) %>% summarise(trt3drps = (sum(drop=="trt3")/n()))
trt4drps <- trial_props %>% group_by(property) %>% summarise(trt4drps = (sum(drop=="trt4")/n()))
trial_drops <- trial_props %>% group_by(property) %>% filter(row_number() == 1)
trial_drops <- merge(trial_drops,stops,by="property") %>%
  merge(trt2drps,by="property") %>%
  merge(trt3drps,by="property") %>%
  merge(trt4drps,by="property") %>%
  dplyr::select(trt_eff_scen,icc,n_per_k,k,stops,trt2drps,trt3drps,trt4drps) %>%
  mutate(k = fct_rev(factor(k))) %>% 
  rename(`Stop for futility` = stops,
         `Arm 2 dropped` = trt2drps,
         `Arm 3 dropped` = trt3drps,
         `Arm 4 dropped` = trt4drps,
         ICC = icc,
         Scenario = trt_eff_scen) %>%
  mutate(Scenario = case_when(Scenario == 1 ~ "Strong effect",
                              Scenario == 2 ~ "Mid/Weak effect",
                              Scenario == 3 ~ "Null effect"))
trial_drops$Scenario <- factor(trial_drops$Scenario, levels = c("Strong effect","Mid/Weak effect","Null effect"))


png(filename=here("Output","trialprobprob_drops.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = trial_drops, 
                 x = "n_per_k", steps = c("ICC","k","Scenario"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion",
                 spu_x_shift = 50,
                 steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
                 hline_intercept = c(0), 
                 ylim = c(-0.75,1),
                 y_breaks = c(0,0.25,0.5,0.75,1),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 5))))+
  scale_colour_manual(values=c("#48157F","#29AF7F","#f7cb48f9","#a65c85ff"))+
  labs(color="Trial outcome",shape="Trial outcome",linetype="Trial outcome",size="Trial outcome")

dev.off()
