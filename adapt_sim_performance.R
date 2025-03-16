#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 SIMULATION PERFORMANCE
#PURPOSE: PLOTTING THE POWER AND ERROR OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,cowplot)

#read in the data
outsim2 <- readRDS(here("Data","prob_outsim.RDS"))
nonadapt_out <- readRDS(here("Data","nonadapt_out.RDS")) %>% mutate(nblock=1)
properties2 <- readRDS(here("Data","properties2.RDS"))
sim_pow <- bind_rows(outsim2,nonadapt_out)

#Power for these sims-----------------------------------------------
#Power defined as number of sims that trt4 is max eff more than 85% of time
trial_success <- sim_pow %>% 
  filter(variable %in% c("pp_trt2","pp_trt3","pp_trt4")) %>%
  pivot_wider(id_cols=c(property,sim,trt_eff_scen,k,icc,n_per_k,nblock,
                        interim1_arm2,interim1_arm3,interim1_arm4,
                        interim2_arm2,interim2_arm3,interim2_arm4,
                        stop_int1,stop_int2),names_from=variable,values_from=mean) %>%
  mutate(type1 = case_when(trt_eff_scen == 3 & is.na(interim1_arm2) & (pp_trt2 >= 0.85 | pp_trt3 >= 0.85 | pp_trt4 >= 0.85) ~ 1, #fixed trials
                   trt_eff_scen == 3 & !is.na(interim1_arm2) & is.na(interim2_arm2) & 
                     ((pp_trt2 >= 0.85 & interim1_arm2 > 0) | (pp_trt3 >= 0.85 & interim1_arm3 > 0) | (pp_trt4 >= 0.85 & interim1_arm4 > 0)) ~ 1, #1 interim trials 
                   trt_eff_scen == 3 & !is.na(interim1_arm2) & !is.na(interim2_arm2) &
                     ((pp_trt2 >= 0.85 & interim2_arm2 > 0) | (pp_trt3 >= 0.85 & interim2_arm3 > 0) | (pp_trt4 >= 0.85 & interim2_arm4 > 0)) ~ 1, #2 interim trials
                   trt_eff_scen == 2 ~ NA,#effect scenario
                   .default = 0), 
         power = case_when(trt_eff_scen == 3 ~ NA,
                           trt_eff_scen == 2 & is.na(interim1_arm2) & pp_trt4 >= 0.85 ~ 1, #fixed trials
                           trt_eff_scen == 2 & !is.na(interim1_arm2) & (pp_trt4 >= 0.85 & interim1_arm4 > 0) ~ 1, #1 interim trials 
                           trt_eff_scen == 2 & !is.na(interim2_arm2) & (pp_trt4 >= 0.85 & interim2_arm4 > 0) ~ 1, #2 interim trials
                           .default = 0),
         correct = ifelse(type1 %in% c(0) | power %in% c(1),1,0)) 


#saveRDS(trial_success,here("Data","trial_success.RDS"))
correct <- trial_success %>% group_by(k,n_per_k,icc,nblock) %>%
  summarise(correct = sum(correct)/n(),
            type1 = sum(type1,na.rm=TRUE)/(n()/2),
            power = sum(power,na.rm=TRUE)/(n()/2)) %>%
  mutate(#k2 = ifelse(k == 15,"1_15","0_25"),
         k3 = factor(k, levels=c("25","15")),
         nint = as.factor(nblock-1)) %>%
  ungroup() %>%
  dplyr::select(-k,-nblock) %>%
  rename(`Mean correct` = correct,
         `Type 1 error` = type1,
         `1 - Type 2 error` = power,
         `N clusters per arm` = k3,
         ICC = icc,
         `N interims` = nint) %>%
  mutate(#`N interims` = forcats::fct_rev(`N interims`),
         `N clusters per arm` = forcats::fct_rev(`N clusters per arm`)) %>%
  dplyr::select(-`Mean correct`)

#Differences in power/type 1 error by trial properties
source(here("Programs","correct_differences.R"))

#Plot type 1 and 2 error
png(filename=here("Output","correct_loopfig1.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = correct, 
                 x = "n_per_k", steps = c("N interims","N clusters per arm","ICC"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 #ylim = c(-0.75,0.05),
                 steps_values_annotate = TRUE, steps_annotation_size = 3, 
                 #steps_annotation_nudge = -0.1,
                 hline_intercept = c(0,1), 
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 8)))) +
  scale_colour_manual(values=c("#9E702B","#B69DE9")) + #"#003331",
  labs(color="Inference",shape="Inference",linetype="Inference",size="Inference")
dev.off()

#Frequency of stops for futility and arm drops
stops <- outsim2 %>% group_by(property,nblock,icc,k,n_per_k,trt_eff_scen) %>% 
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
         `Stopped total` = stop) %>%
  mutate(Scenario = forcats::fct_rev(Scenario),
         `N clusters per arm` = forcats::fct_rev(`N clusters per arm`))


#plot arm dtops and stops for futility
png(filename=here("Output","stops_loopfig2.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = stops, 
                 x = "n_per_k", steps = c("N clusters per arm","ICC","Scenario","N interims"),
                 steps_y_base = -0.04, steps_y_height = 0.03, steps_y_shift = 0.06,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 3, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.35,0.45),
                 y_breaks = c(0,0.1,0.2,0.3,0.4),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 8)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Inference",shape="Inference",linetype="Inference",size="Inference")
dev.off()

#which arms are dropped in which interim
narmdrop <- outsim2 %>% group_by(icc,k,n_per_k,trt_eff_scen,nblock) %>%
  summarise(int1_arm2 = mean(interim1_arm2 == 0 & stop_int1 == 0),
            int1_arm3 = mean(interim1_arm3 == 0 & stop_int1 == 0), 
            int1_arm4 = mean(interim1_arm4 == 0 & stop_int1 == 0),
            int2_arm2 = mean(interim2_arm2 == 0 & stop_int2 == 0),
            int2_arm3 = mean(interim2_arm3 == 0 & stop_int2 == 0), 
            int2_arm4 = mean(interim2_arm4 == 0 & stop_int2 == 0),
            int2_armworst = mean(interim1_arm2 == 0 & interim2_arm3 == 0 & stop_int1 == 0 & stop_int2 == 0)) %>% #this is the best choice for 2 interims effect, removing worst arm
  mutate(nint = nblock-1,
         k2 = factor(k, levels=c("25","15"))) %>%
  ungroup() %>%
  dplyr::select(-nblock,-k) %>%
  pivot_longer(cols = c("int1_arm2","int1_arm3","int1_arm4","int2_arm2","int2_arm3","int2_arm4","int2_armworst"),
               names_to = c("int","arm"),names_sep = "_", values_to = c("prop")) %>%
  pivot_wider(id_cols = c("icc","n_per_k","k2","trt_eff_scen","nint","int"),names_from = "arm",values_from = "prop") %>%
  mutate(trt_eff_scen = ifelse(trt_eff_scen == "2","Effect","Null"),
         int = ifelse(int == "int1","1","2")) %>%
  rename(`Arm 2` = arm2,
         `Arm 3` = arm3,
         `Arm 4` = arm4,
         ICC = icc,
         `Scenario` = trt_eff_scen,
         `N clusters per arm` = k2,
         `Interim` = int,
         `Number of interims` = nint) %>%
  mutate(`N clusters per arm` = forcats::fct_rev(`N clusters per arm`))

#Differences by trial properties and designs
source(here("Programs","arm_drop_differences.R"))


#LOOP PLOT FOR EFFECT SCENARIO
narmdrop_i1 <- narmdrop %>% 
  filter(`Number of interims` == 1,`Interim` == "1",Scenario=="Effect") %>% 
  dplyr::select(-`Number of interims`,-`Interim`,-Scenario,-armworst)
narmdrop_i2 <- narmdrop %>% filter(`Number of interims` == 2,Scenario=="Effect") %>%
  dplyr::select(-`Number of interims`,-Scenario,-armworst)

p1 <- nested_loop_plot(resdf = narmdrop_i1, 
                 x = "n_per_k", steps = c("N clusters per arm","ICC"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                 x_name = " ", y_name = " ",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 ylim = c(-0.5,1),
                 y_breaks = c(0,0.25,0.5,0.75,1),
                 steps_values_annotate = TRUE, steps_annotation_size = 3, 
                 hline_intercept = c(0),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 6)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

p2 <- nested_loop_plot(resdf = narmdrop_i2, 
                       x = "n_per_k", steps = c("N clusters per arm","ICC","Interim"),
                       steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                       x_name = "Sample size per cluster", y_name = " ",
                       spu_x_shift = 30,
                       line_alpha = 0.6,
                       point_alpha = 0.8,
                       ylim = c(-0.75,1),
                       y_breaks = c(0,0.25,0.5,0.75,1),
                       steps_values_annotate = TRUE, steps_annotation_size = 2.6, 
                       hline_intercept = c(0),
                       post_processing = list(
                         add_custom_theme = list(
                           axis.text.x = element_text(angle = -90, 
                                                      vjust = 0.5, 
                                                      size = 6)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")
p2 <- p2 + theme(legend.position="none")

left <- plot_grid(p1,p2,nrow=2)
plots1 <- plot_grid(left,legend,rel_widths = c(1,0.2))

png(filename=here("Output","pdrop_loopfig3.png"),width=8,height=8,res=300,units="in")
plots1 +
  draw_label("Proportion of trials that drop treatment arm", x=0, y=0.5, vjust= 1.5, angle=90,size=10) +
  draw_label("Two interim trials",x=0.75,y=0.48,size=10) +
  draw_label("One interim trials",x=0.75,y=0.98,size=10)
dev.off()

#LOOP PLOT FOR NULL SCENARIO
narmdrop_i1 <- narmdrop %>% 
  filter(`Number of interims` == 1,`Interim` == "1",Scenario=="Null") %>% 
  dplyr::select(-`Number of interims`,-`Interim`,-Scenario,-armworst)
narmdrop_i2 <- narmdrop %>% filter(`Number of interims` == 2,Scenario=="Null") %>%
  dplyr::select(-`Number of interims`,-Scenario,-armworst)

p1 <- nested_loop_plot(resdf = narmdrop_i1, 
                       x = "n_per_k", steps = c("N clusters per arm","ICC"),
                       steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                       x_name = " ", y_name = " ",
                       spu_x_shift = 30,
                       line_alpha = 0.6,
                       point_alpha = 0.8,
                       ylim = c(-0.5,0.5),
                       y_breaks = c(0,0.25,0.5,0.75,1),
                       steps_values_annotate = TRUE, steps_annotation_size = 3, 
                       hline_intercept = c(0),
                       post_processing = list(
                         add_custom_theme = list(
                           axis.text.x = element_text(angle = -90, 
                                                      vjust = 0.5, 
                                                      size = 6)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

p2 <- nested_loop_plot(resdf = narmdrop_i2, 
                       x = "n_per_k", steps = c("N clusters per arm","ICC","Interim"),
                       steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.15,
                       x_name = "Sample size per cluster", y_name = " ",
                       spu_x_shift = 30,
                       line_alpha = 0.6,
                       point_alpha = 0.8,
                       ylim = c(-0.75,0.5),
                       y_breaks = c(0,0.25,0.5,0.75,1),
                       steps_values_annotate = TRUE, steps_annotation_size = 2.6, 
                       hline_intercept = c(0),
                       post_processing = list(
                         add_custom_theme = list(
                           axis.text.x = element_text(angle = -90, 
                                                      vjust = 0.5, 
                                                      size = 6)))) +
  scale_colour_manual(values=c("#B69DE9","#9E702B","#003331")) +
  labs(color="Treatment arm",shape="Treatment arm",linetype="Treatment arm",size="Treatment arm")
p2 <- p2 + theme(legend.position="none")

left <- plot_grid(p1,p2,nrow=2)
plots1 <- plot_grid(left,legend,rel_widths = c(1,0.2))

png(filename=here("Output","pdrop_null_loopfig4.png"),width=8,height=8,res=300,units="in")
plots1 +
  draw_label("Proportion of trials that drop treatment arm", x=0, y=0.5, vjust= 1.5, angle=90,size=10) +
  draw_label("Two interim trials",x=0.75,y=0.48,size=10) +
  draw_label("One interim trials",x=0.75,y=0.98,size=10)
dev.off()
