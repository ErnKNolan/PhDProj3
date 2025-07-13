#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 SIMULATION PERFORMANCE
#PURPOSE: PLOTTING THE POWER AND ERROR OF THE SIMULATIONS
# install.packages("devtools")
#devtools::install_github("matherealize/looplot")

pacman::p_load(here,future.apply,tictoc,car,ggforce,dplyr,cmdstanr,looplot,forcats,tidyr,cowplot)

#read in the data
outsim2 <- readRDS(here("Data","prob_outsim_newscope.RDS"))
outsim2_k25 <- readRDS(here("Data","prob_outsim_newscope_k25.RDS"))
#outsim_newscen <- readRDS(here("Data","prob_outsim_newscen.RDS"))
nonadapt_out <- readRDS(here("Data","nonadapt_out_newscope.RDS")) %>% mutate(nblock=1)
nonadapt_outk25 <- readRDS(here("Data","nonadapt_out_newscopek25.RDS")) %>% mutate(nblock=1)

#primary scenarios
sim_pow <- bind_rows(outsim2,outsim2_k25,nonadapt_out,nonadapt_outk25)
#primary adaptive scenarios
adapt_dat <- bind_rows(outsim2,outsim2_k25)
#to do the new scenarios
#sim_pow <- bind_rows(outsim2,outsim_newscen,nonadapt_out) %>% filter(trt_eff_scen != 4,trt_eff_scen != 2)

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
                   trt_eff_scen %in% c(2,5) ~ NA,#effect scenario
                   .default = 0), 
         power = case_when(trt_eff_scen == 3 ~ NA,
                           trt_eff_scen %in% c(2,5,6) & is.na(interim1_arm2) & pp_trt4 >= 0.85 ~ 1, #fixed trials
                           trt_eff_scen %in% c(2,5,6) & !is.na(interim1_arm2) & (pp_trt4 >= 0.85 & interim1_arm4 > 0) ~ 1, #1 interim trials 
                           trt_eff_scen %in% c(2,5,6) & !is.na(interim2_arm2) & (pp_trt4 >= 0.85 & interim2_arm4 > 0) ~ 1, #2 interim trials
                           .default = 0),
         correct = ifelse(type1 %in% c(0) | power %in% c(1),1,0)) 

#select variables for scaling the power
scale_type1 <- trial_success %>% filter(trt_eff_scen ==3) %>% group_by(k,icc,n_per_k,nblock) %>% filter(row_number()==1) %>%
  dplyr::select(k,icc,n_per_k,nblock)
#get the pp_trts we need to scale
prop1 <- trial_success %>% filter(trt_eff_scen ==3) %>% group_by(k,icc,n_per_k,nblock)%>% 
  mutate(pp_trt2 = case_when(interim1_arm2 %in% c(0) ~ NA, #dont want to include cutpoints for dropped arms
                             interim2_arm2 %in% c(0) ~ NA,
                             .default=pp_trt2),
         pp_trt3 = case_when(interim1_arm3 %in% c(0) ~ NA,
                             interim2_arm3 %in% c(0) ~ NA,
                             .default=pp_trt3),
         pp_trt4 = case_when(interim1_arm4 %in% c(0) ~ NA,
                             interim2_arm4 %in% c(0) ~ NA,
                             .default=pp_trt4)) %>%
  group_split()
#loop function to find the threshold related to a type 1 error of 0.05
scale <- list()
for(i in 1:36){ #1:the number of rows in prop1
  pp_trt2 <- prop1[[i]]$pp_trt2
  pp_trt3 <- prop1[[i]]$pp_trt3
  pp_trt4 <- prop1[[i]]$pp_trt4
  x <- matrix(nrow=1001,ncol=2)
  ind <- 1
  for(y in seq(0,1,by=0.001)){
    x[ind,1] <- y
    x[ind,2] <- sum(pp_trt2 >= y | pp_trt3 >= y | pp_trt4 >= y,na.rm=T)/2500#2500 = number of sims here
    ind <- ind+1
  }
  #finding the relevant quantile that has the type 1 error closest to 0.05
  quantile <- data.frame(prop=x[which.min(abs(x[,2] - 0.05)),1],
                         icc=prop1[[i]][1,]$icc,
                         k=prop1[[i]][1,]$k,
                         n_per_k=prop1[[i]][1,]$n_per_k,
                         nblock=prop1[[i]][1,]$nblock)
  scale[[i]] <- quantile
}
#add in all the thresholds
scale_type1 <- bind_rows(scale)
#determine the scaled power from the thresholds
trial_success2 <- merge(trial_success,scale_type1,by=c("k","icc","n_per_k","nblock")) %>%
  mutate(scale_power = case_when(trt_eff_scen == 3 ~ NA,
                                 trt_eff_scen %in% c(2,5,6) & is.na(interim1_arm2) & pp_trt4 >= prop ~ 1, #fixed trials
                                 trt_eff_scen %in% c(2,5,6) & !is.na(interim1_arm2) & (pp_trt4 >= prop & interim1_arm4 > 0) ~ 1, #1 interim trials 
                                 trt_eff_scen %in% c(2,5,6) & !is.na(interim2_arm2) & (pp_trt4 >= prop & interim2_arm4 > 0) ~ 1, #2 interim trials
                                 .default = 0),
         scale_type1 = case_when(trt_eff_scen == 3 & is.na(interim1_arm2) & (pp_trt2 >= prop | pp_trt3 >= prop | pp_trt4 >= prop) ~ 1, #fixed trials
                                 trt_eff_scen == 3 & !is.na(interim1_arm2) & is.na(interim2_arm2) & 
                                   ((pp_trt2 >= prop & interim1_arm2 > 0) | (pp_trt3 >= prop & interim1_arm3 > 0) | (pp_trt4 >= prop & interim1_arm4 > 0)) ~ 1, #1 interim trials 
                                 trt_eff_scen == 3 & !is.na(interim1_arm2) & !is.na(interim2_arm2) &
                                   ((pp_trt2 >= prop & interim2_arm2 > 0) | (pp_trt3 >= prop & interim2_arm3 > 0) | (pp_trt4 >= prop & interim2_arm4 > 0)) ~ 1, #2 interim trials
                                 trt_eff_scen %in% c(2,5,6) ~ NA,#effect scenario
                                 .default = 0))

#saveRDS(trial_success,here("Data","trial_success.RDS"))
correct <- trial_success2 %>% 
  group_by(k,n_per_k,icc,nblock) %>%
  summarise(correct = sum(correct)/n(),
            type1 = sum(type1,na.rm=TRUE)/(n()/2),
            power = sum(power,na.rm=TRUE)/(n()/2),
            scale_power = sum(scale_power,na.rm=T)/(n()/2),
            scale_type1 = sum(scale_type1,na.rm=T)/(n()/2)) %>%
  mutate(#k2 = ifelse(k == 15,"1_15","0_25"),
         k3 = factor(k, levels=c("25","15")),
         nint = nblock-1,
         n_per_k = factor(n_per_k)) %>%
  ungroup() %>%
  dplyr::select(-k,-nblock) %>%
  rename(`Mean correct` = correct,
         `Type 1 error` = type1,
         `1 - Type 2 error` = power,
         `N clusters per arm` = k3,
         ICC = icc,
         `N interims` = nint,
         `1 - Type 2 error (scaled)` = scale_power) %>%
  mutate(#`N interims` = forcats::fct_rev(`N interims`),
         `N clusters per arm` = forcats::fct_rev(`N clusters per arm`)) %>%
  dplyr::select(-`Mean correct`) %>%
  #monte carlo standard error
  mutate(MCSE_type1min = `Type 1 error`-sqrt((`Type 1 error`*(1-`Type 1 error`))/2500),
         MCSE_type1max = `Type 1 error`+sqrt((`Type 1 error`*(1-`Type 1 error`))/2500),
         MCSE_type2min = `1 - Type 2 error`-sqrt((`1 - Type 2 error`*(1-`1 - Type 2 error`))/2500),
         MCSE_type2max = `1 - Type 2 error`+sqrt((`1 - Type 2 error`*(1-`1 - Type 2 error`))/2500),
         MCSE_type1min2 = scale_type1-sqrt((scale_type1*(1-scale_type1))/2500),
         MCSE_type1max2 = scale_type1+sqrt((scale_type1*(1-scale_type1))/2500),
         MCSE_type2min2 = `1 - Type 2 error (scaled)`-sqrt((`1 - Type 2 error (scaled)`*(1-`1 - Type 2 error (scaled)`))/2500),
         MCSE_type2max2 = `1 - Type 2 error (scaled)`+sqrt((`1 - Type 2 error (scaled)`*(1-`1 - Type 2 error (scaled)`))/2500)) %>%
  dplyr::select(-scale_type1,-MCSE_type1min2,-MCSE_type1max2) %>%
  rename(`N participants per cluster` = n_per_k)

#make the table of the results and MCSE
table_correct <- correct %>%
  mutate(`Type 1 error` = paste0(round(`Type 1 error`,3)," (",round(MCSE_type1min,3),", ",round(MCSE_type1max,3),")"),
         `1 - Type 2 error` = paste0(round(`1 - Type 2 error`,3)," (",round(MCSE_type2min,3),", ",round(MCSE_type2max,3),")"),
         `1 - Type 2 error (scaled)` = paste0(round(`1 - Type 2 error (scaled)`,3)," (",round(MCSE_type2min2,3),", ",round(MCSE_type2max2,3),")"),
         `N participants per arm` = as.numeric(as.character(`N clusters per arm`))*as.numeric(as.character(`N participants per cluster`))) %>%
  dplyr::select(`N interims`,ICC,`N clusters per arm`,`N participants per cluster`,`Type 1 error`,`1 - Type 2 error`,`1 - Type 2 error (scaled)`) %>%
  arrange(`N participants per cluster`,`N clusters per arm`,ICC,`N interims`)
#save table of results
write.csv(table_correct,here("Output","table_correct.csv"))

#Differences in power/type 1 error by trial properties
source(here("Programs","correct_differences.R"))

#Plot type 1 and 2 error
png(filename=here("Output","correct_loopfig1_difforder.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = correct, 
                 pass_through = c("MCSE_type1min","MCSE_type1max","MCSE_type2min","MCSE_type2max",
                                  "MCSE_type2min2","MCSE_type2max2"),
                 x = "N interims", steps = c("N participants per cluster","N clusters per arm","ICC"),
                 steps_y_base = -0.1, steps_y_height = 0.1, steps_y_shift = 0.1,
                 x_name = "Number of interims", y_name = "Proportion of trials",
                 spu_x_shift = 2,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 #ylim = c(-0.75,0.05),
                 steps_values_annotate = TRUE, steps_annotation_size = 3, 
                 #steps_annotation_nudge = -0.1,
                 hline_intercept = c(0,0.05,1), 
                 y_breaks = c(0,0.05,0.25,0.5,0.75,1),
                 post_processing = list(
                   add_custom_theme = list(
                     axis.text.x = element_text(angle = -90, 
                                                vjust = 0.5, 
                                                size = 8)),
                   add_geom_at_position = list(
                     geom = geom_errorbar(
                       aes(ymin=MCSE_type1min,ymax=MCSE_type1max),color="#9E702B",width=0.2,linewidth=0.5
                     )),
                   add_geom_at_position = list(
                     geom = geom_errorbar(
                       aes(ymin=MCSE_type2min,ymax=MCSE_type2max),color="#B69DE9",width=0.2,linewidth=0.5
                     )),
                   add_geom_at_position = list(
                     geom = geom_errorbar(
                       aes(ymin=MCSE_type2min2,ymax=MCSE_type2max2),color="#33638D",width=0.2,linewidth=0.5
                     )))) +
  guides(linetype="none",shape="none",size="none")+
  scale_colour_manual(values=c("#9E702B","#B69DE9","#33638D"),labels=c("Type 1 error","1 - Type 2 error","1 - Type 2 error (scaled)")) + #"#003331",
  labs(color="Inference",shape="Inference",linetype="Inference",size="Inference")
dev.off()

#Frequency of stops for futility and arm drops
#if doing new scen:
#outsim2 <- bind_rows(outsim2,outsim_newscen) %>% filter(trt_eff_scen != 4)


stops <- adapt_dat  %>% group_by(property,nblock,icc,k,n_per_k,trt_eff_scen) %>% 
  filter(nblock != 1) %>%
  summarise(stop1 = mean(stop_int1 > 0, na.rm=TRUE),
            stop2 = mean(stop_int2 > 0 & stop_int1 == 0, na.rm=TRUE)) %>%
  mutate(stop2 = ifelse(nblock < 3, NA, stop2),
         stop = ifelse(!is.na(stop2),stop1+stop2,stop1),
         Scenario = case_when(trt_eff_scen %in% c(2,6) ~ "Effect",
                              trt_eff_scen == 3 ~ "Null",
                              trt_eff_scen == 5 ~ "Unclear opt"),
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

table_stops <- stops %>% rename(`N participants per cluster` = n_per_k) %>%
  dplyr::select(Scenario,`N interims`,ICC,`N clusters per arm`,`N participants per cluster`,`Stopped interim 1`,
                `Stopped interim 2`,`Stopped total`) %>%
  arrange(Scenario,`N participants per cluster`,`N clusters per arm`,ICC,`N interims`)
#Save the table of stops
write.csv(table_stops,here("Output","stops.csv"))


#plot arm drops and stops for futility
png(filename=here("Output","stops_loopfig2.png"),width=10,height=6,res=300,units="in")
nested_loop_plot(resdf = stops, 
                 x = "n_per_k", steps = c("N interims","N clusters per arm","ICC","Scenario"),
                 steps_y_base = -0.04, steps_y_height = 0.03, steps_y_shift = 0.06,
                 x_name = "Sample size per cluster", y_name = "Proportion of trials",
                 spu_x_shift = 30,
                 line_alpha = 0.6,
                 point_alpha = 0.8,
                 steps_values_annotate = TRUE, steps_annotation_size = 3, 
                 hline_intercept = c(0,1), 
                 ylim = c(-0.35,0.30),
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
narmdrop <- adapt_dat %>% group_by(icc,k,n_per_k,trt_eff_scen,nblock) %>%
  summarise(int1_arm2 = mean(interim1_arm2 == 0 & stop_int1 == 0),
            int1_arm3 = mean(interim1_arm3 == 0 & stop_int1 == 0), 
            int1_arm4 = mean(interim1_arm4 == 0 & stop_int1 == 0),
            int2_arm2 = mean(interim2_arm2 == 0 & stop_int2 == 0),
            int2_arm3 = mean(interim2_arm3 == 0 & stop_int2 == 0), 
            int2_arm4 = mean(interim2_arm4 == 0 & stop_int2 == 0)) %>% 
  mutate(nint = nblock-1,
         k2 = factor(k, levels=c("25","15"))) %>%
  ungroup() %>%
  dplyr::select(-nblock,-k) %>%
  pivot_longer(cols = c("int1_arm2","int1_arm3","int1_arm4","int2_arm2","int2_arm3","int2_arm4"),
               names_to = c("int","arm"),names_sep = "_", values_to = c("prop")) %>%
  pivot_wider(id_cols = c("icc","n_per_k","k2","trt_eff_scen","nint","int"),names_from = "arm",values_from = "prop") %>%
  mutate(trt_eff_scen = case_when(trt_eff_scen == "2" ~ "Effect",
                                  trt_eff_scen == "3" ~ "Null",
                                  trt_eff_scen == "5" ~ "Unclear opt",
                                  trt_eff_scen == "6" ~ "Effect"),
         int = ifelse(int == "int1","1","2")) %>%
  rename(`Arm 2` = arm2,
         `Arm 3` = arm3,
         `Arm 4` = arm4,
         ICC = icc,
         `Scenario` = trt_eff_scen,
         `N clusters per arm` = k2,
         `Interim` = int,
         `N interims` = nint,
         `N participants per cluster` = n_per_k) %>%
  mutate(`N clusters per arm` = forcats::fct_rev(`N clusters per arm`))

narmdrop_table <- narmdrop %>%
  dplyr::select(Scenario,`N interims`,ICC,`N clusters per arm`,`N participants per cluster`,
                Interim,`Arm 2`,`Arm 3`,`Arm 4`) %>%
  arrange(Scenario,`N participants per cluster`,`N clusters per arm`,ICC,`N interims`,Interim)

#print this as a table
write.csv(narmdrop_table,(here("Output","narmdrop.csv")))

#Differences by trial properties and designs
source(here("Programs","arm_drop_differences.R"))


#LOOP PLOT FOR EFFECT SCENARIO
narmdrop_i1 <- narmdrop %>% 
  filter(`N interims` == 1,`Interim` == "1",Scenario=="Effect") %>% 
  dplyr::select(-`N interims`,-`Interim`,-Scenario)
narmdrop_i2 <- narmdrop %>% filter(`N interims` == 2,Scenario=="Effect") %>%
  dplyr::select(-`N interims`,-Scenario)

p1 <- nested_loop_plot(resdf = narmdrop_i1, 
                 x = "N participants per cluster", steps = c("N clusters per arm","ICC"),
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
                       x = "N participants per cluster", steps = c("N clusters per arm","ICC","Interim"),
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
  filter(`N interims` == 1,`Interim` == "1",Scenario=="Null") %>% 
  dplyr::select(-`N interims`,-`Interim`,-Scenario)
narmdrop_i2 <- narmdrop %>% filter(`N interims` == 2,Scenario=="Null") %>%
  dplyr::select(-`N interims`,-Scenario)

p1 <- nested_loop_plot(resdf = narmdrop_i1, 
                       x = "N participants per cluster", steps = c("N clusters per arm","ICC"),
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
                       x = "N participants per cluster", steps = c("N clusters per arm","ICC","Interim"),
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
