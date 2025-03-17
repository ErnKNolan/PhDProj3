
#differences by trial property in results section
#this is referenced in adapt_sim_performance.R

#Number of participants per cluster
#Arm 2
narm_nperk_arm2 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(`N clusters per arm`,ICC,`Number of interims`,Interim),
              values_from = `Arm 2`,
              names_from = n_per_k) %>%
  mutate(diff75_50 = `75` - `50`, diff50_25 = `50` - `25`)
#Arm 3
narm_nperk_arm3 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(`N clusters per arm`,ICC,`Number of interims`,Interim),
              values_from = `Arm 3`,
              names_from = n_per_k) %>%
  mutate(diff75_50 = `75` - `50`, diff50_25 = `50` - `25`)
#Arm 4
narm_nperk_arm4 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(`N clusters per arm`,ICC,`Number of interims`,Interim),
              values_from = `Arm 4`,
              names_from = n_per_k) %>%
  mutate(diff75_50 = `75` - `50`, diff50_25 = `50` - `25`)

#Dropping worst arm for the 2 interim trial
narm_nperk_worst <- narmdrop %>% filter(Scenario=="Effect", 
                                  (Interim == 2 & `Number of interims` == 2)) %>%
  pivot_wider(id_cols=c(`N clusters per arm`,ICC,`Number of interims`,Interim),
              values_from = armworst,
              names_from = n_per_k) %>%
  mutate(diff75_50 = `75` - `50`, diff50_25 = `50` - `25`)

#Number of clusters
#Arm 2
narm_k_arm2 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(n_per_k,ICC,`Number of interims`,Interim),
              values_from = `Arm 2`,
              names_from = `N clusters per arm`) %>%
  mutate(diff25_15 = `25` - `15`)

#Arm 3
narm_k_arm3 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(n_per_k,ICC,`Number of interims`,Interim),
              values_from = `Arm 3`,
              names_from = `N clusters per arm`) %>%
  mutate(diff25_15 = `25` - `15`)

#Arm 4
narm_k_arm4 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(n_per_k,ICC,`Number of interims`,Interim),
              values_from = `Arm 4`,
              names_from = `N clusters per arm`) %>%
  mutate(diff25_15 = `25` - `15`)

#Dropping worst arm for 2 interim trial
narm_k_worst <- narmdrop %>% filter(Scenario=="Effect", 
                              (Interim == 2 & `Number of interims` == 2)) %>%
  pivot_wider(id_cols=c(n_per_k,ICC,`Number of interims`,Interim),
              values_from = armworst,
              names_from = `N clusters per arm`) %>%
  mutate(diff25_15 = `25` - `15`)

#ICC
#Arm 2
narm_icc_arm2 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,`Number of interims`,Interim),
              values_from = `Arm 2`,
              names_from = ICC) %>%
  mutate(diff02_005 = `0.2` - `0.05`)

#Arm 3
narm_icc_arm3 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,`Number of interims`,Interim),
              values_from = `Arm 3`,
              names_from = ICC) %>%
  mutate(diff02_005 = `0.2` - `0.05`)

#Arm 4
narm_icc_arm4 <- narmdrop %>% filter(Scenario=="Effect") %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,`Number of interims`,Interim),
              values_from = `Arm 4`,
              names_from = ICC) %>%
  mutate(diff02_005 = `0.2` - `0.05`)

#Dropping worst arm for 2 interim trial
narm_icc_worst <- narmdrop %>% filter(Scenario=="Effect", 
                                    (Interim == 2 & `Number of interims` == 2)) %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,`Number of interims`,Interim),
              values_from = armworst,
              names_from = ICC) %>%
  mutate(diff02_005 = `0.2` - `0.05`)

#Number of interims
#Arm 2
narm_nint_arm2 <- narmdrop %>% filter(Scenario=="Effect", Interim == `Number of interims`) %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,ICC),
              values_from = `Arm 2`,
              names_from = `Number of interims`) %>%
  mutate(diff2_1 = `2` - `1`)

#Arm 3
narm_nint_arm3 <- narmdrop %>% filter(Scenario=="Effect", Interim == `Number of interims`) %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,ICC),
              values_from = `Arm 3`,
              names_from = `Number of interims`) %>%
  mutate(diff2_1 = `2` - `1`)

#Arm 4
narm_nint_arm4 <- narmdrop %>% filter(Scenario=="Effect", Interim == `Number of interims`) %>%
  pivot_wider(id_cols=c(n_per_k,`N clusters per arm`,ICC),
              values_from = `Arm 4`,
              names_from = `Number of interims`) %>%
  mutate(diff2_1 = `2` - `1`)
