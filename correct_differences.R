#Differences in the power / type 1 error by trial property
#This is referenced in adapt_sim_performance.R

#Number of participants per cluster
#Power
power_nperk <- correct %>% 
  pivot_wider(id_cols = c(ICC,`N interims`,`N clusters per arm`),
              values_from=`1 - Type 2 error`,
              names_from=n_per_k) %>% 
  mutate(diff75_50 = `75` - `50`,diff50_25 = `50` - `25`)

#Type 1 error
type1_nperk <- correct %>% 
  pivot_wider(id_cols = c(ICC,`N interims`,`N clusters per arm`),
              values_from=`Type 1 error`,
              names_from=n_per_k) %>% 
  mutate(diff75_50 = `75` - `50`,diff50_25 = `50` - `25`)


#Number of clusters
#Power
power_k <- correct %>% 
  pivot_wider(id_cols = c(n_per_k,`N interims`,ICC),
              values_from=`1 - Type 2 error`,
              names_from=`N clusters per arm`) %>% 
  mutate(diff25_15 = `25` - `15`)

#Type 1 error
type1_k <- correct %>% 
  pivot_wider(id_cols = c(n_per_k,`N interims`,ICC),
              values_from=`Type 1 error`,
              names_from=`N clusters per arm`) %>% 
  mutate(diff25_15 = `25` - `15`)


#ICC
#Power
power_icc <- correct %>% 
  pivot_wider(id_cols = c(n_per_k,`N interims`,`N clusters per arm`),
              values_from=`1 - Type 2 error`,
              names_from=ICC) %>% 
  mutate(diff02 = `0.05` - `0.2`)

#Type 1 error
type1_icc <- correct %>% 
  pivot_wider(id_cols = c(n_per_k,`N interims`,`N clusters per arm`),
              values_from=`Type 1 error`,
              names_from=ICC) %>% 
  mutate(diff02 = `0.05` - `0.2`)

#Number of interim analyses
#Power
power_nint <- correct %>% 
  pivot_wider(id_cols = c(n_per_k,ICC,`N clusters per arm`),
              values_from=`1 - Type 2 error`,
              names_from=`N interims`) %>% 
  mutate(diff12 = `1` - `2`, diff01 = `0`-`1`,diff02 = `0`-`2`)

#Type 1 error
type1_nint <- correct %>% 
  pivot_wider(id_cols = c(n_per_k,ICC,`N clusters per arm`),
              values_from=`Type 1 error`,
              names_from=`N interims`) %>% 
  mutate(diff12 = `1` - `2`, diff01 = `0`-`1`,diff02 = `0`-`2`)

