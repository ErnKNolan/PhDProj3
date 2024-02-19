makeDecision <- function(properties = properties, interim_res = interim_res, j=j,drop_cut=drop_cut,stop_cut=stop_cut,ties=ties){
  
  #Determining what (if any) treatment group to drop
  int_drop <- interim_res %>% 
    filter(grepl("pp_|pred_prob",variable)) %>% 
    dplyr::select(variable,mean) %>%
    pivot_wider(names_from = variable, values_from = mean) %>%
    #I want the columns in a specific order to keep the old code the same
    dplyr::select(!contains("[1]")) %>% #get rid of control pred prob
    dplyr::select(contains(c("pp_","pred_prob_")))
  
  #finding the smallest pred prob
  int_drop$drop <- apply(int_drop[,grepl("pp_",names(int_drop))], 1, min, na.rm = TRUE)
  #determine if stopping rule met
  int_drop$stop <- ifelse(all(int_drop[,grepl("pp_",names(int_drop))] < stop_cut), 1, 0)
  
  #  int_drop$stop <- ifelse(int_drop$pp_trt2 < stop_cut & int_drop$pp_trt3 < stop_cut & int_drop$pp_trt4 < stop_cut, 1, 0)
  
  #dropping only if pp less than drop_cut
  #Optimisation option - if tie then drop higher constraint
  
  if(ties == "ties_opt"){
    int_drop <- int_drop %>%
      mutate(drop = ifelse(drop > drop_cut, NA, drop),
             droptrt = ifelse(!is.na(drop), 
                              tail(names(int_drop[,grepl("pp_",names(int_drop))])[which(int_drop[grepl("pp_",names(int_drop))]==int_drop$drop,arr.ind=T)[,"col"]],n=1),
                              "none"))
  } else{
    
    int_drop <- int_drop %>%
      mutate(drop = ifelse(drop > drop_cut, NA, drop),
             total = rowSums(int_drop[grep("pp_", names(int_drop))] == drop,na.rm=TRUE),
             droptrt = ifelse(total <= 1, 
                              ifelse(!is.na(drop), 
                                     names(int_drop[,grepl("pp_",names(int_drop))])[which(int_drop[grepl("pp_",names(int_drop))]==int_drop$drop,arr.ind=T)[,"col"]],
                                     "none"),
                              names(int_drop[,grepl("pp_",names(int_drop))])[which(int_drop[grepl("pred_prob",names(int_drop))]==min(int_drop[,grepl("pred_prob",names(int_drop))]),arr.ind=T)[,"col"]])) %>%
      #get rid of the pp_ section if it has it
      mutate(droptrt = gsub(".*_","",droptrt))
  }
  
  properties <- properties[j,]
  properties$drop <- int_drop$droptrt
  properties$stop <- int_drop$stop
  stop <- int_drop$stop

  
  return(properties)
  
}