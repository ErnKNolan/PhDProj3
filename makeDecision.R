makeDecision <- function(properties = properties, interim_res = interim_res, j=j,drop_cut=drop_cut,stop_cut=stop_cut,ties=ties,drops=drops,i=i){
  
  #Determining what (if any) treatment group to drop
  int_drop <- interim_res %>% 
    filter(grepl("pp_|pred_prob",variable)) %>% 
    dplyr::select(variable,mean) %>%
    pivot_wider(names_from = variable, values_from = mean) %>%
    #I want the columns in a specific order to keep the old code the same
    dplyr::select(!contains("[1]")) %>% #get rid of control pred prob
    dplyr::select(contains(c("pp_","pred_prob_")))
  
  #finding the smallest pred prob on non-dropped arms
  x1 <- int_drop[,grepl("pp_",names(int_drop))]
  pred1 <- int_drop[grepl("pred_prob",names(int_drop))]
  
  if(sum(drops[,i])==length(drops[,i])){
    x1 <- x1
    pred1 <- pred1
  } else {
    x1 <- x1[,-which(drops[,i]==0)]
    pred1 <- pred1[,-which(drops[,i]==0)]
  }
  int_drop$drop <- apply(x1, 1, min, na.rm = TRUE)
  #determine if stopping rule met
  int_drop$stop <- ifelse(all(x1 < stop_cut), 1, 0)

  #dropping only if pp less than drop_cut
    int_drop <- int_drop %>%
      mutate(drop = ifelse(drop > drop_cut, NA, drop),
             total = rowSums(x1 == drop,na.rm=TRUE),
             droptrt = ifelse(total <= 1, 
                              ifelse(!is.na(drop), 
                                     names(x1)[which(x1==int_drop$drop,arr.ind=T)[,"col"]],
                                     "none"),
                              names(x1)[which(pred1==min(pred1),arr.ind=T)[,"col"]])) %>%
      #get rid of the pp_ section if it has it
      mutate(droptrt = gsub(".*_","",droptrt))
  
  properties <- properties[j,]
  properties$drop <- int_drop$droptrt
  properties$stop <- int_drop$stop
  stop <- int_drop$stop

  
  return(properties)
  
}