#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 RUN ADAPTIVE TRIAL
#PURPOSE: RUN INTERIMS AND FINAL MODEL AND RETURN RESULTS

runSimTrial <- function(properties = properties, mod = mod, outdir=outdir, j=j,t=t, adaption=adaption,drop_cut=drop_cut,stop_cut=stop_cut,nblock=nblock,draws=draws) {
  
  #The interim clusters are full clusters divided by block
  properties$kt2 <- floor(properties$k[j]/(properties$nblock[j]))
  properties$kt3 <- floor(properties$k[j]/(properties$nblock[j]))
  properties$kt4 <- floor(properties$k[j]/(properties$nblock[j]))
  
  #number of interims
  nint <- properties$nblock[j] - 1
  #making the cluster and drop matrices
  mat <- matrix(0,nrow = t-1,ncol=properties$nblock[j]-1)
  drops <- matrix(1,nrow=t-1,ncol=properties$nblock[j]-1)
  #interim values
  properties_int <- properties[j,]
  ints <- 0
  
  #run the interim X number of times
  for(i in 1:nint){
    
    #either the initial interim clusters or the updated value
    if(i==1){
      clusters <- data.frame(interim = properties_int$interim, kt2 = properties_int$kt2, kt3 = properties_int$kt3, kt4 = properties_int$kt4)
    } else {
      clusters <- data.frame(interim = properties_int$interim*i, kt2 = sum(mat[1,])+properties_int$interim, kt3 = sum(mat[2,])+properties_int$interim, kt4 = sum(mat[3,])+properties_int$interim)
    }
    #making clusters for the interim dataset
    intclusters <- makeClusters(t=4,nid=properties$n_per_k[j],
                                t1=clusters$interim,t2=clusters$kt2,
                                t3=clusters$kt3,t4=clusters$kt4)
    if(i == 1){
      #analysing the interim data
      interim <- testInterim(expdat=intclusters,t=4,mod=mod,outdir=outdir,
                             rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j],draws=draws)
    } else{
      #analysing the interim data but using this function to read in the already existing data as well
      interim <- testFull(expdat=intclusters,t=4,mod=mod,outdir=outdir,int_dat=resp_clus,
                          rho=properties_int$icc,t1=properties_int$t1,t2=properties_int$t2,t3=properties_int$t3,t4=properties_int$t4,draws=draws)
    }
    #save the interim results and sim
    resp <- unlist(interim[["resp"]]) 
    interim_res <- interim[[1]]
    
    #add the simulated results into the interim dataset
    resp_clus <- cbind(intclusters,resp)
    
    #making the decision and getting the clusters
    properties_int <- makeDecision(properties = properties, interim_res = interim_res, j=j,
                                   drop_cut = drop_cut, stop_cut = stop_cut, ties = ties, drops = drops,i=i)
    #adding in what arm to drop based on decision
    if(properties_int$drop != "none"){
      drops[as.numeric(gsub(".*?([0-9]+).*", "\\1", properties_int$drop))-1,i:(properties$nblock[j]-1)] <- 0
    } 
    if(properties_int$stop == 1){
      drops[,i:(properties$nblock[j]-1)] <- 0
    }
    #assigning clusters based on whats dropped
    mat[,i] <- assignCluster(mat=mat, k=properties$k[j]-floor(properties$k[j]/(properties$nblock[j])), nblock=properties$nblock[j]-1, drops=drops, i=i)
    #rewriting the properties to read into data generation
    for(k in 1:t-1){
      x <- paste0("kt",k+1)
      properties_int[,x] <- mat[k,i]
    }
    #if we have a stop, we leave the loop as no more interims
    if(properties_int$stop == 1){
      break
    }
    ints <- ints+1
  }
  #adding in the updated cluster assignment
  clusters <- data.frame(k = properties_int$k, kt2 = sum(mat[1,])+properties_int$interim, kt3 = sum(mat[2,])+properties_int$interim, kt4 = sum(mat[3,])+properties_int$interim)
  fullclusters <- makeClusters(t=4,nid=properties_int$n_per_k,
                               t1=clusters$k,
                               t2=clusters$kt2,
                               t3=clusters$kt3,
                               t4=clusters$kt4)
  #save how many interims were run
  properties_int$ints <- ints
  
  #need to append the interim data and new data together
  #then run the full trial model
  if(properties_int$stop==0){
    full <- testFull(expdat=fullclusters,t=4,mod=mod,outdir=outdir,int_dat=resp_clus,
                     rho=properties_int$icc,t1=properties_int$t1,t2=properties_int$t2,t3=properties_int$t3,t4=properties_int$t4,draws=draws)
    
    res <- list(interim=interim[[1]],full=full[[1]],properties=properties_int,drops=drops,mat=mat)
    return(res)
    #if there's a stop then the last interim run is the final
  } else if(properties_int$stop==1){
    res <- list(interim=interim[[1]],full=interim[[1]],properties=properties_int,drops=drops,mat=mat)
  }
}
