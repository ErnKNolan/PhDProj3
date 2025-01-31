
#run the non-adaptive trials 

nonAdaptTrial <- function(properties = properties, mod = mod, outdir=outdir, j=j,t=t,draws=draws) {
  
  fullclusters <- makeClusters(t=4,nid=properties$n_per_k[j],
                               t1=properties$k[j],
                               t2=properties$k[j],
                               t3=properties$k[j],
                               t4=properties$k[j])
  
  full <- testInterim(expdat=fullclusters,t=4,mod=mod,outdir=outdir,
                      rho=properties$icc[j],t1=properties$t1[j],t2=properties$t2[j],t3=properties$t3[j],t4=properties$t4[j],draws=draws)
  
  res <- list(full=full[[1]])
  return(res)
}
