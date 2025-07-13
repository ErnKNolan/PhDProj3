#AUTHOR: Erin Nolan
#TITLE: PROJECT 3 ASSIGN CLUSTERS
#PURPOSE: ASSIGN CLUSTERS BASED ON WHATS DROPPED

#assign out integers
int_split <- function(n, p) n %/% p + (sequence(p) - 1 < n %% p)
#base r sample function if X is variable
resample <- function(x, ...) x[sample.int(length(x), ...)]

assignClusterLate <- function(mat = mat, k = k, nblock = nblock, drops = drops, i=i){
  #matrix if no drops occur
  fullmat <- matrix(nrow = t-1,ncol=2) #forcing this to be 2 columns
  #assign clusters to each arm as if no drops
  for(g in 1:t-1){
    fullmat[g,] <- rev(int_split(k,nblock))
  } 
  #assign out number of clusters to remaining arm
  diffs <- resample(int_split(sum(fullmat[,i]),sum(drops[,i])))
  drops2 <- drops[,i]
  #assign only to arms that arent dropped
  drops2[drops2 != 0] <- diffs
  mat[,i] <- mat[,i] + drops2
  
  return(mat[,i])
}
