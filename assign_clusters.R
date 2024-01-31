
nblock <- 4
k <- 20
t <- 5
mat <- matrix(nrow = t-1,ncol=nblock)
drops <- matrix(nrow=t-1,ncol=nblock)
drops[,1] <- 1
#first interim
drops[,2] <- c(1,1,1,1)
#second interim
drops[,3] <- c(1,1,0,1)
#third interim
drops[,4] <- c(1,0,0,1)

int_split <- function(n, p) n %/% p + (sequence(p) - 1 < n %% p)

for(i in 1:nblock){
if(i==1){
  mat[,i] <- floor(k/nblock)
} else if(i != 1){
  leftovers <- k - (floor(k/nblock)*nblock) #this is for leftover clusters when using the floor method for each block
                                            #e.g. 35 clusters over 3 blocks will have 2 leftovers per group
  for(j in 1:t-1){
    mat[j,i] <- ifelse(drops[j,i]==0,0, #dont assign clusters to dropped arms
                       floor(floor(k/nblock)*(t-1)/sum(drops[,i]))) #clusters for that block divided out to remaining groups
                                        #assign one left over to each group starting from the last block until all used up
  }
  diff <- floor(k/nblock)*(t-1) - sum(mat[,i]) #need to divy out the remaining clusters from the floor(floor())
  if(diff > 0){
    if(sum(drops[,i]) > 1){ #if more than one arm left, then divy out the remaining
      diffs <- sample(c(rep(1,times=diff),rep(0,times=sum(drops[,i])-diff)),sum(drops[,i]))
      drops2 <- drops[,i]
      drops2[drops2 != 0] <- diffs
      mat[,i] <- mat[,i] + drops2
    } else if(sum(drops[,i])==1){ #if only 1 arm, give it all the remaining clusters
      drops2 <- drops[,i]
      drops2[drops2 != 0] <- diffs
      mat[,i] <- mat[,i] + drops2
    }
  }
  if(sum(drops[,i]) > 1){ #if an arm is dropped then divy out its leftovers to the other arms
    #if its time to start assigning left over clusters to each group starting from last block, then take the number of treatment groups divided
    #by the number of treatment groups still included
    l <- (i >= nblock-(leftovers-1))
    s <- sum(drops[,i])
    add <- sample(int_split(l*(t-1),s))#take away the sample() if we want it non-random
    drops2 <- drops[,i]
    drops2[drops2 != 0] <- add
    mat[,i] <- mat[,i] + drops2
  } else if(sum(drops[,i])==1){ #if only 1 arm then give it all the leftovers
    drops2 <- drops[,i]
    drops2[drops2 !=0] <- (l*(t-1))
    mat[,i] <- mat[,i] + drops2
  } else{# I dont think this should occur, if sum drops = 0 then there's no arms
    mat[,i] <- mat[,i]
    print("Warning! No arms assigned")
  }

}
}



