

# Read in the obtained tables/files/lists
  load("listObservedSignals.Rdata"); load("listObservedHotSpots.Rdata"); load("listCandidateWards.Rdata")
  load("listTimeAndPhones.Rdata")

# Initialize obsProbs, sigMeans, sigSDs
  dummy <- rowSums(sapply(1:8,function(x) colSums(listObservedHotSpots[[x]])))/sum(sapply(1:8,function(x) nrow(listObservedHotSpots[[x]])))
  p <- matrix(dummy,ncol(listCandidateWards[[2]]),ncol(listObservedHotSpots[[1]]),byrow = T)
  dummy <- rowMeans(sapply(1:8,function(x) colMeans(listObservedSignals[[x]],na.rm=T)),na.rm = T)
  mu <- matrix(dummy,ncol(listCandidateWards[[2]]),ncol(listObservedHotSpots[[1]]),byrow = T)
  sigma <- matrix(15,ncol(listCandidateWards[[2]]),ncol(listObservedHotSpots[[1]]))
  rm(dummy)
  
# Initialize weight parameters for features (proper starting values!)
  nu <- rep(1,ncol(listCandidateWards[[2]]))
  lambda <- matrix(-2,ncol(listCandidateWards[[2]]),ncol(listCandidateWards[[2]]))
  diag(lambda) <- 0
  
# Initialize membership probabilities matrices
  wardP <- vector("list",length(listTimeAndPhones))
  for (i in 1:length(wardP)){  
    wardP[[i]] <- matrix(0,nrow(listCandidateWards[[i]]),nrow(p))
    print(ncol(listCandidateWards[[i]]))
  }
  wardsP <- vector("list",length(listTimeAndPhones))
  for (i in 1:length(wardP)){  
    wardsP[[i]] <- vector("list",nrow(listCandidateWards[[i]]))
    for (j in 1:length(wardsP[[i]])){
      wardsP[[i]][[j]]<-matrix(0,nrow(p),nrow(p))
    }
  }  
  
# Need to avoid NA entries in best signals... turn them in -1000 (will be ignored in algorithm)
  for (i in 1:length(listObservedSignals)){
    listObservedSignals[[i]][is.na(listObservedSignals[[i]])] <- -1000
  }; rm(i,j)
  
# Initialize a dummy vector of length ncol(bestSignals)
  dummy <- rep(0,ncol(listObservedSignals[[1]]))
  
  
#***********  
# FUNCTIONS 
#***********
  
# Observation functions
  b <- function(phone,t,w){
    idx <- as.logical(listObservedHotSpots[[phone]][t,])
    if( sum(idx)==0 ) return(prod(1 - p[w,]))
    dummy[!idx] <- (1 - p[w,!idx])
    dummy[idx] <- p[w,idx] * dnorm( (as.numeric(listObservedSignals[[phone]][t,idx]) - mu[w,idx]) / sigma[w,idx] ) / sigma[w,idx]
    return(prod(dummy))
  }
  bVector <- function(phone,t){
    probs <- sapply(1:nrow(p), function(x) b(phone,t,x))
    return(probs)
  }
  
# Entropy Classifier Transition Matrix (proportional!), this is gonna be sloooooow
  transitionM <- function(t,phone){
    return(exp(lambda+ t(t(nu)) %*% listCandidateWards[[phone]][t,] ))
  }

# Functions to extract MLE estimators of parameters for each duplet, given membership probabilities 
  getEstimators <- function(k,j){
    memberSum <- 0; memberObsSum <- 0; memberObsSignal <- 0; memberObsDev <- 0
    for (phoneDummy in 1:length(listObservedHotSpots)){
      idx <- which(listObservedHotSpots[[phoneDummy]][,j] == 1)
      memberSum <- memberSum + sum(wardP[[phoneDummy]][,k])
      memberObsSum <- memberObsSum + sum(wardP[[phoneDummy]][idx,k])
      memberObsSignal <- memberObsSignal + sum(wardP[[phoneDummy]][idx,k]*listObservedSignals[[phoneDummy]][idx,j])
    }
    obsProb <- memberObsSum / memberSum
    if (memberObsSum==0) {mean <- mu[k,j]} else {mean <- memberObsSignal / memberObsSum}
    for (phoneDummy in 1:length(listObservedHotSpots)){
      idx <- which((listObservedHotSpots[[phoneDummy]][,j] == 1))
      memberObsDev <- memberObsDev + sum(wardP[[phoneDummy]][idx,k]*(listObservedSignals[[phoneDummy]][idx,j]-mean)^2)
    } 
    if (memberObsSum==0) {stDev <- sigma[k,j]} else {stDev <- sqrt( memberObsDev / memberObsSum )}
    return(c(obsProb,mean,stDev))
  }
  wardEstimators <- function(k){
    return(t(sapply(1:ncol(p), function(x) getEstimators(k,x))))
  }  
  
# Function to extract MLE estimator of weights in transition model!
  scoreLambdasInner <- function(lambdasK,nuK,k,phone,t){
    lambdasK[k] <- 0 # Always... thats the pivot in multinomial logit!
    eprops <- exp(lambdasK+nuK*listCandidateWards[[phone]][t,])
    return(sum(wardsP[[phone]][[t-1]][k,] * log( eprops/sum(eprops))))
  }
  sumTime <- function(lambdasK,nuK,k,phone){
    # Here is when we ignore those entries where the phone was not active!
    return(sum( sapply( which(listTimeAndPhones[[phone]][,4]=="yes")[-1], function(x) scoreLambdasInner(lambdasK,nuK,k,phone,x) ) ))
  }
  sumPhone <- function(lambdasK,nuK,k){
    return(sum( parSapply(c1, 1:length(listTimeAndPhones), function(x) sumTime(lambdasK,nuK,k,x) ) ))
  }
  optimLambdasNuK <- function(x,k) {
    lambdasK <- x[1:nrow(p)]
    nuK<-x[nrow(p)+1]
    return(-sumPhone(lambdasK,nuK,k)) #Negative since optim minimizes!
  }

# Function to compute log-likelihood proportional bit (to assess whether convergence was reached)
  sumJ <- function(phone,t,w){
    iepe <- rep(0,ncol(p))
    idx <- as.logical(listObservedHotSpots[[phone]][t,])
    iepe[!idx] <- (1-p[w,!idx])
    iepe[idx] <- (p[w,idx] * dnorm( (as.numeric(listObservedSignals[[phone]][t,idx]) - mu[w,idx]) / sigma[w,idx] ) / sigma[w,idx])
    iepe <- wardP[[phone]][t,w] * sum( pmax(-100,log(iepe)) ) # We need to bound it... to be able to do parameter roundings later!
    return( iepe )
  }
  sumW <- function(phone,t){
    iep <- sapply(1:nrow(p), function(x) sumJ(phone,t,x))
    return(sum(iep))
  }
  sumT <- function(phone){
    iep <- sapply(1:nrow(wardP[[phone]]), function(x) sumW(phone,x))
    return(sum(iep))    
  }
  
# Full pseudo log-likelihood accounting for membership probabilities!
  logLik <- function(){
    iep <- parSapply(c1,1:length(wardP), sumT)
    ueee <- 0
    for (i in 1:nrow(p)) ueee <- ueee + sumPhone(lambda[i,],nu[i],i)
    return(sum(iep)+ueee)
  }


#****************************************************
# EXPECTATION MAXIMIZATION FORWARD/BACKWARD PROCEDURE
#****************************************************

  library(snow); library(parallel)
  c1 <- makeCluster(7,type="SOCK") #This machine has got 8 cores, check parallel.detectCores()

# We declare the alpha and beta lists
  alpha <- vector(mode="list",length=8) # Keep things just up to proportionality on w
  beta <- vector(mode="list",length=8) # Keep things just up to proportionality on w
  for (phone in 1:8) alpha[[phone]] <- matrix(0,nrow(listCandidateWards[[phone]]),ncol(wardP[[1]]))
  for (phone in 1:8) beta[[phone]] <- matrix(0,nrow(listCandidateWards[[phone]]),ncol(wardP[[1]]))
  
# We define functions for the forward/backward passes; so that we can easily paralellize
  forwardPass <- function(phone){
    alpha[[phone]][1,9] <- 1 # This is because in first entry it must be in Coord Room
    logVectorActPhone <- listTimeAndPhones[[phone]][,4]=="yes"
    for (t in 2:nrow(listCandidateWards[[phone]])){
      if (!logVectorActPhone[t]){
        alpha[[phone]][t,9] <- 1 # Recall that if we know its in coord room we need to account for it!
        # Probability of observation being made from anywhere else than coord room is 0!
      } else {
        alpha[[phone]][t,] <- bVector(phone,t) * alpha[[phone]][t-1,] %*% transitionM(t,phone)
        alpha[[phone]][t,] <- alpha[[phone]][t,] / sum(alpha[[phone]][t,])
      }
    }
    return(alpha[[phone]])
  } 
  backwardPass <- function(phone){
    beta[[phone]][nrow(beta[[phone]]),] <- 1/nrow(p) # We just keep info up to proportionality on ward likelihoods
    logVectorActPhone <- listTimeAndPhones[[phone]][,4]=="yes"
    for (t in seq(nrow(beta[[phone]])-1,1,-1) ){
      if (!logVectorActPhone[t+1]){
        beta[[phone]][t,] <- transitionM(t+1,phone)[,9]  # If at next time we know it must be in coord room, then this is proportional to transitions
        # to that cluster
        beta[[phone]][t,] <- beta[[phone]][t,] / sum(beta[[phone]][t,])
      } else {
        beta[[phone]][t,] <- transitionM(t+1,phone) %*% (bVector(phone,t+1) * beta[[phone]][t+1,])
        beta[[phone]][t,] <- beta[[phone]][t,] / sum(beta[[phone]][t,])
      }
    }
    return(beta[[phone]])
  } 

# We also define a function for the computations of dual memberships, to also be able to parallelize
  dualMemberships <- function(phone){
    logVectorActPhone <- listTimeAndPhones[[phone]][,4]=="yes"
    for(t in 2:length(wardsP[[phone]])){
      if (!logVectorActPhone[t-1] && !logVectorActPhone[t]){
        wardsP[[phone]][[t-1]][,] <- 0
        wardsP[[phone]][[t-1]][9,9] <- 1
      } else if (logVectorActPhone[t-1] && !logVectorActPhone[t]){
        wardsP[[phone]][[t-1]][,] <- 0
        wardsP[[phone]][[t-1]][,9] <- alpha[[phone]][t-1,] * transitionM(t,phone)[,9]
        wardsP[[phone]][[t-1]] <- wardsP[[phone]][[t-1]] / sum(wardsP[[phone]][[t-1]]) # Can do after
        #In this case it can only go to coord room straight in next transition!
      } else {
        wardsP[[phone]][[t-1]] <- t(t(alpha[[phone]][t-1,])) %*% (beta[[phone]][t,] * bVector(phone,t)) 
        wardsP[[phone]][[t-1]] <- wardsP[[phone]][[t-1]] * transitionM(t,phone) # Can do after
        wardsP[[phone]][[t-1]] <- wardsP[[phone]][[t-1]] / sum(wardsP[[phone]][[t-1]]) # Can do after
      }
    }
  return(wardsP[[phone]])
  }
  
  
#****************
# FIRST ITERATION
#****************
  
# Export info to Cluster
  clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                        "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM"))  

# The first iteration is simpler on backward pass so it goes on its own... (its uninformative)
  alpha <- parSapply(c1,1:length(wardP), forwardPass)
  for (phone in 1:8) beta[[phone]] <- matrix(1/nrow(p),nrow(listCandidateWards[[phone]]),nrow(p))
  # Memberships
    for (phone in 1:8) {
      wardP[[phone]] <- alpha[[phone]] * beta[[phone]]
      wardP[[phone]] <- wardP[[phone]] / rowSums(wardP[[phone]])
    }
  # Update Cluster, don't bother... just plug everything!
    clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                          "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                          "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  
  # Dual Memberships
    wardsP <- parSapply(c1,1:length(wardP), dualMemberships) 
  # Update Cluster, don't bother... just plug everything!
    clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                          "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                          "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  
  # Update p,mu,sigma,lambdas and nus
    iepe <- optim(c(lambda[1,],nu[1]),function(x) optimLambdasNuK(x,1),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[1,]<-iepe[1:9]; nu[1]<-iepe[10]
    iepe <- optim(c(lambda[2,],nu[2]),function(x) optimLambdasNuK(x,2),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[2,]<-iepe[1:9]; nu[2]<-iepe[10]
    iepe <- optim(c(lambda[3,],nu[3]),function(x) optimLambdasNuK(x,3),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[3,]<-iepe[1:9]; nu[3]<-iepe[10]
    iepe <- optim(c(lambda[4,],nu[4]),function(x) optimLambdasNuK(x,4),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[4,]<-iepe[1:9]; nu[4]<-iepe[10]
    iepe <- optim(c(lambda[5,],nu[5]),function(x) optimLambdasNuK(x,5),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[5,]<-iepe[1:9]; nu[5]<-iepe[10]
    iepe <- optim(c(lambda[6,],nu[6]),function(x) optimLambdasNuK(x,6),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[6,]<-iepe[1:9]; nu[6]<-iepe[10]
    iepe <- optim(c(lambda[7,],nu[7]),function(x) optimLambdasNuK(x,7),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[7,]<-iepe[1:9]; nu[7]<-iepe[10]
    iepe <- optim(c(lambda[8,],nu[8]),function(x) optimLambdasNuK(x,8),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[8,]<-iepe[1:9]; nu[8]<-iepe[10]
    iepe <- optim(c(lambda[9,],nu[9]),function(x) optimLambdasNuK(x,9),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[9,]<-iepe[1:9]; nu[9]<-iepe[10]
    #NO IDEA WHY I COULD NOT FORCE A LOOP WITHOUT HAVING TO EXPORT TO CLUSTER THE INDEX   
    rm(iepe);
    for (i in 1:nrow(p)) {
      dummy2 <- round(wardEstimators(i),10) # round in order to save memory
      p[i,] <- pmin(1,pmax(0,dummy2[,1])) # allow extreme probs or "other" ward will not get built up!
      mu[i,] <- dummy2[,2]
      sigma[i,] <- pmax(2,dummy2[,3]) # No way we allow the standard deviation to be below 2
    }

# Update Cluster, for likelihood evaluation
  clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                        "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))   
  print(logLik()); #-1711098

  
#**********************
# ADDITIONAL ITERATIONS
#**********************
    
  for(iter in 1:7){
  
  ### FORWARD ###
  alpha <- parSapply(c1,1:length(wardP), forwardPass) 
  #for(i in 1:8) {plot(colMeans(alpha[[i]][,1:6]),ylim=c(0,0.25),type="h")}; alpha[[1]][600,]

  ### BACKWARD ###
  beta <- parSapply(c1,1:length(wardP),backwardPass)
  #for(i in 1:8) {plot(colMeans(beta[[i]][,1:6]),ylim=c(0,0.25),type="h")}; beta[[1]][600,]
  
  ### MEMBERSHIP PROBABILITIES ###
  for (phone in 1:8) {
    wardP[[phone]] <- alpha[[phone]] * beta[[phone]]
    wardP[[phone]] <- wardP[[phone]] / rowSums(wardP[[phone]])
  }
  #wardP[[1]][600,]; wardP[[1]][601,]; wardP[[1]][602,];
  #summary(apply(wardP[[1]][wardP[[1]][,7]<1,],1,max)); # EXTREME PROBABILITIES ALREADY IN ITERATION 2!

  ### DUAL MEMBERSHIP PROBABILITIES ###
  clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                        "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))   
  wardsP <- parSapply(c1,1:length(wardP), dualMemberships) 
  #wardsP[[1]][[600]]

  ### UPDATE P,MU,SIGMA,LAMBDAs and NU ###
  clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                        "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  
  iepe <- optim(c(lambda[1,],nu[1]),function(x) optimLambdasNuK(x,1),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[1,]<-iepe[1:9]; nu[1]<-iepe[10]
  iepe <- optim(c(lambda[2,],nu[2]),function(x) optimLambdasNuK(x,2),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[2,]<-iepe[1:9]; nu[2]<-iepe[10]
  iepe <- optim(c(lambda[3,],nu[3]),function(x) optimLambdasNuK(x,3),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[3,]<-iepe[1:9]; nu[3]<-iepe[10]
  iepe <- optim(c(lambda[4,],nu[4]),function(x) optimLambdasNuK(x,4),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[4,]<-iepe[1:9]; nu[4]<-iepe[10]
  iepe <- optim(c(lambda[5,],nu[5]),function(x) optimLambdasNuK(x,5),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[5,]<-iepe[1:9]; nu[5]<-iepe[10]
  iepe <- optim(c(lambda[6,],nu[6]),function(x) optimLambdasNuK(x,6),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[6,]<-iepe[1:9]; nu[6]<-iepe[10]
  iepe <- optim(c(lambda[7,],nu[7]),function(x) optimLambdasNuK(x,7),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[7,]<-iepe[1:9]; nu[7]<-iepe[10]
  iepe <- optim(c(lambda[8,],nu[8]),function(x) optimLambdasNuK(x,8),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[8,]<-iepe[1:9]; nu[8]<-iepe[10]
  iepe <- optim(c(lambda[9,],nu[9]),function(x) optimLambdasNuK(x,9),method="L-BFGS-B",lower=c(rep(-Inf,9),0),upper=c(rep(0,9),Inf))$par; lambda[9,]<-iepe[1:9]; nu[9]<-iepe[10]
  #lambda; nu; transitionM(600,1)
  rm(iepe);
  for (i in 1:nrow(p)) {
    dummy2 <- round(wardEstimators(i),10) # round in order to save memory
    p[i,] <- pmin(1,pmax(0,dummy2[,1])) # allow extreme probs or "other" ward will not get built up!
    mu[i,] <- dummy2[,2]
    sigma[i,] <- pmax(2,dummy2[,3]) # No way we allow the standard deviation to be below 2
  }
  #for (i in 1:7) plot(p[i,],type="h",main=colnames(listCandidateWards[[1]])[i],ylim=c(0,1))

  ### EVALUATE LIKELIHOOD ###
  clusterExport(c1,list("wardP","sumW","sumT","sumJ","sumTime","sumPhone","p","listObservedHotSpots","listObservedSignals","listCandidateWards",
                        "listTimeAndPhones","nu","wardsP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner")) 
  #print(logLik()); #-1711098, ... ,-2530599
  print(iter)
  }
  
  stopCluster(c1); rm(c1)

  
  
  save(lambda, file = "outputFiles/lambdas.Rdata")
  save(nu, file = "outputFiles/nus.Rdata")  
  save(mu, file = "outputFiles/obsMeans.Rdata")
  save(sigma, file = "outputFiles/obsSDs.Rdata")
  save(p, file = "outputFiles/obsProbs.Rdata")
  save(wardP,file="outputFiles/probabilityWards.Rdata")
  
  for (i in 1:9) plot(p[i,],type="h",main=colnames(listCandidateWards[[1]])[i],ylim=c(0,1))
  