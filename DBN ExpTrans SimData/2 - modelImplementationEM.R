

# Read in the obtained tables/files/lists
  load("candidateTables.Rdata"); load("observedAP.Rdata"); load("observedRS.Rdata")

# Initialize obsProbs, sigMeans, sigSDs
  p <- matrix(colMeans(observedAP),ncol(candidateTables),ncol(observedAP),byrow = T)
  mu <- matrix(colMeans(observedRS),ncol(candidateTables),ncol(observedAP),byrow = T)
  sigma <- matrix(apply(observedRS,2,sd),ncol(candidateTables),ncol(observedAP),byrow = T)

# Initialize weight parameters for features (proper starting values!)
  nu <- rep(1,ncol(candidateTables))
  lambda <- matrix(-2,ncol(candidateTables),ncol(candidateTables))
  diag(lambda) <- 0
  
# Initialize membership probabilities matrices
  tableP <- matrix(0,nrow(candidateTables),nrow(p))
  tablesP <- vector("list",nrow(candidateTables))
  for (j in 1:length(tablesP)){
    tablesP[[j]]<-matrix(0,nrow(p),nrow(p))
  }

# Initialize a dummy vector of length ncol(bestSignals)
  dummy <- rep(0,ncol(observedAP))
  
  
#***********  
# FUNCTIONS 
#***********
  
# Observation functions
  b <- function(t,w){
    idx <- as.logical(observedAP[t,])
    if( sum(idx)==0 ) return(prod(1 - p[w,]))
    dummy[!idx] <- (1 - p[w,!idx])
    dummy[idx] <- p[w,idx] * dnorm( (as.numeric(observedRS[t,idx]) - mu[w,idx]) / sigma[w,idx] ) / sigma[w,idx]
    return(prod(dummy))
  }
  bVector <- function(t){
    probs <- sapply(1:nrow(p), function(x) b(t,x))
    return(probs)
  }
  
# Entropy Classifier Transition Matrix (proportional!)
  transitionM <- function(t){
    return(exp(lambda+ t(t(nu)) %*% candidateTables[t,] ))
  }

# Functions to extract MLE estimators of parameters for each duplet, given membership probabilities 
  getEstimators <- function(k,j){
    idx <- which(observedAP[,j] == 1)
    memberSum <- sum(tableP[,k])
    memberObsSum <- sum(tableP[idx,k])
    memberObsSignal <- sum(tableP[idx,k]*observedRS[idx,j])
    obsProb <- memberObsSum / memberSum
    if (memberObsSum==0) {mean <- mu[k,j]} else {mean <- memberObsSignal / memberObsSum}
    idx <- which((observedAP[,j] == 1))
    memberObsDev <- sum(tableP[idx,k]*(observedRS[idx,j]-mean)^2)
    if (memberObsSum==0) {stDev <- sigma[k,j]} else {stDev <- sqrt( memberObsDev / memberObsSum )}
    return(c(obsProb,mean,stDev))
  }
  tableEstimators <- function(k){
    return(t(sapply(1:ncol(p), function(x) getEstimators(k,x))))
  }  
  
# Function to extract MLE estimator of weights in transition model!
  scoreLambdasInner <- function(lambdasK,nuK,k,t){
    lambdasK[k] <- 0 # Always... thats the pivot in multinomial logit!
    eprops <- exp(lambdasK+nuK*candidateTables[t,])
    return(sum(tablesP[[t-1]][k,] * log( eprops/sum(eprops))))
  }
  sumTime <- function(lambdasK,nuK,k){
    # Here is when we ignore the change in shifts (that's no transition to 4th table)
    return(sum( parSapply(c1, (1:nrow(candidateTables))[-seq(1,nrow(candidateTables),120)]  , function(x) scoreLambdasInner(lambdasK,nuK,k,x) ) ))
  }
  optimLambdasNuK <- function(x,k) {
    lambdasK <- x[1:nrow(p)]
    nuK<-x[nrow(p)+1]
    return(-sumTime(lambdasK,nuK,k)) #Negative since optim minimizes!
  }

# Function to compute log-likelihood proportional bit (to assess whether convergence was reached)
  sumJ <- function(t,w){
    iepe <- rep(0,ncol(p))
    idx <- as.logical(observedAP[t,])
    iepe[!idx] <- (1-p[w,!idx])
    iepe[idx] <- (p[w,idx] * dnorm( (as.numeric(observedRS[t,idx]) - mu[w,idx]) / sigma[w,idx] ) / sigma[w,idx])
    iepe <- tableP[t,w] * sum( pmax(-100,log(iepe)) ) # We need to bound it... to be able to do parameter roundings later!
    return( iepe )
  }
  sumW <- function(t){
    iep <- sapply(1:nrow(p), function(x) sumJ(t,x))
    return(sum(iep))
  }

# Full pseudo log-likelihood accounting for membership probabilities!
  logLik <- function(){
    iep <- parSapply(c1,1:nrow(candidateTables), sumW)
    ueee <- 0
    for (i in 1:nrow(p)) ueee <- ueee + sumTime(lambda[i,],nu[i],i)
    return(sum(iep)+ueee)
  }


#****************************************************
# EXPECTATION MAXIMIZATION FORWARD/BACKWARD PROCEDURE
#****************************************************

  library(snow); library(parallel)
  c1 <- makeCluster(7,type="SOCK") #This machine has got 8 cores, check parallel.detectCores()

# We declare the alpha and beta lists
  alpha <- matrix(0,nrow(candidateTables),ncol(tableP))
  beta <- matrix(0,nrow(candidateTables),ncol(tableP))
  
# We define functions for the forward/backward passes; so that we can easily paralellize
  forwardPass <- function(){
    alpha[1,4] <- 1 # This is because in first entry it must be in 4th table!
    for (t in 2:nrow(candidateTables)){
      if ( (t-1)%%120 ==0 ){
        alpha[t,4] <- 1 # Recall that if we know its in 4th Table
        # Probability of observation being made from anywhere else is 0!
      } else {
        alpha[t,] <- bVector(t) * alpha[t-1,] %*% transitionM(t)
        alpha[t,] <- alpha[t,] / sum(alpha[t,])
      }
    }
    return(alpha)
  } 
  backwardPass <- function(){
    beta[nrow(beta),] <- 1/nrow(p) # We just keep info up to proportionality on ward likelihoods
    for (t in seq(nrow(beta)-1,1,-1) ){
      if ( (t)%%120 ==0 ){
        beta[t,] <- rep(1/6,6)  # RESET to 4th table start, uninformative
      } else {
        beta[t,] <- transitionM(t+1) %*% (bVector(t+1) * beta[t+1,])
        beta[t,] <- beta[t,] / sum(beta[t,])
      }
    }
    return(beta)
  } 

# We also define a function for the computations of dual memberships, to also be able to parallelize
  dualMemberships <- function(){
    for(t in 2:length(tablesP)){
      if ( (t-1)%%120 ==0 ){
        tablesP[[t-1]][,] <- -Inf # I'll put something crazy, just to make sure I'm ignoring it later (check if crash)
        #tablesP[[t-1]][,4] <- alpha[t-1,] * transitionM(t)[,4]
        #tablesP[[t-1]] <- tablesP[[t-1]] / sum(tablesP[[t-1]]) # Can do after
        #In this case it can only go to coord room straight in next transition!
      } else {
        tablesP[[t-1]] <- t(t(alpha[t-1,])) %*% (beta[t,] * bVector(t)) 
        tablesP[[t-1]] <- tablesP[[t-1]] * transitionM(t) # Can do after
        tablesP[[t-1]] <- tablesP[[t-1]] / sum(tablesP[[t-1]]) # Can do after
      }
    }
  return(tablesP)
  }
  
  
#****************
# FIRST ITERATION
#****************
  
# The first iteration is simpler on backward pass so it goes on its own... (its uninformative)
  alpha <- forwardPass()
  beta <- matrix(1/6,nrow(candidateTables),nrow(p))
  # Memberships
    tableP <- alpha * beta
    tableP <- tableP / rowSums(tableP)
  # Dual Memberships
    tablesP <- dualMemberships()

  # Export info to Cluster
    clusterExport(c1,list("table","sumW","sumJ","sumTime","p","observedAP","observedRS","candidateTables",
                          "nu","tablesP","tableP","alpha","beta","forwardPass","dualMemberships",
                          "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  

  # Update p,mu,sigma,lambdas and nus
    iepe <- optim(c(lambda[1,],nu[1]),function(x) optimLambdasNuK(x,1),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[1,]<-iepe[1:6]; nu[1]<-iepe[7]
    iepe <- optim(c(lambda[2,],nu[2]),function(x) optimLambdasNuK(x,2),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[2,]<-iepe[1:6]; nu[2]<-iepe[7]
    iepe <- optim(c(lambda[3,],nu[3]),function(x) optimLambdasNuK(x,3),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[3,]<-iepe[1:6]; nu[3]<-iepe[7]
    iepe <- optim(c(lambda[4,],nu[4]),function(x) optimLambdasNuK(x,4),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[4,]<-iepe[1:6]; nu[4]<-iepe[7]
    iepe <- optim(c(lambda[5,],nu[5]),function(x) optimLambdasNuK(x,5),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[5,]<-iepe[1:6]; nu[5]<-iepe[7]
    iepe <- optim(c(lambda[6,],nu[6]),function(x) optimLambdasNuK(x,6),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[6,]<-iepe[1:6]; nu[6]<-iepe[7]
      #NO IDEA WHY I COULD NOT FORCE A LOOP WITHOUT HAVING TO EXPORT TO CLUSTER THE INDEX   
    rm(iepe);
    for (i in 1:nrow(p)) {
      dummy2 <- round(tableEstimators(i),10) # round in order to save memory
      p[i,] <- pmin(1,pmax(0,dummy2[,1])) # allow extreme probs or "other" ward will not get built up!
      mu[i,] <- dummy2[,2]
      sigma[i,] <- pmax(2,dummy2[,3]) # No way we allow the standard deviation to be below 2
    }

# Update Cluster, for likelihood evaluation
  clusterExport(c1,list("table","sumW","sumJ","sumTime","p","observedAP","observedRS","candidateTables",
                        "nu","tablesP","tableP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  
  print(logLik()); #-262307.7
    

#**********************
# ADDITIONAL ITERATIONS
#**********************
    
  for(iter in 1:6){
  
  ### FORWARD ###
  alpha <- forwardPass()
  #plot(colMeans(alpha),type="h"); 

  ### BACKWARD ###
  beta <- backwardPass()
  #plot(colMeans(beta),type="h"); 
  
  ### MEMBERSHIP PROBABILITIES ###
    tableP <- alpha * beta
    tableP <- tableP / rowSums(tableP)
    #tableP[779:781,];  #should be clusters 2 12 3 respectively
    #summary(apply(tableP,1,max)); # EXTREME PROBABILITIES ALREADY IN ITERATION 2!

  ### DUAL MEMBERSHIP PROBABILITIES ###
  tablesP <- dualMemberships()
  #tablesP[[779]]; 

  ### UPDATE P,MU,SIGMA,LAMBDAs and NU ###
  clusterExport(c1,list("table","sumW","sumJ","sumTime","p","observedAP","observedRS","candidateTables",
                        "nu","tablesP","tableP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  

  iepe <- optim(c(lambda[1,],nu[1]),function(x) optimLambdasNuK(x,1),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[1,]<-iepe[1:6]; nu[1]<-iepe[7]
  iepe <- optim(c(lambda[2,],nu[2]),function(x) optimLambdasNuK(x,2),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[2,]<-iepe[1:6]; nu[2]<-iepe[7]
  iepe <- optim(c(lambda[3,],nu[3]),function(x) optimLambdasNuK(x,3),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[3,]<-iepe[1:6]; nu[3]<-iepe[7]
  iepe <- optim(c(lambda[4,],nu[4]),function(x) optimLambdasNuK(x,4),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[4,]<-iepe[1:6]; nu[4]<-iepe[7]
  iepe <- optim(c(lambda[5,],nu[5]),function(x) optimLambdasNuK(x,5),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[5,]<-iepe[1:6]; nu[5]<-iepe[7]
  iepe <- optim(c(lambda[6,],nu[6]),function(x) optimLambdasNuK(x,6),method="L-BFGS-B",lower=c(rep(-Inf,6),0),upper=c(rep(0,6),Inf))$par; lambda[6,]<-iepe[1:6]; nu[6]<-iepe[7]
  #lambda; nu; transitionM(779)
  rm(iepe);
  for (i in 1:nrow(p)) {
    dummy2 <- round(tableEstimators(i),10) # round in order to save memory
    p[i,] <- pmin(1,pmax(0,dummy2[,1])) # allow extreme probs or "other" ward will not get built up!
    mu[i,] <- dummy2[,2]
    sigma[i,] <- pmax(1,dummy2[,3]) # No way we allow the standard deviation to be below 1
  }

  ### EVALUATE LIKELIHOOD ###
  clusterExport(c1,list("table","sumW","sumJ","sumTime","p","observedAP","observedRS","candidateTables",
                        "nu","tablesP","tableP","alpha","beta","forwardPass","dualMemberships",
                        "mu","sigma","lambda","b","bVector","dummy","backwardPass","transitionM","scoreLambdasInner"))  
  print(logLik()); #-262307.7,-196329.6,-170750.5, -168493.9, -168096.9, -168046.1 , -168036.5
  
  }
  
  stopCluster(c1); rm(c1)

  
  
  
  
  
  save(lambda, file = "lambdas.Rdata")
  save(nu, file = "nus.Rdata")  
  save(mu, file = "obsMeans.Rdata")
  save(sigma, file = "obsSDs.Rdata")
  save(p, file = "obsProbs.Rdata")
  save(tableP,file="probabilityTables.Rdata")
 