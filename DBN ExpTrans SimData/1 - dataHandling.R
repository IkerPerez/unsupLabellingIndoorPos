
# We begin creating an open space environment... with tables, 30x40 meters
  layout <- matrix(0,120,160) #Walkable space or wall? binary
  layout[c(1,120),]<-1; layout[,c(1,160)]<-1; layout[96:100,1]<-0 #Door
  # Add activity tables in greish (2x3 meter ones)
  layout[4:11,146:157] <- 0.35;  layout[110:117,146:157] <- 0.35; layout[4:11,4:15] <- 0.35; layout[110:117,4:15] <- 0.35;
  layout[c(43:48),c(95:102)] <- 0.35; layout[c(43:50),c(103:108)] <- 0.35; #10 and 15 meters off centre
  layout[c(70:77),c(52:57)] <- 0.35; layout[c(72:77),c(58:65)] <- 0.35;

  # We draw the layout
  image(1:120,1:160,layout,zlim=c(0,1),col = grey(seq(1, 0, length = 120)),xaxt='n',yaxt='n',ann=F,axes=FALSE)#,useRaster = T)
  points(c(95,101),c(0,0),pch=15,cex=0.5)
  
# We place the access points...
  points(c(1,55,120),c(100,160,40),pch=19,col="red",cex=2)
  
# Name locations... numerically to put transitions on a list
  text(c(7.5,113.5,7.5,113.5,45.5,74.5),c(151.5,151.5,9.5,9.5,102,58),labels=c("1","2","3","4","5","6"))

# Task clusters and Transition routes, 15s from cluster to cluster
  TaskC <- vector(mode="list",length = 12)
  TaskC[[1]] <- expand.grid(2:13,144:159); points(TaskC[[1]],col="red")
  TaskC[[2]] <- expand.grid(108:119,144:159); points(TaskC[[2]],col="red")
  TaskC[[3]] <- expand.grid(2:13,2:17); points(TaskC[[3]],col="red")
  TaskC[[4]] <- expand.grid(108:119,2:17); points(TaskC[[4]],col="red")
  TaskC[[5]] <- expand.grid(41:52,93:110); points(TaskC[[5]],col="red")
  TaskC[[6]] <- expand.grid(68:79,50:67); points(TaskC[[6]],col="red")
  TaskC[[7]] <- expand.grid(50:70,140:159); points(TaskC[[7]],col="yellow"); text(60,150,label="7")
  TaskC[[8]] <- expand.grid(100:119,70:90); points(TaskC[[8]],col="yellow");text(100,70,label="8")
  TaskC[[9]] <- expand.grid(50:70,2:20); points(TaskC[[9]],col="yellow");text(50,2,label="9")
  TaskC[[10]] <- expand.grid(2:20,70:90); points(TaskC[[10]],col="yellow");text(2,70,label="10")
  TaskC[[11]] <- rbind(expand.grid(30:50,40:60),expand.grid(70:90,100:120)); points(TaskC[[11]],col="yellow");text(c(30,70),c(40,100),labels=c("11","11"))
  TaskC[[12]] <- expand.grid(50:70,70:90); points(TaskC[[12]],col="yellow");text(50,70,label="12")
  
# Trans Matrix index
  transMat <- matrix(c(0,7,10,11,1,5,7,0,12,8,2,2,10,12,0,9,3,3,11,8,9,0,6,4,5,5,5,6,0,5,5,6,6,6,6,0),6,6)

# Signal model (diag 50 meters, then 30x40)
  getSignal <- function(dist){
    return(-30-16*log(1+dist) + rnorm(1,0,2))
  } # If signal below whatever... we do unobserved! use bottom %5 quantile
  getSignals <- function(location){
    location <- loc
    toNorth <- sqrt(sum(((c(55,160)-location)/4)^2))
    toWest <- sqrt(sum(((c(1,100)-location)/4)^2))
    toEast <- sqrt(sum(((c(120,40)-location)/4)^2))
    return(sapply(c(toNorth,toWest,toEast),getSignal))
  }

# Start simulating data, scans will happen every 15s, and transitions from cluster to cluster (if need to move)
  shift <- rep(NA,24000) #100 hours of Scans 100*60*4
  candidateTables<-matrix(0,24000,6); colnames(candidateTables)<-as.character(1:6)
  observedAP<-matrix(1,24000,3); colnames(observedAP)<-c("North","West","East")
    #This one we will populate using logical comparisons on signals, then remove signals
  observedRS<-matrix(NA,24000,3); colnames(observedRS)<-c("North","West","East")
  
# Each shift 30m... 30*4=120 scans... easy to represent on plots
  idx <- 1
  for (sft in 1:200){
    #First iteration
      shift[idx] <- sft
      locCode <- 4
      loc <- TaskC[[locCode]][sample(1:nrow(TaskC[[locCode]]),1),] # Everyone starts at D, they finish where their last task was!
      observedRS[idx,] <- getSignals(loc)
      candidateTables[idx,]<-rep(0,6) # Reset Tasks
      #Task comes up?
        if(runif(1)<0.75) candidateTables[idx+1,sample(1:6,1)] <- 1 
      idx <- idx +1
      for (iter in 1:119) {
        shift[idx] <- sft
        #Choose next location based on tasks, first: any task to finish? if not, move
        if (locCode>6){ # Just get him to whetever he was going
          locCode <- toGo
          loc <- TaskC[[locCode]][sample(1:nrow(TaskC[[locCode]]),1),] 
          observedRS[idx,] <- getSignals(loc)
          candidateTables[idx+1,] <- candidateTables[idx,]
          idx <- idx+1
        } else if (candidateTables[idx,locCode]>0 || sum(candidateTables[idx,])==0){
          loc <- TaskC[[locCode]][sample(1:nrow(TaskC[[locCode]]),1),] 
          observedRS[idx,] <- getSignals(loc)
          candidateTables[idx+1,] <- candidateTables[idx,]
          #Task over? if any
          if(candidateTables[idx,locCode]>0 & runif(1)< (iter<100)*0.2+(iter>=100)*0.5) candidateTables[idx+1,locCode] <- candidateTables[idx+1,locCode] - 1           
          #Task assigned?
          if(runif(1)<0.25 && sum(candidateTables[idx+1,])<5 && iter<100) {
            probs<-rep(1,6); probs[locCode]<-5
            candidateTables[idx+1,] <- candidateTables[idx+1,] + rmultinom(1,1,prob=probs) 
          } #more likelihood to get it assigned where you are
          idx <- idx+1
        } else {
          toGo <- which.max(candidateTables[idx,]) # A transition, not task adding or completing when transitioning
          locCode <- transMat[locCode,toGo]
          loc <- TaskC[[locCode]][sample(1:nrow(TaskC[[locCode]]),1),] 
          observedRS[idx,] <- getSignals(loc)
          candidateTables[idx+1,] <- candidateTables[idx,]
          idx <- idx+1
        }
      }  
  }
  
  #Prepare to save things (and create observed AP)
  cut<-quantile(observedRS,0.05)
  observedAP[which(observedRS<cut)]<-0
  
  save(candidateTables,file="candidateTables.Rdata")
  save(observedRS,file="observedRS.Rdata")
  save(observedAP,file="observedAP.Rdata")
  
  
  plot(density(observedRS[candidateTables[,2]>1,1]))
  colMeans(observedAP[candidateTables[,2]>1,])
  
  
  