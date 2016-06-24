
# Read off all the relevant files we have so far
  load("listObservedSignals.Rdata")
  load("listObservedHotSpots.Rdata")
  load("listCandidateWards.Rdata")
  load("listTimeAndPhones.Rdata")
  
  # Do not know why I had the idea to try with an "other" ward... with such small data set ignore it right away!
  for (i in 1:length(listCandidateWards)) listCandidateWards[[i]] <- listCandidateWards[[i]][,-10]
  
  library(RMySQL)
  myDB <- dbConnect(MySQL(), user="root",
                    password="Anz5Ur8TUPVuh",
                    dbname="wayward",
                    host="localhost")
  
  dbListTables(myDB)
  hanData <- ( fetch( dbSendQuery(myDB, "select * from han") , n=-1) ) # Data from the task info
  drShiftTimes <- ( fetch( dbSendQuery(myDB, "select * from drshift") , n=-1) ) # Relates phones at use with each doctor/shift
  shifts <- read.csv("shiftsTimmings.csv",header=T)    # Jesse's is wrong for some reason
  dbDisconnect(myDB); rm(myDB)
  
  load("outputFiles/obsProbs.Rdata")
  load("outputFiles/obsMeans.Rdata")
  load("outputFiles/obsSDs.Rdata")
  load("outputFiles/probabilityWards.Rdata")
  load("outputFiles/lambdas.Rdata")
  load("outputFiles/nus.Rdata")
  

  ##################
  # ANALYSE STUFF! #  
  ##################
  
# We begin with transition matrix under the influence of no tasks at all!
  library(lattice)
  require(gridExtra)
  A <- exp(lambda+ t(t(nu)) %*% rep(0,9) )
  for (i in 1:9) A[i,] <- A[i,]/sum(A[i,])
  
  plot1 <- levelplot(t(A),column.values = seq(9,1,-1),cuts=15,colorkey=list(height=.5, space="right", at=seq(0, 1, 1/15), cuts=15),
                    main="Transition Matrix 1", xlab="", ylab="",col.regions = gray(seq(15,0,-1)/15),
                    scales=list(y=list(at=seq(9,1,-1), labels=colnames(listCandidateWards[[1]]), cex=1),
                                x=list(at=1:9, labels=colnames(listCandidateWards[[1]]),rot=90,cex=1 )) )
           
  B<-A; diag(B)<-0; for(i in 1:7){B[i,]<-B[i,]/sum(B[i,])}
  plot2 <- levelplot(t(B),column.values = seq(9,1,-1),cuts=15,colorkey=list(height=.5, space="right", at=seq(0, 1, 1/15), cuts=15),
                    main="Transition Matrix 2", xlab="", ylab="",col.regions = gray(seq(15,0,-1)/15),
                    scales=list(y=list(at=seq(9,1,-1), labels=colnames(listCandidateWards[[1]]), cex=1),
                                x=list(at=1:9, labels=colnames(listCandidateWards[[1]]),rot=90,cex=1 )) )
  grid.arrange(plot1,plot2,ncol=2)
  sort(table(hanData$ward))
  

# Overview of signal reception probability of HS from major areas (Also use wards in other data set!)
  mainWardsNames <- c("NorthWest", "Central", "SouthWest", "NurseCoord")
  idx <- sapply(mainWardsNames, function(x) which(colnames(listCandidateWards[[1]])==x))
  par(mfrow = c(2,2),
      oma = c(1,0,0,0) + 0.1, # oma will push graphs to make space for title
      mar = c(1,1,1,-1) + 2) 
  
  for (i in idx) {
    plot(1:ncol(p),as.numeric(p[i,]),type="h",xlim=c(1,ncol(p)),ylim=c(0,1),
                        main=c("North",colnames(listCandidateWards[[1]])[2:9])[i],bty="n",xaxt="n",xlab=" ",ylab=" ",col=gray.colors(1,0.2))
    axis(1, at=seq(0,ncol(p),15),labels=as.character(seq(0,ncol(p),15)),
         col.axis="black", las=1, cex.axis=1, tck=-.05) 
    
    }
  
  rm(mainWardsNames,idx,i)

# Overview of received signals from major HS for South
  idx <-which(colnames(listCandidateWards[[1]])=="South"); colnames(listCandidateWards[[1]])[idx]
  HSidx <- order(p[idx,],decreasing = T)[1:4]
  signals <- c()
  for(iepe in 1:8){
    idx2 <- which(wardP[[iepe]][,idx]>0.95)
    if(length(idx2)>0) signals <- rbind(signals,listObservedSignals[[iepe]][idx2,HSidx])
  }
  par(mfrow = c(2,2),
      oma = c(1,0,0,0) + 0.1, # oma will push graphs to make space for title
      mar = c(1,1,1,-1) + 2) 
  for(i in 1:4){
    dummy<-density(signals[!is.na(signals[,i]),i],adjust=2)
    dummy$y <- dummy$y * p[idx,HSidx[i]]
    plot(dummy,xlim=c(-100,-30),ylim=c(0,0.06),main=paste("AP:",colnames(listObservedSignals[[iepe]])[HSidx[i]]),
         lwd=1.5,bty="n",xlab=" ",ylab=" ");
    lines(seq(-120,0,0.01), 
          p[idx,HSidx[i]]*dnorm(seq(-120,0,0.01),mu[idx,HSidx[i]],sigma[idx,HSidx[i]]),
          col="red",lty=2,lwd=1.5)
  }
  rm(dummy, HSidx, i, idx, idx2)

  
# Now try and filter the position of the phone over a whole shift were it is being used  
  shift <- drShiftTimes[2,]  
  starttime <- shifts$starttime[shifts$shiftname==shift$shift]
  endtime <- shifts$endtime[shifts$shiftname==shift$shift]
  
  ue<-hanData[hanData$timeof_complete-hanData$duration_max_accept*60 >= starttime-60*60 & hanData$timeof_complete <= endtime+60*60 & hanData$name==shift$hash,] # Just checking...
  ue<-ue[order(ue$timeof_complete-ue$duration_max_accept*60),]
  cbind(acceptTime=(ue$timeof_complete-ue$duration_max_accept*60-starttime)/60,
        completeTime=(ue$timeof_complete-starttime)/60,ue[,c(8,9,12)])
  
  listIdx<-0; for (i in 1:8) {if (listTimeAndPhones[[i]][1,2] == shift$tattooNo) listIdx <- i}
  idx <- which(listTimeAndPhones[[listIdx]]$time >= starttime-60*60 & listTimeAndPhones[[listIdx]]$time <= endtime+60*60
               & listTimeAndPhones[[listIdx]]$doctor==shift$hash)
  
  timeAndPhone1 <- listTimeAndPhones[[listIdx]][idx,]
  candWards1 <- listCandidateWards[[listIdx]][idx,]
  probabilityWards1 <- wardP[[listIdx]][idx,]

# We will NOT use a heatmap, it's awful
  par(mar=c(5, 5, 4, 2) + 0.5, oma = c(0,2,0,0) + 0.1, mfrow=c(2,1) )
  plot(c(1,nrow(candWards1)), c(20,20), type="l",col=gray(1),lwd=6,
       xlim=c(1,nrow(candWards1)),ylim=c(0.8,ncol(candWards1)+0.2),yaxt="n",xaxt="n",xlab="Minutes",ylab="",bty="n",
       main="Proportional Activity Plot")
  axis(1, at=seq(0,nrow(candWards1),30),labels=as.character(seq(0,nrow(candWards1),30)),
       col.axis="black", las=1, cex.axis=0.9, tck=-.01)  
  axis(2, at=1:ncol(candWards1),labels=colnames(candWards1),
       col.axis="black", las=2, cex.axis=0.9, tck=-.01)
  for(j in 1:ncol(candWards1)){
    for(i in 1:(nrow(candWards1))){
      points(c(i),c(j),col=gray(   1-(3^candWards1[i,j]/sum(3^candWards1[i,]))  ),pch=15,cex=1)
    }
  }
  plot(c(1,nrow(candWards1)), c(20,20), type="l",col=gray(1),lwd=6,
       xlim=c(1,nrow(candWards1)),ylim=c(0.8,ncol(candWards1)+0.2),yaxt="n",xaxt="n",xlab="Minutes",ylab="",bty="n",
       main="Labelled Location")
  axis(1, at=seq(0,nrow(candWards1),30),labels=as.character(seq(0,nrow(candWards1),30)),
       col.axis="black", las=1, cex.axis=0.9, tck=-.01)  
  axis(2, at=1:ncol(candWards1),labels=colnames(candWards1),
       col.axis="black", las=2, cex.axis=0.9, tck=-.01)
  for(j in 1:ncol(candWards1)){
    for(i in 1:(nrow(candWards1))){
      points(c(i),c(j),col=gray(   1-probabilityWards1[i,j]*0.8  ),pch=15,cex=1)
    }
  }
 
# Now we check the scans james and I got!
  scanFiles <- paste("../citySurvey260116/",list.files(path = "../citySurvey260116"),sep="")
  
  dummy <- rep(0,ncol(listObservedHotSpots[[1]]))
  b <- function(obs,sig,w){
    if( sum(obs)==0 ) return(prod(1 - p[w,]))
    dummy[!obs] <- (1 - p[w,!obs])
    dummy[obs] <- p[w,obs] * dnorm( (sig[obs] - mu[w,obs]) / sigma[w,obs] ) / sigma[w,obs]
    return(prod(dummy))
  }
  bVector <- function(obs,sig){
    probs <- sapply(1:nrow(p), function(x) b(obs,sig,x))
    return(probs/sum(probs))
  }
  for (file in scanFiles[2:(length(scanFiles)-2)]){
    dummyScans <- read.csv(file,header=F)[,c(2,6)]
    dummyScans[,1] <- as.character(dummyScans[,1]); dummyScans[,2] <- as.character(dummyScans[,2])
    for (i in 1:nrow(dummyScans)){
      dummyScans[i,1] <- paste(strsplit(dummyScans[i,1],split=":")[[1]][-6],collapse=".")
      if (!is.na(as.numeric(strsplit(dummyScans[i,1],split="\\.")[[1]][1]))){
        dummyScans[i,1] <- paste(c("X",dummyScans[i,1]),collapse="")
      }
      dummyScans[i,2] <- strsplit(dummyScans[i,2],split=" ")[[1]][1]
    }
    dummyScans <- dummyScans[order(dummyScans[,2]),]
    dummyScans <- dummyScans[!duplicated(dummyScans[,1]),]
    obs <- colnames(listObservedHotSpots[[1]]) %in% dummyScans[,1]
    sig <- rep(0,ncol(listObservedHotSpots[[1]]))
    for (i in which(obs)){
      sig[i] <- dummyScans[which(dummyScans[,1]==colnames(listObservedHotSpots[[1]])[i]),2]
    }
    sig<-as.numeric(sig)
    likel <- bVector(obs,sig)
    print(list(file,sig,rbind(colnames(listCandidateWards[[1]]),round(likel,3))))
    
  }
  