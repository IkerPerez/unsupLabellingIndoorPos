
# Read off all the relevant files we have so far
  load("candidateTables.Rdata"); load("observedAP.Rdata"); load("observedRS.Rdata")

  load("outputFiles/obsProbs.Rdata")
  load("outputFiles/obsMeans.Rdata")
  load("outputFiles/obsSDs.Rdata")
  load("outputFiles/probabilityTables.Rdata")
  load("outputFiles/lambdas.Rdata")
  load("outputFiles/nus.Rdata")
  
  ################
  # Recreate Map #
  ################
  
  library(raster)

  # We begin creating an open space environment... with tables, 30x40 meters
  red <- matrix(1,120,160); green<- red; blue<-red #Walkable space or wall? binary
    #Walls
    red[c(1,120),]<-0; red[,c(1,160)]<-0; red[96:100,1]<-1; green <-red; blue <-red
  
    # Add activity tables in brown (2x3 meter ones)
    r <- 0.6; g <- 0.5; b <- 0.45
    red[4:11,146:157] <- r;  red[110:117,146:157] <- r; red[4:11,4:15] <- r; red[110:117,4:15] <- r;
    red[c(43:48),c(95:102)] <- r; red[c(43:50),c(103:108)] <- r; #10 and 15 meters off centre
    red[c(70:77),c(52:57)] <- r; red[c(72:77),c(58:65)] <- r;
    green[4:11,146:157] <- g;  green[110:117,146:157] <- g; green[4:11,4:15] <- g; green[110:117,4:15] <- g;
    green[c(43:48),c(95:102)] <- g; green[c(43:50),c(103:108)] <- g; #10 and 15 meters off centre
    green[c(70:77),c(52:57)] <- g; green[c(72:77),c(58:65)] <- g;
    blue[4:11,146:157] <- b;  blue[110:117,146:157] <- b; blue[4:11,4:15] <- b; blue[110:117,4:15] <- b;
    blue[c(43:48),c(95:102)] <- b; blue[c(43:50),c(103:108)] <- b; #10 and 15 meters off centre
    blue[c(70:77),c(52:57)] <- b; blue[c(72:77),c(58:65)] <- b;  
    
  #Create a map
  layout <- rgb(red[,seq(160,1,-1)],green[,seq(160,1,-1)],blue[,seq(160,1,-1)]) #red matrix, green matrix, blue matrix
  dim(layout) <- dim(red); layout <- as.raster((t(layout)))
  
  #Plot it
  par(mfrow = c(1,2),
      oma = c(0,0,4,0) + 0.1, # oma will push graphs to make space for title
      mar = c(2,0,0,0) + 2) 
  plot(layout,interpolate=F,ylim=c(0,160), xlim=c(0,120))
  points(c(95,101),c(0,0),pch=15,cex=0.5)
  
  # We place the access points...
  points(c(1,55,119),c(100,159,40),pch=19,col="red",cex=1.2)
  
  # Name locations... numerically to put transitions on a list
  text(c(7.5,113.5,7.5,113.5,45.5,74.5),c(151.5,151.5,9.5,9.5,102,58),labels=c("1","2","3","4","5","6"),font=2,cex=0.85)
  
  # EXAMPLE OF TASK DATA
  plot(c(3,12),c(1,1),type="l",xlim=c(0,30),ylim=c(0.5,10),xlab="",ylab=" ",yaxt="n",bty="n",lty=1,lwd=3,col=1); lines(c(3,3),c(1,9),lty=3); points(c(3,12),c(1,1),pch=19)
  lines(c(5,18),c(2,2),type="l",xlim=c(0,30),lty=1,lwd=3,col=1); lines(c(5,5),c(2,9),lty=3);points(c(5,18),c(2,2),pch=19)
  lines(c(8,19),c(3,3),type="l",xlim=c(0,30),lty=1,lwd=3,col=1);lines(c(8,8),c(3,9),lty=3);points(c(8,19),c(3,3),pch=19)
  lines(c(14,20),c(4,4),type="l",xlim=c(0,30),lty=1,lwd=3,col=1);lines(c(14,14),c(4,9),lty=3);points(c(14,20),c(4,4),pch=19)
  lines(c(22,26),c(5,5),type="l",xlim=c(0,30),lty=1,lwd=3,col=1);lines(c(22,22),c(5,9),lty=3);points(c(22,26),c(5,5),pch=19)
  lines(c(24,27),c(6,6),type="l",xlim=c(0,30),lty=1,lwd=3,col=1);lines(c(24,24),c(6,9),lty=3);points(c(24,27),c(6,6),pch=19)
  text(c(3,5,8,14,22,24),rep(9.5,6),labels=c("1","3","1","4","6","3"),font=1); text(14,10.2,"Table"); text(28,0.4,"Time")

  mtext("Location Map and Activity Data", outer = TRUE, cex = 1.5)

  
  ##############
  # Classifier #
  ##############  
  
  # Now we need four matrices each defining the probability for every cluster at each pixel
  b <- function(w,scans){
    dummy<- rep(0,3)
    observedAP<-as.numeric(scans>-91)
    idx <- as.logical(observedAP)
    if( sum(idx)==0 ) return(prod(1 - p[w,]))
    dummy[!idx] <- (1 - p[w,!idx])
    dummy[idx] <- p[w,idx] * dnorm( (as.numeric(scans[idx]) - mu[w,idx]) / sigma[w,idx] ) / sigma[w,idx]
    return(prod(dummy))
  }
  bVector <- function(scans){
    probs <- sapply(1:nrow(p), function(x) b(x,scans))
    return(probs/sum(probs))
  }
  #Build list of matrices
  probsTable <- vector("list",6)
  for (i in 1:6) probsTable[[i]]<-matrix(1,120,160)
  for (i in 1:120){
    for (j in 1:160){
      signals<- c(sqrt(sum(((c(55,160)-c(i,j))/4)^2)),sqrt(sum(((c(1,100)-c(i,j))/4)^2)),sqrt(sum(((c(120,40)-c(i,j))/4)^2)))
      signals<- -30-16*log(1+signals)
      probs <- bVector(signals)
      for (k in 1:6) probsTable[[k]][i,j] <- probs[k]
    }
  }
  
  #Assign one color to each... and begin the plot again
  reds<-c(69,0,131,153,205,205)/255; greens<-c(139,139,139,50,149,51)/255; blues<-c(0,139,139,204,12,51)/255
  
  # We begin creating an open space environment... with tables, 30x40 meters
  red <- matrix(1,120,160); green<- red; blue<-red #Walkable space or wall? binary
  
  # Colors reflecting probabilities
  for (i in 1:120){
    for (j in 1:160){
      #ue<-rep(0,6)
      #ue[which.max(sapply(1:6,function(x) probsTable[[x]][i,j]))]<-1
      ue<-sapply(1:6,function(x) probsTable[[x]][i,j])
      red[i,j] <- sum(reds * ue)
      green[i,j] <- sum(greens * ue)
      blue[i,j] <- sum(blues * ue)
    }
  }

  #Walls
  red[c(1,120),]<-0; red[,c(1,160)]<-0; red[96:100,1]<-1;
  green[c(1,120),]<-0; green[,c(1,160)]<-0; green[96:100,1]<-1;
  blue[c(1,120),]<-0; blue[,c(1,160)]<-0; blue[96:100,1]<-1;
    
  # Add activity tables in brown (2x3 meter ones)
  r <- 0.6; g <- 0.5; b <- 0.45
  red[4:11,146:157] <- r;  red[110:117,146:157] <- r; red[4:11,4:15] <- r; red[110:117,4:15] <- r;
  red[c(43:48),c(95:102)] <- r; red[c(43:50),c(103:108)] <- r; #10 and 15 meters off centre
  red[c(70:77),c(52:57)] <- r; red[c(72:77),c(58:65)] <- r;
  green[4:11,146:157] <- g;  green[110:117,146:157] <- g; green[4:11,4:15] <- g; green[110:117,4:15] <- g;
  green[c(43:48),c(95:102)] <- g; green[c(43:50),c(103:108)] <- g; #10 and 15 meters off centre
  green[c(70:77),c(52:57)] <- g; green[c(72:77),c(58:65)] <- g;
  blue[4:11,146:157] <- b;  blue[110:117,146:157] <- b; blue[4:11,4:15] <- b; blue[110:117,4:15] <- b;
  blue[c(43:48),c(95:102)] <- b; blue[c(43:50),c(103:108)] <- b; #10 and 15 meters off centre
  blue[c(70:77),c(52:57)] <- b; blue[c(72:77),c(58:65)] <- b;  
  
  #Grey all out
  red <- red*0.8 + c(131,139,139)/255*0.2; green <- green*0.8 + c(131,139,139)/255*0.2; blue <- blue*0.8 + c(131,139,139)/255*0.2;  
  
  #Create a map
  layout <- rgb(red[,seq(160,1,-1)],green[,seq(160,1,-1)],blue[,seq(160,1,-1)]) #red matrix, green matrix, blue matrix
  dim(layout) <- dim(red); 
  
  #image(matrix(1:19200,ncol=160,nrow=120),col=layout,xaxt="n",yaxt="n")
  layout <- as.raster((t(layout)))
  
  #Plot it
  par(mfrow = c(1,1),
      oma = c(0,0,0,0) + 0.1, # oma will push graphs to make space for title
      mar = c(2,2,0,0) + 2) 
  plot(layout,interpolate=F,ylim=c(0,160), xlim=c(0,120))
  points(c(95,101),c(0,0),pch=15,cex=0.5)
  
  # We place the access points...
  points(c(1,55,119),c(100,159,40),pch=19,col="red",cex=1.2)
  
  # Name locations... numerically to put transitions on a list
  text(c(7.5,113.5,7.5,113.5,45.5,74.5),c(151.5,151.5,9.5,9.5,102,58),labels=c("1","2","3","4","5","6"),font=2,cex=0.85)
  
  # draw an axis on the left
  axis(2, at=seq(0,160,40),labels=seq(0,40,10), las=2,tck=-.01,cex.axis=0.7)
  
  # draw an axis on the right, with smaller text and ticks
  axis(1, at=seq(0,120,40),labels=seq(0,30,10), las=1, tck=-.01, pos = c(-7,0),cex.axis=0.7)

  # Add coordinates
  for(i in 1:39){
    lines(c(1,119),c(i*4,i*4),lty=3,lwd=0.2,col=gray.colors(1,0.4))
  }
  for(i in 1:29){
    lines(c(i*4,i*4),c(1,159),lty=3,lwd=0.2,col=gray.colors(1,0.4))
  }

  ##################
  # ANALYSE STUFF! #  
  ##################
  
# We begin with transition matrix under the influence of no tasks at all!
  library(lattice)
  require(gridExtra)
  A <- exp(lambda+ t(t(nu)) %*% rep(0,6) )
  for (i in 1:6) A[i,] <- A[i,]/sum(A[i,])
  
  plot1 <- levelplot(t(A),column.values = seq(6,1,-1),cuts=15,colorkey=list(height=.5, space="right", at=seq(0, 1, 1/15), cuts=15),
                    main="Transition Matrix 1", xlab="", ylab="",col.regions = gray(seq(15,0,-1)/15),
                    scales=list(y=list(at=seq(6,1,-1), labels=colnames(candidateTables), cex=1),
                                x=list(at=1:6, labels=colnames(candidateTables),rot=90,cex=1 )) )
           
  B<-A; diag(B)<-0; for(i in 1:6){B[i,]<-B[i,]/sum(B[i,])}
  plot2 <- levelplot(t(B),column.values = seq(6,1,-1),cuts=15,colorkey=list(height=.5, space="right", at=seq(0, 1, 1/15), cuts=15),
                     main="Transition Matrix", xlab="", ylab="",col.regions = gray(seq(15,0,-1)/15),
                     scales=list(y=list(at=seq(6,1,-1), labels=colnames(candidateTables), cex=1),
                                 x=list(at=1:6, labels=colnames(candidateTables),rot=90,cex=1 )) )
  grid.arrange(plot1,plot2,ncol=2) # Screw this matrix... people walk so quick and two scans could be anywhere!

  
# Overview of signal reception probability from tables
  mainWardsNames <- as.character(seq(1:6))
  idx <- sapply(mainWardsNames, function(x) which(colnames(candidateTables)==x))
  par(mfrow=c(2,3))
  for (i in idx) {plot(1:ncol(p),as.numeric(p[i,]),type="h",xlim=c(1,ncol(p)),ylim=c(0,1),
                        main=colnames(candidateTables)[i],xlab="AP Number",ylab="Probability")}
  
  rm(mainWardsNames,idx,i)

# Overview of received signals from major HS for ward Southwell
  idx <-which(colnames(candidateTables)==c("1")); 
  HSidx <- order(p[idx,],decreasing = T)[1:3]
  signals <- c()
    idx2 <- which(tableP[,idx]>0.00001)
    if(length(idx2)>0) signals <- rbind(signals,observedRS[idx2,HSidx])
  par(mfrow=c(1,3))
  for(i in 1:3){
    dummy<-density(signals[!is.na(signals[,i]),i],adjust=2)
    dummy$y <- dummy$y * p[idx,HSidx[i]]
    plot(dummy,xlim=c(-100,-30),ylim=c(0,0.2),main=colnames(observedRS)[HSidx[i]]);
    lines(seq(-120,0,0.01), 
          p[idx,HSidx[i]]*dnorm(seq(-120,0,0.01),mu[idx,HSidx[i]],sigma[idx,HSidx[i]]),
          col="red",lty=2)
  }
  rm(dummy, HSidx, i, idx, idx2)
  
# Now try and filter the position of the phone over a whole shift were it is being used  
  shifts <- rep(1:200,each=120) #100 hours of Scans 100*60*4
  shift <- which(shifts==1)

  candWards1 <- candidateTables[shift,]
  probabilityWards1 <- tableP[shift,]

# We will NOT use a heatmap, it's awful
  par(mar=c(5, 5, 4, 2) + 0.5, oma = c(0,0,0,0) + 0.1, mfrow=c(2,1) )
  plot(c(1,nrow(candWards1)), c(20,20), type="l",col=gray(1),lwd=6,
       xlim=c(1,nrow(candWards1)),ylim=c(0.8,ncol(candWards1)+0.2),yaxt="n",xaxt="n",xlab="Minutes",ylab="Location",bty="n",
       main="Proportional Activity Plot")
  axis(1, at=seq(0,120,20),labels=as.character(seq(0,30,5)),
       col.axis="black", las=1, cex.axis=0.9, tck=-.01)  
  axis(2, at=1:ncol(candWards1),labels=colnames(candWards1),
       col.axis="black", las=2, cex.axis=0.9, tck=-.01)
  for(j in 1:ncol(candWards1)){
    for(i in 1:(nrow(candWards1))){
        points(c(i),c(j),col=gray(   1-(3^candWards1[i,j]/sum(3^candWards1[i,]))  ),pch=15,cex=1)
    }
  }
  plot(c(1,nrow(candWards1)), c(20,20), type="l",col=gray(1),lwd=6,
       xlim=c(1,nrow(candWards1)),ylim=c(0.8,ncol(candWards1)+0.2),yaxt="n",xaxt="n",xlab="Minutes",ylab="Location",bty="n",
       main="Labelled Location")
  axis(1, at=seq(0,120,20),labels=as.character(seq(0,30,5)),
       col.axis="black", las=1, cex.axis=0.9, tck=-.01)  
  axis(2, at=1:ncol(candWards1),labels=colnames(candWards1),
       col.axis="black", las=2, cex.axis=0.9, tck=-.01)
  for(j in 1:ncol(candWards1)){
    for(i in 1:(nrow(candWards1))){
      points(c(i),c(j),col=gray(   1-probabilityWards1[i,j]*0.8  ),pch=15,cex=1)
    }
  }
  
