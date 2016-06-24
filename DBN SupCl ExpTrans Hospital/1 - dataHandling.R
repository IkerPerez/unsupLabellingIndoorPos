
library(RMySQL)

# Create database connection object
  myDB <- dbConnect(MySQL(), user="...",
                    password="...",
                    dbname="...",
                    host="...")
  
  dbListTables(myDB)
  scansData <- ( fetch( dbSendQuery(myDB, "select * from scans") , n=-1) ) # Data from the phone Scans
  hanData <- ( fetch( dbSendQuery(myDB, "select * from han") , n=-1) ) # Data from the task info
  drShiftTimes <- ( fetch( dbSendQuery(myDB, "select * from drshift") , n=-1) ) # Relates phones at use with each doctor/shift
  phoneNames <- ( fetch( dbSendQuery(myDB, "select * from phones") , n=-1) ) # Relates the phone name with an IMEI identity
  #shifts <- ( fetch( dbSendQuery(myDB, "select * from shifts") , n=-1) ) # Start and End time for shifts :), not that doctors respect these times :P
  shifts <- read.csv("shiftsTimmings.csv",header=T)    # Jesse's is wrong for some reason
  
# Now we disconnect!
  dbDisconnect(myDB); rm(myDB)
  
# Han Data description (rest self-explanatory)
#
#   - timeof_submit:        Integer, Time a task was submitted at (in seconds, rounded to closest minute)
#   - timeof_accept:        Integer, Time a task submitted was accepted at (in seconds, rounded to closest minute)
#   - timeof_complete:      Integer, Time a task accepted was marked complete at (in seconds, rounded to closest minute)
#   - duration_submit:      Integer, minute duration from submission to completion (Feels includes many outliers)
#   - duration_min_accept:  Integer, minute duration from accepting (by First Doctor) to completion (Feels includes many outliers)
#   - duration_max_accept:  Integer, minute duration from accepting (by Doctor who completed) to completion (Feels includes some outliers)
#   - priority:             Integer (1-5), priority assigned by main nurse (can ignore 2 and 4)  
#   - ward:                 String, name of the ward
#   - site:                 String, name of the hospital
#   - category:             String, type of task (30 different, including "None of the above")
#   - close_reason:         String, reason for completion
#   - staff_group:          String, type of staff dealing with task
#   - name:                 String, anonimous name of doctor (20 different, one of them rather ignore since single observation)
#
#   ********************************************************************************************  
  
# Quick fix to priority level 4
  hanData$priority[hanData$priority == 4] <- hanData$priority[hanData$priority == 4] + 1

# Get the time in seconds when the doctor dealing with the task accepted it
  hanData <- data.frame( hanData[,1:2] , timeof_max_accept = hanData$timeof_complete - hanData$duration_max_accept*60 , hanData[,3:13] )

# Filter out unecessary data, by SSID
  citySSID <- tail(names(sort(table(scansData$ssid))),7) # Only the ones deployed within hospital
  toKeep <- scansData$ssid %in% citySSID; rm(citySSID)
  scansData <- scansData[toKeep,]; rm(toKeep) 
  length(table(scansData$ap_id)) # 720 hotspots... hopefully reduced as I merge with han data :)
  table(scansData$imei) # Phone 16 seems to have arrived in the hospital much earlier!

# Filter out unecessary data, by Time window
  lowerEnd <- min(hanData$timeof_submit); upperEnd <- max(hanData$timeof_complete)
  (upperEnd - lowerEnd)/60/60/24; hist(hanData$timeof_complete,50) # Need remove the testing week
  cutPoint <- 1.394e+09
  hanData <- hanData[hanData$timeof_submit > cutPoint, ]; rm(cutPoint)
  lowerEnd <- min(hanData$timeof_submit); upperEnd <- max(hanData$timeof_complete)
  (upperEnd - lowerEnd)/60/60/24; hist(hanData$timeof_complete,50) 
  toKeep <- (scansData$time >= lowerEnd & scansData$time <= upperEnd)
  scansData <- scansData[toKeep,]; rm(toKeep,lowerEnd,upperEnd)
  
# Note the following activation and turning off times for phones
  for (i in unique(scansData$imei)){
      print(i); print(c(as.POSIXct(min(scansData$time[scansData$imei==i]),origin="1970-01-01"),
                        as.POSIXct(max(scansData$time[scansData$imei==i]),origin="1970-01-01")))
  }

# Store that table to help define the time/phones boundaries
  phoneActivationTimes<-data.frame(phone=unique(scansData$imei),turnOn=rep(0,8),turnOff=rep(0,8))
  for (i in 1:8){
    phoneActivationTimes$turnOn[i] <- min(scansData$time[scansData$imei==phoneActivationTimes$phone[i]])
    phoneActivationTimes$turnOff[i] <- max(scansData$time[scansData$imei==phoneActivationTimes$phone[i]])
  } 
  #write.csv(phoneActivationTimes,"phoneActivationTimes.csv",row.names=F)
  
# First phone arrived at the hospital at 1398629672, Unix timestamp
#
#   That is Sun, 27 Apr 2014 20:14:32 GMT.
#   The shifts had been running since very early that day... that's why we have lots of Han data
#     without scan we can associate them to :(.
#   Most phones yet arrived for the Sunday Night shift... so we use 1398636000 as starting point
#     for all data sets.
#   Also, phones were picked up on a Monday Night, but were last taken by doctors on the Sunday
#     day shift. Probably resting at coord room. since. We use end point 1399280400.
#   Note that battery for two phones dies out before end-of-week... need to ignore those phones
#     for those points onwards or I'll get false negatives!
#
#************************************************************************************************
  
# Cap All data using starting value...
  lowerEnd <- 1398636000; upperEnd <- 1399280400
  toKeep <- (hanData$timeof_max_accept >= lowerEnd & hanData$timeof_complete <= upperEnd) # The max accept time is what matters, he's the one who dealt with issue
  hanData <- hanData[toKeep,]; 
  hanData <- hanData[,-c(4,6)] # We only care about the times for the doctor who completed task
  
# Now remove data that with tasks that did not lead to completion
  hanData <- hanData[hanData$close_reason == "Completed",]
  
# Same with scans
  toKeep <- (scansData$time>=lowerEnd & scansData$time <= upperEnd)
  scansData <- scansData[toKeep,];
  rm(toKeep,lowerEnd,upperEnd,i)
  
# And so we start building the tables we need
#
#   These are to be fed to a forward/backward EM algorithm. We will create one table per phone one entry per minute.
#   Then we will try associate these entries to doctors and possible locations.
#
#****************************************************************************************************

# Start sorting han and scan tables by time
  hanData <- hanData[order(hanData$timeof_max_accept,decreasing=FALSE),]
  scansData <- scansData[order(scansData$time,decreasing=FALSE),]
  table(tail(scansData$time,nrow(scansData)-1)-head(scansData$time,nrow(scansData)-1)) # Seems phones were not synced :(

# If you check individually, you will note one phone's scan frequency was higher (phone 16)
#
#   Frequency is every 1 minute for all of them except this one :)
#
#******************************************************************************************

# We build a data set of phones and times, we will split it later!
  timeWindow <- seq(1398636000,1399280400,60)
  timeAndPhones <- data.frame(matrix(NA,length(timeWindow)*8,2))
  names(timeAndPhones) <- c("time","Phone")
  timeAndPhones$time<-rep(timeWindow,each=8); timeAndPhones$Phone<-rep(unique(scansData$imei),length(timeWindow))  
  
# Remove time/phone entries where devices were off
  phoneOn <- function(index){
    dummy <- timeAndPhones[index,]
    if (dummy$time >= phoneActivationTimes$turnOn[phoneActivationTimes$phone==dummy$Phone] &
        dummy$time <= phoneActivationTimes$turnOff[phoneActivationTimes$phone==dummy$Phone]) {
      return(T)
    } else {
      return(F)
    }
  }
  toKeep<- sapply(1:nrow(timeAndPhones), phoneOn)
  timeAndPhones <- timeAndPhones[toKeep,]
  rm(toKeep,timeWindow,phoneOn)
  
# Change IMEI codes to phone names
  returnName <- function(imei) {
    return(as.numeric(phoneNames$phone_name[phoneNames$phone_imei==imei]))
  }
  timeAndPhones$Phone <- sapply(timeAndPhones$Phone, returnName)
  timeAndPhones <- timeAndPhones[order(timeAndPhones$time,timeAndPhones$Phone),]; rm(returnName)
  
  
# Now the hard part starts!!!
#
#   Associate those phones and times to doctors or coordinator room, or try. :)
#   Then extract the possible wards.
#   Finally, add the received signal strengths and obersvation binaries for EACH hotspot.
#
#****************************************************************************************
  
  timeAndPhones$doctor <- "none"
  timeAndPhones$approxShiftActive <- "no" # If the phone was taken at some time close... then we account for the fact
                                          # there is a chance it could be anywhere!
  
# Start with phones to tasks
  drShiftTimes <- data.frame(drShiftTimes, startApprx = rep(0,nrow(drShiftTimes)), endApprx = rep(0,nrow(drShiftTimes)))
  for (i in 1:nrow(drShiftTimes)){
    drShiftTimes$startApprx[i] <- shifts$starttime[shifts$shiftname == drShiftTimes$shift[i] ]
    drShiftTimes$endApprx[i] <- shifts$endtime[shifts$shiftname == drShiftTimes$shift[i] ]
  }
  hanData <- data.frame( hanData , phone = rep("none",nrow(hanData)))
  hanData$phone <- "none"
  for (i in 1:nrow(hanData)){
    shiftInfo <- drShiftTimes[ drShiftTimes$hash==hanData$name[i] &
                                 !(hanData$timeof_complete[i]<drShiftTimes$startApprx) &
                                 !(hanData$timeof_max_accept[i]>drShiftTimes$endApprx) ,]
    if (nrow(shiftInfo) == 1){
      hanData$phone[i] <- shiftInfo$tattooNo
    } else if (nrow(shiftInfo) > 1){
      print("beware!")
    }
  }
  rm(i)  

# If they had no phone on them, their task info is irrelevant
  hanData <- hanData[ hanData$phone!="none",]

# So we add to timeAndPhone info on who should have had the devices
  findDoc <- function(time,phone){
    iep <- drShiftTimes[ drShiftTimes$tattooNo == phone & 
                           time > drShiftTimes$startApprx & 
                           time <= drShiftTimes$endApprx , ]
    if(nrow(iep) == 1){
      return(iep$hash)
    } else if (nrow(iep) > 1){
      print("beware!")
      print(time); print(phone)
    } else {
      return("none")
    }
  }
  timeAndPhones$doctor <- apply(timeAndPhones[,1:2],1,function(x) findDoc(x[1],x[2]))
  table(timeAndPhones$doctor) # Around 83 of scans should be from nurse coord room
  findShift <- function(time,phone){
    iep <- drShiftTimes[ drShiftTimes$tattooNo == phone & 
                           time > drShiftTimes$startApprx - 60*60 & # One hour less
                           time <= drShiftTimes$endApprx + 60*60 , ] # One hour more
    if(nrow(iep) >= 1){
      return("yes")
    } else {
      return("no")
    }
  }
  timeAndPhones$approxShiftActive <- apply(timeAndPhones[,1:2],1,function(x) findShift(x[1],x[2]))
  table(timeAndPhones$approxShiftActive) 
  
# Split in a list by phone and save the object
  listTimeAndPhones <- vector("list", length(unique(timeAndPhones$Phone)))
  for (i in 1:length(listTimeAndPhones)){
    listTimeAndPhones[[i]] <- timeAndPhones[timeAndPhones$Phone==unique(timeAndPhones$Phone)[i],]
    print(dim(listTimeAndPhones[[i]]))
  }
  save(listTimeAndPhones, file = "listTimeAndPhones.Rdata")
  #load("listTimeAndPhones.Rdata")
  rm(findDoc,findShift,i)

# And... record amount of tasks associated to each ward for each time and phone
  wardNames<-names(table(hanData$ward))
  candidateWards <- matrix(0,nrow(timeAndPhones),length(wardNames)+2)
  colnames(candidateWards) <- c(wardNames,"Other","Nurse Coord.")
  findWards <- function(index){
    dummy <- timeAndPhones[index,]
    if(dummy$doctor == "none"){
      return( rep(0,28) ) # No tasks anywhere... if approxShift is not active then we'll have to force Coordinator room
    } else {
      #Check the doctor or phone tasks whose range includes the time in dummy
      dummy2 <- hanData[hanData$phone == dummy$Phone & hanData$name == dummy$doctor
                        & hanData$timeof_max_accept <= dummy$time & dummy$time <= hanData$timeof_complete,]
      if (nrow(dummy2) == 0){
        return(rep(0,28)) # No tasks anywhere
      } else {
        wardTasks <- table(dummy2$ward)
        dummyReturn <- rep(0,28)
        for (j in 1:length(wardTasks)){
          dummyReturn[which(c(wardNames)==names(wardTasks)[j])] <- as.numeric(wardTasks[j])
        }
        return( dummyReturn )
      }
    }
  }
  candidateWards[,]<-t(sapply(1:nrow(timeAndPhones),findWards)); rm(wardNames,findWards)
  sort(colMeans(candidateWards)) # The CSSU task was one of those (accept and complete in a milisecond, doctors annoying!)
  sort(apply(candidateWards,2,max))
  sort(table(hanData$ward))

# DEFINE SUPERCLUSTERS
  NurseCoord <- c("Nurse Coord.")
  NorthWest <- c("Loxley","Lister 1")
  NorthEast <- c("Winifred 2","Nightingale","Patience 1")
  Linden <- c("Linden Lodge")
  central <- c("CSSU","Morris","ACU")
  westGreen <- c("Hogarth","Carrell","Bramley","Gervis Pearson","Fraser","Papplewick")
  TogFletch <- c("Toghill","Fletcher")
  BerFleBee <- c("Beeston","Berman 1","Berman 2","Fleming")
  SeaNewSouSRU <- c("Newell","Seacole","SRU","Southwell")

  # Dock is not anywhere to be found so we ignore it!

# New Candidate Wards matrix based on Clusters... keep old for AP selection...
  candidateWards2 <- matrix(0, nrow(candidateWards),9)
  colnames(candidateWards2) <- c("NorthWest","NorthEast","Linden","Central","SouthWest","South","SouthEast 1","SouthEast 2","NurseCoord")
    candidateWards2[,1] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% NorthWest)])
    candidateWards2[,2] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% NorthEast)])
    candidateWards2[,3] <- candidateWards[,which(colnames(candidateWards) %in% Linden)]
    candidateWards2[,4] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% central)])
    candidateWards2[,5] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% westGreen)])
    candidateWards2[,6] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% TogFletch)])
    candidateWards2[,7] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% BerFleBee)])
    candidateWards2[,8] <- rowSums(candidateWards[,which(colnames(candidateWards) %in% SeaNewSouSRU)])
    candidateWards2[,9] <- candidateWards[,which(colnames(candidateWards) %in% NurseCoord)]

# Split in a list by phone and save the object
  listCandidateWards <- vector("list", length(unique(timeAndPhones$Phone)))
  for (i in 1:length(listCandidateWards)){
    listCandidateWards[[i]] <- candidateWards2[timeAndPhones$Phone==unique(timeAndPhones$Phone)[i],]
    print(dim(listCandidateWards[[i]]))
    print(dim(listTimeAndPhones[[i]]))
  }
  save(listCandidateWards, file = "listCandidateWards.Rdata")
  #load("listCandidateWards.Rdata")

# Now we aim do something similar for all hotSpots (binary variables, observed YES/NO) and Signals!
  hotSpots <- unique(scansData$ap_id)
  observedHotSpots <- matrix(NA,nrow(timeAndPhones),length(hotSpots))
  observedSignals <- matrix(NA,nrow(timeAndPhones),length(hotSpots))
  colnames(observedHotSpots) <- hotSpots
  colnames(observedSignals) <- hotSpots
  
# First add phone names to simplify things later...
  phoneNames$phone_name<-as.numeric(phoneNames$phone_name); phoneNames$phone_imei<-as.numeric(phoneNames$phone_imei)
  scansData$phone<-rep(NA,nrow(scansData)) # We add phone names to scans Data to speed the things that will follow!
  replaceIMEI <- function(index){
    print(index)
    return(phoneNames$phone_name[phoneNames$phone_imei==scansData$imei[index]])
  }
  scansData$phone <- sapply(1:nrow(scansData), replaceIMEI) 

# Now get the observed hotspots and signals
  findHotSpots <- function(index){
    dummy <- timeAndPhones[index,]
    #Check all scans for index phone in a [time-30,time+30] range of the index time
    dummy2 <- scansData[ scansData$time <= dummy$time + 30 & scansData$time > dummy$time-30 
                         & scansData$phone == dummy$Phone,]
    # scansData is sorted by time and RSSI... we need to take either best/random RSSI for each time/hotspot (Try BEST RSSI)
    dummy2 <- dummy2[!duplicated(dummy2$ap_id),] # Remove duplicates and retain best RSSI
    hotSpotDummy <- dummy2$ap_id # Save names of observed hotSpots
    obsBinary <- as.numeric(hotSpots %in% hotSpotDummy) # Get the binary in the table, 1 if observed
    hotSpotDummyIndex <- sapply(hotSpotDummy, function(x) which(hotSpots == x)) # Get the indexes in the table for the observed hotSpots
    obsSignal <- rep(NA,length(hotSpots))
    if (length(hotSpotDummyIndex)>0){
      obsSignal[hotSpotDummyIndex] <- dummy2$rssi # Save the signals in the precise indexes within the table
    }  
    toReturn<-cbind(obsBinary,obsSignal)
    print(index)
    return(toReturn)
  }  
  
# ScansData is soooo big that this will take ages! aaaaages! so let's try snow library...
  library(snow)
  c1<-makeCluster(6,type="SOCK") #This machine has got 8 cores, check parallel.detectCores()
  clusterExport(c1,list("timeAndPhones","scansData","hotSpots"))
  a<-proc.time()
  hugeDummyTable<-t(parSapply(c1,1:nrow(timeAndPhones),findHotSpots)) # This will take an hour at least! get ready
  proc.time() - a
  stopCluster(c1)
  observedHotSpots[,] <- hugeDummyTable[,1:439]
  observedSignals[,] <- hugeDummyTable[,440:878]  
  sort(colSums(observedHotSpots))
  #write.csv(observedHotSpots,"observedHotSpotsFull.csv",row.names=F) #Just keep it full, just in case
  #write.csv(observedSignals,"bestObservedSignalsFull.csv",row.names=F)

  observedSignals <- read.csv("bestObservedSignalsFull.csv",header=T)
  observedHotSpots <- read.csv("observedHotSpotsFull.csv",header=T)  

  
# EACH WARD WILL NOMINATE most relevant %... of APs! 
  dummy3 <- c()
  for (i in 1:(ncol(candidateWards)-2)){
    print(i)
    idddx <- which(sort(colSums(observedHotSpots[candidateWards[,i]>0,]),decreasing=T)[1:5]>50)
    dummy3<-c(dummy3,names(sort(colSums(observedHotSpots[candidateWards[,i]>0,]),decreasing=T)[idddx]))
    print(mean(sort(colSums(observedHotSpots[candidateWards[,i]>0,]),decreasing=T)[idddx]))
    idddx <- which(sort(colSums(observedHotSpots[candidateWards[,i]>1,]),decreasing=T)[1:10]>100)
    dummy3<-c(dummy3,names(sort(colSums(observedHotSpots[candidateWards[,i]>1,]),decreasing=T)[idddx]))
    print(mean(sort(colSums(observedHotSpots[candidateWards[,i]>0,]),decreasing=T)[idddx])) 
  }
  dummy3<-unique(dummy3)

  #Check out test run..
  iepe<-read.csv("../citySurvey260116/2016-01-26-153254.csv",header=F)
  iepe<-unique(as.character(iepe$V2))
  for (i in 1:length(iepe)) iepe[i]<-paste("X",iepe[i],sep="")
  for (i in 1:length(iepe)) iepe[i]<-paste(strsplit(iepe[i],split = ":")[[1]][-6],collapse = ".")
  iepe<-unique(iepe)
  iepe%in%dummy3

  #Add based on our scans... pinchin and iker, since we went one whole year later!
  scanFiles <- paste("../citySurvey260116/",list.files(path = "../citySurvey260116"),sep="")
  dummy3new <- c()
  for (file in scanFiles[2:(length(scanFiles)-2)]){
    dummyScans <- read.csv(file,header=F)[,c(2)]
    dummyScans<- as.character(dummyScans)
    for (i in 1:length(dummyScans)){
      dummyScans[i] <- paste(strsplit(dummyScans[i],split=":")[[1]][-6],collapse=".")
      if (!is.na(as.numeric(strsplit(dummyScans[i],split="\\.")[[1]][1]))){
        dummyScans[i] <- paste(c("X",dummyScans[i]),collapse="")
      }
    }
    dummyScans <- unique(dummyScans)
    dummy3new<-c(dummy3new,names(observedHotSpots)[which(names(observedHotSpots) %in% dummyScans)])
  }
  dummy3new<-unique(dummy3new)
  dummy3new<- names(sort(colSums(observedHotSpots[,names(observedHotSpots)%in%dummy3new]),decreasing=T)[1:50])
  
  dummy3 <- c(dummy3,dummy3new)
  dummy3 <- unique(dummy3)
  
  #Continue...
  idx<-which(colnames(observedHotSpots) %in% dummy3)
  observedHotSpots <- observedHotSpots[,idx]
  observedSignals <- observedSignals[,idx]
  sort(colSums(observedHotSpots))  
  
  
  #Split in list and save objects by phone!
  listObservedHotSpots <- vector("list", length(unique(timeAndPhones$Phone)))
  listObservedSignals <- vector("list", length(unique(timeAndPhones$Phone)))
  for (i in 1:length(listObservedHotSpots)){
    listObservedHotSpots[[i]] <- observedHotSpots[timeAndPhones$Phone==unique(timeAndPhones$Phone)[i],]
    listObservedSignals[[i]] <- observedSignals[timeAndPhones$Phone==unique(timeAndPhones$Phone)[i],]
    print(dim(listObservedHotSpots[[i]]))
    print(dim(listTimeAndPhones[[i]]))
  }
  save(listObservedHotSpots, file = "listObservedHotSpots.Rdata")
  save(listObservedSignals, file = "listObservedSignals.Rdata")

  
  #################################
  
  inCoordRoom <- timeAndPhones$approxShiftActive=="no"
  plot(colMeans(observedHotSpots[inCoordRoom,]),type="h",ylim=c(0,1))
  lines(1:ncol(observedHotSpots)+0.15,colMeans(observedHotSpots[candidateWards2[,2]>0,]),type="h",col="blue") 
  lines(1:ncol(observedHotSpots)+0.30,colMeans(observedHotSpots[candidateWards2[,4]>0,]),type="h",col="red") 
  lines(1:ncol(observedHotSpots)+0.45,colMeans(observedHotSpots[candidateWards2[,5]>0,]),type="h",col="brown") 
  
  iep <- observedSignals[inCoordRoom,1] 
  hist(iep,30,freq=FALSE)

  