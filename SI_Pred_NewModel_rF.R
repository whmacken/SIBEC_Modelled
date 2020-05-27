####This script contains the working steps for estimating SI values. This includes modelling aSMR from rSMR,
####building polynomial models of standardised SI over edatopic grids, and predicting max SI based on temperature data
###Kiri Daust, July 2018

###Model aSMR <-> rSMR
library(tcltk)
library(reshape2)
library(foreach)
library(ggplot2)
library(dplyr)
library(magrittr)
library(data.table)
library(gridExtra)
library(mgcv)
library(rgl)
library(deldir)
library(openxlsx)
library(stringr)

changeNames <- function(x, old, new){
  result <- vector("numeric", length(x))
  for (i in 1:length(x)){
    code <- x[i]
    index <- match(code, old)
    result[i] <- new[index]
  }
  return(result)
}

###########Create full SIBEC Predict Data set from BART model###########
setwd("C:/Users/Kiri Daust/Desktop/SI-Prediction")
sibecOrig <- read.csv("SIBECDatRaw.csv")
climDat <- fread("BECv11_750Pt_FiveYr_2011_2015MSY.csv", data.table = F)
climDat <- climDat[climDat$Eref_sm > -1000,]
climDat <- climDat[,c("ID2","CMD","FFP","Eref_sm","MCMT","SHM","TD","PAS","DD5_sp","MSP")]
colnames(climDat)[1] <- "BGC"
climAve <- climDat %>%
  group_by(BGC) %>%
  summarise_all(mean, na.rm = T)
CMD <- climAve[,c("BGC","CMD")]

for (i in 1:3){
  CMD[,2+i] <- CMD[,1+i]/2
}
colnames(CMD) <- c("BGC","4","5","6","7")
CMD <- CMD[,c(1,3:5,2)]

###for each drier rSMR, previous CMD + 100
for (i in 1:4){
  CMD[,length(CMD)+1] <- CMD[,length(CMD)] + 100
}
colnames(CMD)[6:9] <- c("3","2","1","0")

CMD <- CMD[,order(colnames(CMD))]## creates full grid of CMD values by BGC by rSMR class
CMD <- as.data.frame(CMD)
CMD <- melt(CMD)
colnames(CMD) <- c("BGC","rSMR","EdaCMD")

eda <- read.csv("Edatopic_v11_7.csv")
eda <- eda[is.na(eda$Codes),]
eda <- eda[,-c(5,6)]
climAve <- merge(eda[,-1], climAve, by.x = "MergedBGC", by.y = "BGC", all = T)
climAve <- climAve[!is.na(climAve$CMD),]
climAve$SNR <- str_sub(climAve$Edatopic, 1,1) %>% str_replace_all(c("A" = "1", "B" = "2", "C" = "3", "D" = "4" ,"E" = "5" ))
climAve$rSMR <-str_sub(climAve$Edatopic, -1,-1)
climAve <- merge(climAve, CMD, by.x = c("MergedBGC","rSMR"), by.y = c("BGC","rSMR"), all = TRUE)

suitTbl <- read.csv("TreeSpp_ESuit_v11_18.csv")
suitTbl <- suitTbl[,c("Unit","Spp","ESuit")]
climAve <- merge(climAve, suitTbl, by.x = "SS_NoSpace", by.y = "Unit", all.x = T, all.y = F)
climAve <- climAve[climAve$Spp != "X",] ###tthis is bad, use to set bounds
climAve <- climAve[climAve$ESuit != 5,]
climAve <- climAve[!is.na(climAve$SS_NoSpace),]

##add in aSMR modelled data
aSMRrSMR="modelled_ALLv11_rSMR_aSMR_grid_HalfStep.csv"
aSMR <-read.csv(aSMRrSMR,stringsAsFactors=FALSE,na.strings=".")
aSMR <- aSMR[-1]
colnames(aSMR) <- c("BGC","0","1","2","3","4","5","6","7")
aSMR.list <- melt(aSMR)
colnames(aSMR.list)[2:3] <- c("rSMR","aSMR")
climAve <- merge(climAve,aSMR.list, by.x = c("MergedBGC","rSMR"),by.y = c("BGC","rSMR"), all.x = TRUE)
climAve <- climAve[!is.na(climAve$SS_NoSpace),]
SSVars <- climAve %>%
  group_by(SS_NoSpace) %>%
  summarise_at(c("aSMR","SNR","EdaCMD"), .funs = list(~min(.),~mean(as.numeric(.),na.rm = T),~max(.))) %>%
  ungroup()
climVars <- climAve[,c("SS_NoSpace","MergedBGC","Spp","ESuit","CMD","FFP","Eref_sm","MCMT","SHM","TD","PAS","DD5_sp","MSP")] %>%
  unique()
climVars <- climVars[!is.na(climVars$Spp),]
finalDat <- merge(climVars,SSVars,by = "SS_NoSpace", all.x = T)
colnames(finalDat) <- c("Unit", "BGC", "Spp", "ESuit", "CMD", "FFP", "Eref_sm", 
                        "MCMT", "SHM", "TD", "PAS", "DD5_sp", "MSP", "aSMR_min", "SNR_num_min", "EdaCMD_min",  
                        "aSMR_mean","SNR_num_mean", "EdaCMD_mean", "aSMR_max", "SNR_num_max","EdaCMD_max")
finalDat[,-(1:3)] <- mutate_all(finalDat[,-(1:3)], as.numeric)
write.csv(finalDat,"BartInputDat.csv", row.names = F)
finalDat <- read.csv("BartInputDat.csv")
######################################################################################################

####Create SI models for each species####
###__________________________________________________________________#####

#####Pull in climatic and site data for each siteseries species in sibecOrig (pull from suitability script)
SSSuit_Data <- finalDat
#Harmonize the species labelling in both data sets. Sx for all interior spruce, Pl, Cw, Fd
SSSuit_Data$Spp <- dplyr::recode(SSSuit_Data$Spp, Cwc="Cw", Fdc = "Fd", Se = "Sx", Sw = "Sx", Plc = "Pl", Bgc = "Bg")
sibecOrig$TreeSpp <- dplyr::recode(sibecOrig$TreeSpp, Cwc="Cw", Fdc = "Fd", Se = "Sx", Sw = "Sx", Plc = "Pl", Bgc = "Bg")

########## create 2 merged data sets. One where site index data is available and one where it is not.
# SS_wSibec <- sibecOrig [(sibecOrig$PlotCountSpp > 0),]
# SS_wSibec <- na.omit (SS_wSibec)

SS_wSibec <- sibecOrig
SS_SI <- SS_wSibec [,c(9,5,7) ]
colnames (SS_SI) [1:3] <- c("Unit", "Spp", "SI")
SS_SI2 <- merge (SS_SI, SSSuit_Data, by = c("Unit","Spp"))
SS_SI2 <- distinct(SS_SI2)
Spp_count2 <- SS_SI2 %>%
  group_by(Spp) %>%
  summarise(count = n())


Spp.list <- Spp_count2$Spp[Spp_count2$count > 10] %>% as.character()###list of species where there are more than 5 values
SS_SI2 <- SS_SI2[!is.na(SS_SI2$SI),]
###No sibec data - SHould be all site series. Where species does not occur then Site Index is 0
# SS_noSibec <- sibecOrig [is.na(sibecOrig$PlotCountSpp),]
# SS_SInew <- SS_noSibec [,c(9,5,7) ]
# colnames (SS_SInew) [1:3] <- c("Unit", "Spp", "SI")
# SS_SInew2 <- merge (SS_SInew, SSSuit_Data, by = c("Unit","Spp"))
# SS_SInew2 <- distinct(SS_SInew2)
# ###ignore warning
# #######################
# ###Then build regression model for each species and predict SI for places where values

library(rJava)
options(java.parameters = "-Xmx10000m")
library(bartMachine)
set_bart_machine_num_cores(5)


BFMods <- foreach(Spp = Spp.list, .combine = c)  %do% {
  SI_Spp <- SS_SI2 [(SS_SI2$Spp %in% Spp), ]
  SI_Spp <- SI_Spp[,-c(4,1,2,5)]
  rownames(SI_Spp) <- NULL
  
  bF.fit <- bartMachine(SI_Spp[,-1],SI_Spp$SI, k = 2, nu = 3, q = 0.99, num_trees = 200)
  list(bF.fit)
}

names(BFMods) <- Spp.list
save(BFMods, file = "BFModels.RData")

testDat <- unique(finalDat)
bartSIPred <- foreach(Spp = Spp.list, .combine = rbind) %do% {
  szUns <- unique(testDat$BGC[testDat$Spp == Spp]) %>% as.character()
  testAll <- testDat[testDat$BGC %in% szUns,-(2:4)] %>% unique()
  mod <- BFMods[[Spp]]
  test <- testAll[,mod$training_data_features]
  out <- data.frame(Unit = testAll$Unit, aSMR = testAll$aSMR_mean, SIPred = predict(mod, new_data = test))
  out$Spp <- Spp
  ##cbind(out[,c("SIPred","Spp")],testAll)
  out
}

bartSIPred <- merge(bartSIPred, sibecOrig[,c(9,5,7)], by.x = c("Unit","Spp"), by.y = c("SS_NoSpace","TreeSpp"), all.x = T)
bartSIPred <- merge(bartSIPred, suitTbl, by = c("Unit","Spp"), all.x = T)
bartSIPred <- unique(bartSIPred)
write.csv(bartSIPred[,c("Unit","Spp","SIPred","aSMR")], "BartPredSI.csv", row.names = F)

SI_Spp$PredSI <- predict(rF.fit, SI_Spp[,-c(1)])
SI_Pred <- merge(SI_Units, SI_Spp[,c(1,20)], by = 0)

SI_Sppnew <- SS_SInew2 [(SS_SInew2$Spp %in% Spp), ]
SI_Unitsnew <- SI_Sppnew[, c(4,1,2,5)]
SI_Sppnew <- SI_Sppnew[, -c(4,1,2,5)]
SI_Sppnew$PredSI <- predict(rF.fit, SI_Sppnew[,-c(1)])
SI_Prednew <- merge(SI_Unitsnew, SI_Sppnew[,c(1,20)], by = 0)
SI_Pred2 <- rbind (SI_Pred, SI_Prednew)
SI_Pred2

####model testing
SI_PredAll <- foreach(Spp = Spp.list, .combine = rbind)  %do% {
  options(stringsAsFactors = FALSE)
  SI_Spp <- SS_SI2 [(SS_SI2$Spp %in% Spp), ]
  SI_Spp <- SI_Spp[,-c(4,1,2,5)]
  rownames(SI_Spp) <- NULL
  ind <- sample(rownames(SI_Spp), size = nrow(SI_Spp)/3)
  test <- SI_Spp[ind,]
  train <- SI_Spp[!rownames(SI_Spp) %in% ind,]
  
  # rF.fit <- randomForest(SI ~ .,data=train, nodesize = 2, 
  #                        do.trace = 10, ntree=201, na.action=na.fail, importance=TRUE, proximity=TRUE)
  
  bF.fit <- bartMachine(train[,-1],train$SI, k = 2, nu = 3, q = 0.99, num_trees = 200)
  
  test$BFPred <- predict(bF.fit, test[,-1])
  test$Spp <- Spp
  test
}


##################################################################################################################

  ############test caret ensemble
  
  control <- trainControl(method="cv", number=5, returnResamp = "all",
                          classProbs = TRUE, 
                          search = "random")#, repeats=3
  set.seed (12345)
  #metric <- "ROC"
  #C5.grid <- expand.grid(.cp=0)
  # C5.0
  droplevels(X1.sub$ESuit)
  X1.sub$ESuit <- as.factor(X1.sub$ESuit)
  #fit.c50 <- train(ESuit ~., X1.sub, method="C5.0", metric= "Accuracy", trControl=control)#[,c(1:4)]
  
  # Stochastic Gradient Boosting
  set.seed (12345)
  fit.rf <- train(ESuit ~.,  X1.sub, method='rf', metric="Accuracy", trControl=control, verbose=FALSE,  do.trace = 10, ntree=101)
  varImp(fit.rf)
  
  # summarize results
  #boosting_results <- resamples(list(c5.0=fit.c50, rf=fit.rf))
  #summary(boosting_results)
  #dotplot(boosting_results)
  X1.sub$Pred <- predict(fit.rf, X1.sub[,-c(1)])
  confusionMatrix(X1.sub[,1],X1.sub$Pred)
  ###show confusion matrix of model
  X1.sub2$Pred <- predict(fit.rf, X1.sub2[,-c(1)])
  confusion <- as.matrix (confusionMatrix(X1.sub2[,1],X1.sub2$Pred))
  confusionMatrix(X1.sub2[,1],X1.sub2$Pred)
  #write.csv(confusion, file= paste(Spp,"_ConfusionMatrix.csv", sep=""))
  SppPredict <- merge(X1.unit, X1.sub2, by = 0)
  SppPredict <- na.omit(SppPredict)
  ########output the units that are misclassified
  #compareC5 <- SUsumMatrix[,c(1,length(SUsumMatrix),2)]
  SppPredict$Same <- ifelse(SppPredict$ESuit == SppPredict$Pred,1,0)
  SppPredict <- SppPredict [,c("Unit", "ESuit", "Pred", "Same")]
  SppPredict<- SppPredict[!(SppPredict$ESuit == "E5" & SppPredict$Pred == "E5"),]
  write.csv(SppPredict, file= paste(Spp, "_",State, "_ESuit_C50Suit.csv", sep=""))
  
  DiffSuit <- SppPredict[SppPredict$Same == 0,]
  DiffSuit <- DiffSuit [,c("Unit", "ESuit", "Pred")]
  #write.csv(DiffSuit, file= paste(Spp,"_",State,"_ESuit_C50DiffOnly.csv", sep=""))
  
  ESuitSpp <- ESuit [(ESuit$Spp %in% Spp),]
  ESuitSpp <- merge (ESuitSpp[-4],SppPredict, by = "Unit")
  ESuitSppnew <-""
  
  
  SuitPred <- ESuitSpp ### Use where only the BGCs with previously estimated ESUIT are included
  EstSuit <- rbind (EstSuit, SuitPred)
  EstSuit
  
} # ignore error message here

Date <- Sys.Date()
SuitCompare$Table <- ""
SuitCompare$Table [SuitCompare$Unit %in% ESuit1$Unit] <- "BC_Banner"
SuitCompare$Table [SuitCompare$Unit %in% ESuit2$Unit] <- "USA_Meidinger"
SuitCompare$Table [SuitCompare$Unit %in% ESuit3$Unit] <- "Alberta_Kabzims"
SuitCompare$Table [SuitCompare$Unit %in% ESuit4$Unit] <- "Nelson_AddedPred"
#SuitCompare$Table [(SuitCompare$Table == "" )] <- "Nelson_Refguide"
SuitCompare <- SuitCompare %>% distinct()
SuitCompare <- na.omit(SuitCompare)
SuitCompare2 <- merge (SuitCompare, Codes, by = "Unit")
write.csv(SuitCompare2, file= paste("ComparisonESuit_rF_1", Date, ".csv", sep=""), row.names = FALSE)






















#########Old method from Kiri
 modelList <- foreach(Spp = SppList, .combine = c) %do%
  sibec <- sibecOrig[sibecOrig$TreeSpp == Spp,]##choose species
  sibec <- sibec[!is.na(sibec$BGCUnit),]
  numPlots <- aggregate(MeanPlotSiteIndex ~ BGCUnit, sibec, FUN = length)
  if(Spp == "Lw" | Spp == "Py"){
    cutoff <- 2
  }else{
    cutoff <- 4
  }
  BGC <- numPlots$BGCUnit[numPlots$MeanPlotSiteIndex > cutoff]###remove SI values in zones with little data
  sibec <- sibec[sibec$BGCUnit %in% BGC,]
  
  ###Loop to match SI values with edatopes and convert to aSMR
  out <- foreach(i = unique(sibec$BGCUnit), .combine = rbind) %do% {
    SIsub <- sibec[sibec$BGCUnit == i,c(9,7)]
    edaSub <- Eda[Eda$MergedBGC == i, 3:4]
    if(length(edaSub$SS_NoSpace) > 0){
      SIsub <- merge(SIsub, edaSub, by = "SS_NoSpace", all.x = TRUE)
      if(!any(is.na(SIsub$Edatopic))){
        SIeda <- aggregate(MeanPlotSiteIndex ~ Edatopic, data =  SIsub, FUN = mean)
        colnames(SIeda)[2] <- "SI"
        SIeda$rSMR <- gsub("[[:alpha:]]","", SIeda$Edatopic)
        SIeda$SNR <- gsub("[[:digit:]]","", SIeda$Edatopic)
        SIeda$rSMR <- as.numeric(as.character(SIeda$rSMR))
        SMRsub <- SMRCross[SMRCross$BGC == i,-1]
        SIeda <- merge(SIeda, SMRsub, by = "rSMR", all.x = TRUE)
        temp <- SIeda[SIeda$aSMR %% 1 != 0,] ###if 0.5, give to higher and lower whole num
        temp$aSMR <- ceiling(temp$aSMR)
        SIeda$aSMR <- floor(SIeda$aSMR)
        SIeda <- rbind(SIeda,temp)
        SIeda <- SIeda[,c(4,5,3)]
        SIeda$SI <- SIeda$SI/max(SIeda$SI)###Standardise out of max SI value in BGC
        SIeda$BGC <- i
        SIeda
      }
      
    }
    
  }
  
  SI2 <- out
  SI2 <- dcast(SI2, SNR + aSMR ~ BGC, value.var = "SI", fun.aggregate = mean) ##Convert to matrix
  SI2$NumPl <- apply(SI2[,-c(1:2,62)], 1, function(x){length(x[!is.na(x)])}) ###How many values in each edatopic position?
  SI2 <- SI2[SI2$NumPl >= 4,] ##Remove rows with few values (has to be lower for Lw and Py)
  SI2 <- SI2[,-length(SI2)]
  ##len <- length(SI2)/5
  ##SI2[,1:floor(len)][SI2[,1:floor(len)] == "NaN"] <- 0
  SI2[SI2 == "NaN"] <- NA
  SI2$MeanSI <- rowMeans(SI2[,-(1:2)], na.rm = TRUE)##mean for each edatopic cell
  #SI2$SESI <- apply(SI2[,-c(1:2)],1,sd,na.rm = TRUE)
  
  SIcomb <- SI2[,c(1,2,length(SI2))]
  SImat <- dcast(SIcomb, SNR ~ aSMR)
  SImat <- SImat[,-1]
  SImat <- as.matrix(SImat)
  barplot3d(SImat)### plot distribution
  
  SIcomb$SNR <- changeNames(SIcomb$SNR, old = c("A","B","C","D","E"), new = c(1,2,3,4,5))
  SIcomb$SNR <- as.numeric(SIcomb$SNR)
  fit.poly <- lm(MeanSI ~ poly(SNR, aSMR, degree = 4, raw = TRUE), data = SIcomb)###Create polynomial model of degree 4 (create seperate model for each species)
  tempList <- list()
  tempList[[Spp]] <- fit.poly
  tempList
 {}

###Plot polynomial model####
#zfit <- matrix(fitted(fit.bl), ncol = 8)
#surface3d(1:9, 1:5, zfit, col = "purple")
#aspect3d(1,1,1)

####################################################
####Relationship with Temperature#####################
#####Determine equation relating temp and SI####
#############################################

dat <- fread("ALLv11_500Pt_Normal_1961_1990MSY_REDUCED.csv", data.table = FALSE)
#dat <- dat[,c("ID2","PPT_at","PPT_wt","PAS", "Tave_sp", "Tave_sm", "DD5", "MAP", "MAT")]
dat$ID2 <- gsub("[[:space:]]","",dat$ID2)
#mods <- list(Pl = fit.pl, Sx = fit.sx, Fd = fit.fd, Bl = fit.bl, Lw = fit.lw, Cw = fit.cw) ##list of polynomial models
mods <- modelList ##list of polynomial models

edaPos <- list(A = c("C",5),B = c("B",3),C = c("D",6)) ###Which edatopic positions?

###Calculate slopes and intercepts
slopes <- foreach(Spp = SppList, .combine = rbind) %do% {##foreach species
  sibec <- sibecOrig[sibecOrig$TreeSpp == Spp,]
  sibec <- sibec[!is.na(sibec$BGCUnit),]
  
  foreach(eda = edaPos, .combine = rbind) %do% { ##foreach edatopic position
    ddBGC <- aggregate(cbind(DD5,Tave_sm, Tave_sp) ~ ID2, dat, FUN = mean)#MAT,
    colnames(ddBGC)[1] <- "BGC"
    
    SNR.val <- eda[1]
    aSMR.val <- eda[2]
    
    out <- foreach(i = unique(sibec$BGCUnit), .combine = rbind) %do% { ##foreach BGC
      SIsub <- sibec[sibec$BGCUnit == i,c(9,7)]
      edaSub <- Eda[Eda$MergedBGC == i, 3:4]
      if(length(edaSub$SS_NoSpace) > 0){
        SIsub <- merge(SIsub, edaSub, by = "SS_NoSpace", all.x = TRUE)
        SIsub <- SIsub[!is.na(SIsub$MeanPlotSiteIndex),]
        if(!any(is.na(SIsub$Edatopic))){
          SIeda <- aggregate(MeanPlotSiteIndex ~ Edatopic, data =  SIsub, FUN = mean)
          colnames(SIeda)[2] <- "SI"
          SIeda$rSMR <- gsub("[[:alpha:]]","", SIeda$Edatopic)
          SIeda$SNR <- gsub("[[:digit:]]","", SIeda$Edatopic)
          SIeda$rSMR <- as.numeric(as.character(SIeda$rSMR))
          SMRsub <- SMRCross[SMRCross$BGC == i,-1]
          SIeda <- merge(SIeda, SMRsub, by = "rSMR", all.x = TRUE)
          temp <- SIeda[SIeda$aSMR %% 1 != 0,]
          temp$aSMR <- ceiling(temp$aSMR)
          SIeda$aSMR <- floor(SIeda$aSMR)
          SIeda <- rbind(SIeda,temp)
          SIeda <- SIeda[,c(4,5,3)]
          out <- mean(SIeda[SIeda$SNR == SNR.val & SIeda$aSMR == aSMR.val,3]) ##mean for each cell
          out <- as.data.frame(out)
          out$BGC <- i
          out
        }
        
      }
      
    }
    
    out$BGC <- gsub("[[:space:]]","",out$BGC)
    ddBGC$BGC <- gsub("[[:space:]]","",ddBGC$BGC)
    ddBGC <- merge(ddBGC, out, by = "BGC", all.y = TRUE)
    ddBGC <- ddBGC[!ddBGC$BGC %in% c("CWHvh2", "CWHvm1","CWHvh1","CWHws1"),] ###These coastal units are weird with Pl
    ddBGC <- ddBGC[!is.na(ddBGC$out),]
    
    if(nrow(ddBGC) > 3){
      plot(out ~ Tave_sm, data = ddBGC)
      abline(lm(out ~ Tave_sm, data = ddBGC))
      fit <- lm(out ~ Tave_sm, data = ddBGC)##linear model
      temp <- data.frame(Species = Spp, Eda = paste(eda[1],eda[2],sep = "/"), Intercept = fit$coefficients[1], ##extract slope and intercept for each species
                         Slope = fit$coefficients[2], minSI = min(ddBGC$out), maxSI = max(ddBGC$out),
                         minTemp = min(ddBGC$Tave_sm), maxTemp = max(ddBGC$Tave_sm))
      temp
    }
    
  }
}


###############Fill in SI Using polynomial model and slopes##############################
#########################################################################################

#temp <- read.csv(file.choose())###Need to select BGC Units  (currently reading in Fraser TSA predictions)
temp <- fread("FraserTSA_SSpredicted.csv")
missing <- unique(as.character(temp$SS_NoSpace)) 

dat <- fread("ALLv11_500Pt_Normal_1961_1990MSY_REDUCED.csv", data.table = FALSE) ###Climate data
#dat <- dat[,c(2,5,239,230,183,184)]
dat$ID2 <- gsub("[[:space:]]","",dat$ID2)

##this loops calculates SI values for each species in each unit and averages them by unit####
SIFill <- foreach(Spp = SppList, .combine = rbind) %do% {
  sibec <- sibecOrig[sibecOrig$TreeSpp == Spp,]
  sibec <- sibec[!is.na(sibec$TreeSpp),]
  
  climSlope <- slopes ###Or read from csv
  slopeSub <- climSlope[climSlope$Species == Spp & climSlope$Eda == "D/6",] ##Can choose which edatopic position to use
  
  fill <- foreach(SS = missing, .combine = rbind) %do% {
    SIeda <- Eda[Eda$SS_NoSpace == SS,-c(5,6)]
    if(length(SIeda$Edatopic) > 0){
      zone <- as.character(unique(SIeda$MergedBGC))
      zone <- gsub("[[:space:]]","",zone)
      
      SIeda$rSMR <- gsub("[[:alpha:]]","", SIeda$Edatopic)
      SIeda$SNR <- gsub("[[:digit:]]","", SIeda$Edatopic)
      SIeda$rSMR <- as.numeric(as.character(SIeda$rSMR))
      SMRsub <- SMRCross[SMRCross$BGC == zone,-1]
      SIeda <- merge(SIeda, SMRsub, by = "rSMR", all.x = TRUE) ##convert to aSMR
      datNew <- SIeda[,c(6,7)]
      datNew$SNR <- changeNames(datNew$SNR, old = c("A","B","C","D","E"), new = c(1,2,3,4,5))
      if(length(datNew$SNR) < 2){ ##predict doesn't work with only one row
        datNew <- datNew[c(1,1),]
      }
      datNew$MeanSI <- predict(modelList[[Spp]], newdata = datNew) ##Predict
      datNew$MeanSI[datNew$MeanSI > 1] <- 1 ###don't want proportions > 1
      SIprop <- mean(datNew$MeanSI)
      SImin <- min(datNew$MeanSI)
      SImax <- max(datNew$MeanSI)
      SIsub2 <- sibec[sibec$BGCUnit == zone,]
      SIsub2 <- SIsub2[!is.na(SIsub2$TreeSpp),]
      if(length(SIsub2$TreeSpp) > 2){ ##If there's any SI values in that zone, use max of those
        maxSI <- max(SIsub2$MeanPlotSiteIndex)
      }else{###Otherwise calculate max SI from Tave_Sm and slopes
        climSub <- dat[dat$ID2 == zone,]
        tempVar <- mean(climSub$Tave_sm)
        maxSI <- slopeSub$Slope*tempVar + slopeSub$Intercept
        if(length(maxSI) == 0){
          maxSI <- NA
        }
      }
      out <- data.frame(Unit = SS, SImin = maxSI*SImin, SImean = maxSI*SIprop, SImax = SImax*maxSI)##Multiply proportion by maxSI
      out
    }
    
  }
  
  sibec <- sibec[,c(9,7)]
  colnames(sibec)[1] <- "Unit"
  
  fill <- merge(fill, sibec, by = "Unit", all.x = TRUE) ##Merge actual SI values for comparison
  fill$Spp <- Spp
  fill
}

write.csv(SIFill,"PredSI_FraserTSA4Feb.csv", row.names = FALSE)

#####Now same as above but for each edatopic cell########################################
#########################################################################################
sibec <- read.xlsx("SIBEC_for_Portfolio2.xlsx", sheet = 1)
sibec <- sibec[sibec$TreeSpp == "Lw",]
sibec <- sibec[!is.na(sibec$TreeSpp),]
Spp <- "Lw"
SS="BWBSdk/101"
fill <- foreach(SS = missing, .combine = rbind) %do% {
  SIeda <- Eda[Eda$SS_NoSpace == SS,-c(5,6)]
  if(length(SIeda$Edatopic) > 0){
    zone <- as.character(unique(SIeda$MergedBGC))
    zone <- gsub("[[:space:]]","",zone)
    
    SIeda$rSMR <- gsub("[[:alpha:]]","", SIeda$Edatopic)
    SIeda$SNR <- gsub("[[:digit:]]","", SIeda$Edatopic)
    SIeda$rSMR <- as.numeric(as.character(SIeda$rSMR))
    SMRsub <- SMRCross[SMRCross$BGC == zone,-1]
    SIeda <- merge(SIeda, SMRsub, by = "rSMR", all.x = TRUE)
    datNew <- SIeda[,c(6,7)]
    datNew$SNR <- changeNames(datNew$SNR, old = c("A","B","C","D","E"), new = c(1,2,3,4,5))
    datNew <- unique(datNew)
    if(length(datNew$SNR) < 2){
      datNew <- datNew[c(1,1),]
    }
    datNew$MeanSI <- predict(mods[[Spp]], newdata = datNew) ##Predict
    datNew <- unique(datNew)
    SIsub2 <- sibec[sibec$BGCUnit == zone,]
    SIsub2 <- SIsub2[!is.na(SIsub2$TreeSpp),]
    if(length(SIsub2$TreeSpp) > 2){
      maxSI <- max(SIsub2$MeanPlotSiteIndex)
    }else{
      climSub <- dat[dat$ID2 == zone,]
      tempVar <- mean(climSub$Tave_sm)
      maxSI <- slopeSub$Slope*tempVar + slopeSub$Intercept
    }
    datNew$MeanSI <- datNew$MeanSI*maxSI
    datNew$Unit <- SS
    datNew
  }
  
}

######Create edatopic grids with predicted and actual SIs###########
for (BGC in zoneSelect){
  modDat <- fill[grep(BGC, fill$Unit),]
  SMRsub <- SMRCross[SMRCross$BGC == BGC,-1]
  modDat <- merge(modDat, SMRsub, by = "aSMR", all.x = TRUE)
  modDat <- modDat[,c(2,5,3)]
  modDat <- unique(modDat)
  modDat$SNR <- changeNames(modDat$SNR, old = c(1,2,3,4,5), new = c("A","B","C","D","E"))
  modDat$rSMR <- as.numeric(as.character(modDat$rSMR))
  SBsub <- sibec[sibec$BGCUnit == BGC,c(9,3,7)]
  SBsub <- unique(SBsub)
  SBsub <- merge(SBsub, Eda, by = "SS_NoSpace", all.x = TRUE)
  SBsub <- SBsub[,c(2,3,6)]
  SBsub$rSMR <- gsub("[[:alpha:]]","", SBsub$Edatopic)
  SBsub$SNR <- gsub("[[:digit:]]","", SBsub$Edatopic)
  SBsub$rSMR <- as.numeric(as.character(SBsub$rSMR))
  SBsub2 <- aggregate(MeanPlotSiteIndex ~ SNR +rSMR, SBsub, FUN = mean) ##average where SS overlap
  
  assign(paste("Plot",BGC,sep = ""),(ggplot(data = modDat)+
           geom_tile(data = SBsub, aes(x= SNR, y = (rSMR*-1), fill = SiteSeries))+
           geom_text(data = SBsub2, aes(x= SNR, y = (rSMR*-1)+0.15, label = round(MeanPlotSiteIndex, digits = 1)))+
           geom_text(aes(x= SNR, y = (rSMR*-1)-0.15, label = round(MeanSI, digits = 1)))+
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank())+
           coord_fixed()+
           labs(title = BGC))) 
}
library(gridExtra)
ordLay <- rbind(c(1,2),
                c(3,4))
grid.arrange(PlotICHmc2, PlotIDFdk3, PlotSBSmc2, PlotCWHmm1, layout_matrix=ordLay)

##Create edatopic grids just for predicted SI coloured by value##############
for (BGC in zoneSelect){
  modDat <- fill[grep(BGC, fill$Unit),]
  SMRsub <- SMRCross[SMRCross$BGC == BGC,-1]
  modDat <- merge(modDat, SMRsub, by = "aSMR", all.x = TRUE)
  modDat <- modDat[,c(2,5,3)]
  modDat <- unique(modDat)
  modDat$SNR <- changeNames(modDat$SNR, old = c(1,2,3,4,5), new = c("A","B","C","D","E"))
  modDat$rSMR <- as.numeric(as.character(modDat$rSMR))
  
  
  assign(paste("Plot",BGC,sep = ""), ggplot(data = modDat)+
                                       geom_tile(aes(x= SNR, y = (rSMR*-1), fill = MeanSI))+
                                       geom_text(aes(x= SNR, y = (rSMR*-1), label = round(MeanSI, digits = 1)))+
                                       scale_fill_gradient(low = "red", high = "green")+
                                       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank())+
                                       coord_fixed()+
                                       labs(title = BGC))
}
ordLay <- rbind(c(1,2),
                c(3,4))
grid.arrange(PlotICHmc2, PlotIDFdk3, PlotSBSmc2, PlotCWHmm1, layout_matrix=ordLay)

#####----------------------------------------------------------------------------------------###


portBGC <- read.csv(file.choose(), stringsAsFactors = FALSE)###Select units to calculate SI for
suitBC <- read.csv("BCSuit.csv")
zoneSelect <- unique(portBGC$SS_NoSpace)
zoneSelect <- c("ICHmc2","SBSmc2","IDFdk3","CWHmm1")
zoneSelect <- c("IDFdk3")
##missing <- inRef[!(inRef %in% inSI)]
missing <- zoneSelect

missing <- as.character(suitBC$SS_NoSpace[suitBC$ZoneSubzone %in% zoneSelect])##vector of units (with SS)
missing <- as.character(zoneSelect[zoneSelect %in% suitBC$SS_NoSpace])##vector of units (with SS)






####Working Code##################

###Code to test quality of aSMR crosswalk table 
expGrid <- read.csv("ExpertGrid.csv")
colnames(expGrid)[-1] <- paste(colnames(expGrid)[-1],"_E", sep = "")
expGrid$BGC <- gsub("[[:space:]]","",expGrid$BGC)
comp <- merge(expGrid, test, by = "BGC")

for(i in 1:8){
  comp[,length(comp) + 1] <- comp[,i+1] - comp[,i+9]
}

comp <- comp[,-(2:17)]
colnames(comp) <- colnames(test)
hist(comp$SMR4, col = "purple")

for(i in 1:8){
  hist(comp[,i+1], col = "purple", main = colnames(current)[i+1])
}


comp$Zone <- gsub("[[:lower:]]|[[:digit:]]","",comp$BGC)
zoneQual <- aggregate(SMR4 ~ Zone, comp, mean)
barplot(zoneQual$SMR4, names.arg = zoneQual$Zone)

compLong <- melt(comp[,-10])
colnames(compLong) <- c("BGC","rSMR","Diff")
compLong <- compLong[abs(compLong$Diff) >= 1,]
temp <- melt(expGrid)
colnames(temp) <- c("BGC","rSMR","Expert")
temp$rSMR <- gsub("r","",temp$rSMR)
compLong <- merge(compLong, temp, by = c("BGC","rSMR"), all.x = TRUE)

comp2 <- comp
comp2[comp2 == 0.5 | comp2 == -0.5] <- 0

for(j in unique(comp$Zone)){
  pdf(paste(j,".pdf", sep = ""))
  current <- comp[comp$Zone == j,]
  layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9)),widths=c(1,1,1), heights =c(1,1,1), respect=TRUE)
  par(mai = c(0.5, 0.5, 0.2, 0.2)) #speSEfies the margin size in inches
  for(i in 1:8){
    hist(current[,i+1], breaks = 7, col = "purple", main = colnames(current)[i+1])
  }
  dev.off()
}

##More graphs to look at temp###########
out$BGC <- gsub("[[:space:]]","",out$BGC)
ddBGC$BGC <- gsub("[[:space:]]","",ddBGC$BGC)
ddBGC <- merge(ddBGC, out, by = "BGC", all.y = TRUE)
ddBGC <- ddBGC[!ddBGC$BGC %in% c("CWHvh2", "CWHvm1","CWHvh1","CWHws1"),]
#ddBGC$out[is.na(ddBGC$out)] <- 0

plot(out ~ Tave_sm, data = ddBGC)
text(ddBGC$Tave_sm, ddBGC$out, labels = ddBGC$BGC, cex = 0.6, pos = 3)
abline(lm(out ~ Tave_sm, data = ddBGC))
fit <- lm(out ~ Tave_sm, data = ddBGC)

plot(out ~ Tave_sp, data = ddBGC)
text(ddBGC$Tave_sp, ddBGC$out, labels = ddBGC$BGC, cex = 0.6, pos = 3)
plot(out ~ MAT, data = ddBGC)
text(ddBGC$MAT, ddBGC$out, labels = ddBGC$BGC, cex = 0.7, pos = 3)
scatter.smooth(ddBGC$MAT, ddBGC$out)

Suit <- read.csv(file.choose())
suitSub <- Suit[Suit$Spp == "Pl",]
zones <- unique(suitSub$BGC)
zones <- as.character(zones)
zones <- gsub("[[:space:]]","",zones)
ddBGC$RefGuide <- "red"
ddBGC$RefGuide[ddBGC$BGC %in% zones] <- "green"

ggplot(ddBGC, aes(x = Tave_sm, y = out, colour = RefGuide, label = BGC))+
  geom_point()+
  geom_text_repel(size = 2.6)+
  scale_y_continuous(limits = c(-5,30))

####################See aSMR_X_rSMR script for this
###Create rSMR -> aSMR crosswalk######################
wd <- tk_choose.dir(); setwd(wd)

####This part includes removing monthly CMD based on surplus (doesn't make much difference)#####
allDat <- fread("ALLv11_500Pt_Normal_1961_1990MSY_REDUCED.csv", data.table = FALSE) ###Climate data
temp <- allDat[,grep("CMD",colnames(allDat))]
allDat <- allDat[,c("ID2","PPT_at","PPT_wt","PAS", "Tave_sp", "Tave_sm", "DD5", "MAP")]
allDat <- cbind(allDat,temp)
allDat$PPT.dorm <- allDat$PPT_at + allDat$PPT_wt
CMD <- aggregate( . ~ ID2, allDat, mean) ##
CMD$CMDKiri <- ifelse(CMD$PAS > 1000, CMD$CMD06+CMD$CMD07+CMD$CMD08+CMD$CMD09, 
                      CMD$CMD02+CMD$CMD03+CMD$CMD04+CMD$CMD05+CMD$CMD06+CMD$CMD07+CMD$CMD08+CMD$CMD09)
CMD$CMD <- CMD$CMDKiri
CMD <- CMD[,c("ID2","CMD","PPT.dorm")]
CMD$Def <- 275 - CMD$PPT.dorm
CMD$Def[CMD$Def < 0] <- 0
CMD$CMD <- CMD$CMD + CMD$Def
CMD <- CMD[,c("ID2","CMD")]
###_____________________________________________####

###Now just with ppt#######
allDat <- fread("ALLv11_500Pt_Normal_1961_1990MSY_REDUCED.csv", data.table = FALSE)
allDat <- allDat[,c("ID2","PPT_at","PPT_wt","CMD")]
allDat$PPT.dorm <- allDat$PPT_at + allDat$PPT_wt
CMD <- aggregate(cbind(PPT.dorm, CMD) ~ ID2, allDat, mean)###Mean by BGC

CMD$Def <- 250 - CMD$PPT.dorm ###275 seems to work well
CMD$Def[CMD$Def < 0] <- 0
CMD$CMD <- CMD$CMD + CMD$Def
CMD <- CMD[,c("ID2","CMD")]

###for each wetter rSMR, previous CMD is divided by 2
for (i in 1:3){
  CMD[,2+i] <- CMD[,1+i]/2
}
colnames(CMD) <- c("BGC","SMR4","SMR5","SMR6","SMR7")
CMD <- CMD[,c(1,3:5,2)]

###for each drier rSMR, previous CMD + 75
for (i in 1:4){
  CMD[,length(CMD)+1] <- CMD[,length(CMD)] + 75
}
colnames(CMD)[6:9] <- c("SMR3","SMR2","SMR1","SMR0")

CMD <- CMD[,order(colnames(CMD))]

#############################################Now calculate aSMR with ruleset
rules <- read.csv("aSMR_Rules.csv")

aSMRClass <- function(x){
  for(i in 1:length(ruleSelect$CMD)){
    if(x < ruleSelect$CMD[i+1]){
      x <- ruleSelect$aSMR[i]
      break
    }
  }
  return(x)
}

###Calculate values based on rules###
test <- foreach(SMR = colnames(CMD)[-1], .combine = cbind) %do% {
  temp <- CMD[,SMR]
  if(SMR == "SMR7"){
    ruleSelect <- rules[rules$SMRLevel == 7,-1]
  }else if(SMR == "SMR6"){
    ruleSelect <- rules[rules$SMRLevel == 6,-1]
  }else if(SMR == "SMR5"){
    ruleSelect <- rules[rules$SMRLevel == 5,-1]
  }else{
    ruleSelect <- rules[rules$SMRLevel == 0,-1]
  }
  out <- sapply(temp,FUN = aSMRClass)
  out
}


test <- as.data.frame(test)
test <- cbind(CMD$BGC, test)
colnames(test) <- colnames(CMD)
test$BGC <- gsub("[[:space:]]","",test$BGC)
SMRCross <- melt(test) ###aSMR lookup
