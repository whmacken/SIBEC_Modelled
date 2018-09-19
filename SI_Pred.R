####This script contains the working steps for estimating SI values. This includes modelling aSMR from rSMR,
####building polynomial models of standardised SI over edatopic grids, and predicting max SI based on temperature data
###Kiri Daust, July 2018

###Model aSMR <-> rSMR
.libPaths("E:/R packages")
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

changeNames <- function(x, old, new){
  result <- vector("numeric", length(x))
  for (i in 1:length(x)){
    code <- x[i]
    index <- match(code, old)
    result[i] <- new[index]
  }
  return(result)
}

####Function to create 3d barplots from Stack Overflow####################
binplot.3d <- function(x, y, z, alpha=1, topcol="#ff0000", sidecol="#aaaaaa", linecol="#000000")
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  x1 <- c(rep(c(x[1], x[2], x[2], x[1]), 3), rep(x[1], 4), rep(x[2], 4))
  z1 <- c(rep(0, 4), rep(c(0, 0, z, z), 4))
  y1 <- c(y[1], y[1], y[2], y[2], rep(y[1], 4), rep(y[2], 4), rep(c(y[1], y[2], y[2], y[1]), 2))
  x2 <- c(rep(c(x[1], x[1], x[2], x[2]), 2), rep(c(x[1], x[2], rep(x[1], 3), rep(x[2], 3)), 2))
  z2 <- c(rep(c(0, z), 4), rep(0, 8), rep(z, 8))
  y2 <- c(rep(y[1], 4), rep(y[2], 4), rep(c(rep(y[1], 3), rep(y[2], 3), y[1], y[2]), 2))
  rgl.quads(x1, z1, y1, col=rep(sidecol, each=4), alpha=alpha)
  rgl.quads(c(x[1], x[2], x[2], x[1]), rep(z, 4), c(y[1], y[1], y[2], y[2]), col=rep(topcol, each=4), alpha=1) 
  rgl.lines(x2, z2, y2, col=linecol)
}

barplot3d <- function(z, alpha=1, col="#ff0000", scale=1)
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  z <- as.matrix(z)
  xy <- dim(z)
  x <- seq(xy[1])
  y <- seq(xy[2])
  z <- z / max(z, na.rm=TRUE) * max(x, y) * scale
  for (i in x) 
  {
    for (j in y) 
    {
      binplot.3d(c(i, i+1), c(j, j+1), z[i,j], alpha=alpha, topcol=col)
    }
  }
}
wireframe(zfit)

layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9)),widths=c(1,1,1), heights =c(1,1,1), respect=TRUE)
par(mai = c(0.5, 0.5, 0.2, 0.2)) #speSEfies the margin size in inches

###Create rSMR -> aSMR crosswalk######################
wd <- tk_choose.dir(); setwd(wd)

####This part includes removing monthly CMD based on surplus (doesn't make much difference)#####
allDat <- fread("BGCv10_100Pt_Rnd_Dat.csv", data.table = FALSE) ###Climate data
temp <- allDat[,grep("CMD",colnames(allDat))]
allDat <- allDat[,c("ID2","PPT_at","PPT_wt","PAS")]
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
allDat <- fread("BGCv10_100Pt_Rnd_Dat.csv", data.table = FALSE)
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

####Create SI models for each species####
###__________________________________________________________________#####

SppList <- c("Pl","Sx","Bl","Cw","Hw","Fd","Lw","Py")

Eda <- read.csv(file.choose())###Edatopic table

colnames(SMRCross) <- c("BGC", "rSMR", "aSMR")
SMRCross$rSMR <- gsub("[[:alpha:]]","", SMRCross$rSMR)
SMRCross$aSMR <- round(SMRCross$aSMR, digits = 0)

sibecOrig <- read.xlsx("SIBEC_for_Portfolio2.xlsx", sheet = 1) ###import actual SI values

modelList <- foreach(Spp = SppList, .combine = c) %do% {
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
}

###Plot polynomial model####
#zfit <- matrix(fitted(fit.bl), ncol = 8)
#surface3d(1:9, 1:5, zfit, col = "purple")
#aspect3d(1,1,1)

####################################################
####Relationship with Temperature#####################
#####Determine equation relating temp and SI####
#############################################

dat <- fread("BGCv10_100Pt_Rnd_Dat.csv", data.table = FALSE)
dat <- dat[,c(2,5,239,230,183,184)]
dat$ID2 <- gsub("[[:space:]]","",dat$ID2)

edaPos <- list(A = c("C",5),B = c("B",3),C = c("D",6)) ###Which edatopic positions?

###Calculate slopes and intercepts
slopes <- foreach(Spp = SppList, .combine = rbind) %do% {##foreach species
  sibec <- sibecOrig[sibecOrig$TreeSpp == Spp,]
  sibec <- sibec[!is.na(sibec$BGCUnit),]
  
  foreach(eda = edaPos, .combine = rbind) %do% { ##foreach edatopic position
    ddBGC <- aggregate(cbind(DD5,MAT,Tave_sm, Tave_sp) ~ ID2, dat, FUN = mean)
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

temp <- read.csv(file.choose())###Need to select BGC Units  (currenlty reading in BulkleyTSA predictions)
missing <- unique(as.character(temp$SS_NoSpace)) 

dat <- fread("BGCv10_100Pt_Rnd_Dat.csv", data.table = FALSE) ###Climate data
dat <- dat[,c(2,5,239,230,183,184)]
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

write.csv(SIFill,"PredSIforPort_Sept18.csv", row.names = FALSE)

#####Now same as above but for each edatopic cell########################################
#########################################################################################
sibec <- read.xlsx("SIBEC_for_Portfolio2.xlsx", sheet = 1)
sibec <- sibec[sibec$TreeSpp == "Lw",]
sibec <- sibec[!is.na(sibec$TreeSpp),]
Spp <- "Lw"

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
