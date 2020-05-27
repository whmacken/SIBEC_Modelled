library(tidyr)
library(ggplot2)
library(foreach)
library(data.table)
library(dplyr)


###Site Index testing###############

SIdat <- fread("./inputs/BartPredSI.csv", data.table = F)
edaDat <- fread("./inputs/Edatopic_v11_7.csv", data.table = F)
edaDat <- edaDat[edaDat$Source %in% c("BECv10","BECv11"),c("MergedBGC", "SS_NoSpace","Edatopic")]
suitDat <- fread("./inputs/TreeSpp_ESuit_v11_18.csv", data.table = F)
suitDat <- suitDat[suitDat$Unit %in% unique(edaDat$SS_NoSpace),c("Unit","Spp","ESuit")]
SIdat <- SIdat[,-4]

dat <- merge(SIdat, suitDat, by = c("Unit","Spp"), all = T)
dat <- dat[!is.na(dat$SIPred),]
dat <- dat[!is.na(dat$ESuit) | dat$SIPred >= 20,]
dat$ESuit[is.na(dat$ESuit)] <- 5
colnames(edaDat) <- c("Subzone","Unit","Edatopic")
edaDat <- unique(edaDat)
dat <- merge(dat, edaDat, by = "Unit", all.x = T)
dat$SIPred <- round(dat$SIPred,1)
dat <- unique(dat)
dat <- dat[!is.na(dat$Edatopic),]


out <- foreach(BGC = unique(dat$Subzone), .combine = rbind) %do% {
  sub <- dat[dat$Subzone == BGC,]
  gridDat <- foreach(eda = unique(sub$Edatopic), .combine = rbind) %do% {
    sub2 <- sub[sub$Edatopic == eda,]
    temp <- aggregate(cbind(SIPred, ESuit) ~ Spp + Edatopic, sub2, FUN = mean) %>%
      mutate(SIPred = round(SIPred,1), ESuit = round(ESuit, 0))
    temp <- temp[order(-temp$SIPred),]
    temp$Spp <- paste(temp$Spp,temp$SIPred, sep = ":")
    if(any(temp$ESuit == 1)){
      dom <- paste(temp$Spp[temp$ESuit == 1], collapse = ",")
      dom <- paste("*",dom,"*\n", sep = "")
    }else{
      dom <- NULL
    }
    if(any(temp$ESuit == 2)){
      sec <- temp$Spp[temp$ESuit == 2]
      sec <- paste(sec, collapse = ",")
      sec <- paste(sec,"\n")
    }else{
      sec <- NULL
    }
    if(any(temp$ESuit == 3)){
      un <- temp$Spp[temp$ESuit == 3]
      un <- paste(un, collapse = ",")
      un <- paste("(",un,")\n", sep = "")
    }else{
      un <- NULL
    }
    if(any(temp$ESuit %in% c(4,5))){
      x <- temp$Spp[temp$ESuit %in% c(4,5)]
      x <- paste(x,",",sep = "")
      if(length(x) >= 2){
        j <- 0
        for(i in seq(2,length(x), by = 2)){x <- append(x,"\n", after = i+j);j <- j+1}
      }
      if(x[length(x)] == "\n"){x <- x[-length(x)]}
      x <- paste(x, collapse = "")
      x <- paste("[",x,"]", sep = "")
    }else{
      x <- NULL
    }
    lab <- paste(dom, sec,  un,x, sep = "")
    
    data.frame(Edatopic = eda, Label = lab)
  }
  gridDat$Alpha <- gsub("[[:digit:]]","",gridDat$Edatopic)
  gridDat$Numeric <- gsub("[[:upper:]]","",gridDat$Edatopic) %>% as.numeric()
  ##plot
  pdf(file = paste("./outputs/SIGrid_",BGC,".pdf",sep = ""), height = 10.5, paper = "letter")
  print(ggplot(data = gridDat)+
          geom_tile(aes(x= Alpha, y = Numeric), color = "black", fill = "white")+
          geom_text(aes(x = Alpha, y = Numeric, label = Label), size = 2.4)+
          scale_y_discrete(limits = c("8","7","6","5","4","3","2","1","0"))+
          scale_x_discrete(limits = c("A","B","C","D","E"))+
          labs(x = "Relative Soil Nutrient Regime", y = "Relative Soil Moisture Regime", cex = 1.5, title = BGC)+
          theme_bw(base_size = 10)+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          coord_fixed())
  dev.off()
  gridDat
}
