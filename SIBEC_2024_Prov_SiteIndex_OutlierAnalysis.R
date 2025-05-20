################################################################################
#BC4367 - SIBEC Site Index Outliers Script                                     #
#This script will identify site index value outliers for BEC_SS_TreeSpp groups.# 
#The script populates a yes/ no outlier attribute and produces a basic boxplot #
#with outliers labelled (plot is only produced for groups with outliers)       #
#                                                                              #
#                                                                              #
#March 3, 2025 RKITE                                                           #
################################################################################
#clears the workspaces
rm(list=ls(all=T))

#install.packages("RODBC")
library(RODBC)


#set working directory
setwd("I:\\Projects\\BC4367_SIBEC_2020\\GIS_Data\\R_Analysis\\20250319")

#database file path
db<-file.path("I:\\Projects\\BC4367_SIBEC_2020\\GIS_Data\\R_Analysis\\20250319\\15.1_SIBEC_2024_Prov_Analysis_working_RK.accdb")
#connects to access database
ch<- odbcConnectAccess2007("15.1_SIBEC_2024_Prov_Analysis_working_RK.accdb")

#show tables the are accessible
sqlTables(ch)

#read in subset table
site24 <- sqlFetch(ch, "srctbl_SIBEC_SiteIndex24_subset")

#read in raw table
rawSite24<- sqlFetch(ch, "SIBEC_All_2024_Trees") 


###useful commands###
#see first 6 rows. tail() shows last 6
head(site24, n=25)
#see number of rows
nrow(rawSite24)


#set up an empty data frame to hold results
master.df<- NULL


#generates empty folder for plots in the working directory: Plots_todaysdate
dateDIR <- gsub("-","",as.character(Sys.Date()))
outputDIR <- file.path(getwd(),paste("Plots",dateDIR,sep="_"))
if (!dir.exists(outputDIR)) {dir.create(outputDIR)}


#split the subset table by unique BEC_SS_TreeSpp
burst.list <- split(site24, site24$BEC_SS_TreeSpp)

#iterate over the unique BEC_SS_TreeSpp to plot raw site24 values and test for outliers
for(i in 1:length(burst.list)){
  
  #subset unique BEC_SS_TreeSpp for iteration
  df_becspp <- burst.list[[i]]
  print(paste(c("Starting:", unique(df_becspp$BEC_SS_TreeSpp))))
  
  #Add field to record presence of outlier
  df_becspp$Outlier<- NA
  df_becspp$OutlierVal<- NA
  #subset out raw records with matching BEC_SS_TreeSpp
  df_raw<- rawSite24[rawSite24$BEC_SS_TreeSpp==df_becspp$BEC_SS_TreeSpp,]
  
  #check for outliers
  out <- boxplot.stats(df_raw$AnalysisSiteIndex)$out
  
  #Populate outlier attribute with yes or no. 
  df_becspp$Outlier<- ifelse(length(out>0),"Y","N")
  df_becspp$OutlierVal<- ifelse(length(out>0),paste(as.character(out), collapse=" "),"NA")
  
  
  #If outliers exist, produce a simple boxplot with outliers labelled
  if (length(out)>0){
    #set up plot file path
    plotPath<- file.path(outputDIR,paste(gsub("/","_",unique(df_becspp$BEC_SS_TreeSpp)),".jpeg"),sep="")
    
    jpeg(plotPath)
    boxplot(df_raw$AnalysisSiteIndex,ylab = "Site Index",main = "Boxplot of Site Index")
    mtext(paste("Outliers: ", paste(out, collapse = ", ")))
    dev.off()
  }
  master.df<- rbind(master.df,df_becspp)
  print(paste(c("Finishing:", unique(df_becspp$BEC_SS_TreeSpp))))  
  }

#write to a new table in database
sqlSave(ch, master.df, tablename = "srctbl_SIBEC_SiteIndex24_subset_OUTLIERS", fast = TRUE, safer = FALSE, rownames = FALSE, colnames = FALSE)  

#close channel to database
odbcClose(ch)

  