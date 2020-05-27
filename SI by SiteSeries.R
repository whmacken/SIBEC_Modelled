require(plyr)
require(dplyr)
require(tidyr)
require(TeachingDemos)
require (ggplot2)
require (assertr)
#install.packages("TeachingDemos")
#setwd("C:/Users/Alana2012/Sync/BVRC/SORTIE position/Data/SIBEC")

#SI_raw <- fread("./SI_Prediction/SIBEC_SBSmc2_2019_Trees.csv") ##v1 recalcuated SI using Count Age and TotalHtCalc
SI_raw <- fread("./SI_Prediction/SIBEC_SBSmc2_2019_Trees_v2.csv") ##v2 Use the AgeBH where available - gets rid of some outliers. Sx starts looking better 
#SI_explore <- fread("SBSmc2_SI_Cleaned.csv")
table(SI_raw[,SiteSeries])
sort(unique(SI_raw$TreeSpp))
# SI_explore <- SI_explore %>% 
#   mutate(SiteSeries_rename = ifelse(SiteSeries %in% c('01a', '01b', '01c'), 1, 
#                      ifelse(SiteSeries %in% c('09a', '09b'), 9,
#                             ifelse (SiteSeries %in% c('10a','10b'), 10, SiteSeries))))

##### Clean up Species Codes
SI_raw$TreeSpp <- SI_raw$TreeSpp %>% recode(
  "Sw" = "Sx", "sw" = "Sx", "Se" = "Sx",
  "EP"= "Ep", "FD" = "Fd", "Fd " = "Fd", "PL" = "Pl")
SI_raw <- SI_raw[!SI_raw$TreeSpp == "",]

sort(unique(SI_raw$TreeSpp)) # check for any other miscodes

SI_info <- SI_raw %>% select(Zone, SubZone, SiteSeries, TreeNum, )
#Tree_info <- SI_raw %>% select(PlotNumber,TreeSpp, TreeNum)%>% dplyr::rename(Species = TreeSpp)
### Reduce the variable set for Site Tools
SI_clean <-  SI_raw %>% select(PlotNumber,TreeSpp, NewAge, TotalHtCalc) %>% dplyr::rename(Species = TreeSpp, `BHAge` = NewAge, `Height` = TotalHtCalc )
write.csv(SI_clean, "./SI_Prediction/SIBEC_SBScm2_2019_ForSiteTools_v2.csv", row.names = FALSE) ####submit this file to Site Tools (no file header)
###merge the results back into the dataframe
SI_clean_interior <- fread("./SI_Prediction/SIBEC_SBSmc2_2019_SiteTools4.1_interior_v2.csv")
SI_clean_interior <- cbind(SI_info, SI_clean_interior)
SI_tidy <- SI_clean_interior
#------------Include if needs both interior and coastal SI for complete set of tables
# SI_clean_coast <- fread("./SI_Prediction/SIBEC_All_2019_SiteTools4.1_coast.csv")
# SI_clean_coast <-  SI_clean_coast %>% select('Site Index m') %>% dplyr::rename(SI_coast = 'Site Index m')
# SI_tidy <- cbind(SI_clean_interior, SI_clean_coast)

SI_tidy <- SI_tidy[!duplicated(SI_tidy),]
SI_tidy <- SI_tidy %>% dplyr::rename(BHAge = 'BH Age yrs', Height = 'Height m', SI_interior = 'Site Index m')
SI_tidy <- SI_tidy[complete.cases(SI_tidy)]
SI_tidy <- SI_tidy[! (SI_tidy$SiteSeries == ""),]
write.csv(SI_tidy, "./SI_Prediction/SIBEC_SBSmc2_2019_v2.csv", row.names = FALSE)


table(SI_tidy[, "SiteSeries"])
SI_spp <- SI_tidy %>% filter(Zone == "SBS", SubZone == "mc2")#, Species == "Lw")
SI_spp <- drop_na (SI_spp)
SI_sppsave <- SI_spp
remove <-c("SBC3727", "SBC3724", "SBC3711", "SBC3712", "SBC3686", "SBC7318", "SBC3109", "SBC3108","SBC3422",
             "SBC6092", "SBC3716", "SBC3714", "SBC6145" ) # these are 12 plots with much higher productivity than the second set (2 phases?)
SI_spp <- SI_spp [!(SI_spp$PlotNumber %in% remove), ]
##box blot of site index by site series for selected zone
SI_spp$SiteSeries <- factor(SI_spp$SiteSeries, levels = c("02", "03", "01c", "01b", "01a", "01", "04", "05", "06", "09", "09b", "09a", "10a", "10", "07", "10b", "12a", "12"))
SI_box <- ggplot(SI_spp, aes(x=SiteSeries, y=SI_interior, fill=Species)) + 
  geom_boxplot()+
  facet_wrap(~Species)
plot(SI_box)

ggsave("SI_boxplot_SBSmc2ver2.jpg", SI_box, device="jpeg", path = "./SI_prediction/outputs/")
### shows the plot number of outliers
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r")
boxplot.with.outlier.label(SI_spp$SI_interior~SI_spp$SiteSeries, SI_spp$PlotNumber)

#####################assertr codee

#SI_clean$SiteSeries <-  as.factor((SI_clean$SiteSeries))
#assertr check of data
SI_clean <- complete.cases(SI_clean)
(SI_clean)SI_clean <- 
not.empty <- function(x) if(x=="") return(FALSE)
SI_clean %>% assert(not.empty, SiteSeries)

check_me <- . %>%
  chain_start %>%
  verify(nrow(SI_clean) > 10) %>%
  verify(mpg > 0) %>%
  insist(within_n_sds(4), mpg) %>%
  assert(in_set(0,1), am, vs)
chain_end 

mtcars %>%
  check_me %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))