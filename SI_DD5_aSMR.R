library(data.table)
library(magrittr)
library(Rcpp)

cppFunction('NumericVector calcASMR(NumericVector rSMR,NumericVector CMD, DataFrame Rules) {
  NumericVector rSMRLev = Rules["SMRLevel"];
  NumericVector CMDctf = Rules["CMD"];
  NumericVector aSMR = Rules["aSMR"];
  int len = CMD.length();
  IntegerVector v;
  NumericVector out(len);
  for(int i = 0; i < len; i++){
    v = Rcpp::seq(0, aSMR.length()-1);
    if(rSMR[i] < 5){
      v = v[rSMRLev == 0];
    }else if(rSMR[i] == 5){
      v = v[rSMRLev == 5];
    }else if(rSMR[i] == 6){
      v = v[rSMRLev == 6];
    }else{
      v = v[rSMRLev == 7];
    }
    int j;
    for(j = v[0]; j < v[v.length()]; j++){
      if(CMD[i] <= CMDctf[j]){
        break;
      }
    }
    out[i] = aSMR[j];
  }
  return(out);
}')

calcCMD <- function(x){
  CMD <- x[1]
  rSMR <- x[2]
  diff <- rSMR - 4
  if(diff > 0){
    return(CMD/2^diff)
  }else if(diff < 0){
    return(CMD + 100*abs(diff))
  }else{
    return(CMD)
  }
}


allEnv <- fread("./inputs/SIBEC_All_2019_Env.csv")
clim <- allEnv[,.(PlotNumber,Latitude,Longitude,Elevation)]
clim$ID2 <- seq_along(clim$PlotNumber)
clim <- clim[,.(PlotNumber,ID2,Latitude,Longitude,Elevation)] %>% set_colnames(c("ID1","ID2","lat","long","el"))
clim <- clim [complete.cases(clim),]
fwrite(clim,"SI_Plots_ClimDat.csv")

climDat <- fread("./inputs/SI_Plots_ClimDat_Normal_1961_1990MSY.csv")
SI <- fread("./inputs/SIBEC2020_Cleaned.csv")
SI <- SI[,.(PlotNumber,Species, SI_interior,SI_coast)] %>% unique()

climDat <- climDat[,.(ID1,DD5,CMD)]
colnames(climDat)[1] <- "PlotNumber"
envSub <- allEnv[,.(PlotNumber,MoistureRegime,NutrientRegime)]
climEnv <- envSub[climDat, on = "PlotNumber"]
climEnv$MoistureRegime <- as.numeric(climEnv$MoistureRegime)
climEnv <- climEnv[!is.na(MoistureRegime),]
climEnv$EdaCMD <- apply(climEnv[,c("CMD","MoistureRegime")],1,FUN = calcCMD)

rules <- fread("./inputsGit/aSMR_Rules_HalfStep_v11_09Dec2019.csv")
climEnv$aSMR <- calcASMR(rSMR = climEnv$MoistureRegime,CMD = climEnv$CMD,Rules = rules)
climEnv$aSMR <- round(climEnv$aSMR, digits = 0)
climSI <- SI[climEnv, on = "PlotNumber"]
climSI <- climSI[DD5 > 0,]
climSpp <- climSI[Species == "Fd",]
par(mfrow = c(2,3))
#asmr = 4
for(asmr in 2:6){
  temp <- climSpp[aSMR == asmr,]
  fit <- glm(SI_interior ~ DD5, data = temp, family = poisson)
  print(anova(fit))
  summary(fit)
  yweight <- predict(fit,type="response")
  plot(SI_interior ~ DD5, data = temp, main = asmr)
  abline(fit)
}
plot(SI_interior ~ DD5, data = temp[aSMR == 7,])

require(tidyverse)
climSpp2 <- climSpp %>% dplyr::filter(aSMR == 4)
climSpp2 <- climSpp2 %>% mutate(binned = as.factor(cut(DD5, 100))) %>% mutate(binned = as.numeric(binned)) #%>% group_by(bins) %>%  dplyr::summarize(max_SI = max(SI_interior))

climSpp2 [, maxSI := max(SI_interior), by = binned]
clim.plot <- climSpp2 %>% dplyr::select(DD5, maxSI, binned) %>% distinct(binned, .keep_all = TRUE) %>% arrange(DD5)

SI.add <- data.frame("DD5"=c(0, 100, 200, 300, 3000),
                 "maxSI" = c(0, 0, 0, 0, 0),
                 "binned" = c(0, 0, 0, 0, 999))
clim.plot <- rbind(SI.add, clim.plot)%>% arrange(DD5)
                 
model.wt <- SummarizeGrowth(clim.plot$DD5, clim.plot$maxSI)# %>% as.data.frame
plot(model.wt)
p1 <- ggplot(clim.plot, aes(x=DD5,y=maxSI)) + geom_point(alpha=0.5) + theme_bw()+
  geom_text(aes(label = binned))
p1
remove.bins <- c(18, 23, 24, 25, 29, 30, 31, 33, 34, 35, 36, 37, 45, 48, 52, 53, 54, 55, 62, 100)
keep.bins <- c(0, 3, 12, 22, 41, 49, 59)
#remove.bins <- c(1, 3, 5, 7, 12, 13, 18, 22, 30, 32, 35, 37, 45, 48, 50, 51, 52, 53, 54, 55, 57, 62, 63, 69, 100)
clim.plot2 <- clim.plot %>% filter(binned %in% keep.bins) 

model.wt <- SummarizeGrowth(clim.plot2$DD5, clim.plot2$maxSI)# %>% as.data.frame
plot(model.wt)

df.predicted <- data.frame(DD5 = clim.plot2$DD5, pred.wt = predict(model.wt$model))
p2 <- ggplot(clim.plot2, aes(x=DD5,y=maxSI)) + geom_point(alpha=0.5) + theme_bw() +
  geom_text(aes(label = binned))+ xlim(0,2500) + 
  geom_line(data=df.predicted, aes(y=df.predicted$pred.wt), color="red")
  
p2

fit <- glm(maxSI ~ DD5, data = clim.plot, family = poisson)
plot(maxSI ~ DD5 , data = clim.plot) 
curve(predict(fit, newdata = data.frame(DD5 = x)), from = 0, to = 2500, add = TRUE)

# your data
clim.dat <- as.data.frame(climSpp2)
# find & remove outliers
outliers <- boxplot(clim.dat$SI_interior ~ clim.dat$bins )$out
data <- setdiff(clim.dat$SI_interior, outliers)

# fitting a Gaussian
mu <- mean(data)
sigma <- sd(data)

# testing the fit, check the p-value
reference.data <- rnorm(length(data), mu, sigma)
ks.test(reference.data, data)                                      
