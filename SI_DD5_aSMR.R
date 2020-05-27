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

climDat <- fread("SI_Plots_ClimDat_Normal_1961_1990MSY.csv")
SI <- fread("./inputs/SIBEC2020_Cleaned.csv")
SI <- SI[,.(PlotNumber,Species, SI_interior,SI_coast)] %>% unique()

climDat <- climDat[,.(ID1,DD5,CMD)]
colnames(climDat)[1] <- "PlotNumber"
envSub <- allEnv[,.(PlotNumber,MoistureRegime,NutrientRegime)]
climEnv <- envSub[climDat, on = "PlotNumber"]
climEnv$MoistureRegime <- as.numeric(climEnv$MoistureRegime)
climEnv <- climEnv[!is.na(MoistureRegime),]
climEnv$EdaCMD <- apply(climEnv[,c("CMD","MoistureRegime")],1,FUN = calcCMD)

rules <- fread("./inputs/aSMR_Rules_HalfStep_v11_09Dec2019.csv")
climEnv$aSMR <- calcASMR(rSMR = climEnv$MoistureRegime,CMD = climEnv$CMD,Rules = rules)
climEnv$aSMR <- round(climEnv$aSMR, digits = 0)
climSI <- SI[climEnv, on = "PlotNumber"]
climSI <- climSI[DD5 > 0,]
climSpp <- climSI[Species == "Sx",]
par(mfrow = c(2,3))
for(asmr in 2:7){
  temp <- climSpp[aSMR == asmr,]
  fit <- lm(SI_interior ~ DD5, data = temp)
  print(anova(fit))
  plot(SI_interior ~ DD5, data = temp, main = asmr)
  abline(fit)
}
plot(SI_interior ~ DD5, data = temp[aSMR == 7,])
