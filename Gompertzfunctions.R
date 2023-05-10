
##from http://www.statsathome.com/2017/06/07/fitting-non-linear-groth-curves-in-r/
library(tidyverse)
require(purrr)
set.seed(4)
fit.gompertz <- function(data, time){
  d <- data.frame(y=data, t=time)
  
  # Must have at least 3 datapoints at different times
  if (length(unique(d$t)) < 3) stop("too few data points to fit curve")
  
  # Pick starting values ###
  i <- which.max(diff(d$y))
  starting.values <- c(a=max(d$y), 
                       mu=max(diff(d$y))/(d[i+1,"t"]-d[i, "t"]), 
                       lambda=i)
  print("Starting Values for Optimization: ")
  print(starting.values)
  ##########################
  
  formula.gompertz <- "y~a*exp(-exp(mu*exp(1)/a*(lambda-t)+1))"
  nls(formula.gompertz, d, starting.values)
}

gompertz <- function(time, a, mu, lambda){
  y <- a*exp(-exp(mu*exp(1)/a*(lambda-time)+1))
  return(data.frame(time=time, y=y))
}

d <- gompertz(1:100, 10, 2, 30)
plot(d)

for(i in 1:nrow(d)) d[i,2] <- rnorm(1, d[i,2], 1)

(fit <- fit.gompertz(d$y, d$time))
plot(d, ylab="microbial abundance")
lines(d$time, predict(fit))
safe.fit.gompertz <- safely(fit.gompertz)

safe.fit.gompertz(c(1,2), c(19, 19))

AIC(fit)
