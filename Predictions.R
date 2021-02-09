#Librarys Needed
library(rjags)
library(coda)
library(deSolve)
source('new_functions.R')
source('TPCurves.R')

#load in Data and reformat
temp.vary <- read.csv('Data/TemperatureData.csv')
vary.od <- read.csv('FinalOutput/Vary_OD.csv')
vary.od <- vary.od[,-1]
colnames(vary.od) <- c('TREATMENT', 'day', 'density')
time <- ((temp.vary$hr) - (temp.vary$hr[1]))/24
dla <- vary.od[vary.od$TREATMENT == 'Dry_Low_Air', ]
whf <- vary.od[vary.od$TREATMENT == 'Wet_High_Frog', ]
wlf <- vary.od[vary.od$TREATMENT == 'Wet_Low_Frog', ]


#Load in HM MCMC Output
HM_rat83 <- readRDS('FinalOutput/HM_rat83.RDS')
HM_b2 <- readRDS('FinalOutput/HM_Briere2.RDS')
HM_log <- readRDS('FinalOutput/HM_Logan.RDS')
HM_ikemoto <- readRDS('FinalOutput/HM_Ikemoto.RDS')
HM_stin <- readRDS('FinalOutput/HM_stinner.RDS')

#Load in DIC
rat83_dic <- readRDS('FinalOutput/dic_rat83.RDS')
b2_dic <- readRDS('FinalOutput/dic_briere2.RDS')
log_dic <- readRDS('FinalOutput/dic_Logan.RDS')
ikemoto_dic <- readRDS('FinalOutput/dic_ikemoto.RDS')
stinner_dic <- readRDS('FinalOutput/dic_Stinner.RDS')

#Look at DIC Values
rat83_dic
b2_dic
log_dic
ikemoto_dic
stinner_dic

#Sampling MCMC Output
samps <- function(MCMC){
  output <- data.frame(rbind(MCMC[[1]][seq(1,nrow(MCMC[[1]]),50),],
                             MCMC[[2]][seq(1,nrow(MCMC[[2]]),50),],
                             MCMC[[3]][seq(1,nrow(MCMC[[3]]),50),],
                             MCMC[[4]][seq(1,nrow(MCMC[[4]]),50),],
                             MCMC[[5]][seq(1,nrow(MCMC[[5]]),50),]))
  return(output)
 }

b2.samps <- samps(HM_b2)
rat83.samps <- samps(HM_rat83)
logan.samps <- samps(HM_log)
stinner.samps <- samps(HM_stin)
ikemoto.samps <- samps(HM_ikemoto)

#####################################################################
### Run Predictions #################################################

##  DLA
output.b2.dla <- multi.logistic.ode(t = seq(0,14,by = .001), samps = b2.samps, 
                                    curve.mod = briere2, temp.mod = temperature.dla, 
                                    Kelvin = FALSE)

output.logan.dla <- multi.logistic.ode(t = seq(0,14,by = .001), samps = logan.samps, 
                                       curve.mod = logan10, temp.mod = temperature.dla,
                                       Kelvin = FALSE)

output.stinner.dla <- multi.logistic.ode(t = seq(0,14,by = .001), samps = stinner.samps, 
                                         curve.mod = stinner, temp.mod = temperature.dla,
                                         Kelvin = FALSE)

output.rat83.dla <- multi.logistic.ode(t = seq(0,14,by = .001), samps = rat83.samps, 
                                       curve.mod = rat83, temp.mod = temperature.dla,
                                       Kelvin = TRUE)

output.ikemoto.dla <- multi.logistic.ode(t = seq(0,14,by = .001), samps = ikemoto.samps, 
                                         curve.mod = ikemoto, temp.mod = temperature.dla,
                                         Kelvin = TRUE)




b2.quant.dla <- quant.lines(output.b2.dla, time = seq(0,14,by = .001), quantiles = c(.025,.975))
logan.quant.dla <- quant.lines(output.logan.dla, time = seq(0,14,by = .001), quantiles = c(.025,.975))
rat.quant.dla <- quant.lines(output.rat83.dla, time = seq(0,14,by = .001), quantiles = c(.025,.975))
stinner.quant.dla <- quant.lines(output.stinner.dla, time = seq(0,14,by = .001), quantiles = c(.025,.975))
ikemoto.quant.dla <- quant.lines(output.ikemoto.dla, time = seq(0,14,by = .001), quantiles = c(.025,.975))



#FIT LOGISTIC MODELS WHF

output.b2.whf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = b2.samps, 
                                    curve.mod = briere2, temp.mod = temperature.whf, 
                                    Kelvin = FALSE)

output.logan.whf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = logan.samps, 
                                       curve.mod = logan10, temp.mod = temperature.whf,
                                       Kelvin = FALSE)

output.stinner.whf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = stinner.samps, 
                                         curve.mod = stinner, temp.mod = temperature.whf,
                                         Kelvin = FALSE)

output.rat83.whf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = rat83.samps, 
                                       curve.mod = rat83, temp.mod = temperature.whf,
                                       Kelvin = TRUE)

output.ikemoto.whf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = ikemoto.samps, 
                                         curve.mod = ikemoto, temp.mod = temperature.whf,
                                         Kelvin = TRUE)


b2.quant.whf <- quant.lines(output.b2.whf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
logan.quant.whf <- quant.lines(output.logan.whf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
rat.quant.whf <- quant.lines(output.rat83.whf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
stinner.quant.whf <- quant.lines(output.stinner.whf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
ikemoto.quant.whf <- quant.lines(output.ikemoto.whf, time = seq(0,14,by = .001), quantiles = c(.025,.975))


#FIT LOGISTIC MODELS WLF

output.b2.wlf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = b2.samps, 
                                    curve.mod = briere2, temp.mod = temperature.wlf, 
                                    Kelvin = FALSE)

output.logan.wlf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = logan.samps, 
                                       curve.mod = logan10, temp.mod = temperature.wlf,
                                       Kelvin = FALSE)

output.stinner.wlf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = stinner.samps, 
                                         curve.mod = stinner, temp.mod = temperature.wlf,
                                         Kelvin = FALSE)

output.rat83.wlf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = rat83.samps, 
                                       curve.mod = rat83, temp.mod = temperature.wlf,
                                       Kelvin = TRUE)

output.ikemoto.wlf <- multi.logistic.ode(t = seq(0,14,by = .001), samps = ikemoto.samps, 
                                         curve.mod = ikemoto, temp.mod = temperature.wlf,
                                         Kelvin = TRUE)


b2.quant.wlf <- quant.lines(output.b2.wlf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
logan.quant.wlf <- quant.lines(output.logan.wlf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
rat.quant.wlf <- quant.lines(output.rat83.wlf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
stinner.quant.wlf <- quant.lines(output.stinner.wlf, time = seq(0,14,by = .001), quantiles = c(.025,.975))
ikemoto.quant.wlf <- quant.lines(output.ikemoto.wlf, time = seq(0,14,by = .001), quantiles = c(.025,.975))

##### Rough Figure to look at output
plot(dla$day, dla$density, ylab = 'Optical Density', xlab = 'Time (Days)', main = 'DLA', ylim = c(0,.2), xlim = c(0,14))
draw.quants(b2.quant.dla, time = seq(0,14, by = .001), col = 'red')
draw.quants(logan.quant.dla, time = seq(0,14, by = .001), col = 'blue')
draw.quants(rat.quant.dla, time = seq(0,14, by = .001), col = 'green')

plot(wlf$day, wlf$density, ylab = 'Optical Density', xlab = 'Time (Days)', main = 'WLF', ylim = c(0,.2), xlim = c(0,14))
draw.quants(b2.quant.wlf, time = seq(0,14, by = .001), col = 'red')
draw.quants(logan.quant.wlf, time = seq(0,14, by = .001), col = 'blue')
draw.quants(rat.quant.wlf, time = seq(0,14, by = .001), col = 'green')

plot(whf$day, whf$density, ylab = 'Optical Density', xlab = 'Time (Days)', main = 'WHF', ylim = c(0,.2), xlim = c(0,14))
draw.quants(b2.quant.whf, time = seq(0,14, by = .001), col = 'red')
draw.quants(logan.quant.whf, time = seq(0,14, by = .001), col = 'blue')
draw.quants(rat.quant.whf, time = seq(0,14, by = .001), col = 'green')


## find median and hdi for parameters and adjusted Tmax and Tmin values
t <- seq(0,60, by = .01)

median(apply(rat83.curve, 2, max))
hdi(apply(rat83.curve, 2, max))

median(apply(b2.curve, 2, max))
hdi(apply(b2.curve, 2, max))

median(apply(stin.curve, 2, max))
hdi(apply(stin.curve, 2, max))

median(apply(ikemoto.curve, 2, max))
hdi(apply(ikemoto.curve, 2, max))

median(apply(log.curve, 2, max))
hdi(apply(log.curve, 2, max))


rat83.topt <- NA
for(i in 1:ncol(rat83.curve)){
  rat83.topt[i] <- which(rat83.curve[,i] == max(rat83.curve[,i]))
}
median(t[rat83.topt])
hdi(t[rat83.topt])

log.topt <- NA #Compare this to the equation solution
for(i in 1:ncol(log.curve)){
  log.topt[i] <- which(log.curve[,i] == max(log.curve[,i]))
}
median(t[log.topt])
hdi(t[log.topt])


ikemoto.topt <- NA
for(i in 1:ncol(ikemoto.curve)){
  ikemoto.topt[i] <- which(ikemoto.curve[,i] == max(ikemoto.curve[,i]))
}
median(t[ikemoto.topt])
hdi(t[ikemoto.topt])


stin.topt <- NA
for(i in 1:ncol(stin.curve)){
  stin.topt[i] <- which(stin.curve[,i] == max(stin.curve[,i]))
}
median(t[stin.topt])
hdi(t[stin.topt])

b2.topt <- NA #Check estimate with equation
for(i in 1:ncol(b2.curve)){
  b2.topt[i] <- which(b2.curve[,i] == max(b2.curve[,i]))
}
median(t[b2.topt])
hdi(t[b2.topt])


##### Adjusted Tmin at 0.01
adj.tim.b2 <- NA
for(i in 1:ncol(b2.curve)){
  adj.tim.b2[i] <- t[min(which(b2.curve[,i] > 0.01)-1)]
}
median(adj.tim.b2)
hdi(adj.tim.b2)

adj.tim.rat83 <- NA
for(i in 1:ncol(rat83.curve)){
  adj.tim.rat83[i] <- t[min(which(rat83.curve[,i] > 0.01)-1)]
}
median(adj.tim.rat83)
hdi(adj.tim.rat83)

adj.tim.ik <- NA
for(i in 1:ncol(ikemoto.curve)){
  adj.tim.ik[i] <- t[min(which(ikemoto.curve[,i] > 0.01)-1)]
}
median(adj.tim.ik)
hdi(adj.tim.ik)

adj.tim.stin <- NA
for(i in 1:ncol(stin.curve)){
  adj.tim.stin[i] <- t[min(which(stin.curve[,i] > 0.01)-1)]
}
median(adj.tim.stin)
hdi(adj.tim.stin)

adj.tim.log <- NA
for(i in 1:ncol(log.curve)){
  adj.tim.log[i] <- t[min(which(log.curve[,i] > 0.01)-1)]
}
median(adj.tim.log)
hdi(adj.tim.log)

#### Adj Tmax 0.01 ####

adj.tmax.b2 <- NA
for(i in 1:ncol(b2.curve)){
  adj.tmax.b2[i] <- t[max(which(b2.curve[,i] > 0.01)+1)]
}
median(adj.tmax.b2)
hdi(adj.tmax.b2)

adj.tmax.rat83 <- NA
for(i in 1:ncol(rat83.curve)){
  adj.tmax.rat83[i] <- t[max(which(rat83.curve[,i] > 0.01)+1)]
}
median(adj.tmax.rat83)
hdi(adj.tmax.rat83)

adj.tmax.stin <- NA
for(i in 1:ncol(stin.curve)){
  adj.tmax.stin[i] <- t[max(which(stin.curve[,i] > 0.01)+1)]
}
median(adj.tmax.stin)
hdi(adj.tmax.stin)

adj.tmax.ik <- NA
for(i in 1:ncol(ikemoto.curve)){
  adj.tmax.ik[i] <- t[max(which(ikemoto.curve[,i] > 0.01)+1)]
}
median(adj.tmax.ik)
hdi(adj.tmax.ik)

adj.tmax.log <- NA
for(i in 1:ncol(log.curve)){
  adj.tmax.log[i] <- t[max(which(log.curve[,i] > 0.01)+1)]
}
median(adj.tmax.log)
hdi(adj.tmax.log)
