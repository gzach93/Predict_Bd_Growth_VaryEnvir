#load in data
qld <- read.csv('Data/optical_QLD.csv')

#Varying Optical Density
vary.od <- read.csv('FinalOutput/Vary_OD.csv')
vary.od <- vary.od[,-1]
colnames(vary.od) <- c('TREATMENT', 'day', 'density')
dla <- vary.od[vary.od$TREATMENT == 'Dry_Low_Air', ]
whf <- vary.od[vary.od$TREATMENT == 'Wet_High_Frog', ]
wlf <- vary.od[vary.od$TREATMENT == 'Wet_Low_Frog', ]

#Packages needed
library(coda)
library(rjags)#you have to have JAGS installed on your computer to use this

qld$Treatment <- as.character(qld$Treatment)
qld$Treatment <-substr(qld$Treatment,1,nchar(qld$Treatment)-1)
qld$Treatment <- as.numeric(qld$Treatment)
qld <- qld[order(qld$Treatment),]
trt <- unique(qld$Treatment)

#Briere2
n.chains<-5 #The number of chains the model will set up
n.iter <- 10000 #number of iterations 
model<-jags.model('JAGS/HierarchicalModel/HierarchicalModel_Briere2.bug', #tells rjags where the model/bug file is
                  data = list('Y' = qld$adjusted, 't' = qld$Day, #Y = optical density and t = days
                              'N'= nrow(qld), 'temp' = qld$Treatment), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, K = .1, d = 1,
                             Tmax = 28, t0 = 0, p = .01, bb = 2), n.adapt = 10000) #starting values for the model parameters

OD.Briere2<-coda.samples(model = model, #this function coda.samples is sampling the MCMC
                         variable.names = c('Y0','d','K','sigma', 'p' , 't0', 'Tmax', 'bb'), n.iter = n.iter)
Briere2_dic <- dic.samples(model, n.iter = 10000)
saveRDS(OD.Briere2, 'FinalOutput/HM_Briere2.RDS')
saveRDS(Briere2_dic, 'FinalOutput/dic_briere2.RDS')

##Logan 10
n.chains <- 5 #The number of chains the model will set up
n.iter <- 10000 #number of iterations 
model<-jags.model('JAGS/HierarchicalModel/HierarchicalModel_Logan10.bug', #tells rjags where the model/bug file is
                  data = list('Y' = qld$adjusted, 't' = qld$Day, #Y = optical density and t = days
                              'N'= nrow(qld), 'x' = qld$Treatment), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, K =.1, d = 1,
                             alpha = 1, cc = 70, Tmax = 28, bb = .2, deltaT = 1), n.adapt = 10000) #starting values for the model parameters

OD.Logan<-coda.samples(model = model, #this function coda.samples is sampling the MCMC
                       variable.names = c('Y0','d','K','sigma', 'deltaT' , 'Tmax', 'alpha', 'bb', 'cc'), n.iter = n.iter)
Logan_dic <- dic.samples(model, n.iter = 10000)
saveRDS(Logan_dic, 'FinalOutput/dic_Logan.RDS')
saveRDS(OD.Logan, 'FinalOutput/HM_Logan.RDS')

##Stinner
n.chains<- 5 #The number of chains the model will set up
n.iter <- 1000000 #number of iterations 
model<-jags.model('JAGS/HierarchicalModel/HierarchicalModel_Stinner.bug', #tells rjags where the model/bug file is
                  data = list('Y' = qld$adjusted, 't' = qld$Day, #Y = optical density and t = days
                              'N'= nrow(qld), 'x' = qld$Treatment), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, K =.1, d = 1,
                             k2 = -2, k1 = 32, C = .86, Topt = 21), n.adapt = 1000000) #starting values for the model parameters

OD.stinner<-coda.samples(model = model, #this function coda.samples is sampling the MCMC
                         variable.names = c('Y0','d','K','sigma', 'k1' , 'k2', 'C', 'Topt'), n.iter = n.iter, thin = 100)
Stinner_dic <- dic.samples(model, n.iter = n.iter, thin = 100)
saveRDS(Stinner_dic, 'FinalOutput/dic_Stinner.RDS')
saveRDS(OD.stinner, 'FinalOutput/HM_stinner.RDS')


##Ikemoto
n.chains<-5  #The number of chains the model will set up
n.iter <- 1000000 #number of iterations 
x <- qld$Treatment + 273.15
model<-jags.model('JAGS/HierarchicalModel/HierarchicalModel_Ikemoto.bug', #tells rjags where the model/bug file is
                  data = list('Y' = qld$adjusted, 't' = qld$Day, #Y = optical density and t = days
                              'N'= nrow(qld), 'x' = x), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, d = 1, K = .1,
                             phi = 1,  deltaHA = 16651, deltaHL = -72500, deltaHH = 67500, 
                             TD = 294, TL = 285, TH = 306), n.adapt = 1000000) #starting values for the model parameters

OD.Ikemoto<-coda.samples(model = model, #this function coda.samples is sampling the MCMC
                         variable.names = c('Y0','d', 'K','sigma','phi', 'deltaHH','deltaHA','deltaHL',
                                            'TD','TL','TH'), n.iter = n.iter,thin = 100)
Ikemoto_dic <- dic.samples(model, n.iter = 1000000, thin = 100)
saveRDS(Ikemoto_dic, 'FinalOutput/dic_ikemoto.RDS')
saveRDS(OD.Ikemoto, 'FinalOutput/HM_Ikemoto.RDS')

##Rat83
n.chains <- 5 #The number of chains the model will set up
n.iter <- 10000 #number of iterations 
x <- qld$Treatment + 273.15
model<-jags.model('JAGS/HierarchicalModel/HierarchicalModel_rat83.bug', #tells rjags where the model/bug file is
                  data = list('Y' = qld$adjusted, 't' = qld$Day, #Y = optical density and t = days
                              'N'= nrow(qld), 'x' = x), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, d = 1, K = .1,
                             Tmax = 301,  Tmin = 273, c = 1, l = 1), n.adapt = 10000) #starting values for the model parameters

OD.rat83 <- coda.samples(model = model, #this function coda.samples is sampling the MCMC
                         variable.names = c('Y0','d', 'K','sigma','Tmax', 'Tmin', 'c','l'), n.iter = n.iter)

rat83_dic <- dic.samples(model, n.iter = 10000)
saveRDS(rat83_dic, 'FinalOutput/dic_rat83.RDS')
saveRDS(OD.rat83, 'FinalOutput/HM_rat83.RDS')

########### Run Logistic Model for WHF, WLF, and DLA

#DLA
###Run the Model
n.chains <- 5 #The number of chains the model will set up
n.iter <- 10000 #number of iterations 
model<-jags.model('JAGS/RegularModel/jags-logistic.bug', #tells rjags where the model/bug file is
                  data = list('Y' = dla$density, 't' = dla$day, #Y = optical density and t = days
                              'N'= nrow(dla)), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, d = 1, K = .1), n.adapt = 10000) #starting values for the model parameters

OD.dla <- coda.samples(model = model, #this function coda.samples is sampling the MCMC
                         variable.names = c('Y0','d','r', 'K','sigma'), n.iter = n.iter)

#WHF
###Run the Model
n.chains <- 5 #The number of chains the model will set up
n.iter <- 10000 #number of iterations 
model<-jags.model('JAGS/RegularModel/jags-logistic.bug', #tells rjags where the model/bug file is
                  data = list('Y' = whf$density, 't' = whf$day, #Y = optical density and t = days
                              'N'= nrow(whf)), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, d = 1, K = .1), n.adapt = 10000) #starting values for the model parameters

OD.whf <- coda.samples(model = model, #this function coda.samples is sampling the MCMC
                       variable.names = c('Y0','d','r', 'K','sigma'), n.iter = n.iter)


#WLF
###Run the Model
n.chains <- 5 #The number of chains the model will set up
n.iter <- 10000 #number of iterations 
model<-jags.model('JAGS/RegularModel/jags-logistic.bug', #tells rjags where the model/bug file is
                  data = list('Y' = wlf$density, 't' = wlf$day, #Y = optical density and t = days
                              'N'= nrow(wlf)), n.chains = n.chains, #N = number of data points for for loop
                  inits=list(Y0 = .01, d = 1, K = .1), n.adapt = 10000) #starting values for the model parameters

OD.wlf <- coda.samples(model = model, #this function coda.samples is sampling the MCMC
                       variable.names = c('Y0','d','r', 'K','sigma'), n.iter = n.iter)
