stinner <- function(t, params){
  C <- as.numeric(params['C'])
  k1 <- as.numeric(params['k1'])
  k2 <- as.numeric(params['k2'])
  Topt <- as.numeric(params['Topt'])
output <- ifelse(t<Topt, C/(1+exp(k1+k2*t)), C/(1 + exp(k1+k2*(2*Topt-t))))
return(output)
}

logan10 <- function(t, params){
  cc <- as.numeric(params['cc'])
  bb <- as.numeric(params['bb'])
  deltaT <- as.numeric(params['deltaT'])
  Tmax <- as.numeric(params['Tmax'])
  alpha <- as.numeric(params['alpha'])
  output <- ifelse(Tmax > t,alpha * (1/(1 + cc * exp(- bb * t)) - exp(-((Tmax - t)/deltaT))),0)
  return(output)
}

briere2 <- function(params, t){
  p <- as.numeric(params['p']) 
  t0 <- as.numeric(params['t0'])
  Tmax <- as.numeric(params['Tmax'])
  bb <- as.numeric(params['bb'])
  output<-ifelse(Tmax < t | t0 > t, 0 , p*t*(t-t0)*(abs(Tmax-t))^(1/bb))
  return(output)
}

ikemoto <- function(t, params){
  R <- 1.987
  phi <- as.numeric(params['phi'])
  TD <- as.numeric(params['TD'])
  TL <- as.numeric(params['TL'])
  TH <- as.numeric(params['TH'])
  deltaHH <- as.numeric(params['deltaHH'])
  deltaHA <- as.numeric(params['deltaHA'])
  deltaHL <- as.numeric(params['deltaHL'])
  output <- (phi*(t/TD)*exp((deltaHA/R)*((1/TD)-(1/t))))/(1 + exp((deltaHL/R)*((1/TL)-(1/t))) + exp((deltaHH/R)*((1/TH)-(1/t))))
  return(output)
}

rat83 <- function(t, params){
  Tmax <- as.numeric(params['Tmax'])
  Tmin <- as.numeric(params['Tmin'])
  c <- as.numeric(params['c'])
  l <- as.numeric(params['l'])
  output <- ((c * (t - Tmin) * (1 - exp(l * (t - Tmax))))^2) * (t < Tmax) * (t > Tmin) + 0
  return(output)
}

