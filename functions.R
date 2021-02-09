sampler<-function(results,size){
  n.chains<-length(results)
  comb<-data.frame()
  for(j in 1:n.chains){
    chain<-results[[j]]
    divide<-size/n.chains
    samps<-chain[sample(nrow(chain), size = divide, replace = TRUE),]
    comb<-rbind(samps,comb)
  }
  return(comb)
}



post.solver <- function(t, post.samps, model){
output <- data.frame(matrix(nrow = length(t), ncol = nrow(post.samps)))
for(i in 1:nrow(post.samps)){
  output[,i] <- model(params = post.samps[i,],t = t)
}
return(output)
}


post.lines <- function(t, post.solver,...){
  for(i in 1:ncol(post.solver)){
    lines(t, post.solver[,i], ...)
  }
}


model_vary <- function(curve.model, params, time, temp){
  Y0<- as.numeric(params['Y0'])
  K<- as.numeric(params['K'])
  d<-as.numeric(params['d'])
  r<- curve.model(params = params, t = temp)
  model<- ifelse(time<d, Y0,Y0+K*(Y0)/((Y0)+(K-Y0)*exp(-r*(time-d))))
  return(model) 
}


multi_vary <- function(curve.model, params, time, temp){
  output <- data.frame(matrix(nrow = length(time), ncol = nrow(params)))
  for(i in 1:nrow(params)){
output[,i] <- model_vary(curve.model, params[i,], time, temp)    
  }
  return(output)
}


quant.lines <- function(multivary, time, quantiles){
  output <- data.frame(matrix(nrow = nrow(multivary), ncol = 3))
  for(i in 1:nrow(multivary)){
    output[i,] <- quantile(x = multivary[i,], probs = quantiles)
  }
  output[,3] <- apply(multivary,1,median)
  return(output)
}

draw.quants <- function(output, time, ...){
  lines(time, output[,3], ...)
  lines(time, output[,1], lty = 2, ...)
  lines(time, output[,2], lty =2 , ...)
}


temperature.dla <- function(t, Kelvin = FALSE){
  day <- trunc(t)
  a <- ifelse(Kelvin == TRUE, 273.16, 0)
  t <- t - day
  ifelse(t < 0.03848604, temp <- 18.5 + a,
         ifelse(t > 0.03848604 & t < 0.2530747 ,temp <- (18.858696 + -9.320158 * t)+a,
                ifelse(t > 0.2530747 & t < 0.3321566, temp <- 16.5 + a, 
                       ifelse(t > 0.3321566 & t < 0.6361821, temp <- (11.58362 + 14.80139 * t)+a,
                              ifelse(t > 0.6361821 & t <0.8402777, temp <- 21+a,
                                     ifelse(t > 0.8402777 & t < 0.8844696, temp <- (68.53571 + -56.57143 * t) + a,
                                            ifelse(t > 0.8844696 & t <= 1, temp <- 18.5+a, temp <- NA)))))))
  
}

temperature.whf <-function(t, Kelvin){
  day <- trunc(t)
  a <- ifelse(Kelvin == TRUE, 273.16, 0)
  t <- t - day
  ifelse(t < 0.1273148, temp <- 20.75 + a,
         ifelse(t > 0.1273148 & t < 0.1921296, temp <- (21.732143 +  -7.714286 * t) + a,
                ifelse(t >0.1921296 & t < 0.3402778, temp <- 20.25 + a,
                       ifelse(t > 0.3402778 & t < 0.3472222, temp <- (32.5 +  -36 * t) + a,
                              ifelse(t > 0.3472222 & t < .37499999, temp <- 20 + a,
                                     ifelse(t > .37499999 & t < 0.3819444, temp <- (6.5 +  36 * t) + a,
                                            ifelse(t > 0.3819444 & t < 0.5277778, temp <- 20.25 + a,
                                                   ifelse(t > 0.5277778 & t < 0.5486111, temp <- (8.875  +   21.600 * t) + a,
                                                          ifelse(t > 0.5486111 & t < 0.6805556, temp <- 20.75 + a,
                                                                 ifelse(t > 0.6805556 & t < 0.6875, temp <- (-28.25  +   72.00 * t) + a,
                                                                        ifelse(t > 0.6875 & t <0.8611111, temp <- 21.25 + a, 
                                                                               ifelse(t > 0.8611111 & t < 0.8819444, temp <- (39.825  + -21.600 * t) + a,
                                                                                      ifelse(t >0.8819444 & t <= 1, temp <- 20.75 + a, NA)))))))))))))
}


temperature.wlf <- function(t, Kelvin){
  day <- trunc(t)
  a <- ifelse(Kelvin == TRUE, 273.16, 0)
  t <- t - day
  ifelse(t <  0.1875, temp <- 25 + a,
         ifelse(t > 0.1875 & t < 0.2083333, temp <- (48.625 - 126 * t) + a,
                ifelse(t > 0.2083333 & t < 0.3333333, temp <- 22.5 + a,
                       ifelse(t > 0.3333333 & t < 0.3472222, temp <- (-43.625 + 198 * t) + a,
                              ifelse(t > 0.3472222 & t < 0.5416667, temp <- 24.75 + a,
                                     ifelse(t > 0.5416667 & t < 0.5625, temp <- (-71.075 + 176.4 * t) + a,
                                            ifelse(t > 0.5625 & t < 0.8541667, temp <- 28 + a,
                                                   ifelse(t > 0.8541667 & t < 0.8611111, temp <- (458.5 - 504 * t) + a,
                                                          ifelse(t > 0.8611111 & t <= 1, 25, NA)))))))))
}


logistic.ode <- function(t, state, samps, curve.mod, temp.mod, Kelvin = F){
  with(
    as.list(c(state, samps, curve.mod, temp.mod,Kelvin)),{
      r <- curve.mod(temp.mod(t,Kelvin), params = samps)
      dx <- r*x * (1-x/K) * (t>d) 
      return(list(dx))
    }
  )
}

multi.logistic.ode <- function(t, samps, curve.mod, temp.mod, Kelvin){
  out <- data.frame(matrix(nrow = length(seq(0,14,by = .001)), ncol = nrow(samps)))
  pb = txtProgressBar(min = 0, max = nrow(samps), initial = 0, style = 3) 
  for(i in 1:nrow(samps)){
    state <- c(x = samps[i,'Y0'])
    out[,i] <- ode(y = state, t = t, func = logistic.ode, 
                   parms = samps[i,], curve.mod = curve.mod, temp.mod = temp.mod,
                   Kelvin = Kelvin)[,2]
    setTxtProgressBar(pb,i)
    }
  return(out)
}


