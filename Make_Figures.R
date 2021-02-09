qld <- read.csv('Data/optical_QLD.csv')

###### intro thermal performance curve figure 

params <- c('t0' = 6, 'Tmax' = 30, 'p' = .001)
h.fig <- briere1(params = params, t = seq(0,40, by = .1))
test.dat <- data.frame('r' = h.fig, 't' = seq(0,40, by = .1))

png('examplefigure.png', width = 1399, height = 750)
ggplot(test.dat, aes(x = t, y = r)) + geom_line(size = 1.5) +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 36),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_x_continuous(name= expression('Temperature ('*~degree*C*')'), limits=c(0, 35)) +
  scale_y_continuous(name="Performance", limits=c(0, 1.2)) +
  geom_vline(xintercept = 6, size = 1.5, linetype = 2) +
  geom_vline(xintercept = 30, size = 1.5, linetype = 2) + 
  geom_vline(xintercept = 24.7, size = 1.5, linetype = 2) +
  geom_text(aes(x=5, label="Tmin", y=1.13), angle=90,size=10) +
  geom_text(aes(x=23.7, label="Topt", y=1.13), angle=90, size=10) +
  geom_text(aes(x=29, label="Tmax", y=1.13), angle=90, size=10) 
dev.off()

#Look at temperature treatments for qld
qld$Treatment <- as.character(qld$Treatment)
qld$Treatment <-substr(qld$Treatment,1,nchar(qld$Treatment)-1)
qld$Treatment <- as.numeric(qld$Treatment)
qld <- qld[order(qld$Treatment),]

trt <- unique(qld$Treatment)

par(mfrow = c(2,5))
for(i in 1:10){
  plot(qld[qld$Treatment == trt[i], 'Day'], qld[qld$Treatment == trt[i], 'adjusted'], 
       xlim = c(0,15), ylim = c(0,.3), ylab = 'Optical Density', xlab = 'Days', main = trt[i])
}

qld$Treatment <- paste(qld$Treatment,'C')



####### Temperature and predictions ########
library(viridis)
par(mfrow = c(1,1))
t <- seq(0, 16, by = .01)
plot(t, temperature.dla(t), type = 'l',
     main =  'DLA', xlab = 'Time', ylab = 'Temperature')
plot(t, temperature.wlf(t, Kelvin = F), type = 'l',
     main =  'WLF', xlab = 'Time', ylab = 'Temperature')
plot(t, temperature.whf(t, Kelvin = F), type = 'l',
     main =  'WLF', xlab = 'Time', ylab = 'Temperature')

plot(t, temperature.dla(t), type = 'l',
     ylim = c(15, 30), xlab = 'Time', ylab = 'Temperature',
     lwd = 3)
lines(t, temperature.wlf(t, Kelvin = F), lwd = 3, col = 'green')
lines(t, temperature.whf(t, Kelvin = F), lwd = 3, col = 'blue')
legend('topleft', legend = c('DLA','WHF', 'WLF'), 
       text.col = c('black', 'blue', 'green'))

par(mfrow = c(1,3))
plot(dla$day, dla$density, xlab = 'Time (Days)', 
     ylab = 'Optical Density', main = 'DLA',
     ylim = c(0, .2))
draw.quants(b2.quant.dla, seq(0,16,by = .001), 
            col = 'red', lwd = 3)
draw.quants(logan.quant.dla, seq(0,14,by = .001), 
            col = 'green', lwd = 3)
draw.quants(rat.quant.dla, seq(0,14,by = .001),
            col = 'black', lwd = 3)

plot(whf$day, whf$density, xlab = 'Time (Days)', 
     ylab = 'Optical Density', main = 'WHF',
     ylim = c(0, .2))
draw.quants(b2.quant.whf, seq(0,14,by = .001), 
            col = 'red', lwd = 3)
draw.quants(logan.quant.whf, seq(0,14,by = .001), 
            col = 'green', lwd = 3)
draw.quants(rat.quant.whf, seq(0,14,by = .001), 
            col = 'black', lwd = 3)

plot(wlf$day, wlf$density, xlab = 'Time (Days)', 
     ylab = 'Optical Density', main = 'WLF', 
     ylim = c(0,.2))
draw.quants(b2.quant.wlf, seq(0,14,by = .001), 
            col = 'red', lwd = 3)
draw.quants(logan.quant.wlf, seq(0,14,by = .001), 
            col = 'green', lwd = 3)
draw.quants(rat.quant.wlf, seq(0,14,by = .001), 
            col = 'black', lwd = 3)
draw.quants(stinner.quant.wlf, seq(0,14,by = .001), 
            col = 'purple', lwd = 3)

par(mfrow = c(1,1))
temp <- seq(0, 30, by = .1)
plot(temp, rat83(temp+ 273.16, params = HM_rat83[[1]][1,]), 
     type ='l', lwd = 3, ylim = c(0,1.5))
lines(temp, briere2(HM_Briere2[[1]][1,],  temp), col = 'red',
      lwd = 3)
lines(temp, logan10(temp,HM_logan10[[1]][1,]), col = 'green',
      lwd = 3)

#ggplot version
library(ggplot2)
time <- seq(0,14,by = .001)
t <- seq(0,1, by = .01)

temp.dla <- temperature.dla(seq(0,1, by = .01), Kelvin = FALSE)
temp.wlf <- temperature.wlf(seq(0,1, by = .01), Kelvin = FALSE)
temp.whf <- temperature.whf(seq(0,1, by = .01), Kelvin = FALSE)

plot.dla.temp <- ggplot(data.frame(temp.dla), aes(t, temp.dla)) + geom_line() +
  scale_x_continuous(name="Time (Days)", limits = c(0,1), breaks = c(0,1)) +
  scale_y_continuous(name="Temperature", limits = c(15,30)) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

plot.wlf.temp <- ggplot(data.frame(temp.wlf), aes(t, temp.wlf)) + geom_line() +
  scale_x_continuous(name="Time (Days)", limits = c(0,1), breaks = c(0,1)) +
  scale_y_continuous(name="Temperature", limits = c(15,30)) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

plot.whf.temp <- ggplot(data.frame(temp.whf), aes(t, temp.whf)) + geom_line() +
  scale_x_continuous(name="Time (Days)", limits = c(0,1), breaks = c(0,1)) +
  scale_y_continuous(name="Temperature", limits = c(15,30)) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))


wlf.plot <- ggplot(wlf, aes(x = day, y = density)) + geom_point() +
  theme_classic() +
  geom_line(data = rat.quant.wlf, aes(x = time, y = X3, color = '#238A8DFF'),
            size = 1.5) +
  geom_line(data = rat.quant.wlf, aes(x = time, y = X1, color = '#238A8DFF'),
            linetype = 2, size = 1.5) +
  geom_line(data = rat.quant.wlf, aes(x = time, y = X2, color = '#238A8DFF'),
            linetype = 2, size = 1.5) +
  geom_line(data = b2.quant.wlf, aes(x = time, y = X3, colour = '#DCE319FF'),
            size = 1.5) +
  geom_line(data = b2.quant.wlf, aes(x = time, y = X1, colour = '#DCE319FF'),
            linetype = 2, size = 1.5) +
  geom_line(data = b2.quant.wlf, aes(x = time, y = X2, colour = '#DCE319FF'),
            linetype = 2, size = 1.5) +
  geom_line(data = logan.quant.wlf, aes(x = time, y = X3, colour = '#55C667FF'),
            size = 1.5) +
  geom_line(data = logan.quant.wlf, aes(x = time, y = X1, colour = '#55C667FF'),
            linetype = 2, size = 1.5) +
  geom_line(data = logan.quant.wlf, aes(x = time, y = X2, colour = '#55C667FF'),
            linetype = 2, size = 1.5) +
    geom_line(data = stinner.quant.wlf, aes(x = time, y = X3, colour = '#481567FF'),
            size = 1.5) +
  geom_line(data = stinner.quant.wlf, aes(x = time, y = X1, colour = '#481567FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = stinner.quant.wlf, aes(x = time, y = X2, colour = '#481567FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = ikemoto.quant.wlf, aes(x = time, y = X3, colour = '#39568CFF'),
            size = 1.5) +
  geom_line(data = ikemoto.quant.wlf, aes(x = time, y = X1, colour = '#39568CFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = ikemoto.quant.wlf, aes(x = time, y = X2, colour = '#39568CFF'), 
            linetype = 2, size = 1.5) +
   scale_fill_identity(name = '', guide = 'legend', 
                       labels = c('Logan10','Briere2', 'Ratkowsky', 'Stinner')) +
 scale_colour_manual(name = '', 
                     values =c('#55C667FF'='#55C667FF','#DCE319FF'='#DCE319FF', '#238A8DFF'='#238A8DFF', 
                               '#481567FF' = '#481567FF', '#39568CFF' = '#39568CFF'), 
                     labels = c('Ratkowsky', 'Ikemoto', 'Stinner', 'Logan10', 'Briere2')) +
   theme(legend.position = 'Bottom',
         plot.title = element_text(hjust = 0.5, size = 16),
         axis.title = element_text(size = 22),
         axis.text = element_text(size = 20),
         legend.text = element_text(size = 20)) +
   scale_x_continuous(name="", limits=c(0, 14)) +
   scale_y_continuous(name="Optical Density", limits=c(0, .2)) +
  labs(title = 'WLF')



whf.plot <- ggplot(whf, aes(x = day, y = density)) + geom_point() +
  theme_classic() + 
  geom_line(data = rat.quant.whf, aes(x = time, y = X3, color = '#238A8DFF'),
            size = 1.5) +
  geom_line(data = rat.quant.whf, aes(x = time, y = X1, color = '#238A8DFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = rat.quant.whf, aes(x = time, y = X2, color = '#238A8DFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = b2.quant.whf, aes(x = time, y = X3, colour = '#DCE319FF'),
            size = 1.5) +
  geom_line(data = b2.quant.whf, aes(x = time, y = X1, colour = '#DCE319FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = b2.quant.whf, aes(x = time, y = X2, colour = '#DCE319FF'), 
            linetype = 2, size = 1.5) + 
  geom_line(data = logan.quant.whf, aes(x = time, y = X3, colour = '#55C667FF'),
            size = 1.5) +
  geom_line(data = logan.quant.whf, aes(x = time, y = X1, colour = '#55C667FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = logan.quant.whf, aes(x = time, y = X2, colour = '#55C667FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = stinner.quant.whf, aes(x = time, y = X3, colour = '#481567FF'),
            size = 1.5) +
  geom_line(data = stinner.quant.whf, aes(x = time, y = X1, colour = '#481567FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = stinner.quant.whf, aes(x = time, y = X2, colour = '#481567FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = ikemoto.quant.whf, aes(x = time, y = X3, colour = '#39568CFF'),
            size = 1.5) +
  geom_line(data = ikemoto.quant.whf, aes(x = time, y = X1, colour = '#39568CFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = ikemoto.quant.whf, aes(x = time, y = X2, colour = '#39568CFF'), 
            linetype = 2, size = 1.5) +
  scale_fill_identity(name = '', guide = 'legend', 
                      labels = c('Logan10','Briere2', 'Ratkowsky', 'Stinner')) +
  scale_colour_manual(name = '', 
                      values =c('#55C667FF'='#55C667FF','#DCE319FF'='#DCE319FF', '#238A8DFF'='#238A8DFF', 
                                '#481567FF' = '#481567FF', '#39568CFF' = '#39568CFF'), 
                      labels = c('Ratkowsky', 'Ikemoto', 'Stinner', 'Logan10', 'Briere2'))  +
  theme(legend.position = 'none',axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20)) + 
  scale_x_continuous(name="Time (Days)", limits=c(0, 14)) +
  scale_y_continuous(name="Optical Density", limits=c(0, .2)) +
  labs(title = 'WHF')


dla.plot <- ggplot(dla, aes(x = day, y = density)) + geom_point() +
  theme_classic() + 
  geom_line(data = rat.quant.dla, aes(x = time, y = X3, color = '#238A8DFF'),
            size = 1.5) +
  geom_line(data = rat.quant.dla, aes(x = time, y = X1, color = '#238A8DFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = rat.quant.dla, aes(x = time, y = X2, color = '#238A8DFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = b2.quant.dla, aes(x = time, y = X3, colour = '#DCE319FF'),
            size = 1.5) +
  geom_line(data = b2.quant.dla, aes(x = time, y = X1, colour = '#DCE319FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = b2.quant.dla, aes(x = time, y = X2, colour = '#DCE319FF'), 
            linetype = 2, size = 1.5) + 
  geom_line(data = logan.quant.dla, aes(x = time, y = X3, colour = '#55C667FF'),
            size = 1.5) +
  geom_line(data = logan.quant.dla, aes(x = time, y = X1, colour = '#55C667FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = logan.quant.dla, aes(x = time, y = X2, colour = '#55C667FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = stinner.quant.dla, aes(x = time, y = X3, colour = '#481567FF'),
            size = 1.5) +
  geom_line(data = stinner.quant.dla, aes(x = time, y = X1, colour = '#481567FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = stinner.quant.dla, aes(x = time, y = X2, colour = '#481567FF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = ikemoto.quant.dla, aes(x = time, y = X3, colour = '#39568CFF'),
            size = 1.5) +
  geom_line(data = ikemoto.quant.dla, aes(x = time, y = X1, colour = '#39568CFF'), 
            linetype = 2, size = 1.5) +
  geom_line(data = ikemoto.quant.dla, aes(x = time, y = X2, colour = '#39568CFF'), 
            linetype = 2, size = 1.5) +
  scale_fill_identity(name = '', guide = 'legend', 
                      labels = c('Logan10','Briere2', 'Ratkowsky', 'Stinner')) +
  scale_colour_manual(name = '', 
                      values =c('#55C667FF'='#55C667FF','#DCE319FF'='#DCE319FF', '#238A8DFF'='#238A8DFF', 
                                '#481567FF' = '#481567FF', '#39568CFF' = '#39568CFF'), 
                      labels = c('Ratkowsky', 'Ikemoto', 'Stinner', 'Logan10', 'Briere2'))  +
  theme(legend.position = 'none', axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 24)) + 
  scale_x_continuous(name="", limits=c(0, 14)) +
  scale_y_continuous(name="Optical Density", limits=c(0, .2)) +
  labs(title = 'DLA')

#plot_grid(plot.dla.temp, dla.plot, ncol = 1)

library(cowplot)
legend_b <- get_legend(dla.plot + theme(legend.position="bottom"))



#library(gridExtra)
#grid.arrange(wlf.plot, whf.plot, dla.plot, nrow = 1)


library(viridis)
final.wlf <- ggdraw() +
  draw_plot(wlf.plot + theme(legend.justification = "top"), 0, 0, 1, 1) +
  draw_plot(plot.wlf.temp + scale_color_viridis(discrete = TRUE) + 
              theme(legend.justification = "top"), x = 0.15, 
            y = .5, width = 0.4, height = 0.4) 

final.whf <- ggdraw() +
  draw_plot(whf.plot + theme(legend.justification = "top"), 0, 0, 1, 1) +
  draw_plot(plot.whf.temp + scale_color_viridis(discrete = TRUE) + 
              theme(legend.justification = "top"), x = .01, 
            y = .5, width = 0.4, height = 0.4) 

final.dla <- ggdraw() +
  draw_plot(dla.plot + theme(legend.justification = "top"), 0, 0, 1, 1) +
  draw_plot(plot.dla.temp + scale_color_viridis(discrete = TRUE) + 
              theme(legend.justification = "top"), x = .01, 
            y = .5, width = 0.4, height = 0.4) 

prow <- plot_grid(final.wlf, final.whf, final.dla, align = 'vh',
                  labels = 'AUTO',
                  hjust = -1,
                  nrow = 1)
png('FinalPred_11Mar20.png', width = 1399, height = 750)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .2))
dev.off()

######## New Temperature Figure #############
source('new_functions.R')
temp.vary <- read.csv('Data/TemperatureData.csv')

data <- temp.vary[,c(5,6,8)]

par(mfrow= c(1,3))
for(i in 1:6){
  plot((temp.vary$hr-12.36667)/24, data[,i], type = 'l', xlab = 'Days', ylab = 'Temperature',
       main = colnames(data)[i], ylim = c(15,28))
}

library(reshape2)
n.data <- melt(data)
n.data$time <- rep(((temp.vary$hr-temp.vary$hr[1])/24),times = 3)

n.data$variable <- as.character(n.data$variable)
#n.data[which(n.data$variable == 'dlf'),'variable'] <- 'Dry Low Frog'
#n.data[which(n.data$variable == 'dha'),'variable'] <- 'Dry High Air'
n.data[which(n.data$variable == 'dla'),'variable'] <- 'Dry Low Air'
#n.data[which(n.data$variable == 'dlw'),'variable'] <- 'Dry Low Water'
n.data[which(n.data$variable == 'whf'),'variable'] <- 'Wet High Frog'
n.data[which(n.data$variable == 'wlf'),'variable'] <- 'Wet Low Frog'

n.data$variable <- as.factor(n.data$variable)

temp.time <- seq(0,14.1, by = .01)
y.dla <- temperature.dla(temp.time, Kelvin = FALSE)
y.wlf <- temperature.wlf(temp.time, Kelvin = FALSE)
y.whf <- temperature.whf(temp.time, Kelvin = FALSE)

temp.piece <- melt(cbind(y.dla,y.wlf,y.whf))
temp.piece$Var1 <- temp.time
colnames(temp.piece) <- c('time', 'variable', 'value')
temp.piece$variable <- as.character(temp.piece$variable)
temp.piece[which(temp.piece$variable == 'y.wlf'),'variable'] <- 'Wet Low Frog'
temp.piece[which(temp.piece$variable == 'y.whf'),'variable'] <- 'Wet High Frog'
temp.piece[which(temp.piece$variable == 'y.dla'),'variable'] <- 'Dry Low Air'
temp.piece$variable <- as.factor(temp.piece$variable)


ggplot(n.data, aes(time, value))+ geom_line(size = 1.5)+
  geom_line(data = temp.piece, aes(time, value), color = 'grey', linetype = 'longdash')+
  labs(x = 'Time (Days)', y = expression('Temperature ('*~degree*C*')')) + 
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 24)) +
  scale_fill_identity(name = 'Data', guide = 'legend', 
                      labels = c('Incubator Temperature','Piece-Wise Function')) +
  scale_colour_manual(name = 'Data', 
                      values =c('black'='black','grey'='grey'), 
                      labels = c('Incubator Temperature', 'Piece-Wise Function')) +
  facet_grid(. ~ variable) + theme_classic()



####################### 
#####MCMC PLOTS

pdf('MCMC_Plots_14Aug19.pdf')
plot(HM_briere2, sub = 'Briere 2')
plot(HM_ratkowsky, sub = 'Ratkowsky')
plot(HM_logan, sub = 'Logan 10')
plot(HM_ikemoto, sub = 'Ikemoto')
plot(HM_stinner, sub = 'Stinner')
dev.off()

############## Base Logistic Growth ###########
source('functions.R')
yi <- logistic(params = c('Y0'=.001,'K' = 1, 'd' = 3, 'r' =1), t= seq(0,20, by = .1))
plot(seq(0,20,by=.1), yi, type = 'l')
base.log <- data.frame('Time' = seq(0,20,by=.1), 'Y' = yi)
ggplot(base.log, aes(Time, Y)) + geom_line(size = 1.5) +
  xlab('Time (Days)') + ylab('Optical Density')
