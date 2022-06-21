################################################
#### ALTERNATIVE KF AND ENKF IMPLEMENTATION ####
################################################




library(dlm)
library(forecast)
library(ggplot2)
library(zoo)
library(stats)
library(readr)
library(lmtest)
library(TSA)
library(tseries)


######################
### Experimental data 
######################

X2022_03_10_17_20_36 <- read_csv("2022-03-14-15-41-21.csv")
CE <- X2022_03_10_17_20_36[ 60: 240, ]
CE <- cbind( CE$`657 [m/s]`, CE$`076 [m/s]`)
colnames(CE) <- c("velE", "velC")
CE <- as.data.frame(CE)
x <- index(CE)
g.dlm <- ggplot(data =CE, aes(x= x, y= as.numeric(velE))) + geom_line(color= "navy")+
  geom_line(aes( y=as.numeric(velC)), color = "coral")
g.dlm
velC <- as.numeric(CE$velC)
velE <- as.numeric(CE$velE)

auto.arima(velE)  # WN !! :)
auto.arima(velC)  # WN !! :)
sd(velE)
mean(velE)
sd(velC)

### VelE WN:
acf(velE, main=NULL)
pacf(velE)
eacf(velE)
adf.test(velE)
# p valor petit llavors no té arrel unitària i podria ser estacionari
modWN <- auto.arima(velE)
checkresiduals(modWN)
Box.test(residuals(modWN), type="Ljung-Box") #H0 soroll blanc
ks.test(residuals(modWN), "pnorm")



X2022_03_10_17_20_36 <- read_csv("2022-03-14-15-27-49.csv")
AE <- X2022_03_10_17_20_36[ 120: 300, ]
AE <- cbind( AE$`657 [m/s]`, AE$`076 [m/s]`)
colnames(AE) <- c("velE", "velA")
AE <- as.data.frame(AE)
x <- index(AE)
g.dlm <- ggplot(data =AE, aes(x= x, y= as.numeric(velE)), color= "navy") + geom_line()+
  geom_line(aes( y=as.numeric(velA)), color = "coral")
g.dlm
velA <- as.numeric(AE$velA)
auto.arima(velA)

exp <- as.data.frame(cbind(c(1, 57, 133), c(velE[120], velA[120], velC[120])))



###### EXPERIMENTAL PLOTS

experimentalFilter  <- data.frame(x = index(velA[1:179]), 
                                  values = c(velA[1:179], 
                                             velC[1:179], 
                                             velE[1:179]),
                                  #mmA$mitjana[59:178], PTu$velC06[59:178]), 
                                  Data = c(rep("Point B", 179), 
                                           rep( "Point D",179), 
                                           rep( "Point E", 179)
                                           ))
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c("navy", "darkred", "darkgreen")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c( "solid","solid","solid"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")

auto.arima(velE)
sd(velE)
white <- rnorm(179, mean= 1.8820, sd=sd(velE))
length(white)

experimentalFilter  <- data.frame(x = index(velA[1:179]), 
                                  values = c(white, 
                                             velE[1:179]),
                                  #mmA$mitjana[59:178], PTu$velC06[59:178]), 
                                  Data = c(rep( "WN simulated",179), 
                                           rep( "Experimental velocity E", 179)
                                  ))
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c("green", "darkgreen")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c( "dashed","solid"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")

errorWNE <- (white - velE[1:179])^2
mean(errorWNE)

#mod <- arima(order=c(0,0,0), mean= 1.8820)
#pred <- predict(auto.arima(velE))
#length(pred)

######################
### Simulated data ###
######################

DataL <- read.csv("Long Laminar  - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
Laminar <- as.data.frame(cbind(DataL$X.2[3:153]*10^(-16), DataL$X.5[3:153]*10^(-16)) )
colnames(Laminar) <- c("vel04", "vel035")
DataL <- read.csv("Dades longitudinals - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
VelL <- as.data.frame(DataL$X.3*10^(-16))
colnames(VelL) <- "vel"

DataL <- read.csv("Dades longitudinals - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
VelL <- as.data.frame(DataL$X.3*10^(-16))
colnames(VelL) <- "vel"
DataL <- read.csv("Long Turbulent 045 - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
Turbulent045 <- as.data.frame(DataL$X.3*10^(-16))
colnames(Turbulent045) <- "vel"
Turbulents <- as.data.frame(cbind(VelL, Turbulent045$vel[2:152]))
colnames(Turbulents) <- c("vel04", "vel045")

LongDiam <- read.csv("LongDiam - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
LD <- as.data.frame(cbind(LongDiam$X.2[3:153], LongDiam$X.4[3:153], LongDiam$X.6[3:153], LongDiam$X.8[3:153], LongDiam$X.10[3:153], LongDiam$X.12[3:153], LongDiam$X.14[3:153]))
colnames(LD) <- c("vel039904", "vel039906", "vel0399", "vel0399035", "vel0406", "vel04065", "vel04")

mm <- as.data.frame(cbind((LD$vel04065+LD$vel04)/2, (LD$vel04+ Turbulents$vel045)/2, ((LD$vel04065+LD$vel04)/2 + Turbulents$vel045)/2 ))
colnames(mm) <- c("mitjana1", "mitjana045", "mitjanaX2")

mmL <- as.data.frame((Laminar$vel04 + Laminar$vel035)/2)
colnames(mmL) <- "mitjana"

LLdata <- read.csv("Puntuals 180s - Full 2.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
LL <- as.data.frame(cbind(LLdata$X.1[2:152], LLdata$X.3[2:152], LLdata$X.5[2:152]))
colnames(LL) <- c("vel04", "vel036", "vel037")


###########################
### Punctual simulated data
###########################

PTdata <- read.csv("Puntuals 180s - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
PT <- as.data.frame(cbind(PTdata$X.1[4:362]*10^(-16), PTdata$X.3[4:362]*10^(-16), PTdata$X.5[4:362]*10^(-16), PTdata$X.7[4:362]*10^(-16)))
colnames(PT) <- c("velA05", "velC05", "velC06", "velA06")
row_odd <- seq_len(nrow(PT)) %% 2              # Create row indicator
row_odd  
PTu <- as.data.frame(PT[row_odd == 0, ] ) 

mmA <- as.data.frame((PTu$velA05+ PTu$velA06)/2)
colnames(mmA) <- "mitjana"

## Laminar

PLdata <- read.csv("Puntuals 180s - Full 3.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
PL <- as.data.frame(cbind(PLdata$X[2:359], PLdata$X.2[2:359]))
colnames(PL) <- c("velA037", "velC037")
row_odd <- seq_len(nrow(PL)) %% 2              # Create row indicator
row_odd  
PLa <- as.data.frame(PL[row_odd == 0, ] ) 


###
### Només aplicar DA amb les dades de A i predim C
###

########
## KF ##
########

# Turbulent directament amb MA(2), ja que hme vist que és el model més òptim 
## Però amb el KF no podem tenir la distribució de tota la canonada, no? 
## Preguntar Alejandra, 
## En teoria el KF ens forecastea la distribuciño en els punts dels quals tenim la distribució, no?


auto.arima(mmA$mitjana[59:178]) ## MA(2) with mean 
## log likelihood=457.83
mod1 <- arima(mmA$mitjana[59:178], order = c(2,0,0), include.mean = TRUE)
mod1
##log likelihood = 400.62

## MA(2)
m1.dlm <- dlmModARMA( ma= c(0.8504, 0.2416))
m1.dlm
model.filteredA3 <- dlmFilter(velA-mean(velA), m1.dlm)
A3mean <-  model.filteredA3$f + mean(velA)

model.filteredCA <- dlmFilter(velC-mean(velC), m1.dlm)
CAmean <-  model.filteredCA$f + mean(velC)

experimentalFilter  <- data.frame(x = index(velA[60:179]), 
                                  values = c(velA[58:177], 
                                             #model.filteredA3.wmean[59:178], 
                                             velC[58:177], 
                                             #model.filteredC3.wmean[59:178],
                                             CAmean[59:178],
                                             A3mean[59:178]),
                                             #mmA$mitjana[59:178], PTu$velC06[59:178]), 
                                  Data = c(rep("True velocity B", 120), 
                                           #rep( "Filtered state A MA(2) w intercept", 120),
                                           rep( "True velocity D",120), 
                                           #rep( "Filtered State C AR(2) w intercept", 120),
                                           rep( "FS D MA(2) w/mean", 120),
                                           rep( "FS B MA(2) w/mean", 120)))
                                           #rep("Simulated velocity B", 120), 
                                           #rep("Simulated velocity D", 120)) )
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "blue",  "red", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c( "solid","solid","dotdash","dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")

#################
## KF D Turbulent per predir A

auto.arima(PTu$velC06[59:178]) ## MA(1) with mean
## log likelihood=457.83
mod3 <- arima(PTu$velC05[59:178],  order = c(2,0,0), include.mean = TRUE)
mod3
m4.dlm <- dlmModARMA( ar= c(0.8807,-0.3596))
m4.dlm
model.filteredC3 <- dlmFilter(velC - mean(velC), m4.dlm)
C3mean <-  model.filteredC3$f + mean(velC)

model.filteredAC <- dlmFilter(velA-mean(velA), m4.dlm)
ACmean <-  model.filteredAC$f + mean(velA)


experimentalFilter  <- data.frame(x = index(velA[60:179]), 
                                  values = c(velA[58:177], 
                                             #model.filteredA3.wmean[59:178], 
                                             velC[58:177], 
                                             #model.filteredC3.wmean[59:178],
                                             CAmean[59:178],
                                             A3mean[59:178],
                                             ACmean[59:178],
                                             C3mean[59:178]),
                                  #mmA$mitjana[59:178], PTu$velC06[59:178]), 
                                  Data = c(rep("True velocity B", 120), 
                                           #rep( "Filtered state A MA(2) w intercept", 120),
                                           rep( "True velocity D",120), 
                                           #rep( "Filtered State C AR(2) w intercept", 120),
                                           rep( "FS D MA(2) w/mean", 120),
                                           rep( "FS B MA(2) w/mean", 120),
                                           rep( "FS B AR(2) w/mean", 120),
                                           rep( "FS D AR(2) w/mean", 120)
                                           ))
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "blue", "deepskyblue",  "red",  "coral",  "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c("dotdash", "solid","solid", "dotdash", "dashed","dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")


### Però realment no és massa significatiu, ja que seguim fent servir les dades experimentals en el filtre!!



##########
## EnKF ##
##########

library(mvtnorm)
source("http://www.datall-analyse.nl/R/EnKF.R")
EnKF

## EnKF per a WN laminar en both cases

## Only A measurement

laminarV <- LL$vel037[50:139]
ns <- 120 #number of time itereations
nc <- 140-50 #number of nodes
## Specify noises (variances)
wk <- 5e-5 #process noise
#vk <- 5e-1 #measurement noise
vk <- sd(velA)
##Construct the state transition matrix
A <- matrix(0, ncol=nc, nrow=nc)
#Entries of the interior nodes in the transition matrix
# identity for a white noise
for(i in 2:(nc-1)) {
  A[i, i] <- 1
}

np <- 2 #number of measurement points
### Simulated points
x <- matrix(mean(laminarV), ncol=nc, nrow=ns)
x[1,] <- laminarV
set.seed(12)
for (i in 2:ns) {
  x[i, ] <-  (A%*%x[i-1, ])+ rnorm(nc, 0, sqrt(wk))}
x
dataXT <- as.data.frame(t(x))
dataXT <- dataXT[2:nc-1,]
g.dlm <- ggplot(data =dataXT, aes(x= index(dataXT[,1]), y=dataXT[,1] )) + geom_line(color= "navy")+
  geom_line(aes( y=dataXT[,30]), color = "coral")+
  geom_line(aes( y=dataXT[,60]), color = "turquoise")+
  geom_line(aes( y=dataXT[,90]), color = "green")+
  geom_line(aes( y=dataXT[,120]), color = "pink")+
  ggtitle("Simulated velocity for different times")+
  xlab("Node (cm)") + ylab("Velocity (m/s)")
g.dlm

ex1 <- list(m0=laminarV , #initial state estimates
            #error covariances of the initial state estimates:
            C0=diag(rep(0.1, nc)),
            #measurement noise
            V=diag(rep(vk, np)),
            #process noise
            W=diag(rep(wk, nc)))

#State transition function
GGfunction <- function (x, k){
  (A%*%x)}
#Observation/measurement function
FFfunction <- function (x, k){
  # cbind(x[57-50], x[133-50])}
  #cbind(x[7], x[83])
  cbind(x[7], x[83])
}

dataEx2 <- as.data.frame(cbind(velA, velC))

enkf1 <- EnKF(y=dataEx2, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx2, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx2, mod=ex1, size=50,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx2, mod=ex1, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)

plot(1:nc, laminarV, type="l", col= "black", #col=c(gray(level=.5)),
     ylim=range(c(enkf1$m[121,], enkf2$m[121,], enkf3$m[121,], enkf4$m[121,])),
     xlab="Node, i", ylab="Velocity (i)", main="t=2")
lines(1:nc, enkf1$m[121,], lty=2, col="blue", lwd=1)
lines(1:nc, enkf2$m[121,], lty=2, col="red", lwd=1)
lines(1:nc, enkf3$m[121,], lty=2, col="darkgreen", lwd=1)
lines(1:nc, enkf4$m[121,], lty=2, col="turquoise", lwd=1)
points(c(57-50, 133-50), dataEx2[120,], col="coral", cex=1, pch = 20)
legend("bottomleft", lty=c(1, 2, 2, 2),
       col=c("black", "blue", "red", "darkgreen", "turquoise"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 25 members", "EnKf with 50 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

y <- as.data.frame(laminarV)

#Error plot
par(mfrow=c(1, 1))
e1 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf3$m[i+1,])^2))
e4 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf4$m[i+1,])^2))

rangeError <- range(cbind(e1, e2, e3, e4))
plot(1:(nrow(y)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(y)-10), e2, lty=1, col="red")
lines(1:(nrow(y)-10), e3, lty=1, col="darkgreen")
lines(1:(nrow(y)-10), e4, lty=1, col="turquoise")
legend("topleft", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen", "turquoise"),
       legend=c("EnKf with 10 members",
                "EnKf with 25 members","EnKf with 50 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

enkfWN <- enkf3

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf2$m[62:182,57-50], 
                                             enkf2$m[62:182,133-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 25 members B",121), 
                                           rep( "EnKf with 25 members D", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","deepskyblue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

## És horrible
## EnKF no funciona bé amb un sol punt de mesura (en algun moment peta l'algorsime tal i com està fet)
## Si intento afegir les dades que tenim de la velE i velA dóna uns resultats horribles per al punt D 

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkfWN$m[62:182,57-50],
                                             enkfWN$m[62:182,133-50]
                                             ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("deepskyblue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")


## Provo fent servir dos cops la mateixa distribució experimental com a dos punts de mesura diferents consecutius

FFfunction <- function (x, k){
  # cbind(x[57-50], x[133-50])}
  #cbind(x[7], x[83])
  cbind(x[7], x[8])
}

dataEx4 <- as.data.frame(cbind(velA, velA))

enkf1 <- EnKF(y=dataEx4, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx4, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx4, mod=ex1, size=50,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx4, mod=ex1, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)

y <- as.data.frame(laminarV)

#Error plot
par(mfrow=c(1, 1))
e1 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf3$m[i+1,])^2))
e4 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf4$m[i+1,])^2))

rangeError <- range(cbind(e1, e2, e3, e4))
plot(1:(nrow(y)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(y)-10), e2, lty=1, col="red")
lines(1:(nrow(y)-10), e3, lty=1, col="darkgreen")
lines(1:(nrow(y)-10), e4, lty=1, col="turquoise")
legend("topleft", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen", "turquoise"),
       legend=c("EnKf with 10 members",
                "EnKf with 25 members","EnKf with 50 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

enkfWN2 <- enkf3


experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf2$m[62:182,57-50], 
                                             enkf2$m[62:182,133-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 25 members B",121), 
                                           rep( "EnKf with 25 members D", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","deepskyblue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

## Molt aleatori, resulta que el de size 50 ho fa força bé 

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkfWN$m[62:182,57-50], 
                                             enkfWN$m[62:182,133-50],
                                             enkfWN2$m[62:182,57-50],
                                             enkfWN2$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 50 members B",121), 
                                           rep( "EnKf with 50 members D", 121),
                                           rep( "EnKf 2 with 50 members B", 121),
                                           rep( "EnKf 2 with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","blue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")



#### Suposem que tenim ara les dades experimentals al punt C i volem estimar què passa en A
#### Fem laminar WN també
#### Provarem els dos casos tal i com s'ha fet prèviament


laminarV <- LL$vel037[50:139]
ns <- 120 #number of time itereations
nc <- 140-50 #number of nodes
## Specify noises (variances)
wk <- 5e-5 #process noise
#vk <- 5e-1 #measurement noise
vk <- sd(velC)
##Construct the state transition matrix
A <- matrix(0, ncol=nc, nrow=nc)
#Entries of the interior nodes in the transition matrix
# identity for a white noise
for(i in 2:(nc-1)) {
  A[i, i] <- 1
}

np <- 2 #number of measurement points
### Simulated points
x <- matrix(mean(laminarV), ncol=nc, nrow=ns)
x[1,] <- laminarV
set.seed(12)
for (i in 2:ns) {
  x[i, ] <-  (A%*%x[i-1, ])+ rnorm(nc, 0, sqrt(wk))}
x
dataXT <- as.data.frame(t(x))
dataXT <- dataXT[2:nc-1,]
g.dlm <- ggplot(data =dataXT, aes(x= index(dataXT[,1]), y=dataXT[,1] )) + geom_line(color= "navy")+
  geom_line(aes( y=dataXT[,30]), color = "coral")+
  geom_line(aes( y=dataXT[,60]), color = "turquoise")+
  geom_line(aes( y=dataXT[,90]), color = "green")+
  geom_line(aes( y=dataXT[,120]), color = "pink")+
  ggtitle("Simulated velocity for different times")+
  xlab("Node (cm)") + ylab("Velocity (m/s)")
g.dlm

ex1 <- list(m0=laminarV , #initial state estimates
            #error covariances of the initial state estimates:
            C0=diag(rep(0.1, nc)),
            #measurement noise
            V=diag(rep(vk, np)),
            #process noise
            W=diag(rep(wk, nc)))

#State transition function
GGfunction <- function (x, k){
  (A%*%x)}
#Observation/measurement function
FFfunction <- function (x, k){
  # cbind(x[57-50], x[133-50])}
  #cbind(x[7], x[83])
  cbind(x[1], x[83])
}

dataEx5 <- as.data.frame(cbind(velE, velC))

enkf1 <- EnKF(y=dataEx5, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx5, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx5, mod=ex1, size=50,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx5, mod=ex1, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)

plot(1:nc, laminarV, type="l", col= "black", #col=c(gray(level=.5)),
     ylim=range(c(enkf1$m[121,], enkf2$m[121,], enkf3$m[121,], enkf4$m[121,])),
     xlab="Node, i", ylab="Velocity (i)", main="t=2")
lines(1:nc, enkf1$m[121,], lty=2, col="blue", lwd=1)
lines(1:nc, enkf2$m[121,], lty=2, col="red", lwd=1)
lines(1:nc, enkf3$m[121,], lty=2, col="darkgreen", lwd=1)
lines(1:nc, enkf4$m[121,], lty=2, col="turquoise", lwd=1)
points(c(57-50, 133-50), dataEx2[120,], col="coral", cex=1, pch = 20)
legend("bottomleft", lty=c(1, 2, 2, 2),
       col=c("black", "blue", "red", "darkgreen", "turquoise"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 25 members", "EnKf with 50 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

y <- as.data.frame(laminarV)

#Error plot
par(mfrow=c(1, 1))
e1 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf3$m[i+1,])^2))
e4 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf4$m[i+1,])^2))

rangeError <- range(cbind(e1, e2, e3, e4))
plot(1:(nrow(y)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(y)-10), e2, lty=1, col="red")
lines(1:(nrow(y)-10), e3, lty=1, col="darkgreen")
lines(1:(nrow(y)-10), e4, lty=1, col="turquoise")
legend("topleft", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen", "turquoise"),
       legend=c("EnKf with 10 members",
                "EnKf with 25 members","EnKf with 50 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

enkfWNC <- enkf3

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf2$m[62:182,57-50], 
                                             enkf2$m[62:182,133-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 25 members B",121), 
                                           rep( "EnKf with 25 members D", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","deepskyblue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

## Mateix resultat horrible

FFfunction <- function (x, k){
  # cbind(x[57-50], x[133-50])}
  #cbind(x[7], x[83])
  cbind(x[83], x[84])
}

dataEx6 <- as.data.frame(cbind(velC, velC))

enkf1 <- EnKF(y=dataEx6, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx6, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx6, mod=ex1, size=50,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx6, mod=ex1, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)

y <- as.data.frame(laminarV)

#Error plot
par(mfrow=c(1, 1))
e1 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf3$m[i+1,])^2))
e4 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf4$m[i+1,])^2))

rangeError <- range(cbind(e1, e2, e3, e4))
plot(1:(nrow(y)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(y)-10), e2, lty=1, col="red")
lines(1:(nrow(y)-10), e3, lty=1, col="darkgreen")
lines(1:(nrow(y)-10), e4, lty=1, col="turquoise")
legend("topleft", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen", "turquoise"),
       legend=c("EnKf with 10 members",
                "EnKf with 25 members","EnKf with 50 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

enkfWNC2 <- enkf3


experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf2$m[62:182,57-50], 
                                             enkf2$m[62:182,133-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 25 members B",121), 
                                           rep( "EnKf with 25 members D", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","deepskyblue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

## Aquest completament horrible xd

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkfWNC$m[62:182,57-50], 
                                             enkfWNC$m[62:182,133-50],
                                             enkfWNC2$m[62:182,57-50],
                                             enkfWNC2$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 50 members B",121), 
                                           rep( "EnKf with 50 members D", 121),
                                           rep( "EnKf 2 with 50 members B", 121),
                                           rep( "EnKf 2 with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","blue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkfWN2$m[62:182,57-50], 
                                             enkfWN2$m[62:182,133-50],
                                             enkfWNC2$m[62:182,57-50],
                                             enkfWNC2$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with B, B",121), 
                                           rep( "EnKf with B, D", 121),
                                           rep( "EnKf with D, B", 121),
                                           rep( "EnKf with D, D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","blue",  "deeppink","navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_linetype_manual(values=c("solid", "dotdash", "dotdash","solid","dashed","dashed"))+
  #ggtitle("Experimental and EnKF filtered velocity") + 
  xlab("t (s)") + ylab("Velocity (m/s)")

