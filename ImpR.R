
############################################
#### REGULAR KF AND ENKF IMPLEMENTATION ####
############################################


library(dlm)
library(forecast)
library(ggplot2)
library(zoo)
library(stats)
library(readr)
library(lmtest)


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
#mean(velC) #0.968232
#mean(velE) #1.881989

#ts.plot(velC)
#ts.plot(velE)
auto.arima(velE)  # WN !! :)
auto.arima(velC)  # WN !! :)

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


## Longitudinal velocity prfile over y=0.05m axis, for 0.0399m and 0.04m diameter
DataL <- read.csv("Dif diametresT - Full 1.csv", dec=",", header=TRUE, stringsAsFactors=FALSE)
#velT <- as.data.frame(DataL$X.5[2:152])
velT <- as.data.frame(DataL$X.6[2:152]*10^(-16))
colnames(velT) <- "LongT"
## results are almost equivalents

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


## longitudinal simulated data
Long <- data.frame(x = index(velT$LongT), 
                   values= c(velT$LongT, Laminar$vel04,
                             Turbulents$vel04), 
                   Data= c(rep("Turbulent mean flow", 151), 
                           rep("Laminar flow 0.04m", 151), 
                           rep("Turbulent flow 0.04m", 151)))
ggplot(Long, aes(x, values, col= Data))+geom_line(size=0.5)+
  scale_color_manual(values= c("coral", "navy", "turquoise")) +
  geom_point(data= exp, aes(x = V1, y=V2), color= "deeppink")+
  #theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Simulated velocity") + xlab("x (m)") + ylab("Velocity (m/s)")

## Long rescaled
LongD22 <- data.frame(x = index(LD$vel039904), 
                      values= c(#LD$vel04065, 
                        LD$vel04,
                        #mm$mitjana1,
                        #Laminar$vel04, 
                        #Laminar$vel035, 
                        #mmL$mitjana,
                        mm$mitjana045, 
                        #mm$mitjanaX2
                        #Turbulents$vel045
                        #LL$vel04,
                        #LL$vel036,
                        LL$vel037
                      ), 
                      Data= c(#rep("Turbulent flow, y=0.065", 151),
                        rep("Turbulent flow", 151),
                        #rep("Turbulent flow mean", 151),
                        #rep("Laminar flow mean", 151),
                        #rep("Laminar flow 0.35m", 151),
                        #rep("Laminar flow mean", 151),
                        rep("Turbulent flow mean", 151),
                        #rep("Turbulent flow mean 0.045 2", 151),
                        #rep("Laminar flow", 151), 
                        #rep("Laminar flow 0.036", 151),
                        rep("Laminar flow", 151)
                      ))
ggplot(LongD22, aes(x, values, col= Data))+geom_line(size=0.5)+
  #scale_color_manual(values= c("darkred", "coral", "orange", "navy", "turquoise")) +
  scale_color_manual(values= c("coral", "navy", "turquoise")) +
  geom_point(data= exp, aes(x = V1, y=V2), color= "deeppink")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Simulated longitudinal flow velocity") + 
  xlab("x (m)") + ylab("Velocity (m/s)")+
  ylim(0, 2)


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



Laminar03Filter <- data.frame(x = index(velA[60:179]), 
                              values= c(velA[59:178], velC[59:178],
                                        #PTu$velA05[59:178], 
                                        PTu$velC05[59:178],
                                        #PTu$velC06[59:178],
                                        #PTu$velA06[59:178],
                                        mmA$mitjana[59:178],
                                        PLa$velA037[59:178],
                                        PLa$velC037[59:178]
                              ),
                              
                              Data= c(rep("True velocity B", 120), 
                                      rep( "True velocity D",120), 
                                      #rep("Simulated velocity B Turbulent", 120),
                                      rep("D Turbulent", 120),
                                      #rep("Simulated velocity B Laminar 0.04", 120),
                                      #rep("Simulated velocity D Laminar 0.04", 120), 
                                      #rep("Simulated velocity D Turbulent shifted", 120),
                                      #rep("Simulated velocity B Turbulent shifted", 120),
                                      #rep("Simulated velocity B Turbulent mean", 120),
                                      rep("B Turbulent", 120),
                                      rep("B Laminar", 120),
                                      rep("D Laminar", 120)))
ggplot(Laminar03Filter, aes(x, values, col= Data))+geom_line(size=0.5)+
  scale_color_manual(values= c("turquoise"," blue", "red",  "coral", "navy", "darkred")) +
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Simulated punctual velocity profile over time") + 
  xlab("t (s)") + ylab("Velocity (m/s)")



######################
### Regular KF #######
######################

###############
## KF A Laminar 
auto.arima(PLa$velA037) ## ARIMA(2,2,0) 

m4 <- arima(PLa$velA037, order=c(3,0,0), include.mean = TRUE)
m4
## log likelihood = 1005.9 (no mean? 1.0087 pm 4.0947 ?? ) 
m2 <- arima(PLa$velA037, order=c(2,1,0), include.mean = TRUE)
m2

##ARIMA(2,1,0)
p1 <- 1.5606
p2 <- -0.5620
m1.dlm <- dlmModARMA( ar= c(p1+1,p2-p1, -p2))
m1.dlm
model.filteredAL <- dlmFilter(velA, m1.dlm)
model.filteredAL$f[4:181]

## ARIMA(3,0,0)
meanA <- mean(PLa$velA037)
mA1 <- arima(PLa$velA037- meanA, order=c(3,0,0), include.mean = FALSE)
mA1
## log likelihood = 1005.91

p1 <- 2.5608
p2 <- -2.1231
p3 <- 0.5623
mA1.dlm <- dlmModARMA( ar= c(p1,p2, p3))
mA1.dlm
model.filteredAL1 <- dlmFilter(velA, mA1.dlm)
model.filteredAL1$f[4:181]

###############
## KF C Laminar
auto.arima(PLa$velC037) ## ARIMA(1,2,1) 
# log likelihood=991.84
mm1 <- arima(PLa$velC037, order=c(2,1,0), include.mean = TRUE)
mm1

## ARIMA(2,1,0)
p1 <- 1.5959
p2 <- -0.5974
m2.dlm <- dlmModARMA( ar= c(p1+1,p2-p1, -p2))
m2.dlm
model.filteredCL <- dlmFilter(velC, m2.dlm)
model.filteredCL$f[4:181]

## ARIMA(3,0,1)
meanC <- mean(PLa$velC037)
auto.arima(PLa$velC037- meanC)
mC1 <- arima(PLa$velC037- meanC, order=c(3,0,1), include.mean = FALSE)
mC1
p1 <- 2.9559
p2 <- -2.9126
p3 <- 0.9567
ma1 <- -0.5664
mC.dlm <- dlmModARMA( ar= c(p1,p2, p3), ma=ma1)
mC.dlm
model.filteredCL1 <- dlmFilter(velC, mC.dlm)


experimentalFilter  <- data.frame(x = index(velA[60:179]), 
                                  values = c(velA[58:177], model.filteredAL1$f[59:178], velC[58:177], 
                                             model.filteredCL$f[59:178],
                                             model.filteredCL1$f[59:178],
                                             #model.filteredCL.wsim[59:178],
                                             #model.filteredCL.wexp[59:178],
                                             #model.filteredAL.wexp[59:178],
                                             PLa$velA037[58:177], PLa$velC037[58:177]
                                  ), 
                                  Data = c(rep("True velocity B", 120), 
                                           rep( "FS B ARIMA(3,0,0)", 120),
                                           rep( "True velocity D",120), 
                                           rep( "FS D ARIMA(2,1,0)", 120),
                                           rep( "FS D ARIMA(3,0,1)", 120),
                                           #rep( "Filtered state A ARIMA(2,1,0) (shifted)", 120),
                                           rep("Simulated velocity B", 120), 
                                           rep("Simulated velocity D", 120)) )
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "deepskyblue","coral", "pink", "blue", "deeppink",  "navy", "darkred")) +
  scale_linetype_manual(values=c("solid","solid", "solid","dotdash",   "dotdash","dashed","dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")


#################
## KF B Turbulent

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


#################
## KF D Turbulent

auto.arima(PTu$velC06[59:178]) ## MA(1) with mean
## log likelihood=457.83
mod3 <- arima(PTu$velC05[59:178],  order = c(2,0,0), include.mean = TRUE)
mod3
m4.dlm <- dlmModARMA( ar= c(0.8807,-0.3596))
m4.dlm
model.filteredC3 <- dlmFilter(velC - mean(velC), m4.dlm)
C3mean <-  model.filteredC3$f + mean(velC)


experimentalFilter  <- data.frame(x = index(velA[60:179]), 
                                  values = c(velA[58:177], 
                                             #model.filteredA3.wmean[59:178], 
                                             velC[58:177], 
                                             #model.filteredC3.wmean[59:178],
                                             C3mean[59:178],
                                             A3mean[59:178],
                                             mmA$mitjana[59:178], PTu$velC06[59:178]), 
                                  Data = c(rep("True velocity B", 120), 
                                           #rep( "Filtered state A MA(2) w intercept", 120),
                                           rep( "True velocity D",120), 
                                           #rep( "Filtered State C AR(2) w intercept", 120),
                                           rep( "FS D AR(2) w/mean", 120),
                                           rep( "FS B MA(2) w/mean", 120),
                                           rep("Simulated velocity B", 120), 
                                           rep("Simulated velocity D", 120)) )
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "blue",  "red", "deepskyblue", "coral", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c( "solid","solid","dotdash","dashed", "dotdash", "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")



############
### ENKF ###
############


library(mvtnorm)
source("http://www.datall-analyse.nl/R/EnKF.R")
EnKF



##########
## Laminar



###
### WN
###

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
  cbind(x[7], x[83])}

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

enkfWN <- enkf2

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


###
### AR(2)
###


laminarV <- LL$vel037[50:139]
auto.arima(laminarV) ## ARIMA(0,0,2) with non-zero mean 
##log likelihood=381.63
mod1 <- arima(laminarV, order=c(2,0,0), include.mean = TRUE )
mod1
##log likelihood = 381.34 AR(2)
##log likelihood = 379.65 ARMA(1,1)

ns <- 120 #number of time itereations
nc <- 140-50 #number of nodes
## Specify noises (variances)
wk <- 5e-5 #process noise
#vk <- 5e-1 #measurement noise
vk <- sd(velC)


## EnKF AR(2)
m1.enkf <- dlmModARMA( ar= c(0.8764,-0.3112))
m1.enkf
inter <- 1.0422

p1 <- 0.8764
p2 <- -0.3112

## AR(2) w intercept
even_rows <- seq_len(nc)%%2
B=matrix(0, ncol=nc, nrow=nc)
B[even_rows==0,even_rows==0] <- p1
A <- diag(B[2,])
for(i in 1:(nc-1)) {
  if(A[i,i] == 0){
    A[i, i+1] <- 1
  }
  else if(A[i,i] == p1){
    A[i, i-1] <- p2
  }
}

x <- matrix(inter, ncol=nc, nrow=ns)
x[1,] <- laminarV
set.seed(12)
for (i in 3:ns) {
  x[i, ] <-  (A%*%x[i-1, ])+ rnorm(nc, 0, sqrt(wk))}
x

dataXT <- as.data.frame(t(x))
dataXT <- dataXT[2:nc-1,]

dataXT2 <- dataXT[2:(nc-5),]

g.dlm <- ggplot(data =dataXT2, aes(x= index(dataXT2[,1]), y=dataXT2[,1] )) + geom_line(color= "navy")+
  geom_line(aes( y=dataXT2[,30]+inter), color = "coral")+
  geom_line(aes( y=dataXT2[,60]+inter), color = "turquoise")+
  geom_line(aes( y=dataXT2[,90]+inter), color = "green")+
  geom_line(aes( y=dataXT2[,120]+inter), color = "pink")+
  ggtitle("Simulated velocity for different times")+
  xlab("Node (cm)") + ylab("Velocity (m/s)")
g.dlm

np <- 2
ex1 <- list(m0=laminarV , #initial state estimates
            #error covariances of the initial state estimates:
            C0=diag(rep(0.1, nc)),
            #measurement noise
            V=diag(rep(vk, np)),
            #process noise
            W=diag(rep(wk, nc)))

#State transition function
GGfunction <- function (x, k){
  #  0.9834*vk+(A%*%x)}
  #0.9834*x+(A%*%x)}
  (A%*% x)}
#Observation/measurement function
FFfunction <- function (x, k){
  cbind(x[7], x[83])}

dataEx2 <- as.data.frame(cbind(velA, velC))

enkf1 <- EnKF(y=dataEx2, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx2, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx2, mod=ex1, size=25,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx2, mod=ex1, size=50,
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

#y <- as.data.frame(Laminar$vel035[50:139])
y <- as.data.frame(laminarV)

#Error plot respecte de les simulades
par(mfrow=c(1, 1))
e1 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf1$m[i+1,]+inter)^2))
e2 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf2$m[i+1,]+inter)^2))
e3 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf3$m[i+1,]+inter)^2))
e4 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,]-enkf4$m[i+1,]+inter)^2))

rangeError <- range(cbind(e1, e2, e3, e4))
plot(1:(nrow(y)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(y)-10), e2, lty=1, col="red")
lines(1:(nrow(y)-10), e3, lty=1, col="darkgreen")
lines(1:(nrow(y)-10), e4, lty=1, col="turquoise")
legend("topleft", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen", "turquoise"),
       legend=c("EnKf with 10 members",
                "EnKf with 20 members","EnKf with 25 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

enkfAR2 <- enkf2

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf2$m[62:182,57-50]+mean(velA[61:181]), 
                                             enkf2$m[62:182,133-50]+mean(velC[61:181]),
                                             enkfWN$m[62:182,57-50],
                                             enkfWN$m[62:182,133-50]
                                             #enkf3$m[62:182,57-50]+mean(velA[61:181]),
                                             #enkf3$m[62:182,133-50]+mean(velC[61:181])
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf B AR(2)",121), 
                                           rep( "EnKf D AR(2)", 121),
                                           rep( "EnKf B WN", 121),
                                           rep( "EnKf D WN", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "blue", "coral",  "red", "navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dotdash", "dashed"))+
  #ggtitle("Experimental and EnKF filtered velocity") + 
  xlab("t (s)") + ylab("Velocity (m/s)")




############
## Turbulent


turb <- mm$mitjana045[50:139]


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
x <- matrix(mean(turb), ncol=nc, nrow=ns)
x[1,] <- turb
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

ex1 <- list(m0=turb , #initial state estimates
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
  cbind(x[7], x[83])}

dataEx2 <- as.data.frame(cbind(velA, velC))

enkf1 <- EnKF(y=dataEx2, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx2, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx2, mod=ex1, size=25,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx2, mod=ex1, size=50,
              GGfunction=GGfunction, FFfunction=FFfunction)

plot(1:nc, turb, type="l", col= "black", #col=c(gray(level=.5)),
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

y <- as.data.frame(turb)

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

enkfTWN <- enkf3



###
### ARIMA(2,0,0)
###

#turb <- LD$vel04[50:139]
turb <- mm$mitjana045[50:139]
auto.arima(turb) ## ARIMA(5,2,0) really ARIMA(2,2,0) 
##log likelihood=652.06
mod1 <- arima(turb, order=c(2,1,0), include.mean = TRUE )
mod1
## log likelihood = 621.2 ARIMA(2,1,0)
##log likelihood = 627.44 ARima(2,2,0)
##log likelihood = 624.29 arima(3,1,0)
##log likelihood = 636.1 ARiMA(3,2,0)
##log likelihood = 648.34 ARIMA(4,2,0)

ns <- 120 #number of time itereations
nc <- 140-50 #number of nodes
## Specify noises (variances)
wk <- 5e-5 #process noise
#vk <- 5e-1 #measurement noise
vk <- sd(velC)


## EnKF AR(2)
m1.enkf <- dlmModARMA( ar= c(0.8764,-0.3112))
m1.enkf
inter <- 1.0422

p1 <- 0.4341
p2 <- 0.5654

## AR(2) w intercept
even_rows <- seq_len(nc)%%2
B=matrix(0, ncol=nc, nrow=nc)
B[even_rows==0,even_rows==0] <- p1
A <- diag(B[2,])
for(i in 1:(nc-1)) {
  if(A[i,i] == 0){
    A[i, i+1] <- 1
  }
  else if(A[i,i] == p1){
    A[i, i-1] <- p2
  }
}

x <- matrix(mean(turb), ncol=nc, nrow=ns)
x[1,] <- turb
set.seed(12)
for (i in 3:ns) {
  x[i, ] <-  (A%*%x[i-1, ])+ rnorm(nc, 0, sqrt(wk))}
x

dataXT <- as.data.frame(t(x))
dataXT <- dataXT[2:nc-1,]

dataXT2 <- dataXT[2:(nc-5),]

g.dlm <- ggplot(data =dataXT2, aes(x= index(dataXT2[,1]), y=dataXT2[,1] )) + geom_line(color= "navy")+
  geom_line(aes( y=dataXT2[,30]), color = "coral")+
  geom_line(aes( y=dataXT2[,60]), color = "turquoise")+
  geom_line(aes( y=dataXT2[,90]), color = "green")+
  geom_line(aes( y=dataXT2[,120]), color = "pink")+
  ggtitle("Simulated velocity for different times")+
  xlab("Node (cm)") + ylab("Velocity (m/s)")
g.dlm

## bastant potato

np <- 2
ex1 <- list(m0=turb , #initial state estimates
            #error covariances of the initial state estimates:
            C0=diag(rep(0.1, nc)),
            #measurement noise
            V=diag(rep(vk, np)),
            #process noise
            W=diag(rep(wk, nc)))

#State transition function
GGfunction <- function (x, k){
  #  0.9834*vk+(A%*%x)}
  #0.9834*x+(A%*%x)}
  (A%*% x)}
#Observation/measurement function
FFfunction <- function (x, k){
  cbind(x[7], x[83])}

dataEx2 <- as.data.frame(cbind(velA, velC))

enkf1 <- EnKF(y=dataEx2, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf2 <- EnKF(y=dataEx2, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf3 <- EnKF(y=dataEx2, mod=ex1, size=25,
              GGfunction=GGfunction, FFfunction=FFfunction)

enkf4 <- EnKF(y=dataEx2, mod=ex1, size=50,
              GGfunction=GGfunction, FFfunction=FFfunction)

plot(1:nc, turb, type="l", col= "black", #col=c(gray(level=.5)),
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

#y <- as.data.frame(Laminar$vel035[50:139])
y <- as.data.frame(turb)

#Error plot respecte de les simulades
par(mfrow=c(1, 1))
e1 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,1]-enkf1$m[i+1,]+inter)^2))
e2 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,1]-enkf2$m[i+1,]+inter)^2))
e3 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,1]-enkf3$m[i+1,]+inter)^2))
e4 <- sapply(1:(nrow(y)-10), function (i) mean((y[i,1]-enkf4$m[i+1,]+inter)^2))

rangeError <- range(cbind(e1, e2, e3, e4))
plot(1:(nrow(y)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(y)-10), e2, lty=1, col="red")
lines(1:(nrow(y)-10), e3, lty=1, col="darkgreen")
lines(1:(nrow(y)-10), e4, lty=1, col="turquoise")
legend("topleft", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen", "turquoise"),
       legend=c("EnKf with 10 members",
                "EnKf with 20 members","EnKf with 25 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

## Error plot respecte de les experimentals! 

enkfTurb <- enkf3


### tots els Enkf turbulents
TubrulentFilter  <- data.frame(x = index(velA[61:181]), 
                               values = c(velA[59:179], velC[59:179], 
                                          enkf2$m[62:182,57-50], 
                                          enkf2$m[62:182,133-50],
                                          enkfTWN$m[62:182,57-50],
                                          enkfTWN$m[62:182,133-50]
                                          #enkf3$m[62:182,57-50],
                                          #enkf3$m[62:182,133-50]
                                          #enkf2$m[62:182,57-50]+inter, 
                                          #enkf2$m[62:182,133-50]+inter,
                                          #enkf3$m[62:182,57-50]+inter,
                                          #enkf3$m[62:182,133-50]+inter
                               ), 
                               Data = c(rep("True velocity B", 121), 
                                        rep( "True velocity D", 121),
                                        rep( "EnKf B AR(2) ",121), 
                                        rep( "EnKf D AR(2)", 121),
                                        rep( "EnKf B WN", 121),
                                        rep( "EnKf D WN", 121)))
ggplot(TubrulentFilter, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise","blue",  "coral", "red","navy", "darkred")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid", "dotdash", "dotdash"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Experimental and EnKF filtered velocity") + 
  xlab("t (s)") + ylab("Velocity (m/s)")


## comparació A tots (seleccionant els millors models)
CompFilterA  <- data.frame(x = index(velA[60:179]), 
                           values = c(velA[58:177], 
                                      enkfTurb$m[60:179,57-50],
                                      #model.filteredAL1$f[59:178], 
                                      A3mean[59:178],
                                      #PLa$velA037[58:177], 
                                      #PTu$velA06[58:177],
                                      enkfTWN$m[60:179,57-50],
                                      enkfWN$m[60:179,57-50],
                                      enkfAR2$m[60:179,57-50]+mean(velA[61:180])
                           ), 
                           Data = c(rep("True velocity B", 120), 
                                    rep("EnKF ARIMA(2,1,0)", 120),
                                    #rep( "Filtered state KF laminar AR(3) with mean", 120),
                                    rep( "KF tur MA(2)", 120),
                                    #rep("Simulated velocity laminar", 120), 
                                    #rep("Simulated velocity turbulent", 120), 
                                    rep("EnKF WN", 120),
                                    rep("EnKF lam WN", 120),
                                    #rep("Filtered state EnKF laminar WN", 120),
                                    rep("EnKF lam AR(2)", 120)
                           ) )
ggplot(CompFilterA, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #scale_color_manual(values= c( "blue", "turquoise", "green", "violet", "darkred")) +
  #scale_color_manual(values= c( "blue", "turquoise", "darkred","red","black")) +
  scale_color_manual(values= c( "coral", "darkred","red", "blue" , "turquoise","navy")) +
  scale_linetype_manual(values=c( "solid", "solid", "solid", "solid", "solid", "dotdash"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity A") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")


## Error tots
errorA1 <- data.frame(x = index(velA[61:180]), 
                      values = c((PLa$velA037[59:178]- velA[58:177])^2,
                                 #(PTu$velA06[58:177]- velA[58:177])^2,
                                 #(model.filteredAL1$f[59:178]-velA[58:177])^2,
                                 (enkfTurb$m[60:179,57-50]- velA[58:177])^2,
                                 (enkfWN$m[60:179,57-50]  - velA[58:177])^2,
                                 (enkfTWN$m[60:179,57-50]  - velA[58:177])^2,
                                 (enkfAR2$m[60:179,57-50] + mean(velA) - velA[58:177])^2,
                                 (A3mean[59:178]- velA[58:177])^2
                      ), 
                      Data = c(rep("Simulated", 120), 
                               #rep("True state turbulent - simulated", 120),
                               rep( "EnKF ARIMA(2,1,0)", 120),
                               rep( "EnKF WN",120),
                               rep( "EnKF turb WN",120),
                               rep( "EnKF MA(2)", 120),
                               rep( "KF MA(2)",120) 
                      ))

ggplot(errorA1, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  #scale_color_manual(values= c( "deeppink",  "blue","coral",  "navy", "darkred")) +
  #scale_color_manual(values= c( "coral", "darkred","red", "blue" ,"navy")) +
  #scale_color_manual(values= c( "turquoise", "darkred", "red", "coral", "blue" ,"navy")) +
  scale_color_manual(values= c( "turquoise", "darkred", "red", "violet", "blue" ,"darkgrey")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash","dotdash", "solid",  "solid"))+
  #scale_linetype_manual(values=c("solid", "solid", "solid", "solid",  "dashed"))+
  #ggtitle(" Mean squared error KF and EnKF filtered velocity A") + 
  #scale_y_log10()+
  xlab("t (s)") + ylab("Velocity (m/s)")

ggplot(errorA1, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  #scale_color_manual(values= c( "turquoise", "darkred", "red", "coral", "blue" ,"navy")) +
  scale_color_manual(values= c( "turquoise", "darkred", "red", "violet", "blue" ,"darkgrey")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash","dotdash", "solid",  "solid"))+
  #scale_linetype_manual(values=c("solid", "solid", "solid", "solid",  "dashed"))+
  #ggtitle(" Mean squared error KF and EnKF filtered velocity A") + 
  scale_y_log10()+
  xlab("t (s)") + ylab("Velocity (m/s)")


mean((enkfWN$m[60:179,57-50]  - velA[58:177])^2)
mean((enkfTWN$m[60:179,57-50]  - velA[58:177])^2)
mean((enkfAR2$m[60:179,57-50] + mean(velA) - velA[58:177])^2)
mean((A3mean[59:178]- velA[58:177])^2)
mean((enkfTurb$m[60:179,57-50]- velA[58:177])^2)

## WN laminar tmb

### Error plots escala logaritmica
ErrorSimA <- log10((PLa$velA037[59:178]- velA[58:177])^2)
ErrorSimA
ErrorARIMA2101 <-  log10((enkfTurb$m[60:179,57-50]- velA[58:177])^2)
ErrorWNA <- log10((enkfWN$m[60:179,57-50]  - velA[58:177])^2)
ErrorAR2A <- log10((enkfAR2$m[60:179,57-50] + mean(velA) - velA[58:177])^2)
ErrorKFA <- log10((A3mean[59:178]- velA[58:177])^2)
ErrorKFA
(A3mean[59:178]- velA[58:177])^2
(PLa$velA037[59:178]- velA[58:177])^2

### Errors A escala logarítmica
errorA2 <- data.frame(x = index(velA[61:180]), 
                      values = c(ErrorSimA,
                                 ErrorARIMA2101,
                                 ErrorWNA, 
                                 ErrorAR2A, 
                                 ErrorKFA
                                 ), 
                      Data = c(rep("Simulated", 120), 
                               rep( "EnKF ARIMA(2,1,0)", 120),
                               rep( "EnKF AR(2)", 120),
                               rep( "EnKF WN",120),
                               rep( "KF AR(2)",120) 
                      ))

ggplot(errorA2, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  #scale_color_manual(values= c( "deeppink",  "blue","coral",  "navy", "darkred")) +
  #scale_color_manual(values= c( "coral", "darkred","red", "blue" ,"navy")) +
  scale_color_manual(values= c(  "darkred", "red", "turquoise", "blue" ,"navy")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash", "solid",  "dashed"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid",  "dashed"))+
  #ggtitle(" Mean squared error KF and EnKF filtered velocity A") + 
  xlab("t (s)") + ylab("Velocity (m/s)")


### Error plot amb escala logarítimica no sembla que aporti massa informació ni que 
### millori la visualització dels resultats


#### Caldria comprovar el plot d'errors del punt D

## comparació C òptims
CompFilterC1  <- data.frame(x = index(velA[60:179]), 
                            values = c(velC[58:177], 
                                       enkfTurb$m[60:179,133-50],
                                       #model.filteredAL1$f[59:178], 
                                       C3mean[59:178],
                                       #PLa$velA037[58:177], 
                                       #PTu$velA06[58:177],
                                       enkfTWN$m[60:179,133-50],
                                       enkfWN$m[60:179,133-50],
                                       enkfAR2$m[60:179,133-50]+mean(velC[61:180])
                            ), 
                            Data = c(rep("True velocity D", 120), 
                                     rep("EnKF turb ARIMA(2,1,0)", 120),
                                     #rep( "Filtered state KF laminar AR(3) with mean", 120),
                                     rep( "KF turb AR(2)", 120),
                                     #rep("Simulated velocity laminar", 120), 
                                     #rep("Simulated velocity turbulent", 120), 
                                     rep("EnKF WN", 120),
                                     rep("EnKF lam WN", 120),
                                     #rep("Filtered state EnKF laminar WN", 120),
                                     rep("EnKF lam AR(2)", 120)
                            ) )
ggplot(CompFilterC1, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #scale_color_manual(values= c( "blue", "navy", "coral", "deeppink", "darkred")) +
  scale_color_manual(values= c( "coral", "darkred", "red", "blue", "turquoise","navy")) +
  #scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash", "dotdash", "solid"))+
  scale_linetype_manual(values=c("solid", "solid", "solid","solid", "solid",  "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity C") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")

ggplot(CompFilterC1, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #scale_color_manual(values= c( "blue", "navy", "coral", "deeppink", "darkred")) +
  scale_color_manual(values= c( "coral", "darkred", "red", "blue", "turquoise","navy")) +
  #scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash", "dotdash", "solid"))+
  scale_linetype_manual(values=c("solid", "solid", "solid","solid", "solid",  "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity C") + 
  xlab("Time (s)") + ylab("Velocity (m/s)")


### Els wn del turbulent i del laminar són iguals (estem aplicant un wn igualment, no cal possar-lo dos cops
### als 2 mins ja no fa tanta diferència la condició inicial. Ens quedem amb el WN laminar per fer l'última)

## La resta de diferències seran equivalents a les del wn anterior

## Error tots
errorC1 <- data.frame(x = index(velA[61:180]), 
                      values = c(#(PLa$velC037[59:178]- velC[58:177])^2,
                        (PTu$velC06[58:177]- velC[58:177])^2,
                        #(model.filteredAL1$f[59:178]-velA[58:177])^2,
                        (enkfTurb$m[60:179,133-50]- velC[58:177])^2,
                        (enkfWN$m[60:179,133-50]  - velC[58:177])^2,
                        (enkfAR2$m[60:179,133-50] + mean(velC) - velC[58:177])^2,
                        (C3mean[59:178]- velC[58:177])^2
                      ), 
                      Data = c(#rep("True state laminar - simulated", 120), 
                        rep("Turbulent simulated", 120),
                        rep( "EnKF ARIMA(2,1,0)", 120),
                        rep( "EnKF AR(2)", 120),
                        rep( "EnKF WN",120),
                        rep( "KF AR(2)",120) 
                      ))

ggplot(errorC1, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  #scale_color_manual(values= c( "deeppink",  "blue","coral", "violet",  "navy")) +
  scale_color_manual(values= c( "darkred", "red", "deeppink",  "blue", "navy")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash", "dotdash", "solid",  "solid"))+
  scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash", "solid", "dashed"))+
  #ggtitle(" Mean squared error KF and EnKF filtered velocity A") + 
  xlab("t (s)") + ylab("Velocity (m/s)")


## Error only filters
errorC2 <- data.frame(x = index(velA[61:180]), 
                      values = c(#(PLa$velA037[59:178]- velA[58:177])^2,
                        #(PTu$velA06[58:177]- velA[58:177])^2,
                        #(model.filteredAL1$f[59:178]-velA[58:177])^2,
                        (enkfTurb$m[60:179,133-50]- velC[58:177])^2,
                        (enkfWN$m[60:179,133-50]  - velC[58:177])^2,
                        (enkfTWN$m[60:179,133-50]  - velC[58:177])^2,
                        (enkfAR2$m[60:179,133-50] + mean(velC) - velC[58:177])^2,
                        (C3mean[59:178]- velC[58:177])^2
                      ), 
                      Data = c(#rep("True state laminar - simulated", 120), 
                        #rep("True state turbulent - simulated", 120),
                        rep( "EnKF ARIMA(2,1,0)", 120),
                        rep( "EnKF AR(2)", 120),
                        rep( "EnKF lam WN",120),
                        rep( "EnKF WN",120),
                        rep( "KF AR(2)",120) 
                      ))

ggplot(errorC2, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c( "darkred","red", "violet", "turquoise", "blue")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_linetype_manual(values=c("dotdash", "dotdash", "dotdash", "dotdash",  "solid"))+
  #ggtitle(" Mean squared error KF and EnKF filtered velocity A") + 
  scale_y_log10()+
  xlab("t (s)") + ylab("Velocity (m/s)")


mean((enkfTurb$m[60:179,133-50]- velC[58:177])^2) #0.0003381692
mean((enkfWN$m[60:179,133-50]  - velC[58:177])^2) #0.0003636779
mean((enkfTWN$m[60:179,133-50]  - velC[58:177])^2) # 0.000363
mean((enkfAR2$m[60:179,133-50] + mean(velC) - velC[58:177])^2) #0.0004689085
mean((C3mean[59:178]- velC[58:177])^2) #4.543115e-05

## DEl gràfic maybe el millor EnKF és el WN laminar


#### Estudi de les barres d'errors del KF

KFplot  <- data.frame(x = index(velA[60:179]), 
                                  values = c(velA[58:177], 
                                             #model.filteredA3.wmean[59:178], 
                                             velC[58:177], 
                                             #model.filteredC3.wmean[59:178],
                                             C3mean[59:178],
                                             A3mean[59:178]
                                             #mmA$mitjana[59:178], PTu$velC06[59:178]
                                             ), 
                                  Data = c(rep("True velocity B", 120), 
                                           #rep( "Filtered state A MA(2) w intercept", 120),
                                           rep( "True velocity D",120), 
                                           #rep( "Filtered State C AR(2) w intercept", 120),
                                           rep( "FS D AR(2) w/mean", 120),
                                           rep( "FS B MA(2) w/mean", 120)
                                           #rep("Simulated velocity B", 120), 
                                           #rep("Simulated velocity D", 120)
                                           ) 
                                           )
ggplot(KFplot, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "blue",  "red", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c( "solid","solid","dotdash","dashed", "dotdash", "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  geom_pointrange(KFplot, aes(ymin=values-sd, ymax=values+sd))+
  xlab("Time (s)") + ylab("Velocity (m/s)")

# "deepskyblue", "coral",

model.filteredA3 
model.filteredC3   

?dlmFilter
#The functions applies Kalman filter to compute filtered values of the state vectors, 
#together with their variance/covariance matrices. By default the function returns an object 
#of class "dlmFiltered". Methods for residuals and tsdiag for objects of class "dlmFiltered" exist.
  
# D.C Together with U.C, it gives the SVD of the variances of the estimation errors. 
#The variance of m_t-θ_t is given by U.C[[t]] %*% diag(D.C[t,]^2) %*% t(U.C[[t]]).
  
# D.R Together with U.R, it gives the SVD of the variances of the prediction errors. 
#The variance of a_t-θ_t is given by U.R[[t]] %*% diag(D.R[t,]^2) %*% t(U.R[[t]]).

# f: Time series (or matrix) of one-step-ahead forecast of the observations.

# a: Time series (or matrix) of predicted values of the state vectors given the observations up and including the previous time unit.
 
# m Time series (or matrix) of filtered values of the state vectors. The series starts one time unit before the first observation.


model.filteredA3$f[59:178]+ mean(velA)
model.filteredC3$f[59:178]+ mean(velC)


model.filteredA3$a[59:178]+ mean(velA)
model.filteredC3$a[59:178]+ mean(velC)

model.filteredA3$m[59:178]+ mean(velA)
model.filteredC3$m[59:178]+ mean(velC)

## plot for the different estimates

KFplot2  <- data.frame(x = index(velA[60:179]), 
                      values = c(velA[58:177], 
                                 #model.filteredA3.wmean[59:178], 
                                 velC[58:177], 
                                 #model.filteredC3.wmean[59:178],
                                 C3mean[59:178],
                                 A3mean[59:178],
                                 model.filteredA3$a[59:178]+ mean(velA),
                                 model.filteredC3$a[59:178]+ mean(velC), 
                                 model.filteredA3$m[59:178]+ mean(velA),
                                 model.filteredC3$m[59:178]+ mean(velC)
                                 #mmA$mitjana[59:178], PTu$velC06[59:178]
                      ), 
                      Data = c(rep("True velocity B", 120), 
                               #rep( "Filtered state A MA(2) w intercept", 120),
                               rep( "True velocity D",120), 
                               #rep( "Filtered State C AR(2) w intercept", 120),
                               rep( "Forecast D AR(2) w/mean", 120),
                               rep( "Forecast B MA(2) w/mean", 120),
                               rep( "Predicted B AR(2) w/mean", 120),
                               rep( "Predicted D MA(2) w/mean", 120), 
                               rep( "Filtered B AR(2) w/mean", 120),
                               rep( "Filtered D MA(2) w/mean", 120)
                               #rep("Simulated velocity B", 120), 
                               #rep("Simulated velocity D", 120)
                      ) 
)
ggplot(KFplot2, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "deepskyblue", "coral","blue",  "red", "green", "violet", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c("dotdash","dashed", "solid","solid","dotdash", "dashed", "dotdash", "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  #geom_pointrange(KFplot, aes(ymin=values-sd, ymax=values+sd))+
  xlab("Time (s)") + ylab("Velocity (m/s)")

## són essencialment el mateix xd

# D.R Together with U.R, it gives the SVD of the variances of the prediction errors. 
#The variance of a_t-θ_t is given by U.R[[t]] %*% diag(D.R[t,]^2) %*% t(U.R[[t]]).
t(model.filteredA3$U.R)

diag(model.filteredA3$U.R)

U.R <- model.filteredA3$U.R
D.R <- model.filteredA3$D.R

model.filteredA3$mod
# FF= 1,0,0
## prenem el primer element de la matriu de covariances definida en cada cas

variances=rep(0, 120)
for(i in (1:120)){
    variances[i] <- U.R[[i]] %*% diag(D.R[i,]^2) %*% t(U.R[[i]])
}
variances


U.R[[30]]
diag(D.R[30,]^2)
t(U.R[[30]])

U.R[[30]] %*% diag(D.R[30,]^2) %*% t(U.R[[30]])
U.R[[60]] %*% diag(D.R[60,]^2) %*% t(U.R[[60]])
U.R[[10]] %*% diag(D.R[10,]^2) %*% t(U.R[[10]])

## Estaciona?


U.C <- model.filteredA3$U.C
D.C <- model.filteredA3$D.C

U.C[[30]] %*% diag(D.C[30,]^2) %*% t(U.C[[30]])
U.C[[60]] %*% diag(D.C[60,]^2) %*% t(U.C[[60]])
U.C[[10]] %*% diag(D.C[10,]^2) %*% t(U.C[[10]])
U.C[[120]] %*% diag(D.C[120,]^2) %*% t(U.C[[120]])

variances=rep(0, 120)
for(i in (2:120)){
  #variances[i] <- sqrt(U.C[[i]] %*% diag(D.C[i,]^2) %*% t(U.C[[i]]))[3,3]
  variances[i] <- (U.C[[i]] %*% diag(D.C[i,]^2) %*% t(U.C[[i]]))[3,3]
}
variances

variances2=rep(2.617918e-15, 120)



KFplot  <- data.frame(x = index(velA[60:179]), 
                      values = c(velA[58:177], 
                                 #model.filteredA3.wmean[59:178], 
                                 velC[58:177], 
                                 #model.filteredC3.wmean[59:178],
                                 C3mean[59:178],
                                 A3mean[59:178]
                                 #mmA$mitjana[59:178], PTu$velC06[59:178]
                      ), 
                      Data = c(rep("True velocity B", 120), 
                               #rep( "Filtered state A MA(2) w intercept", 120),
                               rep( "True velocity D",120), 
                               #rep( "Filtered State C AR(2) w intercept", 120),
                               rep( "FS D AR(2) w/mean", 120),
                               rep( "FS B MA(2) w/mean", 120)
                               #rep("Simulated velocity B", 120), 
                               #rep("Simulated velocity D", 120)
                      ),
                      sd= variances2
)
ggplot(KFplot, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "blue",  "red", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c( "solid","solid","dotdash","dashed", "dotdash", "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  geom_pointrange(aes(ymin=values-sd, ymax=values+sd))+
  xlab("Time (s)") + ylab("Velocity (m/s)")





#compute the confidence intervals for the filtered state estimates
varcovFilteredState <- dlmSvd2var(model.filteredA3$U.C, model.filteredA3$D.C)
#95% ci for x-position
seFX <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[1,1])))
seFX
ciFX <- model.filteredA3$f[59:178]+ mean(velA) + qnorm(.05/2)*seFX[59:178]%o%c(1, -1)
ciFX
ciFX2 <- model.filteredA3$f[59:178]+ mean(velA) + qnorm(.1/2)*seFX[59:178]%o%c(1, -1)
ciFX2


#x-position
t <- index(velA[58:177])
#t <- seq(0.1, 3, .1)
plot(index(velA[58:177]), velA[58:177], type="l", xlab="time", ylab="Velocity (m/s)",
     #xlim=c(0,3), ylim=c(-5, 85), col=gray(level=.7), 
     lwd=1.5)
lines(t, ciFX[,1],  type="l", col=gray(level=.7))
lines(t, ciFX[,2],  type="l", col=gray(level=.7))
lines(t, A3mean[59:178],  type="l", col="blue")
#arrows(c(A3mean,t), ciFX[,1], c(A3mean,t), ciFX[,2],
       #code=3, length=0.05, angle=90, col="blue")


### La variància també sembla que és molt petita en comparació al que s'esperaria pot ser? 
## el valor experimental no entra dins de l'interval de confiança del 95%

## provem el del 90%
## tampoc es difernecien



#### Estudi residuals wtf?
#raw residuals (innovations)
plot(t, residuals(model.filteredA3, sd=FALSE, type="raw")[59:178], type="o",
     ylab="raw residuals")
abline(h=0, lty=2, col="darkgray")

#standardized residuals
sresx <- residuals(model.filteredA3, sd=FALSE, type="standardized")[59:178]
plot(t, sresx, type="o", ylab="standardized residuals")
abline(h=0, lty=2, col="darkgray")

#qq-plot for standardized residuals
qqnorm(sresx)
abline(a=0, b=1, col="darkgray", lty=2)

#acf plot for standardized residuals
acf(residuals(model.filteredA3, sd=FALSE), main="x-position")


