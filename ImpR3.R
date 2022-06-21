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


##########
## EnKF ##
##########

library(mvtnorm)
source("http://www.datall-analyse.nl/R/EnKF.R")
EnKF

## WN Laminar best approximation for both experimental distributions


laminarV <- LL$vel037[50:139]
auto.arima(laminarV)
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
enkfWN4 <-  enkf4

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf2$m[62:182,57-50], 
                                             enkf2$m[62:182,133-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50],
                                             PLa$velA037[62:182]
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
                                           rep( "EnKf with 50 members D", 121),
                                           rep( "Simulated state laminar", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "coral","deepskyblue",  "deeppink", "blue", "navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")


### Punt intermig 

## Tenim les dades dels punts B i D a 57 cm i 133 cm

## A: 
47.1
## B:
47.1+9.1
## que ha estat 57

## C:
47.1+9.1+64.2
## 120
## D:
120+14.7
## que ha estat 133

## Prove de predir a 120cm
## punt C

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf3$m[62:182,120-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50]
                                             #enkf2$m[62:182,57-50]+inter, 
                                             #enkf2$m[62:182,133-50]+inter,
                                             #enkf3$m[62:182,57-50]+inter,
                                             #enkf3$m[62:182,133-50]+inter
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 50 members C", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("turquoise", "violet", "coral", "navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

## una predicció força horrible 

experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             enkf3$m[62:182,120-50],
                                             enkf3$m[62:182,57-50],
                                             enkf3$m[62:182,133-50],
                                             enkf3$m[62:182,100-50],
                                             enkf3$m[62:182,70-50],
                                             enkf3$m[62:182,60-50],
                                             PLa$velC037[62:182],
                                             mmA$mitjana[62:182]
                                     ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           rep( "EnKf with 50 members 120", 121),
                                           rep( "EnKf with 50 members B", 121),
                                           rep( "EnKf with 50 members D", 121),
                                           rep( "EnKf with 50 members 100", 121),
                                           rep( "EnKf with 50 members 70", 121),
                                           rep( "EnKf with 50 members 60", 121),
                                           rep( "Simulated laminar", 121),
                                           rep( "Simulated turbulent", 121)))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("deeppink", "green", "violet", "gold", "turquoise", "coral", "blue", "pink", "navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Experimental and EnKF filtered velocity") + xlab("t (s)") + ylab("Velocity (m/s)")

## Ajusta una mica millor la simulada el de 120
## podria agafar les dades simulades d'equell punt i comparar-les
## inclús les experimentals?

experimentalFilter  <- data.frame(x = index(laminarV), 
                                  values = c(laminarV, enkf1$m[121,], enkf2$m[121,], 
                                             enkf3$m[121,], enkf4$m[121,] ), 
                                  Data = c(rep("True State", 90), 
                                           rep( "EnKf with 10 members", 90),
                                           rep( "EnKf with 25 members",90), 
                                           rep( "EnKf with 50 members", 90),
                                           rep( "EnKf with 100 members", 90) ))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_color_manual(values= c( "coral","deepskyblue",  "green", "violet","navy")) +
  ggtitle("Simulated and filtered velocity") + xlab("x (m)") + ylab("Velocity (m/s)")


### Longitudinal ajusta molt millor el AR 2


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

enkf2 <- EnKF(y=dataEx2, mod=ex1, size=25,
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

exp2 <- as.data.frame(cbind(c( 57-51, 133-51), c( velA[120], velC[120])))
experimentalFilter  <- data.frame(x = index(laminarV[1:(nc-2)]), 
                                  values = c(laminarV[1:(nc-2)], 
                                             enkf1$m[121,1:(nc-2)]+inter,
                                             enkf2$m[121,1:(nc-2)]+inter, 
                                             enkf3$m[121,1:(nc-2)]+inter, 
                                             enkf4$m[121,1:(nc-2)]+inter
                                             #enkf1$m[121,1:(nc-2)]+inter,
                                             #enkf2$m[121,1:(nc-2)]+inter, 
                                             #enkf3$m[121,1:(nc-2)]+inter, 
                                             #enkf4$m[121,1:(nc-2)]+inter
                                  ), 
                                  Data = c(rep("True State", 88), 
                                           rep( "EnKf with 10 members", 88),
                                           rep( "EnKf with 25 members", 88), 
                                           rep( "EnKf with 50 members", 88),
                                           rep( "EnKf with 100 members", 88) ))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c( "coral","deepskyblue",  "green", "violet","navy")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  ggtitle("Simulated and filtered velocity") + xlab("x (m)") + ylab("Velocity (m/s)")

enkfAR2 <- enkf3

experimentalFilter  <- data.frame(x = index(laminarV[1:(nc-2)]), 
                                  values = c(laminarV[1:(nc-2)], 
                                             enkfAR2$m[121,2:(nc-1)]+inter,
                                             #enkfWN$m[121,2:(nc-1)],
                                             enkfWN4$m[121,2:(nc-1)]
                                             ), 
                                  Data = c(rep("True State", 88), 
                                           rep( "EnKf AR(2)", 88),
                                           rep( "EnKf WN", 88)
                                           #rep( "EnKf WN 2", 88)
                                            ))
ggplot(experimentalFilter, aes(x, values, col= Data))+geom_line(size=0.5)+
  geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c( "coral","deepskyblue", "navy")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Simulated and filtered velocity") + 
  xlab("x (m)") + ylab("Velocity (m/s)")



experimentalFilter  <- data.frame(x = index(velA[61:181]), 
                                  values = c(velA[59:179], velC[59:179], 
                                             #enkfWN4$m[62:182,120-50],
                                             enkfWN4$m[62:182,57-50],
                                             enkfWN4$m[62:182,133-50],
                                             enkfWN4$m[62:182,100-50],
                                             enkfWN4$m[62:182,70-50],
                                             enkfWN4$m[62:182,60-50],
                                             PLa$velC037[59:179]
                                             #mmA$mitjana[62:182]
                                  ), 
                                  Data = c(rep("True velocity B", 121), 
                                           rep( "True velocity D", 121),
                                           #rep( "EnKf 120", 121),
                                           rep( "EnKf B", 121),
                                           rep( "EnKf D", 121),
                                           rep( "EnKf 100 cm", 121),
                                           rep( "EnKf 70 cm", 121),
                                           rep( "EnKf 60 cm", 121),
                                           rep( "Simulated laminar", 121)
                                           #rep( "Simulated turbulent", 121)
                                           ))
ggplot(experimentalFilter, aes(x, values, col= Data, linetype= Data))+geom_line(size=0.5)+
  #geom_point(data= exp2, aes(x = V1, y=V2), color= "deeppink")+
  scale_color_manual(values= c("deeppink", "green", "violet",  "turquoise", "coral", "blue", "navy", "darkred")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid","solid", "dotdash", "dashed","dashed"))+
  #ggtitle("Experimental and EnKF filtered velocity") + 
  xlab("t (s)") + ylab("Velocity (m/s)")

