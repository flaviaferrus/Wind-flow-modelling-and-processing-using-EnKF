### KF forecasting



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

X2022_03 <- read_csv("2022-03-14-15-27-49.csv")
AE <- X2022_03[ 120: 300, ]
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


## KF

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

auto.arima(PTu$velC06[59:178]) ## MA(1) with mean
## log likelihood=457.83
mod3 <- arima(PTu$velC05[59:178],  order = c(2,0,0), include.mean = TRUE)
mod3
m4.dlm <- dlmModARMA( ar= c(0.8807,-0.3596))
m4.dlm
model.filteredC3 <- dlmFilter(velC - mean(velC), m4.dlm)
C3mean <-  model.filteredC3$f + mean(velC)


length(velC)
length(velA)



#### Forecasting
length(model.filteredA3$f)
length(model.filteredA3$a)
length(model.filteredA3$m)

## hi ha molts de a

model.filteredA3$a

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

nrow(X2022_03_10_17_20_36)
## hi ha 15 :)( punts on podem fer forecast

vC <- as.numeric(X2022_03_10_17_20_36$`076 [m/s]`[60:254])
vC
length(vC)

vA <- as.numeric(X2022_03$`076 [m/s]`[120:312])
vA <- c(vA, c(1.03, 1.02))
length(vA)

KFplot2  <- data.frame(x = index(vA[58:195]), 
                       values = c(vA[58:195], 
                                  #model.filteredA3.wmean[59:178], 
                                  vC[58:195], 
                                  #model.filteredC3.wmean[59:178],
                                  #C3mean[59:178],
                                  #A3mean[59:178],
                                  model.filteredA3$a[59:196]+ mean(velA),
                                  model.filteredC3$a[59:196]+ mean(velC), 
                                  model.filteredA3$m[59:196]+ mean(velA),
                                  FilteredCWTF
                                  #mmA$mitjana[59:178], PTu$velC06[59:178]
                       ), 
                       Data = c(rep("True velocity B", 138), 
                                #rep( "Filtered state A MA(2) w intercept", 120),
                                rep( "True velocity D",138), 
                                #rep( "Filtered State C AR(2) w intercept", 120),
                                #rep( "Forecast D AR(2) w/mean", 120),
                                #rep( "Forecast B MA(2) w/mean", 120),
                                rep( "Predicted B AR(2) w/mean", 138),
                                rep( "Predicted D MA(2) w/mean", 138), 
                                rep( "Filtered B AR(2) w/mean", 138),
                                rep( "Filtered D MA(2) w/mean", 138)
                                #rep("Simulated velocity B", 120), 
                                #rep("Simulated velocity D", 120)
                       ) 
)
ggplot(KFplot2, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "deepskyblue", "coral","blue",  "red", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c("solid","solid","dotdash", "dotdash", "dashed",  "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  #geom_pointrange(KFplot, aes(ymin=values-sd, ymax=values+sd))+
  xlab("Time (s)") + ylab("Velocity (m/s)")

FilteredCWTF1 <- model.filteredC3$m[49:183]+ mean(velC)
FilteredCWTF2 <- model.filteredC3$m[185:197]+ mean(velC)
FilteredCWTF <- c(FilteredCWTF1, FilteredCWTF2 )
model.filteredC3$m[184]+ mean(velC)

length(FilteredCWTF)


KFplot2  <- data.frame(x = index(vA[48:195]), 
                       values = c(vA[48:195], 
                                  #model.filteredA3.wmean[59:178], 
                                  vC[48:195], 
                                  #model.filteredC3.wmean[59:178],
                                  #C3mean[59:178],
                                  #A3mean[59:178],
                                  model.filteredA3$a[49:196]+ mean(velA),
                                  model.filteredC3$a[49:196]+ mean(velC), 
                                  model.filteredA3$m[49:196]+ mean(velA),
                                  FilteredCWTF
                                  #mmA$mitjana[59:178], PTu$velC06[59:178]
                       ), 
                       Data = c(rep("True velocity B", 148), 
                                #rep( "Filtered state A MA(2) w intercept", 120),
                                rep( "True velocity D",148), 
                                #rep( "Filtered State C AR(2) w intercept", 120),
                                #rep( "Forecast D AR(2) w/mean", 120),
                                #rep( "Forecast B MA(2) w/mean", 120),
                                rep( "Predicted B AR(2)", 148),
                                rep( "Predicted D MA(2)", 148), 
                                rep( "Filtered B AR(2)", 148),
                                rep( "Filtered D MA(2)", 148)
                                #rep("Simulated velocity B", 120), 
                                #rep("Simulated velocity D", 120)
                       ) 
)
ggplot(KFplot2, aes(x, values, col= Data, linetype=Data))+geom_line(size=0.5)+
  scale_color_manual(values= c( "deepskyblue", "coral","blue",  "red", "navy", "darkred")) +
  #scale_linetype_manual(values=c("dotdash", "dashed",  "dotdash","dashed","solid","solid", "solid","solid"))+
  scale_linetype_manual(values=c("solid","solid","dotdash", "dotdash", "dashed",  "dashed"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  #ggtitle("Observed, simulated and filtered punctual velocity") + 
  #geom_pointrange(KFplot, aes(ymin=values-sd, ymax=values+sd))+
  xlab("Time (s)") + ylab("Velocity (m/s)")


errorA <- (model.filteredA3$a[178:196] + mean(velA) - vA[177:195])^2
mean(errorA)
errorC <- (model.filteredC3$a[178:196] + mean(velC) - vC[177:195])^2
mean(errorC)


errorA2 <- (model.filteredA3$a[48:196] + mean(velA) - vA[47:195])^2
mean(errorA2)
errorC2 <- (model.filteredC3$a[48:196] + mean(velC) - vC[47:195])^2
errorC2

