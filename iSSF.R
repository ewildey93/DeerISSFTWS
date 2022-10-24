library(sp)
library(sf)
library(raster)
library(rgdal)
library(adehabitatHR)
library(adehabitatLT)
library(dplyr)
library(data.table)
library(Hmisc)
library(ggplot2)
library(tidyverse)
library(data.table)
library(amt)

filenames <- list.files(path = "C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS/Data",pattern = ".csv",full.names = TRUE)
objectnames <- list.files(path = "C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS/Data",pattern = ".csv",full.names = FALSE)
for(i in 1:length(filenames)) assign(objectnames[i], read.csv(filenames[i]))
Deer <- do.call("list", mget(objectnames))
#Deer <- read.csv("C:/Users/eliwi/OneDrive/Documents/Salida/DeerM4940.csv")
Deer <- lapply(Deer, select, c(4,5,6,7,11,15,16,20))
#colnames(Deer)[4] <- 'DateTime'
colnames <- c('DateTime','Lat','Long','hdop','alt','nsats','vdop','temp')
Deer <- lapply(Deer, setNames, colnames) 
#Deer$DateTime <- as.POSIXct(Deer$DateTime,format= '%Y-%m-%d %H:%M:%S', tz= 'GMT')
listnames <- grep(pattern = ".csv",x = ls(), value = TRUE)
listnames <- gsub(pattern = ".csv", replacement = '',x = listnames)
#listnames <- gsub(pattern = "_3.31", replacement = '',x = listnames)
Deer <- lapply(Deer, function (x) {x$DateTime <- as.POSIXct(x$DateTime,format= '%Y-%m-%d %H:%M:%S', tz="GMT")
                                  x <- x[x$Lat > 0,]
                                  x <- x[x$hdop < 10 & x$vdop < 10 & x$nsats >2,];
                                  x})
Deer <- lapply(Deer, function (x) {attributes(x$DateTime)$tzone <- "America/Denver"
                                    ;x})
names(Deer) <- listnames
Deer <- Map(cbind, Deer, ID = names(Deer))
######################################################################
#str(Deer)
#attributes(Deer$DateTime)$tzone <- "America/Denver"
#attributes(Deer$DateTime)$tzone
#Deer <- subset(Deer, lat > 0)
######################################################################
Deer2 <- rbindlist(Deer)
str(Deer2)
Deer2 <- Deer2[c(-1253:-1595,-1762),]

sf_deer<-st_as_sf(Deer2, coords=c("Long", "Lat"), crs=CRS("+init=epsg:4326"))
sf_deer<-st_transform(sf_deer, CRS("+init=epsg:32613"))
Deer2 <- cbind(Deer2, st_coordinates(sf_deer))
write.csv(Deer2, "C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS/Deer2.csv")
st_coordinates(sf_deer)
Deer2 <- read.csv("C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS/Deer2.csv")
##########################################################################
#sfDeer <- Map(st_as_sf, Deer, coords=c("Long", "Lat"), crs=CRS("+init=epsg:4326"))
###############################################################################
#Converting dates to POSIXct
#sf_cari$Date2<-as.POSIXct(strptime(sf_cari$Date, format="%m/%d/%Y %H:%M"))

#Fit trajectory
traj<-as.ltraj(st_coordinates(sf_deer), date=sf_deer$DateTime, id=sf_deer$ID)
x <- st_coordinates(sf_deer)
# Create burst for large missing gaps 
traj2<-cutltraj(traj, "dt>24*3600")

#### sett0
t <- c('07:03','09:23', '18:20', "03:42",'15:00','18:26','11:21','13:17','15:14',"23:03", "19:08", "02:55","02:44","15:35",'18:35','20:42','03:35',
       '20:52','18:32','19:26','21:15','20:01','21:20','14:56','02:47','07:54','12:02','03:46','14:56','03:51','21:13','06:51','21:24','10:26','02:21','08:02','15:27')
refda<-strptime(t, "%H:%M",tz = 'America/Denver')
dftraj2 <- ld(traj2)

z <- as.data.frame(table(dftraj2$dt))
is.regular(traj2)
traj3 <- redisltraj(traj2, "14400", type="time")
df_traj3 <- ld(traj3)
write.csv(df_traj3, "C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS/dftraj3.csv")
is.regular(traj3)
f <- ld(traj3)
g <- f[f$dt> 10800,]
h <- f[f$dt< 900,]
table(f$dt)
#### setNA
traj4<-setNA(traj3, refda, dt=14400, tol=2880, units="sec")
traj4<-sett0(traj4, refda, dt=14400, tol=2880, units="sec")



traj_df<-ld(traj4) ##Convert trajectory to data.frame 
traj_df[traj_df$dist > 1500,]
#traj_df <- traj_df[traj_df$dist < 20000,]


write.csv(traj_df, "C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS/traj_df.csv")
######## Fit step length and turning angle distribution 
#install.packages("MASS")
library(MASS)
?fitdistr

tt<-rnorm(10000, mean=0, sd=1)
hist(tt)
fitdistr(tt, "normal")

##Start with step lenght
dist<-fitdistr(na.omit(traj_df$dist[traj_df$dist>0]), "gamma")
dist

step<-rgamma(10000,  dist$estimate[1], dist$estimate[2])
#z <- rgamma(10000, .001, .001)
#hist(z, xlim=c(0,100), breaks=100)
#plot(dgamma(x,shape=.001,rate=.001),main="Gamma",type="l", ylim=c(0,0.02))
#x<-seq(from=0,to=6,length.out=100)
par(mfrow=c(1,2))
hist(traj_df$dist, breaks=20, main="Observed")
hist(step, breaks=20, main="Step Distance")
##########################################
par(mfrow=c(1,2))
hist(traj_df$dist, breaks=20, main="Observed Step Distance", xlab= "Distance (m)")
hist(traj_df$rel.angle, breaks=20, main="Observed Relative Step Angle", xlab= "Relative Angle (rad)")

##### Turning angle distribution
#install.packages('circular')
library(circular)

VMdist <- fit_distr(traj_df$rel.angle, dist_name = "vonmises", na.rm = TRUE)

angle<-rvonmises(n = 10000,mu = VMdist[["params"]][["mu"]],kappa = VMdist[["params"]][["kappa"]] )
par(mfrow=c(1,2))
hist(traj_df$rel.angle, breaks=20)
par(mfrow=c(1,1))
hist(angle2$x, breaks=20)
angle2 <- as.data.frame(angle)
rand.dist<-data.frame(dist=step, rel.angle=angle, id="Theory")
head(rand.dist)
###### Generate random steps
library(devtools)
install_github("basille/hab")
library(hab)
?rdSteps

ssf<-rdSteps(traj3, nrs=30, rand.dis=rand.dist)


### Calculate end location
ssf$x_end<-ssf$x+ssf$dx
ssf$y_end<-ssf$y+ssf$dy
write.csv(ssf, "C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/ssfobj.csv")
ssf <- read.csv("C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/ssfobj.csv")
##################
res2 <- rcorr(as.matrix(my_data))
res2
library(geosphere)
?dist2Line
ssf_sf <- st_as_sf(x = ssf,coords = c(16,17),crs=CRS("+init=epsg:32613"))
ssf_sf2 <- st_transform(ssf_sf, crs=CRS('+init=epsg:4326'))
x <- st_coordinates(ssf_sf)
ssf_sp<-SpatialPoints(ssf[,c("x_end", "y_end")], CRS("+init=epsg:32613"))
ssf_sp2 <- spTransform(ssf_sp, CRS("+init=epsg:4326"))
Trails <- readOGR(dsn='C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers', layer='DissolvedTrails')
TrailsSF <- st_as_sf(Trails)
Trails2 <- spTransform(Trails,CRS("+init=epsg:4326"))

#get distance to trail for all points
Dist2Trail <- st_distance(ssf_sf, TrailsSF)
Dist2TrailDF <- as.data.frame(Dist2Trail)
Dist2TrailDF$num <- as.numeric(Dist2TrailDF$Dist2Trail)
ssf <- cbind(ssf, Dist2TrailDF)
write.csv(ssf, 'C:/Users/eliwi/OneDrive/Documents/Salida/ssf2.csv')

###################################################
ssf$hr <- format(as.POSIXct(ssf$date), format = '%H')
head(ssf)
ssf$DOW <- weekdays(as.Date(ssf$date))
table(ssf$day)
ssf$DOW<- as.factor(ssf$DOW)
ssf$DOW <- recode_factor(ssf$DOW,Monday = 'weekday',Tuesday='weekday', Wednesday='weekday', Thursday='weekday', Friday='weekday', Saturday='weekend', Sunday='weekend')
ssf$weekday<-ifelse(ssf$weekday %in% 'weekday', 1, 0) #This is call dummy code 
ssf$weekend<-ifelse(ssf$weekday %in% 'weekend', 1, 0)



traj_df$hr <- format(as.POSIXct(traj_df$date), format = '%H')
DistbyHr <- traj_df%>%na.omit(dist)%>%group_by(hr)%>%summarise("MedianDist"=median(dist), "SumDist"=sum(dist), sd=sd(dist))
traj_df$day <- weekdays(as.Date(traj_df$date))
table(traj_df$day)
traj_df$day<- as.factor(traj_df$day)
traj_df$day <- recode_factor(traj_df$day,Monday = 'weekday',Tuesday='weekday', Wednesday='weekday', Thursday='weekday', Friday='weekday', Saturday='weekend', Sunday='weekend')

ggplot(data=DistbyHr, aes(x=hr, y=MedianDist, group=1)) +
  geom_line(lwd=2, color="green" ) +
  labs(x = "Time of Day (hr)", y = "Median Distance Moved (m)",title = 'Mule Deer Activity over Diel Cycle') +
  scale_x_discrete(breaks=c('00','03','06','09','12','15','18','21'),labels= c('0','3','6','9','12','15','18','21')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
Weekend <- traj_df%>%na.omit(dist)%>%group_by(day,hr)%>%summarise("MedianDist"=median(dist), "SumDist"=sum(dist), sd=sd(dist))
ggplot(data=Weekend, aes(x=hr, y=MedianDist, group=day, color=day)) +
  geom_line(lwd=2) +
  labs(x = "Time of Day (hr)", y = "Median Distance Moved (m)",title = "Mule Deer Activity over Diel Cycle") +
  scale_x_discrete(breaks=c('00','03','06','09','12','15','18','21'),labels= c('0','3','6','9','12','15','18','21')) +
  theme_bw() +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  theme(plot.title = element_text(hjust = 0.5))
##################################################
lc <- raster("C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers/NLCD2019Clip4.6.tif")
projection(lc)
ssf_sf<-st_transform(ssf_sf, projection(lc))
ssf$lc<-extract(lc, ssf_sf)
table(ssf$lc)
ssf$developed<-ifelse(ssf$lc %in% c(21,22,23,24), 1, 0) #This is call dummy code 
ssf$forest<-ifelse(ssf$lc %in% c(41,42), 1, 0)
ssf$shrub<-ifelse(ssf$lc %in% c(52), 1, 0)
ssf$herb <- ifelse(ssf$lc %in% c(71,81,82), 1, 0)
ssf$wetlands <- ifelse(ssf$lc %in% c(90,95), 1, 0)
table(is.na(ssf$lc))
a <- ssf[is.na(ssf$lc),]
write.csv(a, "C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/checkLC.csv")
###########################################
TRI <- raster("C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers/MergeTRI.tif")
projection(TRI)
ssf_sf<-st_transform(ssf_sf, projection(TRI))
ssf$TRI<-extract(TRI, ssf_sf)
table(is.na(ssf$TRI))
b <- ssf[is.na(ssf$TRI),]
write.csv(b, "C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/checkTRI.csv")
###########################################
Aspect <- raster("C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers/Aspect4.6.tif")
projection(Aspect)
ssf_sf<-st_transform(ssf_sf, projection(Aspect))
ssf$Aspect<-extract(Aspect, ssf_sf)
table(is.na(ssf$Aspect))
c <- ssf[is.na(ssf$Aspect),]
write.csv(c, "C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/checkLC.csv")
################################################
library(amt)
x <- Deer2$X
y <- Deer2$Y
ssf$date <- as.POSIXct(ssf$date, format= "%Y-%m-%d %H:%M:%S", tz = "")
t <- Deer2$DateTime
tagvec <- Deer2$ID
trackDeer <- track(x = x, y = y, t=t, TagId=tagvec, crs=CRS("+init=EPSG:32613"))
#time_of_day function to get day, night, dawn, dusk
trackDeer <- time_of_day(x = trackDeer,solar.dep = 6,include.crepuscule = TRUE)
#combine dawn and dusk into crepuscular
library(dplyr)
trackDeer$tod_ <- as.factor(trackDeer$tod_)
trackDeer$tod_ <- recode_factor(trackDeer$tod_,dusk = 'crepuscular',dawn='crepuscular')

table(trackDeer$tod_)

ssf$TOD <- trackDeer$tod_
ssf$day<-ifelse(ssf$TOD %in% "day", 1, 0) #This is call dummy code 
ssf$crepuscular<-ifelse(ssf$TOD %in% "crepuscular", 1, 0)
ssf$night<-ifelse(ssf$TOD %in% "night", 1, 0)
###################################################
ssf <- ssf%>%
        mutate(forest_day = as.numeric(forest == 1 & TOD == "day"),
               forest_night = as.numeric(forest == 1 & TOD == "night"),
               forest_cre = as.numeric(forest == 1 & TOD == "crepuscular"),
               shrub_day = as.numeric(shrub == 1 & TOD == "day"),
               shrub_night = as.numeric(shrub == 1 & TOD == "night"),
               shrub_cre = as.numeric(shrub == 1 & TOD == "crepuscular"),
               herb_day = as.numeric(herb == 1 & TOD == "day"),
              herb_night = as.numeric(herb == 1 & TOD == "night"),
              herb_cre = as.numeric(herb == 1 & TOD == "crepuscular"))
ssf <- ssf%>%
  mutate(TRI_day = as.numeric(forest == 1 & TOD == "day"),
         TRI_night = as.numeric(forest == 1 & TOD == "night"),
         TRI_cre = as.numeric(forest == 1 & TOD == "crepuscular"),
         shrub_day = as.numeric(shrub == 1 & TOD == "day"),
         shrub_night = as.numeric(shrub == 1 & TOD == "night"),
         shrub_cre = as.numeric(shrub == 1 & TOD == "crepuscular"),
         herb_day = as.numeric(herb == 1 & TOD == "day"),
         herb_night = as.numeric(herb == 1 & TOD == "night"),
         herb_cre = as.numeric(herb == 1 & TOD == "crepuscular"))
write.csv(ssf, "C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/3.31/DeerSSF2")

###################################################
ssf$Dist2Trail <- Dist2Trail[,1]
#make speed column in meters/hr
ssf$Speed <- ssf$dist/3
write.csv(ssf, "C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/DeerSSF.csv")
hist(log(ssf$Dist2Trail))
#run ssf
install.packages("survival")
library(survival)
?clogit
str(ssf)
ssf <- read.csv("C:/Users/eliwi/OneDrive/Documents/Salida/Deer Data/3.31/DeerSSF2")
#control#
ssf1<-fit_issf(ssf,
              #response            
              case~
               #movement
               dist + log(dist)+ rel.angle +
              #random effect
                frailty(id) +
               #stratum
               strata(strata), model= TRUE)
ssf1
#habitat model
ssf2 <- fit_issf(ssf,
              #response            
              case~
              #habitat
              developed+forest+shrub+herb+wetlands+scale(TRI)+
              #movement
              dist + log(dist)+ rel.angle +
              #random effect
              frailty(id) +
              #stratum
              strata(strata), model=TRUE)
ssf2
#human
ssf3 <- fit_issf(ssf,
                 #response            
                 case~
                   #human
                   scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), model=TRUE)
ssf3
#global
ssf4 <- fit_issf(ssf,
                 #response            
                 case~
                   #habitat
                   developed+forest+shrub+herb+wetlands+scale(TRI)+
                   #human
                   #human
                   scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), model=TRUE)
ssf4

#hab + TOD
ssf5 <- fit_issf(ssf,
                 #response            
                 case~
                   #habitat + TOD
                   forest_day+forest_night+forest_cre+shrub_day+shrub_night+shrub_cre+herb_day+herb_night+herb_cre+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), model=TRUE)
#significant from global model
ssf6 <- fit_issf(ssf,
  #response            
  case~
    #human
    scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+ scale(TRI) +
    #movement
    dist + log(dist)+ rel.angle +
    #random effect
    frailty(id) +
    #stratum
    strata(strata), model=TRUE)
#TOD* TRI interaction?
ssf7 <- fit_issf(ssf,
  #response            
  case~
    #human
    scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+ scale(TRI)*TOD+
    #movement
    dist + log(dist)+ rel.angle +
    #random effect
    frailty(id) +
    #stratum
    strata(strata), model=TRUE)

#########Model Selection###############
#control#
ssf1A<-clogit(
                #response            
               case~
                 #movement
                 dist + log(dist)+ rel.angle +
                 #random effect
                 frailty(id) +
                 #stratum
                 strata(strata), data=ssf)
ssf1
#habitat model
ssf2A <- clogit(
                 #response            
                 case~
                   #habitat
                   developed+forest+shrub+herb+wetlands+scale(TRI)+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), data=ssf)
ssf2
#human
ssf3A <- clogit(
                 #response            
                 case~
                   #human
                   scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), data=ssf)
ssf3
#global
ssf4A <- clogit(
                 #response            
                 case~
                   #habitat
                   developed+wetlands+shrub+herb+forest+scale(TRI)+
                   #human
                   #human
                   scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), data=ssf)
summary(ssf4)
ssf4B <- clogit(
  #response            
  case~
    #habitat
    developed+wetlands+shrub+herb+forest+scale(TRI)+
    #human
    #human
    scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+
    #movement
    dist + log(dist)+ rel.angle +
    #stratum
    strata(strata), data=ssf)
#hab + TOD
ssf5A <- clogit(
                 #response            
                 case~
                   #habitat + TOD
                   forest_day+forest_night+forest_cre+shrub_day+shrub_night+shrub_cre+herb_day+herb_night+herb_cre+
                   #movement
                   dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), data=ssf)
#significant from global model
ssf6A <- clogit(
                 #response            
                 case~
                   #human
                   scale(Dist2Trail)+scale(Dist2Trail)*scale(Speed)+ scale(TRI)+
                 #movement
                 dist + log(dist)+ rel.angle +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), data=ssf)
#TOD* TRI interaction?
ssf7A <- clogit(
                 #response            
                 case~
                   #human
                   scale(Dist2Trail)+
                   #movement
                   dist + log(dist)+ rel.angle + Dist2Trail*TOD +
                   #random effect
                   frailty(id) +
                   #stratum
                   strata(strata), data=ssf)

install.packages("AICcmodavg")
library(AICcmodavg)

modlist1 <- list(control=ssf1A,
  habitat=ssf2A, human=ssf3A, global=ssf4A, habTOD=ssf5A)
modlist2 <- list(control=ssf1A,
                 habitat=ssf2A, human=ssf3A, global=ssf4A, habTOD=ssf5A, globalsig=ssf6A,TODTRI=ssf7A, ssf4B=ssf4B)
aictab(modlist1)
aictab(modlist2)
summary(ssf4A)
###########Graph#################
install.packages('sjPlot')
library(sjPlot)
library(ggplot2)
p <- plot_model(ssf4B, rm.terms= c('rel.angle','log(dist)','Speed'))
class(p)
p + labs(title='SSF Coefficients',) +
  ylab("Odds Ratios") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- plot_model(ssf4B, type='std')
p2 + labs(title='Beta Coefficients',) +
  theme(plot.title = element_text(hjust = 0.5))
#################################################################
## Figure
# If we want to use this to make a figure, we can pass a sequence
# of values to x1. Remember, x2 must always be 1 row. Let's 
# visualize the RSS for forage vs. mean forage.
x1 <- data.frame(rel.angle = seq(-3.140270, 6.283153, length.out = 100), 
                 scale(Dist2Trail):scale(Speed) = 0, 
                 scale(TRI) = 50, log(dist) = log(50),
                 dist = 0, scale(Dist2Trail) = 200, scale(Speed)=0)

x2 <- data.frame(rel.angle = mean(values(ssf$rel.angle)), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

logRSS <- log_rss(m1, x1, x2, ci = "se", ci_level = 0.95)

# We have a plot method for 'log_rss' objects to make a very basic figure.
plot(logRSS)

# But if we want more control, we can use ggplot with the 'df' data.frame.
ggplot(logRSS$df, aes(x = forage_x1, y = exp(log_rss), 
                      ymin = exp(lwr), ymax = exp(upr))) +
  geom_ribbon(color = "black", fill = "gray80", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  xlab("Forage at x1") +
  ylab("RSS vs. Mean Forage") +
  theme_bw()
