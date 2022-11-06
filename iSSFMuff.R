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
library(zoo)
library(INLA)
library(ggregplot)
library(magrittr)
setwd("C:/Users/eliwi/OneDrive/Documents/R/DeerISSFTWS")

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
Deer[["F49862"]] <- Deer[["F49862"]][!(Deer[["F49862"]]$DateTime %in% F49862out$DateTime),]
Deer <- Deer[-4]
saveRDS(Deer, file="./DeerList.rds")
readRDS("./DeerList.rds")
######################################################################
#str(Deer)
#attributes(Deer$DateTime)$tzone <- "America/Denver"
#attributes(Deer$DateTime)$tzone
#Deer <- subset(Deer, lat > 0)
######################################################################
Deer2 <- rbindlist(Deer)
str(Deer2)
#Deer2 <- Deer2[c(-1253:-1595,-1762),]


sf_deer<-st_as_sf(Deer2, coords=c("Long", "Lat"), crs=CRS("+init=epsg:4326"))
sf_deer<-st_transform(sf_deer, CRS("+init=epsg:32613"))
sf_deer <- cbind(sf_deer, st_coordinates(sf_deer))
Deer2 <- cbind(Deer2, st_coordinates(sf_deer))
Deer2 <- split(x = Deer2, f=Deer2$ID)
saveRDS(Deer2, file="./DeerList2.rds")
readRDS("./DeerList2.rds")
#uneven no. of observations per animal not a problem in Muff et al. 2020 example
#################################################################################
# We can format our data as a 'track_*' object using 'make_track()'.
trkList <- lapply(Deer2, function (x) make_track(x, X, Y, DateTime, id = ID, crs = CRS("+init=EPSG:32613")))
#trkA <- trk %>% group_by(id)%>%nest()

#don't need to break into behavioral states using HMM because of coarse interval
# We can do that with the function 'steps()'.
#trk1 <- trk %>% nest(-'id')
trk1 <- lapply(trkList, function (x) track_resample(x,rate=hours(4), tolerance=minutes(36)))
          #tolerance recommended to be less than 20% of interval https://github.com/ewildey93/Statistical-Methods-Seminar-Series/blob/main/avgar-smith_issa/Q_and_A.md
trk2 <- lapply(trk1, function (x)filter_min_n_burst(x,min_n=3))         
steps <- lapply(trk2, function (x) steps_by_burst(x))
saveRDS(steps, "./DeerSteps.rds")
steps <- readRDS('./DeerSteps.rds')

max <- max(sapply(steps, function(x) as.character(max(x$t2_, na.rm=TRUE))))
min <- min(sapply(steps, function(x) as.character(min(x$t2_, na.rm=TRUE))))
min <- as.POSIXct(min,format= '%Y-%m-%d %H:%M:%S', tz="America/Denver")
max <- as.POSIXct(max,format= '%Y-%m-%d %H:%M:%S', tz="America/Denver")
min <- min-5*60*60
max <- max
min <- as.POSIXct(format(round(min, units="hours"), format="%Y-%m-%d %H:%M"),tz="America/Denver")
max <- as.POSIXct(format(round(max, units="hours"), format="%Y-%m-%d %H:%M"),tz="America/Denver")
min <- "2021-04-05 17:00:00 MDT"
max <- "2022-10-09 11:00:00 MDT"

#####################################################################################
#####################################################################################
# Step 1. Sample ----
RndSteps <- lapply(steps, function (x) random_steps(x, n_control = 20))

#distributions fir using RndSteps4 df so theres one distribution instead of 10
tadist <- fit_distr(StepsDF$ta_ , "vonmises")
gamma <-  fit_distr(StepsDF$sl_, "gamma")
table(RndSteps4$sl_ == 0)
StepsDF <- bind_rows(steps, .id="column_label")
# Since we used the gamma distribution as our tentative step-length
# distribution, we need to include step length and log(step length)
# in our iSSF. We already have a column, 'sl_', but we need to add 
# 'log_sl_'. Likewise, we need the cosine of the turn angle to update
# the von Mises distribution, so we need to add cos_ta_:
RndSteps <- lapply(RndSteps, function (x) {
  mutate(x,
         log_sl_ = log(x$sl_),
         cos_ta_ = cos(x$ta_))})
saveRDS(RndSteps, "./DeerRndSteps.rds")
write.csv(rbindlist(RndSteps), "./RndStepsdf.csv")
RndSteps <- readRDS("./DeerRndSteps.rds")
#########################################################################
# Step 2. Attribute ---                                                 #
#########################################################################
#TOD
RndSteps <- lapply(RndSteps, function (x)time_of_day(x, where = "end",include.crepuscule=TRUE))

#Rnd Steps to sf object:
RndStepsSF <- lapply(RndSteps, function (x) st_as_sf(x, coords= c("x2_","y2_"), crs=CRS("+init=epsg:32613")))
saveRDS(RndStepsSF, "./RndStepsSF.rds")
#Trails:
Trails <- readOGR(dsn='C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers', layer='DissolvedTrails')
TrailsSF <- st_as_sf(Trails)


#get distance to trail for all points
Dist2Trail <- lapply(RndStepsSF, y=TrailsSF, function (x, y)st_distance(x, y))
Dist2TrailDF <- lapply(Dist2Trail, as.numeric)
Dist2TrailDF <- lapply(Dist2Trail, function (x) as.data.frame(x, col.names="Dist2Trail"))
#D2TList <- lapply(Dist2TrailDF, setNames, "Dist2Trail") 
#D2TList <- rbindlist(D2TList,idcol = TRUE)

#RndStepsA <- Map(bind_cols, RndSteps, ID = names(Deer))
#use tidyverse functions isntead of base R to rpeserve as class "steps"
RndSteps2 <- mapply(bind_cols, RndSteps, Dist2TrailDF, SIMPLIFY=FALSE)
#RndSteps2B <- lapply(RndSteps, function (x) left_join(x,D2TList, by=c(''=".id"), all=T))
saveRDS(RndSteps2, "./RndSteps2.rds")
RndSteps2 <- readRDS("./RndSteps2.rds")
max <- max(sapply(RndSteps2, function(x) as.character(max(x$x, na.rm=TRUE))))
min <- min(sapply(RndSteps2, function(x) as.character(min(x$x, na.rm=TRUE))))
#Traffic volume-daily
LRTrail <- read.csv("./CoVs/Little Rainbow Trail RAW.csv", header=FALSE)
SpartanTH <- read.csv("./CoVs/Spartan T.H. Parking Lot RAW.csv", header=FALSE)
colnames(LRTrail)[1:2] <- c("DateTime", "V1")
colnames(SpartanTH)[1:2] <- c("DateTime", "V2")
TrailActList <- list(LRTrail, SpartanTH)
TrailActList <- lapply(TrailActList, function (x) {x$DateTime <- as.POSIXct(x$DateTime,format= '%Y-%m-%d %H:%M:%S', tz="America/Denver");x})
TrailActList <- lapply(TrailActList, function (x) filter(x,DateTime > min & DateTime < max))
TrailActdf <- merge(TrailActList[[1]], TrailActList[[2]], by="DateTime", all=TRUE)
TrailActdf$Total <- TrailActdf$V1 + TrailActdf$V2
TrailActdf$RA <- rollmean(TrailActdf$Total, 4, fill=NA, align="right")
RA <- TrailActdf[,c(1,4,5)]
colnames(RA)[1] <- "RoundDate"
#extract time covariates
RndSteps3 <- lapply(RndSteps2, function (x) mutate (x, RoundDate= as.POSIXct(format(round(t2_, units="hours"), format="%Y-%m-%d %H:%M:%S"),tz="America/Denver")))
RndSteps3 <- lapply(RndSteps3, function (x) left_join(x,RA, by="RoundDate", all.y=F, all.x=T))
saveRDS(RndSteps3, "./RndSteps3.rds")

#landcover
lc <- raster("./CoVs/NLCD_2019_Land_Cover_L48_20210604_JbsuwO6GkIW9V4xHbi6d.tiff")
projection(lc)
RndSteps4 <- bind_rows(RndSteps3, .id="column_label")
RndSteps4SF<-st_as_sf(RndSteps4, coords=c("x2_", "y2_"), crs=CRS("+init=epsg:32613"))
RndSteps4SF<-st_transform(RndSteps4SF, projection(lc))
RndSteps4$lc<-raster::extract(lc, RndSteps4SF)
table(RndSteps4$lc)
RndSteps4$developed<-ifelse(RndSteps4$lc %in% c(21,22,23,24), 1, 0) #This is call dummy code 
RndSteps4$forest<-ifelse(RndSteps4$lc %in% c(41,42), 1, 0)
RndSteps4$shrub<-ifelse(RndSteps4$lc %in% c(52), 1, 0)
RndSteps4$herb <- ifelse(RndSteps4$lc %in% c(71,81,82), 1, 0)
RndSteps4$wetlands <- ifelse(RndSteps4$lc %in% c(90,95), 1, 0)
table(is.na(RndSteps4$lc))
saveRDS(RndSteps4, "./RndSteps4.rds")

#TOD recode
RndSteps4$tod_end_ <- recode_factor(RndSteps4$tod_end_,dusk = 'crepuscular',dawn='crepuscular')
saveRDS(RndSteps4, "./RndSteps4.rds")

#TRI
TRI <- raster("C:/Users/eliwi/OneDrive/Documents/Salida/GeospatialLayers/TRI10.26.tif")
RndSteps4$TRI<-raster::extract(TRI, RndSteps4SF)
saveRDS(RndSteps4, "./RndSteps4.rds")
RndSteps4 <- readRDS("./RndSteps4.rds")

#make combo categories
RndSteps4 <- rename(RndSteps4, 'TOD'='tod_end_')
RndSteps4 <- RndSteps4%>%
  mutate(forest_day = as.numeric(forest == 1 & TOD == "day"),
         forest_night = as.numeric(forest == 1 & TOD == "night"),
         forest_cre = as.numeric(forest == 1 & TOD == "crepuscular"),
         shrub_day = as.numeric(shrub == 1 & TOD == "day"),
         shrub_night = as.numeric(shrub == 1 & TOD == "night"),
         shrub_cre = as.numeric(shrub == 1 & TOD == "crepuscular"),
         herb_day = as.numeric(herb == 1 & TOD == "day"),
         herb_night = as.numeric(herb == 1 & TOD == "night"),
         herb_cre = as.numeric(herb == 1 & TOD == "crepuscular"))
RndSteps4 <- RndSteps4%>%
  mutate(         develop_day = as.numeric(developed == 1 & TOD == "day"),
                  develop_night = as.numeric(developed == 1 & TOD == "night"),
                  develop_cre = as.numeric(developed == 1 & TOD == "crepuscular"))

RndSteps4 <- RndSteps4%>%
  mutate(Stratam=as.factor(paste(column_label, burst_)))
RndSteps4 <- RndSteps4 %>%                               # Replacing values
  mutate(Stratum = as.integer(Stratam))
RndSteps4$Stratum <- as.factor(RndSteps4$Stratum)
class(RndSteps4)
saveRDS(RndSteps4, "./RndSteps4.rds")

#add variable for daily traffic
Daily <- read.csv("./CoVs/TRAFx Daily totals (2021-03-31 to 2022-08-14).csv")
Daily$Daily <- Daily$Little.Rainbow.Trail + Daily$Spartan.T.H..Parking.Lot
Daily2 <- Daily[,c(1,4)]
colnames(Daily2)[1] <- "Date"
Daily2$Date <- as.Date(Daily2$Date)
RndSteps4 <- mutate(RndSteps4, Date= as.Date(RndSteps4$t2_))
RndSteps4 <- left_join(x = RndSteps4, y=Daily2, by="Date", all.x=T)
RndSteps4 <- RndSteps4[-(9787:10122),]
RndSteps4 <- RndSteps4 %>% 
  mutate(sl_ = replace(sl_, sl_==0, 1))
saveRDS(RndSteps4, "./RndSteps4.rds")
RndSteps4 <- readRDS("./RndSteps4.rds")
hist(RndSteps4$log_sl_)
hist(RndSteps4$sl_)
range(RndSteps4$Total)
#correlation test
library(corrplot)
CoVs <- RndSteps4[,c(18,21,28)]
CoVs <- CoVs[-(9787:10122),]
corrplot(cor(CoVs),
         method = "number",
         type = "upper" # show only upper side
)
cor(CoVs$x,CoVs$RA)
cor(CoVs$RA,CoVs$TRI)

#remove data without traffic volume

#RndSteps5 <- RndSteps4[-(9787:10122),]
class(RndSteps5)
table(is.na(RndSteps5$RA))
range(RndSteps5$Total)

#scale continous variables
hist(RndSteps4$cos_ta_)
hist(RndSteps4$ta_)
hist(RndSteps5$log_sl_)
hist(RndSteps5$x)
hist(RndSteps5$RA)
RndSteps5 <- RndSteps4 %>% 
  mutate(sl_ = replace(sl_, sl_==0, 1))
RndSteps5 <- mutate(RndSteps5, log_sl_= log(sl_))
RndSteps5$sl_ <- scale(RndSteps5$sl_)
RndSteps5$RA <- scale(RndSteps5$RA)
RndSteps5$x <- scale(RndSteps5$x)
RndSteps5$log_sl_ <- scale(RndSteps5$log_sl_)
saveRDS(RndSteps5, "./RndSteps5.rds")
RndSteps5G <- readRDS("./RndSteps5.rds")
###########################################################################
####                           Modelling-iSSF         #####################
###########################################################################
#' Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

RndSteps5G$ANIMAL_ID <- as.numeric(as.factor(RndSteps5$column_label))
RndSteps5G$case <- as.integer(as.logical(RndSteps5$case_))
saveRDS(RndSteps5, "./RndSteps5.rds")
#To fit the model with random slopes in INLA, we need to generate new (but identical) variables of individual ID (ID cannot be used multiple times in the model formula):
RndSteps5G$ANIMAL_ID1 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID2 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID3 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID4 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID5 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID6 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID7 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID8 <- RndSteps5$ANIMAL_ID
RndSteps5G$ANIMAL_ID9 <- RndSteps5$ANIMAL_ID
saveRDS(RndSteps5, "./RndSteps5.rds")
#Control model
formula.control <- case ~ -1 + 
                    #movement kernel
                    sl_ + cos_ta_ + log_sl_ + log_sl_*cos_ta_ +
                    f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-4),fixed=T)))

r.inla.control <- inla(formula.control, family ="Poisson", data=RndSteps5G,
                     control.fixed = list(
                       mean = mean.beta,
                       prec = list(default = prec.beta)
                     ),control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE)
)

r.inla.control$summary.fixed
Efxplot(list(r.inla.control))
c(WAIC=r.inla.control$waic$waic, DIC=r.inla.control$dic$dic)

#library(glmmTMB)
#TMBStruc.control = glmmTMB(case ~ sl_ + cos_ta_ +
#                        (1|Stratum), 
#                      family=poisson, data=RndSteps5, doFit=FALSE) 
#TMBStruc.control$parameters$theta[1] = log(1e3) 
#TMBStruc.control$mapArg = list(theta=factor(c(NA)))
#glmm.TMB.control = glmmTMB:::fitTMB(TMBStruc.control) 
#summary(glmm.TMB.control)



#habitat model
formula.habitat <- case ~ -1 +
  developed+shrub+herb+wetlands+scale(TRI)+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_ + log_sl_*cos_ta_ +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID2,shrub,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID4,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))




r.inla.habitat <- inla(formula.habitat, family ="Poisson", data=RndSteps5,
                       control.fixed = list(
                         mean = mean.beta,
                         prec = list(default = prec.beta)
                       ),control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE)
)

r.inla.habitat$summary.fixed
c(WAIC=r.inla.habitat$waic$waic, DIC=r.inla.habitat$dic$dic)
Efxplot(list(r.inla.habitat))

#human model
formula.human <- case ~ -1 +
  x*RA+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_+ log_sl_*cos_ta_ +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,x,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


r.inla.human <- inla(formula.human, family ="Poisson", data=RndSteps5,
                       control.fixed = list(
                         mean = mean.beta,
                         prec = list(default = prec.beta)
                       ),control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE)
)

r.inla.human$summary.fixed
Efxplot(list(r.inla.human))
r.inla.human$summary.hyperpar
c(WAIC=r.inla.human$waic$waic, DIC=r.inla.human$dic$dic)
#global
formula.global <- case ~ -1 +
  #habitat
  developed+forest+herb+wetlands+scale(TRI)+
  #human
  x+
  #movement kernel    
  sl_ *RA+ cos_ta_ + log_sl_*RA +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,x,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID2,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,developed,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID4,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID5,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) 

r.inla.global <- inla(formula.global, family ="Poisson", data=RndSteps5,
                     control.fixed = list(
                       mean = mean.beta,
                       prec = list(default = prec.beta)
                     ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.global))

#global2
formula.global.2G <- case ~ -1 +
  #habitat
  developed+forest*RA+herb+wetlands+scale(TRI)+
  #human
  x+
  #movement kernel    
  sl_ * RA  + log_sl_*RA +cos_ta_ + 
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,x,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID2,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,developed,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID4,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID5,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) 


r.inla.global.2 <- inla(formula.global.2G, family ="Poisson", data=RndSteps5G,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)
                      ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)

Efxplot(list(r.inla.global.2))
r.inla.global.2$summary.fixed
ModelList <- list(control=r.inla.control, habitat=r.inla.habitat,human=r.inla.human, global=r.inla.global,global2=r.inla.global.2)
ModelSelTable <- INLA.model.sel(ModelList)
#################################################################
#####################Model Selection#############################
#################################################################
#see function at bottom
ModelList <- list(control=r.inla.control, habitat=r.inla.habitat,human=r.inla.human, global=r.inla.global,global2=r.inla.global.2)
ModelSelTable <- INLA.model.sel(ModelList)
saveRDS(ModelSelTable, "./ModelSelTable.rds")
write.csv(ModelSelTable, "./ModelSelTable.csv")
#
RndSteps6 <- RndSteps5
RndSteps6$Total <- scale(RndSteps6$Total)
RndSteps6$Daily <- scale(RndSteps6$Daily)


##########################################################
################iSSF by TOD###############################
##########################################################
Day <- filter(RndSteps6, TOD == "day")
Night <- filter(RndSteps6, TOD == "night")
Cre <- filter(RndSteps6, TOD == "crepuscular")
#control day
r.inla.control.day <- inla(formula.control, family ="Poisson", data=Day,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.control.day))
#habitat day
formula.habitat.day <- case ~ -1 +
  developed+shrub+herb+wetlands+scale(TRI)+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_ +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID2,shrub,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID4,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))




r.inla.habitat.day <- inla(formula.habitat, family ="Poisson", data=Day,
                       control.fixed = list(
                         mean = mean.beta,
                         prec = list(default = prec.beta)
                       ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)

r.inla.habitat.day$summary.fixed
c(WAIC=r.inla.habitat.day$waic$waic, DIC=r.inla.habitat.day$dic$dic)
Efxplot(list(r.inla.habitat.day))

#human day
Day2 <- Day[,c(2,4,5,7,30,37,38)]
seq <- rep(g$round, times=9)
#g=41 obs
x <- rep(g$`RndSteps5$x`, times=9)
newobs <- data.frame(sl_=0,log_sl_=0,cos_ta_=0,x=x,UnscaleDist=seq,Stratum=seq(424:792),case=NA,ANIMAL_ID1=rep(1:9, each=41))
newsobs <- newobs[order(newobs$ANIMAL_ID1,newobs$x),]

formula.human.day <- case ~ -1 +
  x+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_+RA*log_sl_
f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,x,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

r.inla.human.day <- inla(formula.human, family ="Poisson", data=Day,
                         control.fixed = list(
                           mean = mean.beta,
                           prec = list(default = prec.beta)
                         ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE),control.predictor = list(compute=TRUE),
)

c(WAIC=r.inla.human.day$waic$waic, DIC=r.inla.human.day$dic$dic)
r.inla.human.day$summary.fixed
Efxplot(list(r.inla.human.day))

#global
formula.global.day <- case ~ -1 +
  #habitat
  developed+shrub+herb+wetlands+scale(TRI)+
  #human
  x+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_ + log_sl_*RA +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,x,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID2,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,shrub,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID4,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID5,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

r.inla.global.day <- inla(formula.global, family ="Poisson", data=Day,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)
                      ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.global.day))
r.inla.global.day$summary.fixed

#global 2 day
r.inla.global.day2 <- inla(formula.global.2, family ="Poisson", data=Day,
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)
                          ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.global.day2))
#Model selection
DayList <- list(control.day=r.inla.control.day,habitat.day=r.inla.habitat.day,human.day=r.inla.human.day,global.day=r.inla.global.day, global.day.2=r.inla.global.day2)
DayTable <- INLA.model.sel(DayList)

##################################################################################

#habitat night
formula.habitat.night <- case ~ -1 +
  developed+shrub+herb+wetlands+scale(TRI)+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_ +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID2,shrub,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID4,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))




r.inla.habitat.night <- inla(formula.habitat, family ="Poisson", data=Night,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(dic = TRUE, waic = TRUE)
)

r.inla.habitat.night$summary.fixed
c(WAIC=r.inla.habitat.night$waic$waic, DIC=r.inla.habitat.night$dic$dic)
Efxplot(list(r.inla.habitat.night))

#control night
r.inla.control.night <- inla(formula.control, family ="Poisson", data=Night,
                             control.fixed = list(
                               mean = mean.beta,
                               prec = list(default = prec.beta)
                             ),control.compute = list(dic = TRUE, waic = TRUE)
)

#human night
r.inla.human.night <- inla(formula.human, family ="Poisson", data=Night,
                             control.fixed = list(
                               mean = mean.beta,
                               prec = list(default = prec.beta)
                             ),control.compute = list(dic = TRUE, waic = TRUE)
)

#global night
formula.global

r.inla.global.night <- inla(formula.global, family ="Poisson", data=Night,
                             control.fixed = list(
                               mean = mean.beta,
                               prec = list(default = prec.beta)
                             ),control.compute = list(dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.global.night))

#global 2 night
formula.global.2
r.inla.global.night2 <- inla(formula.global.2, family ="Poisson", data=Night,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(cpo=TRUE,dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.global.night2))
#Model Selection
NightList <- list(control.night=r.inla.control.night,habitat.night=r.inla.habitat.night,human.night=r.inla.human.night,global.night=r.inla.global.night, global.night2=r.inla.global.night2)
NightTable <- INLA.model.sel(NightList)


#################################################################################
#habitat cre
formula.habitat.cre <- case ~ -1 +
  developed+shrub+herb+wetlands+scale(TRI)+
  #movement kernel    
  sl_ + cos_ta_ + log_sl_ +
  f(Stratum, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,forest,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID2,shrub,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,herb,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ANIMAL_ID4,wetlands,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))




r.inla.habitat.cre <- inla(formula.habitat, family ="Poisson", data=Cre,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(dic = TRUE, waic = TRUE)
)

r.inla.habitat.cre$summary.fixed
c(WAIC=r.inla.habitat.cre$waic$waic, DIC=r.inla.habitat.cre$dic$dic)
Efxplot(list(r.inla.habitat.cre))

#control cre
r.inla.control.cre <- inla(formula.control, family ="Poisson", data=Cre,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(dic = TRUE, waic = TRUE)
)

#human.cre
r.inla.human.cre <- inla(formula.human, family ="Poisson", data=Cre,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(dic = TRUE, waic = TRUE)
)
Efxplot(list(r.inla.human.cre))
#global.cre
r.inla.global.cre <- inla(formula.global, family ="Poisson", data=Cre,
                           control.fixed = list(
                             mean = mean.beta,
                             prec = list(default = prec.beta)
                           ),control.compute = list(dic = TRUE, waic = TRUE)
)

Efxplot(list(r.inla.global.cre))

#global2.cre
r.inla.global2.cre <- inla(formula.global.2, family ="Poisson", data=Cre,
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)
                          ),control.compute = list(dic = TRUE, waic = TRUE)
)

Efxplot(list(r.inla.global2.cre))
#Model Selection
CreList <- list(control.cre=r.inla.control.cre,habitat.cre=r.inla.habitat.cre,human.cre=r.inla.human.cre,global.cre=r.inla.global.cre, global2=r.inla.global2.cre)
CreTable <- INLA.model.sel(CreList)

####################scrap############################
o <- read.csv("C:/Users/eliwi/Downloads/d_otter.csv")
table(o$Loc)
table(is.na(trk2df$ta_))
F49862out <- Deer2[c(1253:1595,1762),]
attr(F49862out$DateTime, "tzone") <- "America/Denver"
Out <- Deer[["F49862"]][(Deer[["F49862"]]$DateTime %in% F49862out$DateTime),]

z <- RndStepsSF[["F46538"]][1:40,]
st_distance(z, TrailsSF)
table(is.na(CoVs$RA))
which(is.na(CoVs$RA))
what <- RndSteps4[9787:10122,]
which(table(RndSteps4$step_id_) != 10 )
table(Day2$ANIMAL_ID1)
hist(Night$sl_)
hist(Night$log_sl_)
hist(Day2$cos_ta_)
hist(RndSteps4$x)
boxplot(Night$sl_)
boxplot(Night$log_sl_)
Outlier <- Night[Night$sl_ > 10,]
scalex <- cbind.data.frame(RndSteps4$x,RndSteps5$x)
scalex$round <- round(scalex$`RndSteps4$x`)
e <- scalex[scalex$round %in% seq,]
table(e$round)
range(RndSteps4$x)
f <- scalex[as.numeric(scalex$`RndSteps4$x`) < 2650,]
g <- e%>%distinct(round, .keep_all = T)
cover = factor("grassland",
               levels = c("grassland", "forest", "wetland"))
NightCase <- Night[Night$case== 1,]
DayCase <- Day[Day$case== 1,]
ggplot(DayCase, aes(x=sl_, y=log_sl_)) + geom_point()
tawtf <- RndSteps4[,c("cos_ta_","ta_")]
###################################################################
#functions#########################################################
###################################################################
INLA.model.sel <- function (x) {
  df <- data.frame(WAIC=numeric(),DIC=numeric(),CPO=numeric(),ML=numeric())
  
  a <- lapply(x, function (x) rbind(df,data.frame(WAIC=x$waic$waic, DIC=x$dic$dic, CPO=-1*sum(log(x$cpo$cpo)), ML=as.numeric(x$mlik[1,1]))))
  b <- rbindlist(a,idcol = T)
  b$DeltaWAIC <- b$WAIC-min(b$WAIC)
  b$DeltaDIC <- b$DIC-min(b$DIC)
  b <- b[order(b$DeltaWAIC,b$DeltaDIC),]
  return(b)
}
ModelList <- list(control=r.inla.control, habitat=r.inla.habitat,human=r.inla.human)
ModelSelTable <- INLA.model.sel(ModelList)
