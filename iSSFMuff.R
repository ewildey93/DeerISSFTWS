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

max <- max(sapply(steps, function(x) as.character(max(x$t2_, na.rm=TRUE))))
min <- min(sapply(steps, function(x) as.character(min(x$t2_, na.rm=TRUE))))
min <- as.POSIXct(min,format= '%Y-%m-%d %H:%M:%S', tz="America/Denver")
max <- as.POSIXct(max,format= '%Y-%m-%d %H:%M:%S', tz="America/Denver")
min <- min-5*60*60
max <- max
min <- as.POSIXct(format(round(min, units="hours"), format="%Y-%m-%d %H:%M"),tz="America/Denver")
max <- as.POSIXct(format(round(max, units="hours"), format="%Y-%m-%d %H:%M"),tz="America/Denver")
#min <- "2021-04-05 17:00:00 MDT"
#max <- "2022-10-09 11:00:00 MDT"

#####################################################################################
#####################################################################################
# Step 1. Sample ----
RndSteps <- lapply(steps, function (x) random_steps(x, n_control = 20))

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

RndSteps2 <- mapply(cbind, RndSteps, Dist2TrailDF, SIMPLIFY=FALSE)
saveRDS(RndSteps2, "./RndSteps2.rds")

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
RndSteps4 <- rbindlist(RndSteps3, idcol=T)
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
####################scrap############################
o <- read.csv("C:/Users/eliwi/Downloads/d_otter.csv")
table(o$Loc)
table(is.na(trk2df$ta_))
F49862out <- Deer2[c(1253:1595,1762),]
attr(F49862out$DateTime, "tzone") <- "America/Denver"
Out <- Deer[["F49862"]][(Deer[["F49862"]]$DateTime %in% F49862out$DateTime),]

z <- RndStepsSF[["F46538"]][1:40,]
st_distance(z, TrailsSF)
