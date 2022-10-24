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

####################scrap############################
o <- read.csv("C:/Users/eliwi/Downloads/d_otter.csv")
table(o$NA_ANIMAL)
table(is.na(trk2df$ta_))
F49862out <- Deer2[c(1253:1595,1762),]
attr(F49862out$DateTime, "tzone") <- "America/Denver"
Out <- Deer[["F49862"]][(Deer[["F49862"]]$DateTime %in% F49862out$DateTime),]
