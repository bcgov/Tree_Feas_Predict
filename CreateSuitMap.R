###Kiri Daust, Nov 2017###
###This code imports a lat/long grid and a tree suitibility grid###
###and creates a map of suitibility by BGC###

.libPaths("E:/R packages")
rm(list=ls())

library(scales)
library(MASS)   
library(stats)
library(rgl)
library(RColorBrewer)
library(FNN)
library(igraph)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(colorRamps)
library(rgeos)
library(rgdal)
library(foreign)
library(randomForest)
require(tcltk)
library(bindr)
library(lattice)
library(plyr)
library(ggplot2)
library(rasterVis)

wd=tk_choose.dir()
setwd(wd)

Original <- read.csv("TreeSppSuit_v10.csv", stringsAsFactors = TRUE)
SuitTable <- read.csv("PredSuit_AllMods_2025.csv", stringsAsFactors = TRUE) ##Read in suitability table from model
SuitTable <- SuitTable[,c(2,3,5)]
BEC <- read.csv("BGCv10_2k_Clipped.csv", stringsAsFactors = TRUE) ##Read in BGC grid
BEC <- BEC[,c(3,4,2)]
BEC$BGC <- gsub(" ","", BEC$ID2) ##remove spaces from BGC labels
BEC <- BEC[,-3]
Original <- Original[, c(1,3,4)]

colnames(SuitTable) <- c("Site.Series", "Spp", "FutureSuit")
SuitTable <- SuitTable[grep("01", SuitTable$Site.Series),] ###Use zonal SS
SuitTable$Site.Series <- gsub(".{3}$", "", SuitTable$Site.Series) ##remove last three characters
SuitTable$Site.Series <- gsub('/','', SuitTable$Site.Series) ### if 101

treeList <- unique(SuitTable$Spp) ##list of tree species 
###SuitTable now has all species but only zonal Site Series#####

###Read shape file outline of BC
CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
BC_out <- readShapePoly("ProvincialOutline.shp", proj4string = CRS.albers)

###loop to create maps for each species in treeList###
for(j in 1:length(treeList)){
SuitShort <- subset(SuitTable, SuitTable$Spp == treeList[j]) ##Subset table for one tree species
Original$Spp <- as.character(Original$Spp)
OrigSuit <- subset(Original, Original$Spp == treeList[j])
SuitShort <- SuitShort[,-2]
OrigSuit <- OrigSuit[,-c(2)]

##Use below line if using min suitability in each subzone instead of 01 site series###
##SuitShort <- SuitShort[SuitShort$FSuit == ave(SuitShort$FSuit, SuitShort$BGC, FUN = min),] ##convert suitibility to minimum for each subzone

###link spatial BGC to suitibility (slow but works)
##BEC <- subset(BEC, BEC$BGC %in% OrigSuit$BGC == TRUE)
for (i in 1:length(BEC$BGC)){
  index <- which(SuitShort$Site.Series == BEC$BGC[i])[1]
  if(is.na(index)){
    BEC$Suit[i] <- 0
  }else{
    BEC$Suit[i] <- SuitShort$FutureSuit[index]
  }
  index <- which(OrigSuit$BGC == BEC$BGC[i])[1]
  if(is.na(index)){
    BEC$OrigSuit[i] <- 0
  }else{
    BEC$OrigSuit[i] <- 0.5
  }
}

BEC$Suit <- ifelse(BEC$Suit == 13, 3.5,
                   ifelse(BEC$Suit ==12, 2.5,
                          ifelse(BEC$Suit == 11, 1.5, BEC$Suit)))
##Make two different sets
TestOrig <- BEC[,-c(3:4)]
TestOrig <- subset(TestOrig, TestOrig$OrigSuit != 0) ##remove areas where species doesn't exist

TestMod <- BEC[,-c(3,5)]
TestMod <- subset(TestMod, TestMod$Suit != 0)

if(nrow(TestMod) != 0){ #just incase a species gets in that doesn't have any data

###Project from Lat Long to BC Albers
coordinates(TestMod) <- c("long","lat")
proj4string(TestMod) <- CRS("+init=epsg:4326")
TestMod <- spTransform(TestMod, CRS.albers)  # standard albers projection for BC gov't data
##Suitibility needs to be a factor
TestMod$Suit <- as.factor(TestMod$Suit)

##for original suitability
coordinates(TestOrig) <- c("long","lat")
proj4string(TestOrig) <- CRS("+init=epsg:4326")
TestOrig <- spTransform(TestOrig, CRS.albers)  # standard albers projection for BC gov't data
TestOrig$Suit <- as.factor(TestOrig$OrigSuit)

##Plot as PDF with darker colours where previously didn't exist
name = paste(treeList[j], "2025", "diff", ".pdf", sep = "")
pdf(name)
plot(BC_out)
plot(TestOrig, col = "grey", pch = 15, cex = 0.1, add = TRUE)
palette(c("Green","Green4", "Yellow", "gold", "orangered", "orangered3"))
plot(TestMod, col = TestMod$Suit, pch = 15, cex = 0.1, add = TRUE)
dev.off()

TestMod$Suit <- ifelse(TestMod$Suit == 3.5, 3,
                   ifelse(TestMod$Suit == 2.5, 2,
                          ifelse(TestMod$Suit == 1.5, 1, TestMod$Suit)))

##Plot as png with three standard colours
name2 = paste(treeList[j],"2025",".png", sep = "")
png(name2, width = 8, height = 8, units="in", res = 1000)
plot(BC_out)
plot(TestOrig, col = "grey", pch = 15, cex = 0.1, add = TRUE)
palette(c("Green","Yellow", "orangered"))
plot(TestMod, col = TestMod$Suit, pch = 15, cex = 0.1, add = TRUE)
dev.off()
}
}

##Old code
ptime <- system.time({
BEC <- merge.data.frame(BEC, SuitShort, by.x = "BGC", by.y = "Site.Series", all = TRUE)
})
BEC <- BEC[complete.cases(BEC),]

BC_cities <- readShapePoints("MajorCitiest.shp",proj4string = CRS.albers)

combineSuit <- function(x,y,z){
  if(is.na(index)){
    result <- 0
  }else{
    result <- z[index]
  }
  return(result)
}

BEC$Suit2 <- combineSuit(BEC$BGC, SuitShort$BGC, SuitShort$FSuit)
BEC$Suit2 <- lapply(BEC$BGC, FUN = combineSuit(BEC$BGC, SuitShort$BGC, SuitShort$FSuit))
