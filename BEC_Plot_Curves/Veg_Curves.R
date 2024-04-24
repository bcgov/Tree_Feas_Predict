#### Fit distribution curves of tree species to climate variable based on BEC plot data
#### Kiri Daust May 2020

library(data.table)
library(tidyverse)
library(sf)
library(raster)
library(ggstance) ##need this for boxplot
library(foreach)

setwd("D:/GitHub/TreeSuitabilityPrediction/")


aspectCalc1 <- function(lat, slope, aspect){##eqn 1 from McCune 2007
  lat <- lat*pi/180
  slope <- slope*pi/180
  aspect <- aspect*pi/180
  HL <- abs(pi - abs(aspect - (5*pi)/4))
  res <- -1.467+1.582*cos(lat)*cos(slope)-1.5*cos(HL)*sin(slope)*sin(lat)-0.262*sin(lat)*sin(slope)+0.607*sin(HL)*sin(slope)
  return(exp(res))
}

aspectCalc2 <- function(lat, slope, aspect){##eqn 2 from McCune 2007
  lat <- lat*pi/180
  slope <- slope*pi/180
  aspect <- aspect*pi/180
  HL <- abs(pi - abs(aspect - (5*pi)/4))
  res <- -1.236+1.35*cos(lat)*cos(slope)-1.376*cos(HL)*sin(slope)*sin(lat)-0.331*sin(lat)*sin(slope)+0.375*sin(HL)*sin(slope)
  return(exp(res))
}

aspectCalc3 <- function(lat, slope, aspect){##eqn 3 from McCune 2007
  lat <- lat*pi/180
  slope <- slope*pi/180
  aspect <- aspect*pi/180
  HL <- abs(pi - abs(aspect - (5*pi)/4))
  res <- 0.339+0.808*cos(lat)*cos(slope)-0.196*sin(lat)*sin(slope)-0.482*cos(HL)*sin(slope)
  return(res)
}

###change the name of aspectCalc in this function to use a different equation
aspectRatio <- function(lat,slope,aspect,neutSlope = 38, neutAsp = 45){
  new <- aspectCalc3(lat,slope,aspect)
  standard <- aspectCalc3(lat = lat, slope = neutSlope, aspect = neutAsp)
  return(new/standard)
}

##modified function from US script
removeOutlier <- function(dat, alpha,numIDvars){
  out <- foreach(curr = unique(as.character(dat$Zone)), .combine = rbind) %do% {
    temp <- dat[dat$Zone == curr,]
    md <- tryCatch(mahalanobis(temp[,-c(1:numIDvars),with = F],
                               center = colMeans(temp[,-c(1:numIDvars), with = F]), 
                               cov = cov(temp[,-c(1:numIDvars), with = F])), error = function(e) e)
    if(!inherits(md,"error")){
      ctf <- qchisq(1-alpha, df = ncol(temp)-1)
      outl <- which(md > ctf)
      cat("Removing", length(outl), "outliers from",curr, "; ")
      if(length(outl) > 0){
        temp <- temp[-outl,]
      }
    }
    temp
    
  }
  return(out)
}

##load data
###USA FIA data#####
us_veg <- fread("./InputData/AllUSandFIA16Apr2020_Veg.csv")
cols <- grep("Cover",colnames(us_veg))[1:4] ##only include cover 1,2,3,4
us_veg <- us_veg[,c(1,2,cols), with = F]
us_veg <- us_veg[,Cover := rowSums(.SD, na.rm = T),.SDcols = -c("PlotNumber","Species")]
us_veg <- us_veg[,.(PlotNumber,Species,Cover)]

envAll <- fread("./InputData/BECMaster19_Env.csv")
us_env <- fread("./InputData/AllUSandFIA16Apr2020_Env.csv")
envAll <- rbind(envAll, us_env)
load("./InputData/VegDat_Clean.RData")
vegData <- as.data.table(vegData)
vegData <- rbind(vegData, us_veg)
climDat <- fread("./InputData/EcoPlot_Locations_Normal_1961_1990MSY.csv")#, select = c("ID1","ID2","Latitude","Longitude", "Elevation"))
#fwrite(climDat ,"EcoPlot_Locations.csv")
##BGC climate data
# bgcDat <- fread("BGC_ClimDat_Veg_2000-2018.csv",select = c("Year", "ID1","ID2","Latitude", "DD5","Tmax_sm","FFP","Tave_sm","Tmin_wt","CMD"))
# bgcDat <- bgcDat[Year %in% 2000:2018,]
# bgcDat <- bgcDat[,lapply(.SD,mean),by = .(ID1,ID2,Latitude), .SDcols = -("Year")]
# fwrite(bgcDat,"BGC_ClimDat_Veg_2000-2018.csv")
bgcDat <- fread("./InputData/BECv11_250Pt_Normal_1961_1990MSY.csv")


##Prepare data
envAll <- envAll[Date > as.Date("2001-01-01"),] ##remove old plots
vegData <- vegData[Cover >= 1,] ### remove entries with low cover
envPlt <- unique(envAll$PlotNumber)
vegData <- vegData[PlotNumber %in% envPlt,]
envAll <- envAll[PlotNumber %in% vegData$PlotNumber,]

plotBGC <- envAll[Aspect >= 0 & Aspect <= 360,
                  c("PlotNumber","Longitude","Latitude","Zone","SubZone","SiteSeries","MoistureRegime","NutrientRegime","SlopeGradient","Aspect")]
vegSub <- merge(vegData, plotBGC, by = "PlotNumber")
vegSub <- vegSub[!is.na(Longitude) & !is.na(SlopeGradient),]
vegSub$Aspect[vegSub$SlopeGradient == 0] <- 180
vegSub <- vegSub[!is.na(Aspect) & Aspect != 999 & SlopeGradient != 999,] ##remove 999
##covert to degrees - part of function?
#vegSub$Slope <- atan(vegSub$SlopeGradient/100)
#vegSub$Slope <- vegSub$Slope*180/pi ##to radians
##calculate aspect ration
vegSub$AspectRatio <- aspectRatio(lat = vegSub$Latitude, slope = vegSub$Slope, aspect = vegSub$Aspect)
climSub <- climDat[,c("ID1","DD5","Tmax_sm","FFP","Tave_sm","Tmin_wt")]##select climate variables
vegSub <- merge(vegSub, climSub, by.x = "PlotNumber",by.y = "ID1")
##apply aspectratio adjustment to climate variables
vegSub[,13:18] <- lapply(vegSub[,13:18],FUN = function(x){x*vegSub$AspectRatio})
vegAll <- vegSub

###Since not all plots have a zone attached, for outlier analysis, 
###it was necessary to overlay with BEC wNA and use those zones where missing
ZoneDat <- fread("./InputData/NA_Plot_Zones.csv")
vegAll <- vegAll[,-("Zone")]
vegAll <- merge(vegAll,ZoneDat, by = "PlotNumber")

##Choose species and variable
spp <- "PINUCON"
var <- "DD5"

###The commented out bit below selects only plots that are in the maximum edatopic position for that species
###But it doesn't work with the US plot data, and doesn't seem to make much difference anyways. So I think ok to not use.

#veg <- vegAll[Species == spp,]
# temp <- veg[NutrientRegime != "",.(CovSum = sum(Cover)),by = .(Zone,MoistureRegime,NutrientRegime)]
# use <- temp[CovSum > 500,]
# vegSub <- merge(vegAll, use, by = c("Zone","MoistureRegime","NutrientRegime"), all = F)
# ##or without removing plots

##Select coastal plots or interior
# coastal <- c("CDF", "CWH", "MH", "MHRF", "CMA", "CMX", "CRF", "CWF")
# vegSub <- vegAll[!Zone %in% coastal,]

vegSub <- vegAll

climSub <- unique(vegSub[,c("PlotNumber","Zone", "AspectRatio", "DD5", "Tmax_sm", "FFP", "Tave_sm", "Tmin_wt")])
climSub <- climSub[Tmax_sm > -20,] ##remove rows with missing data
climSub <- removeOutlier(dat = climSub, alpha = 0.01, numIDvars = 3)
climSub <- climSub[!duplicated(climSub, by = "PlotNumber"),] ##there are a few duplicates, I'm not really sure why. I should look into this

brk <- seq(min(climSub[[var]]),max(climSub[[var]]),length.out = 40)
matHistAll <- hist(climSub[[var]], breaks = brk) ##create histogram of all plots for specified variable
histAll <- matHistAll$counts

###now for specified species
veg <- vegSub[Species == spp & PlotNumber %in% climSub$PlotNumber,]
veg <- veg[Tmax_sm > -10,]
veg <- veg[!duplicated(veg, by = "PlotNumber"),]
veg <- veg[,c(1:13,19,14:18),with = F]
veg <- removeOutlier(dat = veg, alpha = 0.1, numIDvars = 14) ###also remove outliers once it's just plots containing said species

library(weights) ##for weighted histogram
vegHist <- wtd.hist(veg[[var]], plot = T,weight = veg$Cover, breaks = brk)
minNumPlots <- 10 ##remove bins where there are few plots in total
brk <- brk[histAll > minNumPlots]
h <- vegHist$counts[histAll > minNumPlots]/histAll[histAll > minNumPlots] ##ratio of weighted veg hist to total num plots in each bin
if(length(brk) > length(h)){brk <- brk[-1]} ##depending on the variable sometimes breaks don't line up

plot(h ~ brk)
dat <- data.frame(x = brk, r = h) ##create data for nls
muStart <- brk[h == max(h)] ###starting values for optimisation
sdStart <- muStart/3

#nrm <- formula(r ~ k*exp(-1/2*(x-mu)^2/sigma^2))
library(minpack.lm) ##this nls algorithm is better
library(sn)

###non-linear least squares fit for optimal parameters
###note that we have to multiply the dsn function by a constant, otherwise integrates to 1
###Occasionally the nls algorithm gets stuck if the starting values aren't close enough - if this happens, try and make the initial values better
fit <- nlsLM(r ~ k * dsn(x,xi = mu,omega = sigma, alpha = a), start=c(mu=muStart,sigma=sdStart,k=muStart,a = 0) , 
             data = dat)
v <- summary(fit)$parameters[,"Estimate"] ###v contains the 4 parameters for the best fit curve

##set up BGC subzone data
###note that you can adjust the neutral parameters in the aspectRation function
bgDat <- bgcDat[,c("ID2","Latitude",var), with = F]
colnames(bgDat)[3] <- "Var"
t1 <- aspectRatio(lat = bgDat$Latitude, slope = 0, aspect = 180,neutSlope = 30,neutAsp = 90)
bgDat$Warm <- bgDat$Var*t1 ###"neutral" aspect
t1 <- aspectRatio(lat = bgDat$Latitude, slope = 35, aspect = 225,neutSlope = 30,neutAsp = 90)
bgDat$vWarm <- bgDat$Var*t1 ###"warm" aspect

###select units to display
bgDat <- bgDat[ID2 %in% c("BAFAun", "BWBSdk", "ESSFdk1", "SBSdk","IDFdk3", "IDFxc","CWHmm1", "ICHmw1", "PPxh1"),]

library(ggthemes)

###horizonal boxplot
p2 <- ggplot(bgDat)+
  geom_boxploth(aes(x = Warm, y = ID2, group = ID2))+
  geom_boxploth(aes(x = vWarm, y = ID2, group = ID2), colour = "red")+
  geom_boxploth(aes(x = Var, y = ID2, group = ID2), colour = "blue")+
  xlim(c(min(dat$x-10),max(dat$x + 10)))+
  labs(y = "BGC", x = var)+
  theme_few()

##function to plot best fit skewed normal
plotNrm <- function(x) {v[3] * dsn(x, v[1],v[2],v[4])}
mX <- optimize(plotNrm, interval = c(min(dat$x-10),max(dat$x + 10)), maximum = T) ##find maximum
probMx <- psn(mX$maximum,xi = v[1],omega = v[2],alpha = v[4]) ##find area below maximum
##set up quantiles for colouring
quants <- qsn(c(probMx/10,probMx/2,probMx+(1-probMx)*0.5,probMx+(1-probMx)*0.9),xi = v[1],omega = v[2],alpha = v[4])
quants <- c(min(dat$x-10),quants,max(dat$x + 10))
if(quants[1] > quants[2]){quants[1] <- quants[2]-10}
if(quants[6] < quants[5]){quants[6] <- quants[5]+10}

##plot response curve
p1 <- ggplot(dat)+
  geom_point(aes(x = x, y = r),color = "black")+
  stat_function(fun = plotNrm, colour = "black")+
  geom_ribbon(data = data.frame(y = plotNrm(seq(quants[3],quants[4],0.1)),x = seq(quants[3],quants[4],0.1)), aes(ymin = 0,ymax = y,x =x), alpha = 0.5, fill = "green")+
  geom_ribbon(data = data.frame(y = plotNrm(seq(quants[2],quants[3],0.1)),x = seq(quants[2],quants[3],0.1)), aes(ymin = 0,ymax = y,x =x), alpha = 0.5, fill = "yellow")+
  geom_ribbon(data = data.frame(y = plotNrm(seq(quants[4],quants[5],0.1)),x = seq(quants[4],quants[5],0.1)), aes(ymin = 0,ymax = y,x =x), alpha = 0.5, fill = "yellow")+
  geom_ribbon(data = data.frame(y = plotNrm(seq(quants[1],quants[2],0.1)),x = seq(quants[1],quants[2],0.1)), aes(ymin = 0,ymax = y,x =x), alpha = 0.5, fill = "red")+
  geom_ribbon(data = data.frame(y = plotNrm(seq(quants[5],quants[6],0.1)),x = seq(quants[5],quants[6],0.1)), aes(ymin = 0,ymax = y,x =x), alpha = 0.5, fill = "red")+
  xlim(c(min(dat$x-10),max(dat$x + 10)))+
  ylim(c(0,NA))+
  labs(y = "Standardised Frequency", title = paste(spp, var))+
  theme_few()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

###stack vertically using cowplot library
cowplot::plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = c(0.66, 0.34))





##############################################################################################
#        I'm working on this part to create the data tables efficiently - not done yet
##########################################################################################

dat <- foreach(Spp = Spp.list,.combine = rbind) %:%
  foreach(Var = var.list, .combine = rbind) %do% {
    vegSub <- vegAll
    
    climSub <- unique(vegSub[,c("PlotNumber","Zone", "AspectRatio", "DD5", "Tmax_sm", "FFP", "Tave_sm", "Tmin_wt")])
    climSub <- climSub[Tmax_sm > -20,]
    climSub <- removeOutlier(dat = climSub, alpha = 0.01, numIDvars = 3)
    climSub <- climSub[!duplicated(climSub, by = "PlotNumber"),]
    
    brk <- seq(min(climSub[[var]]),max(climSub[[var]]),length.out = 40)
    matHistAll <- hist(climSub[[var]],plot = F, breaks = brk)
    histAll <- matHistAll$counts
    veg <- vegSub[Species == spp & PlotNumber %in% climSub$PlotNumber,]
    veg <- veg[Tmax_sm > -10,]
    veg <- veg[!duplicated(veg, by = "PlotNumber"),]
    veg <- veg[,c(1:13,19,14:18),with = F]
    veg <- removeOutlier(dat = veg, alpha = 0.1, numIDvars = 14)
    
    vegHist <- wtd.hist(veg[[var]], plot = F,weight = veg$Cover, breaks = brk)
    brk <- brk[histAll > 10]
    h <- vegHist$counts[histAll > 10]/histAll[histAll > 10]
    if(length(brk) > length(h)){brk <- brk[-1]}
    dat <- data.frame(x = brk, r = h)
    muStart <- brk[h == max(h)]
    sdStart <- muStart/3
    
    fit <- nlsLM(r ~ k * dsn(x,xi = mu,omega = sigma, alpha = a), start=c(mu=muStart,sigma=sdStart,k=muStart,a = 0), 
                 data = dat)
    v <- summary(fit)$parameters[,"Estimate"] ###v contains the three parameters for the best fit curve
    
    plotNrm <- function(x) {v[3] * dsn(x, v[1],v[2],v[4])}
    mX <- optimize(plotNrm, interval = c(min(dat$x-10),max(dat$x + 10)), maximum = T) ##find maximum
    probMx <- psn(mX$maximum,xi = v[1],omega = v[2],alpha = v[4])
    quants <- qsn(c(probMx/10,probMx/2,probMx+(1-probMx)*0.5,probMx+(1-probMx)*0.9),xi = v[1],omega = v[2],alpha = v[4])
    quants <- c(min(dat$x-10),quants,max(dat$x + 10))
    if(quants[1] > quants[2]){quants[1] <- quants[2]-10}
    if(quants[6] < quants[5]){quants[6] <- quants[5]+10}
    
  }












###parametric fitting of curves -- old###
vegAvg <- setDT(veg)[,.(Cover = mean(Cover)),by = .(DD5)]
plot(Cover ~ DD5, data = vegAvg)
dat <- data.table(x=as.numeric(as.character(vegAvg$DD5)), r=vegAvg$Cover)

if(aspect == "cool"){
  vegSub <- vegSub[(Aspect >= 0 & Aspect <= 90) | Aspect == 360,] ###warm or cool aspect
}else if(aspect == "warm"){
  vegSub <- vegSub[(Aspect >= 180 & Aspect <= 270),] ###warm or cool aspect
}

require(scales)
dat$x <- rescale(dat$x, to = c(-1,1))
plot(r ~ x, data = dat)

nrm <- formula(r ~ k*exp(-1/2*(x-mu)^2/sigma^2))
poly <- formula(r ~ a*x^2 + b*x + c)

fit <- nls(poly, start=c(a=-1,b=0,c=20) , 
           data = dat, nls.control(maxiter = 1000, minFactor = 0.000001))
fit <- nls(nrm, start=c(mu=1500,sigma=1000,k=0.01) , 
           data = dat, nls.control(maxiter = 1000, minFactor = 0.000001))
v <- summary(fit)$parameters[,"Estimate"] ###v contains the three parameters for the best fit curve
plot(r ~ x, data = dat, pch = 4)
plot(function(x) v[1]*x^2 +v[2]*x + v[3],add=T,xlim=range(dat$x), col = "red")
plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),add=T,xlim=range(dat$x), col = "red")

veg$Zone <- gsub("[[:lower:]]|[[:digit:]]","", veg$BGC_noSpace)
vegAvg <- setDT(veg)[,.(Cover = mean(Cover)),by = .(DD5,Zone)]
plot(Cover ~ DD5, data = vegAvg[Zone == "IDF",])


PlotCount <- count (vegData, c("PlotNumber"))
#PlotCount2 <- aggregate(Species ~ PlotNumber, vegData, length)
PlotCount3 <- as.data.frame (unique(vegData$PlotNumber)) 
#vegData <- read.table("BecMaster15VegData.txt", header = TRUE)

codeCross <- read.csv("TreeCodeCrosswalk.csv", stringsAsFactors = FALSE)


#####OLD CODE############

drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user = "postgres", password = "Kiriliny41", host = "localhost", 
                 port = 5432, dbname = "bec_plot_data")
dbGetQuery(con, "SELECT column_name FROM information_schema.columns WHERE table_name = 'env_full'")

### set to get get climate data ###
location <- dbGetQuery(con, "SELECT plotnumber,longitude,latitude from env_full WHERE longitude IS NOT NULL")
location$longitude <- -(location$longitude)
wd <- tk_choose.dir(); setwd(wd)
dem <- raster("WNA_DEM_3005_clipped.tif") ###Read DEM
p <- st_as_sf(location, coords = c("longitude", "latitude"), 
              crs = "+init=epsg:4326", agr = "constant")
p <- st_transform(p, st_crs(dem))
p$el <- raster::extract(dem,p)
p <- st_transform(p,"+init=epsg:4326")
location$el <- p$el
location$ID2 <- seq_along(location$el)
location <- location[,c("plotnumber","ID2", "latitude","longitude","el")]
colnames(location) <- c("ID1","ID2", "lat","long", "el")
fwrite(location, "BECPlot_Location.csv")

##get us climate data##
library(stars)
env <- fread(file.choose())
env <- env[!is.na(Latitude) & !is.na(Longitude),.(PlotNumber,Longitude,Latitude)]
dem <- raster("../WNA_Climate/WNA_DEM_3005_clipped.tif") ###Read DEM
p <- st_as_sf(env, coords = c("Longitude", "Latitude"), 
              crs = "+init=epsg:4326", agr = "constant")
p <- st_transform(p, st_crs(dem))
p <- p[!is.na(st_dimension(p)),]
p$el <- raster::extract(dem,p,na.rm = T)
p <- st_transform(p,"+init=epsg:4326")
p <- cbind(p,st_coordinates(p)) %>% st_drop_geometry()
p$ID2 <- seq_along(p$el)
p <- p[,c("PlotNumber","ID2","Y","X","el")] %>% set_colnames(c("ID1","ID2","lat","long","el"))
fwrite(p, "USFIA_Locations.csv")
###############################

p <- envAll[,.(PlotNumber,Latitude,Longitude,Elevation)]
p <- p[complete.cases(p),]
p$ID2 <- seq_along(p$PlotNumber)
p <- p[,.(PlotNumber,ID2,Latitude,Longitude,Elevation)] %>% set_colnames(c("ID1","ID2","lat","long","el"))
p$long <- -abs(p$long)
fwrite(p,"BEC_NA_Locations.csv")

###join with WNA shapefile to get BGC zones
# library(sf)
# wnaBEC <- st_read(dsn = "../WNA_Climate", layer = "WNA_BGC_v11_31Dec2019")
# plotLoc <- unique(vegAll[,.(PlotNumber,Latitude,Longitude,Zone)])
# plotLoc <- st_as_sf(plotLoc, coords = c("Longitude", "Latitude"), 
#                     crs = "+init=epsg:4326", agr = "constant") %>% st_transform(st_crs(wnaBEC))
# temp <- st_join(plotLoc, wnaBEC, join = st_intersects, left = T)
# temp <- st_drop_geometry(temp)
# temp$Zone.y <- as.character(temp$Zone.y)
# temp$Zone <- ifelse(is.na(temp$Zone.y),temp$Zone.x,temp$Zone.y)
# temp <- temp[temp$Zone != "",]
# ZoneDat <- temp[,c("PlotNumber","Zone")]

# ###
# temp <- climDat[,.(ID1,Elevation)] %>% set_colnames(c("PlotNumber","ElDem"))
# t2 <- envAll[,.(PlotNumber,Elevation)]
# temp <- merge(temp,t2,by = "PlotNumber")
# veg <- merge(veg, temp, by = "PlotNumber", all.y = F)
# veg$ElDiff <- abs(veg$ElDem - veg$Elevation)
# ####

# fit <- nls(nrm2, start=c(mu=muStart,sigma=sdStart,k=sum(dat$r)) , 
#            data = dat, nls.control(maxiter = 1000, minFactor = 0.000001))
# fit <- nlsLM(r ~ k * dnorm(x,mu,sigma), start=c(mu=muStart,sigma=sdStart,k=sum(dat$r)) , 
#            data = dat)