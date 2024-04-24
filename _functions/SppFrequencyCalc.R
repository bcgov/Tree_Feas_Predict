#dat = CCISS_ESuit; eda = edatopic
frequency_bgc <- function(dat, eda){
QRY <- paste0("select bgc,ss_nospace,sppsplit,spp,feasible
              from feasorig where feasible in (1,2,3,4,5)")
#feas <- as.data.table(dbGetQuery(con, QRY))
feas <- as.data.table(dat, QRY)
eda <- as.data.table(eda)
minDist <- feas
#minDist <- feas[,.SD[feasible == min(feasible, na.rm = T)],by = .(bgc,sppsplit)]
###unremark to show where a non-local species has already been added to maps
#added <- minDist[OHR == "OHR"] %>% dplyr::select(bgc, sppsplit) %>% distinct %>% mutate(Freq = "Expanded")
minDist <- minDist[!OHR == "OHR"]
minDist <- minDist[feasible %in% c(1,2,3),]
# abUnits <- minDist[grep("[[:alpha:]] */[[:alpha:]]+$",ss_nospace),]
# noAb <- minDist[!grepl("[[:alpha:]] */[[:alpha:]]+$",ss_nospace),]
# abUnits <- eda[abUnits, on = "ss_nospace"] ##merge
# abUnits <- abUnits[,.(Temp = if(any(grepl("C4",edatopic))) paste0(ss_nospace,"_01") else ss_nospace, feasible = feasible[1]),
#                    by = .(bgc,ss_nospace,sppsplit,spp)]
# abUnits[,ss_nospace := NULL]
# setnames(abUnits,old = "Temp",new = "ss_nospace")
# minDist <- rbind(noAb,abUnits)
toMatch <- c("/01", "/101")
minDist[,ID := if(any(grepl(paste(toMatch,collapse="|"), ss_nospace)) & feasible[1] == 1) T else F, by = .(bgc,sppsplit)]
minDist[,Freq := NA_character_]
minDist[(ID),Freq := "High"]

minDist2 <- minDist[ID == F,]
minDist2[,ID := if(any(grepl(paste(toMatch,collapse="|"), ss_nospace))) T else F, by = .(bgc,sppsplit)]
minDist2[(ID),Freq := "Moderate"]

minDist3 <- minDist2[ID == F,]
minDist3[,Freq := "Low"]
# minEda <- eda[minDist3, on = "ss_nospace"]
# minEda <- minEda[,.(AvgEda = mean(smr)), by = .(bgc,sppsplit,ss_nospace,feasible)]
# minEda[,CentEda := abs(AvgEda - 3.5)]
# minEda <- minEda[,.SD[CentEda == min(CentEda, na.rm = T)], by = .(bgc,sppsplit)]
# lookupTab <- data.table(AvgEda = c(0,2,5,7),Freq = c("Low","Low","Low","Low"))
# temp <- lookupTab[minEda, on = "AvgEda", roll = T]

t1 <- minDist[!is.na(Freq),.(Freq = Freq[1]), by = .(bgc,sppsplit)]
t2 <- minDist2[!is.na(Freq),.(Freq = Freq[1]), by = .(bgc,sppsplit)]
t3 <- minDist3[!is.na(Freq),.(Freq = Freq[1]), by = .(bgc,sppsplit)]
#t3 <- temp[,.(Freq = Freq[1]), by = .(bgc,sppsplit)]
allFreq <- rbind(t1,t2,t3)

# if(nrow(added) > 0){
#   added2 <- anti_join(added, allFreq, by = c("bgc","sppsplit"))
#   allFreq <- rbind(allFreq, added2)
  
#}
return(allFreq)
}
