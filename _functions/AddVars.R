addVars <- function(dat){
  dat$PPT_MJ <- dat$PPT05 + dat$PPT06  # MaY/June precip
  dat$PPT_JAS <- dat$PPT07 + dat$PPT08 + dat$PPT09  # July/Aug/Sept precip
  dat$PPT.dormant <- dat$PPT_at + dat$PPT_wt  # for calculating spring deficit
  dat$CMD.def <- 500 - (dat$PPT.dormant)  # start of growing season deficit original value was 400 but 500 seems better
  dat$CMD.def[dat$CMD.def < 0] <- 0  #negative values set to zero = no deficit
  dat$CMDMax <- dat$CMD07
  dat$CMD.total <- dat$CMD.def + dat$CMD
  dat$CMD.grow <- dat$CMD05 + dat$CMD06 +dat$CMD07 +dat$CMD08 +dat$CMD09
  dat$DD5.grow <- dat$DD5_05 + dat$DD5_06 + dat$DD5_07 + dat$DD5_08 + dat$DD5_09
  dat$CMDMax <- dat$CMD07 # add in so not removed below
  dat$DDgood <- dat$DD5 - dat$DD18
  dat$DDnew <- (dat$DD5_05 + dat$DD5_06 +dat$DD5_07  + dat$DD5_08)  - (dat$DD18_05 + dat$DD18_06 +dat$DD18_07 +dat$DD18_08)
  dat$TmaxJuly <- dat$Tmax07
  dat$tmaxSum20 <-  dat %>%
    dplyr::select( starts_with("Tmax")) %>% -20 %>%
    mutate_all(., list(~ifelse(. < 0, 0 , .))) %>%
    rowSums(na.rm = TRUE)
  
  dat$tmaxSum25 <-  dat %>%
    dplyr::select( starts_with("Tmax")) %>% -25 %>%
    mutate_all(., list(~ifelse(. < 0, 0 , .))) %>%
    rowSums(na.rm = TRUE)
  
  dat$tmaxSum30 <-  dat %>%
    dplyr::select( starts_with("Tmax")) %>% -30 %>%
    mutate_all(., list(~ifelse(. < 0, 0 , .))) %>%
    rowSums(na.rm = TRUE)
  
  dat$tmaxSum35 <-  dat %>%
    dplyr::select( starts_with("Tmax")) %>% -35 %>%
    mutate_all(., list(~ifelse(. < 0, 0 , .))) %>%
    rowSums(na.rm = TRUE)
  
  # remove some redundant variables considered undesireable
  month <- c("01", "02", "03", "04", "05", "06","07", "08", "09", "10", "11", "12")
  dat <- dat  %>% dplyr::select(-ends_with(month)) %>% #removes all monthly variables
    dplyr::select(-starts_with("Rad"))# %>% ##remove other non-biological variables
  #dplyr::select(-starts_with("RH")) %>%
  #dplyr::select (-contains("DD_0")) %>%
  #dplyr::select  (-contains("DD18")) %>%
  #dplyr::select  (-contains("DD_18"))  %>%
  #dplyr::select( -PPT_sp, -PAS_sm, -PPT_at, -PAS_at, - MAP, -TD, -MAT, -FFP)
  return(dat)
}