# load info on herds
load("New_S.Dublin_SICR_ini.Rdata")
rm("model", "NUTS_10M", "sweden")

# mark dairy herds and HR herds
dn <- geo$row[geo$dairy] # node number of dairy herds

HR <- geo$oland | geo$hr.skane  # high risk area (all herds)
HR <- HR[dn]  # keep only dairy herds

# mark Öland
ol <- geo$oland[dn]

# list useful var references
infect <- c("Ic", "Iy", "Ia")
# calves <- c("Sc", "Ic", "Cc", "Rc")
# youngs <- c("Sy", "Iy", "Cy", "Ry")
adults <- c("Sa", "Ia", "Ca", "Ra")
all <- c("Sc", "Ic", "Cc", "Rc", "Sy", "Iy", "Cy", "Ry", "Sa", "Ia", "Ca", "Ra")

# take out active status --> does not change over years
dt <- readRDS("SimulData/Model_17/trj_0001.RDS")
dt$week <- as.numeric(as.factor(dt$time))
dt$tot <- rowSums(dt[,all])
actHerd <- dt$tot > 0

time <- unique(dt$time)

BHP <- c()
BHP_O <- c()
AVG <- dt[1:2]
AVG <- cbind(AVG, as.data.frame(matrix(0, ncol = 12, nrow = nrow(dt))))
names(AVG) <- names(dt[1:14])
IC <- IY <- IA <- c()
CC <- CY <- CA <- c()
RC <- RY <- RA <- c()
CAL <- YOU <- ADU <- c()




# __________ LOOP _____________

# load simulated data
for (i in 1:100) {
  print(paste("Running simulation", i))
  
  # load trj data
  dt <- readRDS(sprintf("SimulData/Model_17/trj_%04i.RDS",i))
  dt$week <- as.numeric(as.factor(dt$time))
  
  dt$inf <- rowSums(dt[,c("Ic","Iy","Ia")])
  dt$adu <- rowSums(dt[,c("Sa","Ia","Ca","Ra")])
  dt$you <- rowSums(dt[,c("Sy","Iy","Cy","Ry")])
  dt$cal <- rowSums(dt[,c("Sc","Ic","Cc","Rc")])
  
  dt$det <- ifelse(dt$adu==0, 0, dt$Ra / dt$adu)
  
  dt$ih <- dt$inf > 0
  dt$dh <- dt$det >= 0.15
  
  
  # extract status matrix
  ########################
  active <- as.data.frame(do.call("cbind", tapply(actHerd, dt$week, function(x) { x == TRUE})))
  infected <- as.data.frame(do.call("cbind", tapply(dt$ih, dt$week, function(x) { x == TRUE})))
  detected <- as.data.frame(do.call("cbind", tapply(dt$dh, dt$week, function(x) { x == TRUE})))
  
  herd_pr <- colSums(infected)/colSums(active)
  herd_pr_O <- colSums(infected[ol,])/colSums(active[ol,])
  
  
  ic <- as.data.frame(matrix(dt$Ic, nrow= length(unique(dt$node)), ncol = 444))
  cc <- as.data.frame(matrix(dt$Cc, nrow= length(unique(dt$node)), ncol = 444))
  rc <- as.data.frame(matrix(dt$Rc, nrow= length(unique(dt$node)), ncol = 444))
  iy <- as.data.frame(matrix(dt$Iy, nrow= length(unique(dt$node)), ncol = 444))
  cy <- as.data.frame(matrix(dt$Cy, nrow= length(unique(dt$node)), ncol = 444))
  ry <- as.data.frame(matrix(dt$Ry, nrow= length(unique(dt$node)), ncol = 444))
  ia <- as.data.frame(matrix(dt$Ia, nrow= length(unique(dt$node)), ncol = 444))
  ca <- as.data.frame(matrix(dt$Ca, nrow= length(unique(dt$node)), ncol = 444))
  ra <- as.data.frame(matrix(dt$Ra, nrow= length(unique(dt$node)), ncol = 444))
  cal <- as.data.frame(matrix(dt$cal, nrow= length(unique(dt$node)), ncol = 444))
  you <- as.data.frame(matrix(dt$you, nrow= length(unique(dt$node)), ncol = 444))
  adu <- as.data.frame(matrix(dt$adu, nrow= length(unique(dt$node)), ncol = 444))
  isInf <- as.data.frame(matrix(dt$ih, nrow= length(unique(dt$node)), ncol = 444))
  
  
  # stack the results of each simulation -----------
  BHP <- rbind(BHP, herd_pr)
  BHP_O <- rbind(BHP_O, herd_pr_O)
  AVG[3:14] <- AVG[3:14] + dt[3:14]
  
  IC <- rbind(IC, colSums(ic*isInf))
  CC <- rbind(CC, colSums(cc*isInf))
  RC <- rbind(RC, colSums(rc*isInf))
  IY <- rbind(IY, colSums(iy*isInf))
  CY <- rbind(CY, colSums(cy*isInf))
  RY <- rbind(RY, colSums(ry*isInf))
  IA <- rbind(IA, colSums(ia*isInf))
  CA <- rbind(CA, colSums(ca*isInf))
  RA <- rbind(RA, colSums(ra*isInf))
  CAL <- rbind(CAL, colSums(cal*isInf))
  YOU <- rbind(YOU, colSums(you*isInf))
  ADU <- rbind(ADU, colSums(adu*isInf))
  
  rm(list=setdiff(ls(), c("infect", "adults", "all", "HR", "dn", "ol",
                          "actHerd",  "time", "BHP","BHP_O", "AVG", 
                          "IC", "CC", "RC", "IY", "CY", "RY", 
                          "IA", "CA", "RA", "CAL", "YOU", "ADU")))
  gc()
  
}

save("HR", "dn", "ol", "BHP","BHP_O", "AVG", "IC", "CC", "RC", "IY", "CY", "RY", 
     "IA", "CA", "RA", "CAL", "YOU", "ADU", file = "DataSummary_M17.RData")

