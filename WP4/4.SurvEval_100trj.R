# load info on herds
load("New_S.Dublin_SICR_ini.Rdata")
rm("model", "NUTS_10M", "sweden")

# mark dairy herds and HR herds
dn <- geo$row[geo$dairy] # node number of dairy herds

HR <- geo$oland | geo$hr.skane  # high risk area (all herds)
HR <- HR[dn]  # keep only dairy herds

# list useful var references
infect <- c("Ic", "Iy", "Ia")
adults <- c("Sa", "Ia", "Ca", "Ra")
samp <- c(7, 20, 33, 46)

ACT <- c()
INF <- c()
ACT_t <- c()
INF_t <- c()
ACT_r <- c()
INF_r <- c()

TR_Q1 <- c()
TR_Q2 <- c()
TR_Q3 <- c()
TR_Q4 <- c()
TR_Y <- c()
RB_Q1 <- c()
RB_Q2 <- c()
RB_Q3 <- c()
RB_Q4 <- c()
RB_Y <- c()


# __________ LOOP _____________

# load simulated data
for (i in 1:100) {
  print(paste("Running simulation", i))
  
  # load trj data
  dt <- readRDS(sprintf("SimulData/Model_17/trj_%04i.RDS",i))
  
  
  # keep only 2013
  dt13 <- subset(dt, time>="2013-01-01")
  dt13$week <- as.numeric(as.factor(dt13$time))
  
  dt13$tot <- rowSums(dt13[,3:14])
  dt13$inf <- rowSums(dt13[,c("Ic", "Iy", "Ia")])
  dt13$adu <- rowSums(dt13[,c("Sa", "Ia", "Ca", "Ra")])
  dt13$det <- ifelse(dt13$adu==0, 0, dt13$Ra / dt13$adu)
  
  dt13$ah <- dt13$tot > 0
  dt13$ih <- dt13$inf > 0
  dt13$dh <- dt13$det >= 0.15
  
  
  # extract status matrix
  ########################
  active <- as.data.frame(do.call("cbind", tapply(dt13$ah, dt13$week, function(x) { x == TRUE})))
  infected <- as.data.frame(do.call("cbind", tapply(dt13$ih, dt13$week, function(x) { x == TRUE})))
  detected <- as.data.frame(do.call("cbind", tapply(dt13$dh, dt13$week, function(x) { x == TRUE})))
  
  
  ############################
  # TRADITIONAL SURVEILLANCE #
  ############################
  
  # Every herd sampled every quarter of the year
  
  # update status 4 weeks after detection --> herd no longer under surveillance (ah=FALSE, ih=FALSE)
  ###########################################################################################
  active_t <- active[which(rowSums(active)>0),]
  infected_t <- infected[which(rowSums(active)>0),]
  detected_t <- detected[which(rowSums(active)>0),]
  
  # update after 1st quarter (week 7) --------------
  # remove TP herds from surveillance
  active_t[infected_t[,7] & detected_t[,7],11:52] <- FALSE
  infected_t[infected_t[,7] & detected_t[,7],11:52] <- FALSE
  detected_t[infected_t[,7] & detected_t[,7],11:52] <- FALSE

  # update after 2nd quarter (week 20)--------------------
  active_t[infected_t[,20] & detected_t[,20],24:52] <- FALSE
  infected_t[infected_t[,20] & detected_t[,20],24:52] <- FALSE
  detected_t[infected_t[,20] & detected_t[,20],24:52] <- FALSE

  # update after 3rd quarter (week 33)
  active_t[infected_t[,33] & detected_t[,33],37:52] <- FALSE
  infected_t[infected_t[,33] & detected_t[,33],37:52] <- FALSE
  detected_t[infected_t[,33] & detected_t[,33],37:52] <- FALSE

  # update after 4th quarter (week 46) 
  active_t[infected_t[,46] & detected_t[,46],50:52] <- FALSE
  infected_t[infected_t[,46] & detected_t[,46],50:52] <- FALSE
  detected_t[infected_t[,46] & detected_t[,46],50:52] <- FALSE

  
  # Summary per quarter --------------------
  for (j in samp) {
    
    nam <- paste("trad_Q", which(samp==j), sep = "")
    
    assign(nam, data.frame(act = sum(active_t[,j]),   # active herds
                           inf = sum(infected_t[,j]), # infected herds
                           samp = sum(active_t[,j]),   # sampled herds
                           det = sum(detected_t[,j]), # detected herds
                           TP = sum(infected_t[,j] & detected_t[,j]), # True Positive herds
                           FP = sum(!infected_t[,j] & detected_t[,j]), # False Positive herds
                           FN = sum(infected_t[,j] & !detected_t[,j]), # False Negative herds
                           TN = sum(active_t[, j] & !infected_t[,j] & !detected_t[,j]), # True Negative herds
                           Prev = round((sum(infected_t[,j]) / sum(active_t[,j]))*100, 3),  # true prevalence in the population
                           TPr = round((sum(infected_t[,j]) / sum(active_t[,j]))*100, 3),  # true prevalence in the sample
                           APr = round((sum(detected_t[,j]) / sum(active_t[,j]))*100, 3))) # apparent prevalence in the sample
    # apparent prevalence
    
  }
  
  
  # Summary per year -------------------------
  trad_Y <- data.frame(act = sum(rowSums(active_t)>0),
                       inf = sum(rowSums(infected_t)>0),  # OBS! herds infected for only short time??
                       # outb = sum(rowSums(infected_t)>1),
                       samp = sum(trad_Q1$samp, trad_Q2$samp, trad_Q3$samp, trad_Q4$samp),
                       det = sum(trad_Q1$det, trad_Q2$det, trad_Q3$det, trad_Q4$det),
                       TP = sum(trad_Q1$TP, trad_Q2$TP, trad_Q3$TP, trad_Q4$TP),
                       FP = sum(trad_Q1$FP, trad_Q2$TP, trad_Q3$TP, trad_Q4$TP),
                       FN = sum(trad_Q1$FN, trad_Q2$FN, trad_Q3$FN, trad_Q4$FN),
                       TN = sum(trad_Q1$TN, trad_Q2$TN, trad_Q3$TN, trad_Q4$TN))
  
  trad_Y$prev <- with(trad_Y, round((inf/act)*100, 2))
  trad_Y$Se <- with(trad_Y, round(TP/(TP+FN), 2))   # Surveillance Se
  trad_Y$Sp <- with(trad_Y, round(TN/(TN+FP), 2))   # Surveillance Sp
  trad_Y$DF <- with(trad_Y, round(TP/inf, 2))       # Detection Fraction 
  trad_Y$cost <- with(trad_Y, round(samp/TP, 0))    # "cost" of detection = nr. of samples to get one true detection
  
  
  
  
  ###########################
  # RISK-BASED SURVEILLANCE #
  ###########################
  
  # Every herd in HR area sampled every quarter of the year, herds in LR area sampled only at Q3
  
  # update status after detection --> herd no longer under surveillance (ah=FALSE, ih=FALSE)
  ###########################################################################################
  active_r <- active[which(rowSums(active)>0),]
  infected_r <- infected[which(rowSums(active)>0),]
  detected_r <- detected[which(rowSums(active)>0),]
  HR13 <- HR[which(rowSums(active)>0)]
  
  # update after 1st quarter (week 7) --------------
  # remove TP herds from surveillance --> OBS only HR under surveillance for this quarter
  active_r[HR13 & infected_r[,7] & detected_r[,7],11:52] <- FALSE
  infected_r[HR13 & infected_r[,7] & detected_r[,7],11:52] <- FALSE
  detected_r[HR13 & infected_r[,7] & detected_r[,7],11:52] <- FALSE

  # update after 1st quarter (week 20) --------------
  # -->  only HR under surveillance for this quarter
  active_r[HR13 & infected_r[,20] & detected_r[,20],24:52] <- FALSE
  infected_r[HR13 & infected_r[,20] & detected_r[,20],24:52] <- FALSE
  detected_r[HR13 & infected_r[,20] & detected_r[,20],24:52] <- FALSE

  # update after 3rd quarter (week 33)
  # --> all herds under surveillance
  active_r[infected_r[,33] & detected_r[,33],37:52] <- FALSE
  infected_r[infected_r[,33] & detected_r[,33],37:52] <- FALSE
  detected_r[infected_r[,33] & detected_r[,33],37:52] <- FALSE

  # update after 4th quarter (week 46)
  # -->  only HR under surveillance for this quarter
  active_r[HR13 & infected_r[,46] & detected_r[,46],50:52] <- FALSE
  infected_r[HR13 & infected_r[,46] & detected_r[,46],50:52] <- FALSE
  detected_r[HR13 & infected_r[,46] & detected_r[,46],50:52] <- FALSE


  # Summary per quarter --------------------
  for (k in samp[-3]) {
    
    nam2 <- paste("risk_Q", which(samp==k), sep = "")
    
    assign(nam2, data.frame(act = sum(active_r[,k]),   # active herds
                            inf = sum(infected_r[,k]), # infected herds
                            samp = sum(active_r[HR13,k]),   # sampled herds  OBS Q3 is different
                            det = sum(detected_r[HR13,k]), # detected herds
                            TP = sum(infected_r[HR13,k] & detected_r[HR13,k]), # True Positive herds
                            FP = sum(!infected_r[HR13,k] & detected_r[HR13,k]), # False Positive herds
                            FN = sum(infected_r[HR13,k] & !detected_r[HR13,k]), # False Negative herds
                            TN = sum(active_r[HR13, k] & !infected_r[HR13,k] & !detected_r[HR13,k]), # True Negative herds
                            Prev = round((sum(infected_r[,k]) / sum(active_r[,k]))*100, 3),  # true prevalence in the population
                            TPr = round((sum(infected_r[HR13,k]) / sum(active_r[HR13,k]))*100, 3),  # true prevalence (in the sample)
                            APr = round((sum(detected_r[HR13,k]) / sum(active_r[HR13,k]))*100, 3))) # apparent prevalence (in the sample)
    
  }
  
  # OBS! in Q3 all herds are sampled
  
  risk_Q3 <- data.frame(act = sum(active_r[,k]),   # active herds
                        inf = sum(infected_r[,k]), # infected herds
                        samp = sum(active_r[,k]),   # sampled herds  OBS Q3 is different
                        det = sum(detected_r[,k]), # detected herds
                        TP = sum(infected_r[,k] & detected_r[,k]), # True Positive herds
                        FP = sum(!infected_r[,k] & detected_r[,k]), # False Positive herds
                        FN = sum(infected_r[,k] & !detected_r[,k]), # False Negative herds
                        TN = sum(active_r[, k] & !infected_r[,k] & !detected_r[,k]), # True Negative herds
                        Prev = round((sum(infected_r[,k]) / sum(active_r[,k]))*100, 3),  # true prevalence in the population
                        TPr = round((sum(infected_r[,k]) / sum(active_r[,k]))*100, 3),  # true prevalence (in the sample)
                        APr = round((sum(detected_r[,k]) / sum(active_r[,k]))*100, 3))  # apparent prevalence (in the sample)
  
  # Summary per year -------------------------
  risk_Y <- data.frame(act = sum(rowSums(active_r)>0),
                       inf = sum(rowSums(infected_r)>0),  # OBS! herds infected for only short time??
                       # outb = sum(rowSums(infected_r)>1),
                       samp = sum(risk_Q1$samp, risk_Q2$samp, risk_Q3$samp, risk_Q4$samp),
                       det = sum(risk_Q1$det, risk_Q2$det, risk_Q3$det, risk_Q4$det),
                       TP = sum(risk_Q1$TP, risk_Q2$TP, risk_Q3$TP, risk_Q4$TP),
                       FP = sum(risk_Q1$FP, risk_Q2$TP, risk_Q3$TP, risk_Q4$TP),
                       FN = sum(risk_Q1$FN, risk_Q2$FN, risk_Q3$FN, risk_Q4$FN),
                       TN = sum(risk_Q1$TN, risk_Q2$TN, risk_Q3$TN, risk_Q4$TN))
  
  risk_Y$prev <- with(risk_Y, round((inf/act)*100, 2))
  risk_Y$Se <- with(risk_Y, round(TP/(TP+FN), 2))   # Surveillance Se
  risk_Y$Sp <- with(risk_Y, round(TN/(TN+FP), 2))   # Surveillance Sp
  risk_Y$DF <- with(risk_Y, round(TP/inf, 2))       # Detection Fraction 
  risk_Y$cost <- with(risk_Y, round(samp/TP, 0))    # "cost" of detection = nr. of samples to get one true detection
  

  
 # stack the results of each simulation -----------
  ACT <- rbind(ACT, colSums(active))
  INF <- rbind(INF, colSums(infected))
  ACT_t <- rbind(ACT_t, colSums(active_t))
  INF_t <- rbind(INF_t, colSums(infected_t))
  ACT_r <- rbind(ACT_r, colSums(active_r))
  INF_r <- rbind(INF_r, colSums(infected_r))
  
  TR_Q1 <- rbind(TR_Q1, trad_Q1)
  TR_Q2 <- rbind(TR_Q2, trad_Q2)
  TR_Q3 <- rbind(TR_Q3, trad_Q3)
  TR_Q4 <- rbind(TR_Q4, trad_Q4)
  TR_Y <- rbind(TR_Y, trad_Y)
  RB_Q1 <- rbind(RB_Q1, risk_Q1)
  RB_Q2 <- rbind(RB_Q2, risk_Q2)
  RB_Q3 <- rbind(RB_Q3, risk_Q3)
  RB_Q4 <- rbind(RB_Q4, risk_Q4)
  RB_Y <- rbind(RB_Y, risk_Y)
  
  rm(list=setdiff(ls(), c("samp", "HR", "HR13", "dn", "adults", "infect", 
                          "ACT", "INF", "ACT_t", "ACT_r", "INF_t", "INF_r",
                          "TR_Q1", "TR_Q2", "TR_Q3", "TR_Q4", "TR_Y", 
                          "RB_Q1", "RB_Q2", "RB_Q3", "RB_Q4", "RB_Y")))
  gc()
  
}

save("samp", "HR", "HR13", "dn", "adults", "infect", 
     "ACT", "INF", "ACT_t", "ACT_r", "INF_t", "INF_r",
     "TR_Q1", "TR_Q2", "TR_Q3", "TR_Q4", "TR_Y", 
     "RB_Q1", "RB_Q2", "RB_Q3", "RB_Q4", "RB_Y", file = "SurvEval2013_M17.RData")

