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

ACT_r2 <- c()
INF_r2 <- c()
ACT_h <- c()
INF_h <- c()

RB2_Q1 <- c()
RB2_Q2 <- c()
RB2_Q3 <- c()
RB2_Q4 <- c()
RB2_Y <- c()

HR_Q1 <- c()
HR_Q2 <- c()
HR_Q3 <- c()
HR_Q4 <- c()
HR_Y <- c()


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
  dt13$inf <- rowSums(dt13[,infect])
  dt13$adu <- rowSums(dt13[,adults])
  dt13$det <- ifelse(dt13$adu==0, 0, dt13$Ra / dt13$adu)
  
  dt13$ah <- dt13$tot > 0
  dt13$ih <- dt13$inf > 0
  dt13$dh <- dt13$det >= 0.15
  
  
  # extract status matrix
  ########################
  active <- as.data.frame(do.call("cbind", tapply(dt13$ah, dt13$week, function(x) { x == TRUE})))
  infected <- as.data.frame(do.call("cbind", tapply(dt13$ih, dt13$week, function(x) { x == TRUE})))
  detected <- as.data.frame(do.call("cbind", tapply(dt13$dh, dt13$week, function(x) { x == TRUE})))
  
 
  ##################################
  # RISK-BASED SURVEILLANCE --> Q4 #
  ##################################
  
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
  # -->  only HR under surveillance for this quarter
  active_r[HR13 & infected_r[,33] & detected_r[,33],37:52] <- FALSE
  infected_r[HR13 & infected_r[,33] & detected_r[,33],37:52] <- FALSE
  detected_r[HR13 & infected_r[,33] & detected_r[,33],37:52] <- FALSE

  # update after 4th quarter (week 46)
  # --> all herds under surveillance
  active_r[infected_r[,46] & detected_r[,46],50:52] <- FALSE
  infected_r[infected_r[,46] & detected_r[,46],50:52] <- FALSE
  detected_r[infected_r[,46] & detected_r[,46],50:52] <- FALSE


  # Summary per quarter --------------------
  for (k in samp[-4]) {
    
    nam2 <- paste("risk_Q", which(samp==k), sep = "")
    
    assign(nam2, data.frame(act = sum(active_r[,k]),   # active herds
                            inf = sum(infected_r[,k]), # infected herds
                            samp = sum(active_r[HR13,k]),   # sampled herds  OBS Q4 is different
                            det = sum(detected_r[HR13,k]), # detected herds
                            TP = sum(infected_r[HR13,k] & detected_r[HR13,k]), # True Positive herds
                            FP = sum(!infected_r[HR13,k] & detected_r[HR13,k]), # False Positive herds
                            FN = sum(infected_r[HR13,k] & !detected_r[HR13,k]), # False Negative herds
                            TN = sum(active_r[HR13, k] & !infected_r[HR13,k] & !detected_r[HR13,k]), # True Negative herds
                            Prev = round((sum(infected_r[,k]) / sum(active_r[,k]))*100, 3),  # true prevalence in the population
                            TPr = round((sum(infected_r[HR13,k]) / sum(active_r[HR13,k]))*100, 3),  # true prevalence (in the sample)
                            APr = round((sum(detected_r[HR13,k]) / sum(active_r[HR13,k]))*100, 3))) # apparent prevalence (in the sample)
    
  }
  
  # OBS! in Q4 all herds are sampled
  
  risk_Q4 <- data.frame(act = sum(active_r[,k]),   # active herds
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
  

  
  #######################
  # ONLY HIGH RISK AREA #
  #######################
  
  # ONly herds in HR area are sampled every quarter of the year
  
  # update status 4 weeks after detection --> herd no longer under surveillance (ah=FALSE, ih=FALSE)
  ###########################################################################################
  active_h <- active[which(rowSums(active)>0),]
  infected_h <- infected[which(rowSums(active)>0),]
  detected_h <- detected[which(rowSums(active)>0),]
  
  # update after 1st quarter (week 7) --------------
  # remove TP herds from surveillance
  active_h[HR13 & infected_h[,7] & detected_h[,7],11:52] <- FALSE
  infected_h[HR13 & infected_h[,7] & detected_h[,7],11:52] <- FALSE
  detected_h[HR13 & infected_h[,7] & detected_h[,7],11:52] <- FALSE

  
  # update after 2nd quarter (week 20)--------------------
  active_h[HR13 & infected_h[,20] & detected_h[,20],24:52] <- FALSE
  infected_h[HR13 & infected_h[,20] & detected_h[,20],24:52] <- FALSE
  detected_h[HR13 & infected_h[,20] & detected_h[,20],24:52] <- FALSE
 
  # update after 3rd quarter (week 33)
  active_h[HR13 & infected_h[,33] & detected_h[,33],37:52] <- FALSE
  infected_h[HR13 & infected_h[,33] & detected_h[,33],37:52] <- FALSE
  detected_h[HR13 & infected_h[,33] & detected_h[,33],37:52] <- FALSE

  # update after 4th quarter (week 46)
  active_h[HR13 & infected_h[,46] & detected_h[,46],50:52] <- FALSE
  infected_h[HR13 & infected_h[,46] & detected_h[,46],50:52] <- FALSE
  detected_h[HR13 & infected_h[,46] & detected_h[,46],50:52] <- FALSE

  
  # Summary per quarter --------------------
  for (j in samp) {
    
    nam <- paste("high_Q", which(samp==j), sep = "")
    
    assign(nam, data.frame(act = sum(active_h[,j]),   # active herds
                           inf = sum(infected_h[,j]), # infected herds
                           samp = sum(active_h[HR13,j]),   # sampled herds
                           det = sum(detected_h[HR13,j]), # detected herds
                           TP = sum(infected_h[HR13,j] & detected_h[HR13,j]), # True Positive herds
                           FP = sum(!infected_h[HR13,j] & detected_h[HR13,j]), # False Positive herds
                           FN = sum(infected_h[HR13,j] & !detected_h[HR13,j]), # False Negative herds
                           TN = sum(active_h[HR13, j] & !infected_h[HR13,j] & !detected_h[HR13,j]), # True Negative herds
                           Prev = round((sum(infected_h[,j]) / sum(active_h[,j]))*100, 3),  # true prevalence in the population
                           TPr = round((sum(infected_h[HR13,j]) / sum(active_h[HR13,j]))*100, 3),  # true prevalence in the sample
                           APr = round((sum(detected_h[HR13,j]) / sum(active_h[HR13,j]))*100, 3))) # apparent prevalence in the sample
    # apparent prevalence
    
  }
  
  
  # Summary per year -------------------------
  high_Y <- data.frame(act = sum(rowSums(active_h)>0),
                       inf = sum(rowSums(infected_h)>0),  # OBS! herds infected for only short time??
                       # outb = sum(rowSums(infected_h)>1),
                       samp = sum(high_Q1$samp, high_Q2$samp, high_Q3$samp, high_Q4$samp),
                       det = sum(high_Q1$det, high_Q2$det, high_Q3$det, high_Q4$det),
                       TP = sum(high_Q1$TP, high_Q2$TP, high_Q3$TP, high_Q4$TP),
                       FP = sum(high_Q1$FP, high_Q2$TP, high_Q3$TP, high_Q4$TP),
                       FN = sum(high_Q1$FN, high_Q2$FN, high_Q3$FN, high_Q4$FN),
                       TN = sum(high_Q1$TN, high_Q2$TN, high_Q3$TN, high_Q4$TN))
  
  high_Y$prev <- with(high_Y, round((inf/act)*100, 2))
  high_Y$Se <- with(high_Y, round(TP/(TP+FN), 2))   # Surveillance Se
  high_Y$Sp <- with(high_Y, round(TN/(TN+FP), 2))   # Surveillance Sp
  high_Y$DF <- with(high_Y, round(TP/inf, 2))       # Detection Fraction 
  high_Y$cost <- with(high_Y, round(samp/TP, 0))    # "cost" of detection = nr. of samples to get one true detection
  
  
  
  
  
  
  
  
 # stack the results of each simulation -----------
  
  ACT_r2 <- rbind(ACT_r2, colSums(active_r))
  INF_r2 <- rbind(INF_r2, colSums(infected_r))
  ACT_h <- rbind(ACT_h, colSums(active_h))
  INF_h <- rbind(INF_h, colSums(infected_h))
  
  
  RB2_Q1 <- rbind(RB2_Q1, risk_Q1)
  RB2_Q2 <- rbind(RB2_Q2, risk_Q2)
  RB2_Q3 <- rbind(RB2_Q3, risk_Q3)
  RB2_Q4 <- rbind(RB2_Q4, risk_Q4)
  RB2_Y <- rbind(RB2_Y, risk_Y)
  HR_Q1 <- rbind(HR_Q1, high_Q1)
  HR_Q2 <- rbind(HR_Q2, high_Q2)
  HR_Q3 <- rbind(HR_Q3, high_Q3)
  HR_Q4 <- rbind(HR_Q4, high_Q4)
  HR_Y <- rbind(HR_Y, high_Y)
  
  rm(list=setdiff(ls(), c("samp", "HR", "HR13", "dn", "adults", "infect", 
                          "ACT_r2", "ACT_h", "INF_r2", "INF_h",
                          "HR_Q1", "HR_Q2", "HR_Q3", "HR_Q4", "HR_Y", 
                          "RB2_Q1", "RB2_Q2", "RB2_Q3", "RB2_Q4", "RB2_Y")))
  gc()
  
}

save("samp", "HR", "HR13", "dn", "adults", "infect", 
     "ACT_r2", "ACT_h", "INF_r2", "INF_h",
     "HR_Q1", "HR_Q2", "HR_Q3", "HR_Q4", "HR_Y", 
     "RB2_Q1", "RB2_Q2", "RB2_Q3", "RB2_Q4", "RB2_Y", file = "SurvEval2013more_M17.RData")

