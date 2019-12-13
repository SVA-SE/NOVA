library(SimInf)
library(SICR)
library(sp)


#load the empty initialized model
load("New_S.Dublin_SICR_ini.Rdata")


#################
# RUN the model #
#################
a <- SimInf::run(model)


# get ID of dairy herds & dairy herds in Öland
dn <- geo$row[geo$dairy] # node number of dairy herds

########################
# EXTRACT ALL THE DATA #
########################

for (i in 1:100) {
  print(paste("running simulation", i))
  a <- SimInf::run(model)
  all <- trajectory(a, node=dn)
  all$time <- as.Date(all$time, origin="2004-12-31")
  saveRDS(all, file=sprintf("SimulData/Model_17/trj_%04i.RDS",i))
}



