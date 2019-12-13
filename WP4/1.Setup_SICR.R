

library(Matrix)

install.packages("SICR", repos=NULL, type="source")
library(SICR)

###########################
#    MODEL PREPARATION    #
###########################

# Load events, distance, season 
# -----------------------------
data(events)
data(distance)
data(season)
data(u0)

time_offset <- as.integer(strftime(as.Date('2005-07-01'), "%j"))  # = 182
events$time <- events$time + time_offset
tspan <- seq(from=min(events$time)-1, to=max(events$time), by=7)



# set up parameters (as in the first manuscript) 
# ----------------------------------------------

model <- SICR(u0 = u0, 
              tspan = seq(182, max(events$time), by = 7), 
              events = events,
              phi       = rep(0, nrow(u0)),
              upsilon_c = 0.0042,   # uptake rate
              upsilon_y = 0.0042/2,
              upsilon_a = 0.0042/3,     
              gamma_c   = 1/17,      # recovery rate
              gamma_y   = 1/10,
              gamma_a   = 1/10,
              chi_c   = 0.025,       # carrier rate
              chi_y   = 0.03,
              chi_a   = 0.03,
              rho_c   = 1/300,       # carrier recovery rate
              rho_y   = 1/300,
              rho_a   = 1/300,
              psi_c = 1/100,         # immunity loss rate
              psi_y = 1/150,
              psi_a = 1/250,
              alpha_c   = 10,        # bacterial shedding rate
              alpha_y   = 5,
              alpha_a   = 1, 
              sigma     = 1/100,     # infectivity of carriers
              beta_t1   = 1.12e-1,   # seasonal decay rate
              beta_t2   = 0.94e-1,
              beta_t3   = 1.18e-1,
              beta_t4   = 1.3e-1,
              end_t1    = season$end_t1,
              end_t2    = season$end_t2,
              end_t3    = season$end_t3,
              end_t4    = season$end_t4,
              distance  = d,
              coupling  = rep(0.1, nrow(u0))) 

# Change spatial coupling for Öland (Öland effect)
model@ldata[1,geo$oland] <- 0.28

model # this can be saved and run subsequently

# save model for future play around
# save(model, file="S.Dublin_SICR.Rdata")


# seed the initial values  (5% infected dairy herds as in manuscript)
# ----------------------------------------------------------------------

u0.s <- seed.u0(ini_bhp=0.05, seed=123, whp_c=0.082, whp_y=0.007, whp_a=0.024) 
# reorder classes to match model order
u0.s <- u0.s[, c("Sc", "Ic", "Cc", "Rc", 
                 "Sy", "Iy", "Cy", "Ry", 
                 "Sa", "Ia", "Ca", "Ra")]

# incorporate inits
inits <- t(as.matrix.data.frame(u0.s))
mode(inits) <- "integer"
model@u0 <- inits

# add map
library(svar)
data(NUTS_10M)
sweden <- NUTS_10M[NUTS_10M@data$NUTS_ID=="SE",] # keep only boundaies of sweden

# drop unused data
rm(events, d, season)
rm(inits, u0)
gc()

# save SEEDED model (+ geo information)
# save(model, geo, sweden, NUTS_10M, file="New_S.Dublin_SICR_ini.Rdata")


# run the model (as it is) and inspect output
# -------------------------------------------

a <- SimInf::run(model)
summary(a)

time <- get.time(a)

# show nr animals in each node at each time point
head(trajectory(a)) 

# show prevalece of infected animals 
head(prevalence(a, Ic+Iy+Ia~., type="pop"))

# show proportion of positive herds (ie. having at least one I animal)
head(prevalence(a, Ic+Iy+Ia~., type="nop"))

# show within herd prevalence (obs! all herds not only positive)
head(prevalence(a, Ic+Iy+Ia~., type="wnp"))



# plot proportion positive dairy herds
plot(prevalence~as.Date(time, origin = "2004-12-31"),prevalence(a, Ic+Iy+Ia~., type="nop", node=which(geo$dairy)), 
     type="l", main="positive dairy herds", ylim=c(0, 0.05))
# plot proportion positive non-dairy herds
plot(prevalence~as.Date(time, origin = "2004-12-31"),prevalence(a, Ic+Iy+Ia~., type="nop", node=which(!geo$dairy)), 
     type="l", main="positive non-dairy herds", ylim=c(0, 0.02))



# plot proportion positive non-dairy herds


plot(Ic+Iy+Ia~time,trajectory(a,node=5), type="l")

bb<-(prevalence(a, Ic+Iy+Ia~., type="wnp", node=5:7))

