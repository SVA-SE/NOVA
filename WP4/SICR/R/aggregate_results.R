# function list:
# - aggregate.result()
# - herd.inspect()
# - get.time()

# from a model U matrix extract number of animals in each compartment at each timestep
aggregate.result <- function(model)
{
  sc <- SimInf:::extract_U(model, "Sc")       # Susceptible Calves (ALL Herds)
  ic <- SimInf:::extract_U(model, "Ic")       # Infected Calves
  cc <- SimInf:::extract_U(model, "Cc")       # Carrier Calves
  rc <- SimInf:::extract_U(model, "Rc")       # Recovered Calves
  sy <- SimInf:::extract_U(model, "Sy")       # Susceptible Youngs
  iy <- SimInf:::extract_U(model, "Iy")       # Infected Youngs
  cy <- SimInf:::extract_U(model, "Cy")       # Carrier Youngs
  ry <- SimInf:::extract_U(model, "Ry")       # Recovered Youngs
  sa <- SimInf:::extract_U(model, "Sa")       # Susceptible Adults
  ia <- SimInf:::extract_U(model, "Ca")       # Infected Adults
  ca <- SimInf:::extract_U(model, "Ia")       # Carrier Adults
  ra <- SimInf:::extract_U(model, "Ra")       # Recovered Adults
  
  sus <- sc + sy + sa                         # Total number of SUSCEPTIBLE in the herd
  inf <- ic + iy + ia                         # Total number of INFECTED in the herd
  car <- cc + cy + ca                         # Total number of CARRIER in the herd
  rec <- rc + ry + ra                         # Total number of RECOVERED in the herd
  
  Nc <- sc + ic + cc + rc                     # Total number of CALVES in the herd
  Ny <- sy + iy + cy + ry                     # Total number of YOUNGS in the herd
  Na <- sa + ia + ca + ra                     # Total number of ADULTS in the herd
  
  N <- Nc + Ny + Na                           # total number of ANIMALS in the herd
  
  data <- list(sc=sc, ic=ic, cc=cc, rc=rc,
               sy=sy, iy=iy, cy=cy, ry=ry,
               sa=sa, ia=ia, ca=ca, ra=ra,
               sus=sus, inf=inf, car=car, rec=rec,
               Nc=Nc, Ny=Ny, Na=Na, N=N)
  
  return(data)   
}


# take the output of a Siminf model check whether each herd containsat least one 
# infected animal, any animal or carrier animal, per time step

herd.inspect <- function(model, what=c("infected", "any", "carrier"))
{ 
  what <- match.arg(what)
  
  if (what=="infected") {
      inf <- SimInf:::extract_U(model, "Ic") +
             SimInf:::extract_U(model, "Iy") +
             SimInf:::extract_U(model, "Ia")
      herd <- inf>0
    }  
    
  if (what=="any")   {
      all <- SimInf:::extract_U(model, "Sc") +
             SimInf:::extract_U(model, "Sy") +
             SimInf:::extract_U(model, "Sa") +
             SimInf:::extract_U(model, "Ic") +
             SimInf:::extract_U(model, "Iy") +
             SimInf:::extract_U(model, "Ia") +
             SimInf:::extract_U(model, "Cc") +
             SimInf:::extract_U(model, "Cy") +
             SimInf:::extract_U(model, "Ca") +
             SimInf:::extract_U(model, "Rc") +
             SimInf:::extract_U(model, "Ry") +
             SimInf:::extract_U(model, "Ra")
      herd <- all>0 
    }  
    
  if (what=="carrier")   {
      car <- SimInf:::extract_U(model, "Cc") +
             SimInf:::extract_U(model, "Cy") +
             SimInf:::extract_U(model, "Ca")
      herd <- car>0
    } 
 
   return(herd)
}


# get vector of time span from model extracted data using previous function
get.time <- function(model, origin="2005-07-01")
{
  time_offset <- as.integer(strftime(as.Date(origin), "%j"))
  time <- as.Date(model@tspan - time_offset, origin=origin)
  return(time)
}
