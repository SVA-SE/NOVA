seed.u0 <- function (ini_bhp=0.05, seed=0, weighted=T, whp_c=0.082, whp_y=0.007, whp_a=0.024) 
{

   # load geo info on herds
   data(geo)
  
   # Add ID (obs! U0 is already sorted by id: 0 - 37220)
      u0$id <- 0:37220

   # Identify dairy herds
      dairy <- u0$id %in% geo$ID[geo$dairy]
      
   # Deterimine herd size and dairy herds with at least one animal at the
   # first day of the simulation 
      herd_size <- u0$Sc + u0$Sy + u0$Sa
      index <- herd_size > 0  # to sample only from active herds
   
   # list all NUTS3 codes and make an extra one for Öland
      nuts3o <- paste(geo$NUTS3,"-", as.numeric(geo$oland), sep="")

   # Proportion of positive farms in BTM national screening 2013 by län
   # add also small p for those positive in 2007 (same order as nuts3)
      ns <- data.frame(nuts3o = sort(unique(nuts3o)), 
                       prev = c(0, 0.012, 0.007, 0.007, 0, 
                                0, 0, 0.005, 0.001, 0.149, 
                                0, 0, 0.012, 0, 0.002,
                                0, 0, 0, 0, 0, 0.004, 0.009)  )
               
    # create vector of weights for each farm by retrieveing the "prev"
    # from BMS for that län(+Öland) 

      weight<- ns[match(nuts3o, ns$nuts3o), 2][match(u0$id, geo$ID)] 
      
    # add additioal weight to increase sampling prob in risk area
    # double the weight in skane risk area
      weight[geo$hr.skane] <- weight[geo$hr.skane]+weight[geo$hr.skane] 
      
    # match prevelance with nuts3o and reorder it according to u0 order

      # fix seed number to get reproducible results (if asked)
      if (seed>0) {     
       set.seed(seed)
        }
      # Start with 1% infected DAIRY herds to meet Estelle's survey in 2013 --> try 5
      # (higher probability goes to larger herds and herds located in Öland
      index[!geo$dairy] <- FALSE  # to sample from dairy herds only
      
      if (weighted==T) {
         i <- sample(which(index == TRUE),
                     round(sum(index==TRUE)*ini_bhp),
                     prob = herd_size[index==T]*weight[index==T]) 
          }
      else {
         i <- sample(which(index == TRUE),
                     round(sum(index==TRUE)*ini_bhp)) }
      
      
      ## Ågren 2013 - prevalence of salmonella in 15 bulk milk positive dairy herds
      ## Animals age    %Culture+ (mean, min-max)
      ## Calves  0-6m   8.2 % (0 - 33.3%)
      ## Youngs  6-48   0.7 % (0 - 3.1%)
      ## Adults  >48m   2.4 % (0 - 17.2%)
      
      ## Age-specific within herd prevalence (infected herds)
#       whp_c <- 0.082 
#       whp_y <- 0.007 
#       whp_a <- 0.024 
      
      ## Distribute infected animals in infected herds
      u0$Ic[i] <- (u0$Sc[i] * whp_c)
      u0$Iy[i] <- (u0$Sy[i] * whp_y)
      u0$Ia[i] <- (u0$Sa[i] * whp_a)
      
      ## Force to have at least one infected animal in case rounding leads
      ## to 0s in all 3 age groups. This is done by setting the larger
      ## value amongst I_calves, I_young and I_adults  to be 1
      infected <- round(u0$Ic[i]) + round(u0$Iy[i]) + round(u0$Ia[i])
      j <- which(infected==0)
      maxinf <- pmax(u0$Ic[i][j], u0$Iy[i][j], u0$Ia[i][j])
      u0$Ic[i][j][which(u0$Ic[i][j] == maxinf)] <- 1
      u0$Iy[i][j][which(u0$Iy[i][j] == maxinf)] <- 1
      u0$Ia[i][j][which(u0$Ia[i][j] == maxinf)] <- 1
      
      u0$Ic[i]<- round(u0$Ic[i])
      u0$Iy[i] <- round(u0$Iy[i])
      u0$Ia[i] <- round(u0$Ia[i])
      
      u0$Sc[i] <- u0$Sc[i] - u0$Ic[i]
      u0$Sy[i] <- u0$Sy[i] - u0$Iy[i]
      u0$Sa[i] <- u0$Sa[i] - u0$Ia[i]
   
   
 # drop ID
   u0$id <- NULL
      
   return(u0)
   
}


