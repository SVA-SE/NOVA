# function list:
# - length.true.check()
# - status.check()


# the function checks wether the length of consecutive TRUE values 
# in a logical vector x is greater than "span". It returns T for 
# every time the length was >= span and the value was T
# and with F otherwise 

# --> OUTBREAK application: it owerwirtes T with F when the 
#     length was too short to be considered an outbreak

length.true.check <- function(x, span=12, compare=c(">=", "<=", "=" ))
{  
   case <- match.arg(compare)
   tmp <- rle(x)
   
   if (case==">=") x <- rep(tmp$lengths>=span & tmp$values==T, times=tmp$lengths)
   if (case=="=") x <- rep(tmp$lengths==span & tmp$values==T, times=tmp$lengths)
   if (case=="<=") x <- rep(tmp$lengths<=span & tmp$values==T, times=tmp$lengths)
   
   return(x)
}


# Check any sequence greater than span, regardless if it is true or false.
# If it is less than span, then assign the same T/F status it had before

status.check <- function(x, span=52)
{  
   tmp <- rle(x)
   l <- length(tmp$lengths)
   y <- rep(tmp$values[1], times=tmp$lengths[1])
   
   if (l > 1) 
   {
      for (i in 2:length(tmp$lengths)) 
      {
         if(tmp$lengths[i]<span)   y <- c(y, rep(tail(y,1), times=tmp$lengths[i]))
         if(tmp$lengths[i]>=span)  y <- c(y, rep(tmp$values[i], times=tmp$lengths[i]))
      }
   }
   
   return(y)
}