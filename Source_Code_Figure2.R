rm(list=ls())   # Clear all variables and functions

library(deSolve)

######## Parameters are :
# k  -- the environment capacity
# r -- growth rate

####################################### THE MODEL

logmod <- function(times,x,parameters){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(x),parameters),{
  dx <- r * x * (1 - x/k)
    # Note: Population size is constant, so don't need to specify dRdt
    list(c(dx))
  })
}

# ------------------------------ THE GRAPH

jpeg("loggrowthplot.jpeg", width = 600, height = 400)
states <- c(x = 1)
parameters <- c(r = 1, k = 50)
time <- seq(0,10,by = 0.05)
resx <-  ode(y = states, times = time, func = logmod, parms = parameters)
plot(resx[,1], resx[,2], xlab="time-t", ylab = "CD4 count", type = "l", main = "Logistic model", 
     lwd = 2, cex=1.5, cex.lab = 1.5)
abline(a = 50, b =0, lwd = 2, lty = 2, col = "red")
text(0,47,labels = "k", col = "red")
abline(a = 25, b =0, lwd = 2, lty = 2, col = "blue")
points(3.9,25, pch = 16, col = "blue")
text(0,22,labels = "k/2", col = "blue")
dev.off()


#######################################################################
# SAVING DATA SOURCES

# folder <- 'C:/Users/Documents/Github/RM_Code_eLife/'
# write.table(resx,file=paste0(folder,'Source_Data_Figure2.txt'),sep='\t',row.names=FALSE,quote=FALSE)




