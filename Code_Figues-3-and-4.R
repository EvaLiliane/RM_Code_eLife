RatioModel <- function(tim, K,Q,r,s,z0){
  (K*z0/Q)*((1 + (Q - 1)*exp(-s*tim))/(1 + (K - 1)*exp(-r*tim)))
}

RatioModel.v2 <- function(tim, K,Q,r,s,z0){
  ((K*z0*(Q-1))/(Q*(K-z0)))*exp((r - s)*tim)*((1 + exp(s*tim)/(Q - 1))/(1 + exp(r*tim)*z0 /(K - z0)))
} 

AsymModel <- function(tim,Asy,R0,c){
  Asy - (Asy - R0)*exp(-c*tim)
} 


agemth <- seq(0,12*17,1)  # simulated TOFU

#========================================== Fixed parameter - No covariate for RM
# Children parameters
pars <- c(K = 3.4, Q = 0.9, r =  0.35, s = 0.017, z0 = 0.18)
finpars <- c(K = 1.03, Q = 0.68, r =  0.27, s = 0.011, z0 = 0.75)
y <- RatioModel(agemth,pars[1], pars[2], pars[3], pars[4], pars[5])
finy <- RatioModel(agemth,finpars[1], finpars[2], finpars[3], finpars[4], finpars[5])

param <- c(Asym = 0.77, R0 = 0.18, c = 0.12)
z <- AsymModel(agemth,param[1], param[2], param[3])

# Adults parameters
pars.a <- c(K = 2.54, Q = 0.38, r =  1.23, s = 0.01, z0 = 0.13)
finpars.a <- c(K = 1.75, Q = 0.498, r =  0.49, s = 0.022, z0 = 0.20)
y.a <- RatioModel(agemth,pars.a[1], pars.a[2], pars.a[3], pars.a[4], pars.a[5])
finy.a <- RatioModel(agemth,finpars.a[1], finpars.a[2], finpars.a[3], finpars.a[4], finpars.a[5])

param.a <- c(Asym = 0.66, R0 = 0.06, c = 0.14)
z.a <- AsymModel(agemth,param.a[1], param.a[2], param.a[3])

pal(tol7qualitative)

par(mfrow=c(1,1))
plot(agemth,z, type = "l", col = "blue", xlab = "TOFU  (months)", ylab = "Scaled CD4", ylim = c(0,1.2), lwd = 2,
     main = "Population simulation of scaled CD4 for children")
lines(agemth,y, col = "green", type = "l", lwd = 2)
lines(agemth,finy, col = "red", type = "l", lwd = 2)
legend("bottomright", legend = c("AM", "RM", "Adj RM"), col = c("blue","green", "red"), lty = 1, lwd = 2)


plot(agemth,z.a, type = "l", col = "blue", xlab = "TOFU  (months)", ylab = "Scaled CD4", ylim = c(0,1) , lwd = 2,
     main = "Data scaled by a constant")
lines(agemth,y.a, col = "green", type = "l", lwd = 2)
lines(agemth,finy.a, col = "red", type = "l", lwd = 2)
legend("bottomright", legend = c("AM", "RM", "Adj RM"), col = c("blue","green", "red"), lty = 1, lwd = 2)

#========================================== Rescaled 
# Adults parameters

pars.a2 <- c(K = 2.5, Q = 0.39, r =  1.24, s = 0.01, z0 = 0.13)
finpars.a2 <- c(K = 1.75, Q = 0.5, r =  0.48, s = 0.02, z0 = 0.25)
y.a2 <- RatioModel(agemth,pars.a[1], pars.a[2], pars.a[3], pars.a[4], pars.a[5])
finy.a2 <- RatioModel(agemth,finpars.a[1], finpars.a[2], finpars.a[3], finpars.a[4], finpars.a[5])

param.a2 <- c(Asym = 0.63, R0 = 0.14, c = 0.06)
z.a2 <- AsymModel(agemth,param.a[1], param.a[2], param.a[3])


plot(agemth,z.a2, type = "l", col = "blue", lwd = 2, xlab = "TOFU  (months)", ylab = "Scaled CD4", ylim = c(0,1),
     main = "Population simulation of scaled CD4 for adults")
lines(agemth,y.a2, col = "green", type = "l", lty = 1, lwd = 2)
lines(agemth,finy.a2, col = "red", type = "l", lty = 1, lwd = 2)
legend("bottomright", legend = c("AM", "RM", "Adj RM"), col = c("blue","green", "red"), lty = 1, lwd = 2)


#######################################################################
# SAVING DATA SOURCES

# #save predicted trajectories
# d.c <- data.frame(agemth,z,y,finy)
# d.a <- data.frame(agemth, z.a2, y.a2, finy.a2)
# colnames(d.c) <- c("TOFU", "AM", "RM", "AdjRM")
# colnames(d.a) <- c("TOFU", "AM", "RM", "AdjRM")
# 
# folder <- 'C:/Users/Documents/Github/RM_Code_eLife/'
# write.table(d.c,file=paste0(folder,'Source_Data_Figure3.txt'),sep='\t',row.names=FALSE,quote=FALSE)
# write.table(d.c,file=paste0(folder,'Source_Data_Figure4.txt'),sep='\t',row.names=FALSE,quote=FALSE)

