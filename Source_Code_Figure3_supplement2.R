model1cpt<-function(psi,mtimart) {
  #dose<-xidep[,1]
  k<-psi[1]
  q<-psi[2]
  z0<-psi[3]
  s <- psi[4]
  r <- psi[5]
  #k<-CL/V
  ypred <- (k*z0/q)*((1 + (q - 1)*exp(-s*mtimart))/(1 + (k - 1)*exp(-r*mtimart)))
  #ypred<-ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

# #######################################################################
# # READING DATA SOURCES

# Change the below tree lines to make fig 3 - supplement 2
pop_par <- c(3.4,0.9,0.18,0.017,0.35)
ylimID <- c(rep(0.8,3),1,1.7,1.2,1,2,rep(0.8,4))
ind_par <- read.delim("~/GitHub/RM_Code_eLife/Figure 3-figure supplement 2-source data 1.txt")
len <- length(ind_par$patient)
tofu <- 1:60

par(mfrow=c(3,4))  
for(i in 1:len){
  plot(tofu,model1cpt(as.numeric(ind_par[i,2:6]),tofu),type="l",ylim=c(0,ylimID[i]),ylab="Scaled CD4", xlab = "TOFU (month)")
  lines(tofu,model1cpt(pop_par,tofu),type="l",col="red")
}