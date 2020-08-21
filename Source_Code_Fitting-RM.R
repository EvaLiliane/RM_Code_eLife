# This code has not been evaluated. A dataset containing patients data is required forits usage. 
# Necessary features: Patients IDs, Sex, Age at ART initiation, Log-transformed RNA viral load, and a binary variable
# that  specify whtehr ornot the patient has suppressed viral load within 12-months of therapy.
# The outcome variable is scaled CD4.

library(saemix)


# ================================================== FUNCTIONS
# Selecting a subgroup of people
useSubset <- T
subselect <- function(addata,n){ 
  num.to.do <- n ## to save computing time when playing with analyses
  pids.to.run1 <- sample(unique(addata$patient), num.to.do)
  addata <- addata[addata$patient %in% pids.to.run1,] 
  return(addata)
}

#################################### Defining the Ratio Model
model1cpt<-function(psi,id,xidep) {
  #psi is the set of parameters
  #id  contains patient's ID (identification)
  #xidep contains the time of follow-up for each patient
  
  mtimart <-xidep[,1]
  k<-psi[id,1]
  q<-psi[id,2]
  z0<-psi[id,3]
  s <- psi[id,4]
  r <- psi[id,5]
  #k<-CL/V
  ypred <- (k*z0/q)*((1 + (q - 1)*exp(-s*mtimart))/(1 + (k - 1)*exp(-r*mtimart)))
  #ypred<-ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

##################################### Defining the asymptotic model
model1asym <- function(psi,id,xidep) {
  #dose<-xidep[,1]
  mtimart <- xidep[,1]
  Asym <- psi[id,1]
  c <- psi[id,2]
  R0 <- psi[id,3]
  
  ypred <-  Asym -(Asym - R0) * exp(-c*mtimart)
  #ypred <- (k*z0/q)*((1 + (q - 1)*exp(-s*mtimart))/(1 + (k - 1)*exp(-r*mtimart)))
  return(ypred)
}

###############################################################################################################
# =======================================================  DATA 

# Read your dataset. Make sure that the directory is correctly specified. 
# Assuming that you call the dataset inDatach, 
subdat <- subselect(inDatach,1000)

# Defining data to fit the model to.
saemix.data <-saemixData(name.data=subdat,header=TRUE,sep=" ",na=NA,
                         name.group=c("patient"),name.predictors=c("mtimart"),
                         name.response=c("scalcd4a"),name.covariates=c("sex","mthage", "logrna","suppress"),
                         units=list(x="months",y="-",covariates=c("-","months","-","log(copies/mL)","-")), name.X="mtimart")


# ============================================================================= EVALUATING ASYMPTOTIC MODEL

# Model specifications
Adsaemix.Asym00 <- saemixModel(model=model1asym,
                               description="Adults Ratio Model for Immune reconstitution. ",
                               psi0=matrix(c(0.6,0.05,0.2),ncol=3, byrow=TRUE,dimnames=list(NULL,
                                             c("Asym","c","R0"))),transform.par=c(1,1,1),
                               fixed.estim=c(1,1,1),
                               covariance.model=matrix(c(1,1,1,1,1,1,1,1,1),ncol=3,byrow=TRUE),
                               omega.init=matrix(c(1,1,1,1,1,1,1,1,1),ncol=3,byrow=TRUE),
                               error.model="proportional")

dateday <- Sys.Date()
dirname <- paste0("AdOutModelAsym00", dateday)
saemix.options <- list(seed=632545,save=TRUE,save.graphs=TRUE,
                       nb.chains=3, nbiter.saemix = c(300, 150),
                       directory = dirname )


#Fitting
Adsaemix.fitAsym00 <- saemix(Adsaemix.Asym00, saemix.data, saemix.options)
save(Adsaemix.fitAsym00, file = "simAdsaemixAsym00.RData")


# ============================================================================= EVALUATING RATIO MODEL

# Model specifications
Adsaemix.model6Full <- saemixModel(model=model1cpt,
                                   description="Adults Ratio Model for Immune reconstitution. ",
                                   psi0=matrix(c(3,1.1,0.2,0.05,0.5),ncol=5, byrow=TRUE,dimnames=list(NULL,
                                                c("k","q","z0", "s", "r"))),transform.par=c(1,1,1,1,1),
                                   covariate.model=matrix(c(1,1,1,1,0,1,1,1,1,0,1,0,1,0,1,1,0,0,0,1),ncol=5,byrow=TRUE),
                                   fixed.estim=c(1,1,1,1,1),
                                   covariance.model=matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                                           ncol=5,byrow=TRUE),
                                   omega.init=matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),ncol=5,byrow=TRUE),
                                   error.model="proportional")

dateday <- Sys.Date()
dirname <- paste0("AdOutModel6Full", dateday)
saemix.options <- list(seed=632545,save=TRUE,save.graphs=TRUE,
                       nb.chains=3, nbiter.saemix = c(300, 150),
                       directory = dirname )


#Fitting
Adsaemix.fit6Full <- saemix(Adsaemix.model6Full, saemix.data, saemix.options)
save(Adsaemix.fit6Full, file = "sim.AdsaemixModel6Full.RData")
