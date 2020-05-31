## ----setup, include = FALSE, cache = FALSE-------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = TRUE
)

## ---- echo = T-----------------------------------------------------------
library(Disequilibrium)
set.seed(1775) # Semper Fi
data(fjdata)

sfjdata = as.data.frame(scale(fjdata))
mod1 = DE(HS ~ T + HL1 + RML2 | T + DK16L1 + DH13L2 + RML1, data = sfjdata, control = list(MaskRho = 0))

## ---- echo = T-----------------------------------------------------------

# Parameter estimates
mod1$par
# Likelihood hessian
round(mod1$hessian,2)

## ---- echo = T-----------------------------------------------------------
# rho fixed to 0 and data not standardized
mod2 = DE(HS ~ T + HL1 + RML2 | T + DK16L1 + DH13L2 + RML1, data = fjdata, control = list(MaskRho = 0))

# rho free and data not standardized (converges but hessian is not invertible)
mod3 = DE(HS ~ T + HL1 + RML2 | T + DK16L1 + DH13L2 + RML1, data = fjdata, control = list(MaskRho = FALSE))
solve(mod3$hessian)

# rho free and data standardized
mod4 = DE(HS ~ T + HL1 + RML2 | T + DK16L1 + DH13L2 + RML1, data = sfjdata, control = list(MaskRho = FALSE))


## ---- echo = T-----------------------------------------------------------
NStartingValues = 10
k = length(mod1$par)
startingvalues = matrix(rnorm(NStartingValues * k,0,.5),NStartingValues,k)
startingvalues[,k - c(0,1)] = abs(startingvalues[,k - c(0,1)]) + 1
modelstorage = vector("list",NStartingValues + 1)
modelstorage[[1]] = mod1

for(i in 1:NStartingValues){
  modelstorage[[i+1]] = DE(HS ~ T + HL1 + RML2 | T + DK16L1 + DH13L2 + RML1, data = sfjdata, par = startingvalues[i,], control = list(MaskRho = 0))
}

llhoodvalues = sapply(modelstorage,function(x){x$value})
convergedIndex = sapply(modelstorage,function(x){x$convergence}) == 0

finalmod = modelstorage[[which(llhoodvalues == min(llhoodvalues[convergedIndex]))]]


## ---- echo =T------------------------------------------------------------
sfinalmod = summary(finalmod)

round(sfinalmod,3)

## ---- echo =T------------------------------------------------------------
pfinalmod = predict(finalmod)
round(head(pfinalmod),3)

## ---- echo = T-----------------------------------------------------------
sfinalmodrob = summary(finalmod, robust = T)

round(sfinalmodrob,3)

## ---- echo = T, error = F------------------------------------------------
#Simulate some fake clusters
clusters = sample(LETTERS,nrow(sfjdata),replace = T)
#Cluster standard errors in the unrestricted covariance space
vcov1 = sandwich::vcovCL(finalmod, cluster = clusters)
round(vcov1,3)

#Cluster standard errors in the natural positive definite space
newpar = GetDeltaMethodParameters(mu = finalmod$par,covmat = vcov1, MaskRho = finalmod$MaskRho)
vcov2 = newpar$covmat
round(vcov2,3)

