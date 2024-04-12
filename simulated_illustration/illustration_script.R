library(here)

# FUNCTIONS

# Package Installation
installer <- function(vector) {
  for(i in vector) {
    if(!require(i, character.only = T)) install.packages(i)
    library(i, character.only = T)
  }
}

# Condensed Formula Object
form <- function(covariates, outcome = NULL, quotes = T) {
  equation <- as.formula(paste0(outcome, " ~ NULL"))
  for(i in covariates) {
    if(substr(i, 1, 1) %in% c(LETTERS, letters) & substr(i, nchar(i), nchar(i)) %in%  c(LETTERS, letters) |
       quotes == T) {
      equation <- update(equation, as.formula(paste0(" ~ . + `", i,"`"))) 
    }
    else{
      equation <- update(equation, as.formula(paste0(" ~ . + ", i))) 
    }
  }  
  equation
}

# Monotone Regression Function
monofunc <- function(covs,
                     spec.use = spec.short,         # Choice of Number of Iterations/Burn-In
                     parm.use = parm.nomono         # Choice of Model Parametrization
) {
  (begin <- Sys.time())
  results <- monosurv(niter=spec.use$niter, burnin=spec.use$burnin, adapt=spec.use$adapt,
                      refresh = 10, thin = 5,
                      rhoa=0.1, rhob=0.1, deltai=0.5, drange=10.0,
                      predict = c(rep(1, covs$nobs), rep(1, nrow(covs$test))),   # Predict on Both Training and Test Sets
                      include = c(rep(1, covs$nobs), rep(0, nrow(covs$test))),
                      casestatus = c(covs$y, rep(0, nrow(covs$test))),
                      axes = rbind(covs$train, covs$test)[,covs$add.id],
                      covariates = rbind(covs$train, covs$test)[,covs$covs.id],
                      birthdeath = parm.use$birthdeath,
                      settozero = parm.use$settozero,
                      package = parm.use$package
  )
  (tdiff <- Sys.time() - begin)
  results$time <- tdiff
  return(results)
}

# BASE PACKAGES
installer(c("tidyverse", "here", "monoreg"))



# DATA PREPARATION
load(here("simulated_illustration", "data", "sim_illustration.RData"))

wide.data <- sim_illustration$`Wide-Format Data`
long.data <- sim_illustration$`Long-Format Data`

# Grid of Doses and Volumes
dose.range <- c(0, 70) # Range of Dose Distribution
dvals <- seq(dose.range[1], dose.range[2], by = 1)  # 71 (dose indices)
gvals <- seq(0, 1, by = 0.02)   # 51 (volume indices)
num.d <- length(dvals) # Number of Dose Indices
num.g <- length(gvals) # Number of Volume Indices


# Stochastic Intervention for Bladder DVHs
# Volume at 40 Gy <= 30%
q_bladder <- 0.3
d.star_bladder <- 40; d.star_bladder <- dvals[findInterval(d.star_bladder, dvals)]
d.star.data_bladder <- long.data[which(long.data$Dose == d.star_bladder),]


# Stochastic Intervention for Skin DVHs
# Volume at 20 Gy <= 20%
q_skin <- 0.2
d.star_skin <- 20; d.star_skin <- dvals[findInterval(d.star_skin, dvals)]
d.star.data_skin <- long.data[which(long.data$Dose == d.star_skin),]


# FIGURE 1: VISUALIZATION OF DVHs
# source(here("figures", "Figure 1", "figure1_code.R"))



# BOOTSTRAP





