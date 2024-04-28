# MARGINAL MODEL for Observed Dataset

# SKIN DVHs & GI TOXICITY
# Estimated Run Time: 20 minutes

# Outcome: GI Toxicity
y.out <- long.data$`GI Toxicity`

# Covariates Design Matrix (without Intercept)
cov.labs <- c("Age65_numeric")    # can include additional covariates
cov.mat <- data.matrix(long.data[,cov.labs])
colnames(cov.mat) <- cov.labs

# Initialization
nobs <- nrow(long.data)
dose.ind <- (long.data$Dose - min(dvals))/(diff(range(dvals)))  # scaled to [0, 1]
vol.ind <- long.data$Volume_Skin
grid.dose.ind <- (dvals - min(dvals))/diff(range(dvals))
grid.vol.ind <- gvals
grid.dat <- expand.grid(grid.dose.ind, grid.vol.ind, apply(cov.mat, 2, unique))
colnames(grid.dat) <- c("Dose", "Volume", colnames(cov.mat))
cov.dstar <- d.star.data_skin$Age65_numeric

# Training and Testing
cov.idx <- 3:(3 + ncol(cov.mat) - 1)
train <- data.matrix(cbind(rep(1, nobs), dose.ind, vol.ind, dose.ind^2, vol.ind^2, dose.ind^3, vol.ind^3, cov.mat))
test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-cov.idx], grid.dat[,-cov.idx]^2, grid.dat[,-cov.idx]^3, grid.dat[,cov.idx]))
alldat <- rbind(train, test)
colnames(train) <- colnames(test) <- colnames(alldat) <- c("Intercept", "Dose", "Volume",
                                                           "Dose-Squared", "Volume-Squared",
                                                           "Dose-Cubed", "Volume-Cubed", colnames(cov.mat))

# Covariates
cov.idx2 <- 8:(8 + ncol(cov.mat) - 1)
covs.additive <- list(`train` = train, `test` = test, y = y.out, nobs = nobs, add.id = 2:3, covs.id = c(1, cov.idx2))
covs.noadditive <- list(`train` = train, `test` = test, y = y.out, nobs = nobs, add.id = 2:3, covs.id = c(1:3, cov.idx2))
covs.noadditive.poly <- list(`train` = train, `test` = test, y = y.out, nobs = nobs, add.id = 2:3, covs.id = c(1:5, cov.idx2))

# Time Saving Options
spec.default <- list(niter = 5000, burnin=2500, adapt=2500)   # 3.37 hours
spec.short <- list(niter = 500, burnin=250, adapt=250)   # 17.9 minutes
spec.minadapt <- list(niter = 2000, burnin=1000, adapt=1000)   # 1.38 hours
spec.long <- list(niter = 15000, burnin=5000, adapt=5000)   # 8+ hours (don't run)

# Model Parametrization
parm.bivariate <- list(birthdeath = 10, settozero = getcmat(2), package=rep(1,3))           # Bivariable Monotone
parm.twoadditive <- list(birthdeath = 10, settozero = getcmat(2)[1:2,], package=c(1,2))     # Two Additive Components
parm.singleadditive <- list(birthdeath = 10, settozero = getcmat(2)[1,], package = 1)       # Single Additive Component
parm.nomono <- list(birthdeath = 0, settozero = getcmat(2)[1:2,], package = 1:2)            # No Additive

# Monotone Regression Function
monofunc <- function(covs,
                     spec.use = spec.short,         # Choice of Number of Iterations/Burn-In
                     parm.use = parm.nomono,        # Choice of Model Parametrization
                     iter = 0, seeding = 0          # Starting & Replicate Seeds for Reproducibility
) {
  (begin <- Sys.time())
  results <- monosurv(niter=spec.use$niter, burnin=spec.use$burnin, adapt=spec.use$adapt,
                      refresh = 10, thin = 5, seed = iter + seeding,
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

# Fit Statistics
brier.fun <- function(p) {mean((p - y.out)^2, na.rm = T)}
manual_DIC <- function(model) {
  devi <- -2*model$loglik
  eff.num.parms <- var(devi)/2    # Gelman (2004) - BDA
  mean.devi <- mean(devi)
  return(eff.num.parms + mean.devi)
}

# Empty Vectors
mod.size <- 4
mod.preds.train <- matrix(NA, nrow = nrow(train), ncol = mod.size)
mod.preds.test <- matrix(NA, nrow = nrow(test), ncol = mod.size)
deviance <- eff.parms <- brier <- dic <- rep(NA, mod.size)
param.names <- c("Linear", "Polynomial", "Additive", "Bivariable")
colnames(mod.preds.train) <- colnames(mod.preds.test) <- param.names 

# Model Fitting under different Functional Forms

# COVARIATE-ADJUSTED MODELS

# MODEL 1: LINEAR (LOGISTIC)
sim.mod1 <- monofunc(covs = covs.noadditive, spec.use = spec.minadapt, parm.use = parm.nomono)

model <- sim.mod1; midx <- 1

# Predictions
mod.preds.train[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train[, midx]

# Fit Statistics
deviance[midx] <- mean(-2*model$loglik)
eff.parms[midx] <- var(-2*model$loglik)/2
brier[midx] <- brier.fun(p.fit)
dic[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod1); rm(model); gc()

# MODEL 2: POLYNOMIAL (LOGISTIC)
sim.mod2 <- monofunc(covs = covs.noadditive.poly, spec.use = spec.minadapt, parm.use = parm.nomono)

model <- sim.mod2; midx <- 2

# Predictions
mod.preds.train[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train[, midx]

# Fit Statistics
deviance[midx] <- mean(-2*model$loglik)
eff.parms[midx] <- var(-2*model$loglik)/2
brier[midx] <- brier.fun(p.fit)
dic[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod2); rm(model); gc()


# MODEL 3: DOUBLY ADDITIVE
sim.mod3 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.twoadditive)

model <- sim.mod3; midx <- 3

# Predictions
mod.preds.train[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train[, midx]

# Fit Statistics
deviance[midx] <- mean(-2*model$loglik)
eff.parms[midx] <- var(-2*model$loglik)/2
brier[midx] <- brier.fun(p.fit)
dic[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod3); rm(model); gc()


# MODEL 4: BIVARIABLE
sim.mod4 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.bivariate)

model <- sim.mod4; midx <- 4

# Predictions
mod.preds.train[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train[, midx]

# Fit Statistics
deviance[midx] <- mean(-2*model$loglik)
eff.parms[midx] <- var(-2*model$loglik)/2
brier[midx] <- brier.fun(p.fit)
dic[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod4); rm(model); gc()




# COVARIATE-UNADJUSTED MODELS

# Training and Testing
train.unadj <- data.matrix(cbind(rep(1, nobs), dose.ind, vol.ind, dose.ind^2, vol.ind^2, dose.ind^3, vol.ind^3))
test.unadj <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-cov.idx], grid.dat[,-cov.idx]^2, grid.dat[,-cov.idx]^3))
alldat.unadj <- rbind(train.unadj, test.unadj)
colnames(train.unadj) <- colnames(test.unadj) <- colnames(alldat.unadj) <- c("Intercept", "Dose", "Volume",
                                                                             "Dose-Squared", "Volume-Squared",
                                                                             "Dose-Cubed", "Volume-Cubed")

# Covariates
covs.additive <- list(`train` = train.unadj, `test` = test.unadj, y = y.out, nobs = nobs, add.id = 2:3, covs.id = c(1))
covs.noadditive <- list(`train` = train.unadj, `test` = test.unadj, y = y.out, nobs = nobs, add.id = 2:3, covs.id = c(1:3))
covs.noadditive.poly <- list(`train` = train.unadj, `test` = test.unadj, y = y.out, nobs = nobs, add.id = 2:3, covs.id = c(1:5))

# Time Saving Options
spec.default <- list(niter = 5000, burnin=2500, adapt=2500)   # 3.37 hours
spec.short <- list(niter = 500, burnin=250, adapt=250)   # 17.9 minutes
spec.minadapt <- list(niter = 2000, burnin=1000, adapt=1000)   # 1.38 hours
spec.long <- list(niter = 15000, burnin=5000, adapt=5000)   # 8+ hours (don't run)

# Model Parametrization
parm.bivariate <- list(birthdeath = 10, settozero = getcmat(2), package=rep(1,3))           # Bivariable Monotone
parm.twoadditive <- list(birthdeath = 10, settozero = getcmat(2)[1:2,], package=c(1,2))     # Two Additive Components
parm.singleadditive <- list(birthdeath = 10, settozero = getcmat(2)[1,], package = 1)       # Single Additive Component
parm.nomono <- list(birthdeath = 0, settozero = getcmat(2)[1:2,], package = 1:2)            # No Additive

# Empty Vectors
mod.size <- 4
mod.preds.train.unadj <- matrix(NA, nrow = nrow(train.unadj), ncol = mod.size)
mod.preds.test.unadj <- matrix(NA, nrow = nrow(test.unadj), ncol = mod.size)
deviance.unadj <- eff.parms.unadj <- brier.unadj <- dic.unadj <- rep(NA, mod.size)
colnames(mod.preds.train.unadj) <- colnames(mod.preds.test.unadj) <- param.names


# Model Fitting

# MODEL 1: LINEAR (LOGISTIC)
sim.mod1 <- monofunc(covs = covs.noadditive, spec.use = spec.minadapt, parm.use = parm.nomono)

model <- sim.mod1; midx <- 1

# Predictions
mod.preds.train.unadj[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test.unadj[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train.unadj[, midx]

# Fit Statistics
deviance.unadj[midx] <- mean(-2*model$loglik)
eff.parms.unadj[midx] <- var(-2*model$loglik)/2
brier.unadj[midx] <- brier.fun(p.fit)
dic.unadj[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod1); rm(model); gc()

# MODEL 2: POLYNOMIAL (LOGISTIC)
sim.mod2 <- monofunc(covs = covs.noadditive.poly, spec.use = spec.minadapt, parm.use = parm.nomono)

model <- sim.mod2; midx <- 2

# Predictions
mod.preds.train.unadj[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test.unadj[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train.unadj[, midx]

# Fit Statistics
deviance.unadj[midx] <- mean(-2*model$loglik)
eff.parms.unadj[midx] <- var(-2*model$loglik)/2
brier.unadj[midx] <- brier.fun(p.fit)
dic.unadj[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod2); rm(model); gc()


# MODEL 3: DOUBLY ADDITIVE
sim.mod3 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.twoadditive)

model <- sim.mod3; midx <- 3

# Predictions
mod.preds.train.unadj[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test.unadj[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train.unadj[, midx]

# Fit Statistics
deviance.unadj[midx] <- mean(-2*model$loglik)
eff.parms.unadj[midx] <- var(-2*model$loglik)/2
brier.unadj[midx] <- brier.fun(p.fit)
dic.unadj[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod3); rm(model); gc()


# MODEL 4: BIVARIABLE
sim.mod4 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.bivariate)

model <- sim.mod4; midx <- 4

# Predictions
mod.preds.train.unadj[, midx] <- colMeans(model$risk[,1:nobs])
mod.preds.test.unadj[, midx] <- colMeans(model$risk[,-(1:nobs)])
p.fit <- mod.preds.train.unadj[, midx]

# Fit Statistics
deviance.unadj[midx] <- mean(-2*model$loglik)
eff.parms.unadj[midx] <- var(-2*model$loglik)/2
brier.unadj[midx] <- brier.fun(p.fit)
dic.unadj[midx] <- manual_DIC(model)

# Remove Variable to Free up RAM
rm(sim.mod4); rm(model); gc()


# Empty Vectors to Collect Estimates
fit.risk.unadj <- fit.risk <- array(NA, dim = c(num.d, num.g, mod.size))
fit.srisk.unadj <- fit.srisk <- lapply(1:mod.size, function(x) { NA })

# P(X = x)
emp.dist <- long.data[,cov.labs] %>% group_by_all() %>% summarize(prop = n()/nrow(long.data))

# f(V_d | X) / f(G_d | X)
# The code below is for the special case of a binary covariate

# Volumes relevant to Stochastic Intervention; F(G_d = q | X)
obs.g.0 <- d.star.data_skin$Volume_Skin[cov.dstar == 0 & d.star.data_skin$Volume_Skin <= q_skin]
obs.g.1 <- d.star.data_skin$Volume_Skin[cov.dstar == 1 & d.star.data_skin$Volume_Skin <= q_skin]

# Importance Sampling Weights
est.isws.0 <- 1/length(obs.g.0)
est.isws.1 <- 1/length(obs.g.1)

# POINT ESTIMATES (Covariate-Adjusted Models - Linear, Polynomial, Additive, Bivariable)
for(z in 1:mod.size) {
  
  # POINTWISE CAUSAL NTCP
  pred <- mod.preds.test[,z]
  test_eval <- as_tibble(cbind(test, pred)) %>% dplyr::select(c("Dose", "Volume", cov.labs, "pred")) %>%
    group_by(Dose, Volume) %>% left_join(emp.dist, by = cov.labs) %>%
    mutate(`cpred` = prop*pred)
  test_eval <- test_eval %>%
    summarise(cpred = sum(cpred)) %>% ungroup
  cpred_mat <- matrix(test_eval$cpred, length(dvals), length(gvals), byrow = T)
  fit.risk[,,z] <- cpred_mat
  
  
  # STOCHASTIC CAUSAL NTCP
  
  # Stochastic Risk (Conditional on X)
  cond.means <- data.frame(
    `g` = long.data$Volume_Skin[which(long.data[,"Dose"] == d.star_skin)],
    `risk`= mod.preds.train[,z][which(long.data[,"Dose"] == d.star_skin)]
  ) %>% cbind(long.data[which(long.data[,"Dose"] == d.star_skin), cov.labs])
  
  # P(Y = 1 | d^*, g, x) relevant to Stochastic Intervention
  # The code below is for the special case of a binary covariate
  est.risks.0 <- cond.means$risk[cov.dstar == 0 & cond.means$g <= q_skin]
  est.risks.1 <- cond.means$risk[cov.dstar == 1 & cond.means$g <= q_skin]
  
  # Estimated
  est.stochast.0 <- sum(est.risks.0*est.isws.0)
  est.stochast.1 <- sum(est.risks.1*est.isws.1)
  
  fit.srisk[[z]] <- mean(c(est.stochast.0, est.stochast.1)[d.star.data_skin[,cov.labs] %>% pull + 1])
}



# POINT ESTIMATES (Covariate-Unadjusted Models - Linear, Polynomial, Additive, Bivariable)
for(z in 1:mod.size) {
  
  # POINTWISE CAUSAL NTCP
  pred <- mod.preds.test.unadj[,z]
  test_eval <- as_tibble(cbind(test.unadj, pred)) %>% dplyr::select(Dose, Volume, pred) %>%
    group_by(Dose, Volume) %>% mutate(`mpred` = pred)
  test_eval <- test_eval %>% summarise(mpred = mean(mpred)) %>% ungroup
  mpred_mat <- matrix(test_eval$mpred, num.d, num.g, byrow = T)
  fit.risk.unadj[,,z] <- mpred_mat
  
  # STOCHASTIC CAUSAL NTCP
  
  # Stochastic Risk (Conditional on X)
  cond.means <- data.frame(
    `g` = long.data$Volume_Skin[which(long.data[,"Dose"] == d.star_skin)],
    `risk`= mod.preds.train.unadj[,z][which(long.data[,"Dose"] == d.star_skin)]
  )
  
  # P(Y = 1 | d^*, g, x) relevant to Stochastic Intervention
  # The code below is for the special case of a binary covariate
  est.risks <- cond.means$risk[cond.means$g <= q_skin]
  
  # Importance Sampling Weight
  est.isws <- 1/length(d.star.data_skin$Volume_Skin[d.star.data_skin$Volume_Skin <= q_skin])
  
  # Estimate
  est.stochast <- sum(est.risks*est.isws)
  
  fit.srisk.unadj[[z]] <- est.stochast
}

# Empirical Mean
obs.risk <- prop.table(table(wide.data$`GI Toxicity`))[2]




# POINTWISE CAUSAL NTCPs (Contour Plots of Point Estimates)

# Covariate-Unadjusted Models - Linear, Polynomial, Additive, Bivariable
op <- par(las = 1, mfrow = c(2, 2), mar = c(4.5, 5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
for(z in 1:length(param.names)) {
  contour(x = dvals, y = gvals, z = fit.risk.unadj[,,z], nlevels = 10, col = "black",
          xlim = range(dvals), lwd = 1.5, xlab = "Dose (Gy)", ylab = "Volume")
  abline(h = c(0, 1), v = range(dvals), col = "gray")
}
par(op)

# Covariate-Adjusted Models - Linear, Polynomial, Additive, Bivariable
op <- par(las = 1, mfrow = c(2, 2), mar = c(4.5, 5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
for(z in 1:length(param.names)) {
  contour(x = dvals, y = gvals, z = fit.risk[,,z], nlevels = 10, col = "black",
          xlim = range(dvals), lwd = 1.5, xlab = "Dose (Gy)", ylab = "Volume")
  abline(h = c(0, 1), v = range(dvals), col = "gray")
}
par(op)

# STOCHASTIC CAUSAL RISK RATIO (Stochastic Intervention)

# COVARIATE-UNADJUSTED MODEL
# Point Estimate (no CI reported)
scrr.unadj <- unlist(fit.srisk.unadj)/obs.risk
names(scrr.unadj) <- param.names; scrr.unadj

# In-Sample Metrics
insample.metrics.unadj <- matrix(c(brier.unadj, deviance.unadj, eff.parms.unadj, dic.unadj), byrow = T,
                                 nrow = 4, ncol = length(param.names))
rownames(insample.metrics.unadj) <- c("Brier Score", "Deviance", "Effective Num. Parameters", "DIC")
colnames(insample.metrics.unadj) <- param.names; round(insample.metrics.unadj, 4)

# COVARIATE-ADJUSTED MODEL
# Point Estimate (no CI reported)
scrr <- unlist(fit.srisk)/obs.risk
names(scrr) <- param.names; scrr

# In-Sample Metrics
insample.metrics <- matrix(c(brier, deviance, eff.parms, dic), byrow = T,
                           nrow = 4, ncol = length(param.names))
rownames(insample.metrics) <- c("Brier Score", "Deviance", "Effective Num. Parameters", "DIC")
colnames(insample.metrics) <- param.names; round(insample.metrics, 4)

