# MARGINAL MODEL for Observed Dataset

# BLADDER DVHs & GU TOXICITY

# Initialization
nobs <- nrow(long.data)
x1 <- (long.data$Dose - min(dvals))/(diff(range(dvals)))  # scaled to [0, 1]
x2 <- long.data$Volume_Bladder
grid.x1 <- (dvals - min(dvals))/diff(range(dvals))
grid.x2 <- gvals
grid.dat <- expand.grid(`Dose` = grid.x1, `Volume` = grid.x2, `Age65` = c(0, 1))
cov.dstar <- d.star.data_bladder
cov.long <- long.data$Age65_numeric

# Training and Testing
train <- data.matrix(cbind(rep(1, nobs), x1, x2, x1^2, x2^2, x1^3, x2^3, cov.long))
test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-3], grid.dat[,-3]^2, grid.dat[,-3]^3, grid.dat[,3]))
alldat <- rbind(train, test)
colnames(train) <- colnames(test) <- colnames(alldat) <- c("Intercept",
                                                           "Dose", "Volume", "Dose-Squared", "Volume-Squared", "Dose-Cubed",
                                                           "Volume-Cubed", "Age65")

# Covariates
covs.additive <- list(`train` = train, `test` = test, y = long.data$`GU Toxicity`, nobs = nobs, add.id = 2:3, covs.id = c(1, 8))
covs.noadditive <- list(`train` = train, `test` = test, y = long.data$`GU Toxicity`, nobs = nobs, add.id = 2:3, covs.id = c(1:3, 8))
covs.noadditive.poly <- list(`train` = train, `test` = test, y = long.data$`GU Toxicity`, nobs = nobs, add.id = 2:3, covs.id = c(1:5, 8))

# Time Saving Options
spec.default <- list(niter = 5000, burnin=2500, adapt=2500)   # 3.37 hours
spec.short <- list(niter = 500, burnin=250, adapt=250)   # 17.9 minutes
spec.minadapt <- list(niter = 2000, burnin=1000, adapt=1000)   # 1.38 hours
spec.long <- list(niter = 15000, burnin=5000, adapt=5000)   # 8+ hours (don't run)

# Model Parametrization
parm.bivariate <- list(birthdeath = 10, settozero = getcmat(2), package=rep(1,3))           # Bivariable Monotone
parm.twoadditive <- list(birthdeath = 10, settozero = getcmat(2)[1:2,], package=c(1,2))     # Two Additive Components
parm.singleadditive <- list(birthdeath = 10, settozero = getcmat(2)[1,], package = 1)       # Single Additive Component (?)
parm.nomono <- list(birthdeath = 0, settozero = getcmat(2)[1:2,], package = 1:2)            # No Additive

# Model Fitting
sim.mod1 <- monofunc(covs = covs.noadditive, spec.use = spec.minadapt, parm.use = parm.nomono)
sim.mod2 <- monofunc(covs = covs.noadditive.poly, spec.use = spec.minadapt, parm.use = parm.nomono)
sim.mod3 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.twoadditive)
sim.mod4 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.bivariate)

p.age65.emp <- prop.table(table(wide.data$age65))

# Summarize in a List
mod.list <- list(
  sim.mod1, sim.mod2, sim.mod3, sim.mod4
)
mod.preds.train <- lapply(mod.list, function(model) {colMeans(model$risk[,1:nobs])})
mod.preds.test <- lapply(mod.list, function(model) {colMeans(model$risk[,-(1:nobs)])})
names(mod.list) <- names(mod.preds.train) <- names(mod.preds.test) <- c(
  "LR", "PLR", "Additive", "Bivariable"
)

# ePMF Functions
fx.0 <- epmf(d.star.data_bladder$Volume[which(cov.dstar == 0)])
fx.1 <- epmf(d.star.data_bladder$Volume[which(cov.dstar == 1)])

# Empirical CDF at g = q
Fq.0 <- mean(d.star.data_bladder$Volume[which(cov.dstar == 0)] <= q_bladder)
Fq.1 <- mean(d.star.data_bladder$Volume[which(cov.dstar == 1)] <= q_bladder)

# Empirical DVHs
obs.g.0 <- d.star.data_bladder$Volume[cov.dstar == 0 & d.star.data_bladder$Volume <= q_bladder]
obs.g.1 <- d.star.data_bladder$Volume[cov.dstar == 1 & d.star.data_bladder$Volume <= q_bladder]

# f(g | x)
est.isws.0 <- 1/length(obs.g.0)
est.isws.1 <- 1/length(obs.g.1)

fit.risk <- array(NA, dim = c(length(dvals), length(gvals), 2 + 1, length(mod.list)))
fit.srisk <- lapply(1:length(mod.list), function(x) { NA })


# Predicted Risk by DVH
for(z in 1:length(mod.list)) {
  pred <- mod.preds.test[[z]]
  test_eval <- as_tibble(cbind(test, pred)) %>% dplyr::select(Dose, Volume, Age65, pred) %>%
    group_by(Dose, Volume) %>%
    mutate(`cpred` = (Age65 == 0)*p.age65.emp[1]*pred + (Age65 == 1)*p.age65.emp[2]*pred)
  test_eval_0 <- test_eval %>% dplyr::filter(Age65 == 0)
  test_eval_1 <- test_eval %>% dplyr::filter(Age65 == 1)
  test_eval <- test_eval %>%
    summarise(cpred = sum(cpred)) %>% ungroup
  pred_mat <- list(
    `x = 0` = matrix(test_eval_0$pred, length(dvals), length(gvals)),
    `x = 1` = matrix(test_eval_1$pred, length(dvals), length(gvals))
  )
  cpred_mat <- matrix(test_eval$cpred, length(dvals), length(gvals), byrow = T)
  
  fit.risk[,,1,z] <- pred_mat$`x = 0`
  fit.risk[,,2,z] <- pred_mat$`x = 1`
  fit.risk[,,3,z] <- cpred_mat
  
  # Stochastic Risk (Conditional on X)
  cond.means <- data.frame(
    `g` = long.data$Volume[which(long.data[,"Dose"] == d.star_bladder)],
    `risk`= mod.preds.train[[z]][which(long.data[,"Dose"] == d.star_bladder)],
    `x` = long.data$age65[which(long.data[,"Dose"] == d.star_bladder)]
  )
  
  # P(Y = 1 | d^*, g, x)
  est.risks.0 <- cond.means$risk[cov.dstar == 0 & cond.means$g <= q_bladder]
  est.risks.1 <- cond.means$risk[cov.dstar == 1 & cond.means$g <= q_bladder]
  
  # Estimated
  est.stochast.0 <- sum(est.risks.0*est.isws.0)
  est.stochast.1 <- sum(est.risks.1*est.isws.1)
  
  fit.srisk[[z]] <- mean(c(est.stochast.0, est.stochast.1)[as.numeric(d.star.data_bladder$age65)])
}
obs.risk <- prop.table(table(wide.data$anyderm))[2]

p.fit <- mod.preds.train
brier.fun <- function(p) {mean((p - long.data$`GU Toxicity`)^2, na.rm = T)}
brier <- sapply(p.fit, brier.fun)

deviance <- sapply(mod.list, function(model) {mean(-2*model$loglik)})
eff.parms <- sapply(mod.list, function(model) {var(-2*model$loglik)/2})
manual_DIC <- function(model) {
  devi <- -2*model$loglik
  eff.num.parms <- var(devi)/2    # Gelman (2004) - BDA
  mean.devi <- mean(devi)
  return(eff.num.parms + mean.devi)
}
dic <- sapply(mod.list, manual_DIC)

# Data Visualization
# source(here(code.dir, "Marginal Model Visualizations.R"))