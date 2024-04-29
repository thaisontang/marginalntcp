library(here)
res.dir <- here("supplementary simulation", "results") # Results Directory
source(here("functions", "functions.R"))

# Packages
installer(c(c("tidyverse", "compiler"),
            c("gam", "akima", "fields", "abind",
              "monoreg", "pheatmap")))

# Simulation of Functional Monotonic Regression
wd <- here("data")
extension <- wd

# Generate Grid of True Risks
source(here("main simulation", "true_risk_weakconf.R"))

meanfunc.sim <- function(num.cores = 8, seeding = 1) {
  
  # Parallel Processing Set-Up
  total <- 104
  # total <- 16 # 8 cores / 2 rounds
  # total <- 88 # 8 cores / 11 rounds
  R <- total %/% num.cores                      # Simulations per Core
  total.sim <- num.cores * R
  
  para.sim <- function(core.no) {
    
    # FUNCTIONS
    installer <- function(vector) {
      for(i in vector) {
        if(!require(i, character.only = T)) install.packages(i)
        library(i, character.only = T)
      }
    }
    
    # Formula Generator
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
    
    expit <- function(q) {1/(1 + exp(-q))}
    logit <- function(p) {log(p/(1 - p))}
    
    # Custom Save Function
    mysave <- function(..., path = NULL) {
      if(!is.null(path)) {
        save(..., file = paste(path, paste0(deparse(substitute(...)), ".RData"), sep = "/"))
      } else {
        save(..., file = paste0(deparse(substitute(...)), ".RData"))
      }
    }
    
    # PACKAGES
    installer(c(c("tidyverse", "compiler"),
                c("scam", "gam", "akima", "fields", "abind", "monoreg",
                  "pheatmap")))
    enableJIT(3)
    
    
    
    # Initial Parameter Values
    bmin <- 1      
    bmax <- 2      # bmax <- 10
    amin <- 35     # amin <- 5
    amax <- 45     # amax <- 76.5
    
    # DVH Parameters
    up.scale <- 8/3
    beta.parms <- c(0, 0.2, 0.5)
    
    p.grade <- 0.4
    n <- 500
    d.star <- 40
    
    # Reference Dose and DVH Values
    num.d <- 26
    num.g <- 26
    resolution <- 0.001
    min.range <- 30; max.range <- 50
    dvals <- seq(min.range, max.range, length.out = num.d + 2)[-c(1, num.d + 2)]
    gvals <- seq(0, 1, length.out = num.g + 2)[-c(1, num.g + 2)]
    
    # Visualization/Test Set Purposes
    grid.num.d <- 12
    grid.num.g <- 12
    dvals.grid <- seq(min.range, max.range, length.out = grid.num.d + 2)[-c(1, grid.num.d + 2)]
    gvals.grid <- seq(0, 1, length.out = grid.num.g + 2)[-c(1, grid.num.g + 2)]
    
    d.star <- dvals[findInterval(d.star, dvals)]
    g.seq <- seq(0, 1, by = 0.005)
    
    # Functional Mean Model
    out.parms <- c(-18, 0.45, 0.5, 0.5)
    
    ## DGM FUNCTIONS ##
    dg.calc <- function(mu, sd, d) {
      dnorm(d, mean = mu, sd = sd)
    }
    g.calc <- function(mu, sd, d) {
      1 - pnorm(d, mean = mu, sd = sd)
    }
    fun.mod <- function(mu, id) { expit(as.numeric(c(1, mu, x1 = wide.data$x1[id], wide.data$x2[id]) %*% out.parms)) }
    nan.recode <- function(x) {na_if(x, "NaN")}
    expit <- function(q) {1/(1 + exp(-q))}
    logit <- function(p) {log(p/(1 - p))}
    
    # Matrix to Long-Format Data
    mat_pairs <- function(mat, xvals, yvals) {
      grid.vals <- expand.grid(xvals, yvals)
      data.fit2 <- data.frame(xvals = grid.vals$Var1, yvals = grid.vals$Var2, `zval` = as.vector(mat))
      data.fit2[complete.cases(data.fit2),]
    }
    
    # Bayesian Monotone Regression Model
    monofunc <- function(covs,
                         spec.use = spec.short,         # Choice of Number of Iterations/Burn-In
                         parm.use = parm.nomono         # Choice of Model Parametrization
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
    
    # Empirical PMF
    eepmf <- function(orig.vec) {
      uniq.vec <- sort(unique(orig.vec))
      return(data_frame(
        `x` = uniq.vec,
        `pmf` = (ecdf(orig.vec)(uniq.vec) - lag(ecdf(orig.vec)(uniq.vec), default = 0))
      ))
    }
    
    # Beta Function
    betafun <- function(x) {
      beta.parms <- beta.parms.x0	
      if(x == 1) {beta.parms <- beta.parms.x1}	
      alpha.val <- beta.parms[1]; beta.val <- beta.parms[2]
      gamma(alpha.val)*gamma(beta.val)/gamma(alpha.val + beta.val)
    }
    
    # EMPTY VECTORS
    param.names <- c("Logistic", "Polynomial Logistic",
                     "Additive (Doubly)", "Bivariable")
    model.names <- c("Adjusted", "Unadjusted")
    all.param.names <- c(apply(expand.grid(param.names, model.names), 1, function(x) paste(x, collapse = " - ")))
    params.output <- array(NA, dim = c(grid.num.d, grid.num.g, length(all.param.names), R))
    deviance.output <- matrix(NA, nrow = R, ncol = length(all.param.names))
    colnames(deviance.output) <- all.param.names
    brier.output <- deviance.output 
    eff.parms.output <- dic.output <- deviance.output
    
    
    from <- seq(1, total.sim, by = R)[core.no]
    to <- seq(R, total.sim, by = R)[core.no]
    
    for(r in 1:R) {
      
      iter <- from + (r - 1)
      set.seed(iter + seeding)
      
      wide.data <- tibble(
        id = 1:n,
        x1 = rbinom(n, 1, p.grade),
        x2 = rnorm(n),
        eta.breg = as.numeric(cbind(1, x1, x2) %*% beta.parms),
        alpha.breg = up.scale*expit(eta.breg),
        beta.breg = up.scale*(1 - expit(eta.breg)),
        mus = amin + (amax - amin)*(rbeta(n, shape1 = alpha.breg, shape2 = beta.breg)),      # Simulated Means
        sds = runif(n, bmin, bmax)                             # Simulated Standard Deviations
      )
      wide.data$y <- rbinom(n, 1, prob = sapply(1:nrow(wide.data), function(id) {
        fun.mod(mu = wide.data$mus[id], id = id)
      }))
      
      long.data <- data.frame(
        id = rep(1:length(wide.data$mus), each = length(dvals)),
        x1 = rep(wide.data$x1, each = length(dvals)),
        x2 = rep(wide.data$x2, each = length(dvals)),
        mus = rep(wide.data$mus, each = length(dvals)),
        sds = rep(wide.data$sds, each = length(dvals)),
        y = rep(wide.data$y, each = length(dvals)),
        dose = rep(dvals, nrow(wide.data)),
        dvh = as.vector(sapply(1:nrow(wide.data), function(id) {
          g.calc(mu = wide.data$mus[id], sd = wide.data$sds[id], d = dvals)
        })),
        ddvh = as.vector(sapply(1:nrow(wide.data), function(id) {
          dg.calc(mu = wide.data$mus[id], sd = wide.data$sds[id], d = dvals)
        }))
      )
      
      ## FITTED MODELS ##
      
      # Initialization
      nobs <- nrow(long.data)
      x1 <- (long.data$dose - min.range)/(max.range - min.range)
      x2 <- long.data$dvh
      grid.x1 <- (dvals.grid - min.range)/(max.range - min.range)
      grid.x2 <- gvals.grid
      grid.dat <- expand.grid(`Dose` = grid.x1, `Volume` = grid.x2, `x1` = c(0, 1), x2 = unique(wide.data$x2))
      
      # ADJUSTED MODEL
      
      # Training and Testing
      train <- data.matrix(cbind(rep(1, nobs), x1, x2, x1^2, x2^2, x1^3, x2^3, long.data$x1, long.data$x2))
      test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-(3:4)], grid.dat[,-(3:4)]^2, grid.dat[,-(3:4)]^3, grid.dat[,3], grid.dat[,4]))
      alldat <- rbind(train, test)
      colnames(train) <- colnames(test) <- colnames(alldat) <- c("Intercept",
                                                                 "Dose", "Volume", "Dose-Squared", "Volume-Squared", "Dose-Cubed",
                                                                 "Volume-Cubed", "x1", "x2")
      
      # Covariates
      covs.additive <- list(`train` = train, `test` = test, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1, 8:9))
      covs.noadditive <- list(`train` = train, `test` = test, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1:3, 8:9))
      covs.noadditive.poly <- list(`train` = train, `test` = test, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1:5, 8:9))
      
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
      
      # Fit Statistics
      brier.fun <- function(p) {mean((p - long.data$y)^2, na.rm = T)}
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
      colnames(mod.preds.train) <- colnames(mod.preds.test) <- c(
        "LR", "PLR", "Additive", "Bivariable"
      )
      
      # Model Fitting
      
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
      
      
      # UNADJUSTED MODEL
      
      # Training and Testing
      train.mispec <- data.matrix(cbind(rep(1, nobs), x1, x2, x1^2, x2^2, x1^3, x2^3))
      test.mispec <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-(3:4)], grid.dat[,-(3:4)]^2, grid.dat[,-(3:4)]^3))
      alldat.mispec <- rbind(train.mispec, test.mispec)
      colnames(train.mispec) <- colnames(test.mispec) <- colnames(alldat.mispec) <- c("Intercept",
                                                                                      "Dose", "Volume", "Dose-Squared", "Volume-Squared", "Dose-Cubed",
                                                                                      "Volume-Cubed")
      
      # Covariates
      covs.additive <- list(`train` = train.mispec, `test` = test.mispec, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1))
      covs.noadditive <- list(`train` = train.mispec, `test` = test.mispec, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1:3))
      covs.noadditive.poly <- list(`train` = train.mispec, `test` = test.mispec, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1:5))
      
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
      
      # Empty Vectors
      mod.size <- 4
      mod.preds.train.mispec <- matrix(NA, nrow = nrow(train.mispec), ncol = mod.size)
      mod.preds.test.mispec <- matrix(NA, nrow = nrow(test.mispec), ncol = mod.size)
      deviance.mispec <- eff.parms.mispec <- brier.mispec <- dic.mispec <- rep(NA, mod.size)
      colnames(mod.preds.train.mispec) <- colnames(mod.preds.test.mispec) <- c(
        "LR", "PLR", "Additive", "Bivariable"
      )
      
      # Model Fitting
      
      # MODEL 1: LINEAR (LOGISTIC)
      sim.mod1 <- monofunc(covs = covs.noadditive, spec.use = spec.minadapt, parm.use = parm.nomono)
      
      model <- sim.mod1; midx <- 1
      
      # Predictions
      mod.preds.train.mispec[, midx] <- colMeans(model$risk[,1:nobs])
      mod.preds.test.mispec[, midx] <- colMeans(model$risk[,-(1:nobs)])
      p.fit <- mod.preds.train.mispec[, midx]
      
      # Fit Statistics
      deviance.mispec[midx] <- mean(-2*model$loglik)
      eff.parms.mispec[midx] <- var(-2*model$loglik)/2
      brier.mispec[midx] <- brier.fun(p.fit)
      dic.mispec[midx] <- manual_DIC(model)
      
      # Remove Variable to Free up RAM
      rm(sim.mod1); rm(model); gc()
      
      # MODEL 2: POLYNOMIAL (LOGISTIC)
      sim.mod2 <- monofunc(covs = covs.noadditive.poly, spec.use = spec.minadapt, parm.use = parm.nomono)
      
      model <- sim.mod2; midx <- 2
      
      # Predictions
      mod.preds.train.mispec[, midx] <- colMeans(model$risk[,1:nobs])
      mod.preds.test.mispec[, midx] <- colMeans(model$risk[,-(1:nobs)])
      p.fit <- mod.preds.train.mispec[, midx]
      
      # Fit Statistics
      deviance.mispec[midx] <- mean(-2*model$loglik)
      eff.parms.mispec[midx] <- var(-2*model$loglik)/2
      brier.mispec[midx] <- brier.fun(p.fit)
      dic.mispec[midx] <- manual_DIC(model)
      
      # Remove Variable to Free up RAM
      rm(sim.mod2); rm(model); gc()
      
      
      # MODEL 3: DOUBLY ADDITIVE
      sim.mod3 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.twoadditive)
      
      model <- sim.mod3; midx <- 3
      
      # Predictions
      mod.preds.train.mispec[, midx] <- colMeans(model$risk[,1:nobs])
      mod.preds.test.mispec[, midx] <- colMeans(model$risk[,-(1:nobs)])
      p.fit <- mod.preds.train.mispec[, midx]
      
      # Fit Statistics
      deviance.mispec[midx] <- mean(-2*model$loglik)
      eff.parms.mispec[midx] <- var(-2*model$loglik)/2
      brier.mispec[midx] <- brier.fun(p.fit)
      dic.mispec[midx] <- manual_DIC(model)
      
      # Remove Variable to Free up RAM
      rm(sim.mod3); rm(model); gc()
      
      
      # MODEL 4: BIVARIABLE
      sim.mod4 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.bivariate)
      
      model <- sim.mod4; midx <- 4
      
      # Predictions
      mod.preds.train.mispec[, midx] <- colMeans(model$risk[,1:nobs])
      mod.preds.test.mispec[, midx] <- colMeans(model$risk[,-(1:nobs)])
      p.fit <- mod.preds.train.mispec[, midx]
      
      # Fit Statistics
      deviance.mispec[midx] <- mean(-2*model$loglik)
      eff.parms.mispec[midx] <- var(-2*model$loglik)/2
      brier.mispec[midx] <- brier.fun(p.fit)
      dic.mispec[midx] <- manual_DIC(model)
      
      # Remove Variable to Free up RAM
      rm(sim.mod4); rm(model); gc()
      
      # PACKAGING FOR EXPORT
      
      # Empty Vectors
      fit.risk.mispec <- fit.risk <- array(NA, dim = c(grid.num.d, grid.num.g, mod.size))
      
      # Predicted Risk by DVH (Adjusted Model)
      for(z in 1:mod.size) {
        pred <- mod.preds.test[,z]
        test_eval <- as_tibble(cbind(test, pred)) %>% dplyr::select(Dose, Volume, x1, x2, pred) %>%
          group_by(Dose, Volume) %>% mutate(`mpred` = pred)
        test_eval <- test_eval %>% summarise(mpred = mean(mpred)) %>% ungroup
        mpred_mat <- matrix(test_eval$mpred, grid.num.d, grid.num.g, byrow = T)
        fit.risk[,,z] <- mpred_mat
      }
      
      # Predicted Risk by DVH (Unadjusted Model)
      for(z in 1:mod.size) {
        pred <- mod.preds.test.mispec[,z]
        test_eval <- as_tibble(cbind(test.mispec, pred)) %>% dplyr::select(Dose, Volume, pred) %>%
          group_by(Dose, Volume) %>% mutate(`mpred` = pred)
        test_eval <- test_eval %>% summarise(mpred = mean(mpred)) %>% ungroup
        mpred_mat <- matrix(test_eval$mpred, grid.num.d, grid.num.g, byrow = T)
        fit.risk.mispec[,,z] <- mpred_mat
      }
      
      # Combine Adjusted + Unadjusted + True Predictions
      mod.preds <- abind(fit.risk, fit.risk.mispec, along = 3)
      
      # Array Manipulation Toy
      # arr <- array(1:18, dim = c(3, 3, 2)); arr[,,1]; arr[,,2]
      
      # Save
      params.output[,,,r] <- mod.preds
      eff.parms.output[r,] <- nan.recode(c(eff.parms, eff.parms.mispec))
      deviance.output[r,] <- nan.recode(c(deviance, deviance.mispec))
      brier.output[r,] <- nan.recode(c(brier, brier.mispec))
      dic.output[r,] <- nan.recode(c(dic, dic.mispec))
      
      # cat(r,"\n")
    }
    
    save(params.output, file = file.path(outpath, paste0("bivar.params.output_", core.no, ".RData")))
    save(brier.output, file = file.path(outpath, paste0("bivar.brier.output_", core.no, ".RData")))
    save(deviance.output, file = file.path(outpath, paste0("bivar.deviance.output_", core.no, ".RData")))
    save(eff.parms.output, file = file.path(outpath, paste0("bivar.eff.parms.output_", core.no, ".RData")))
    save(dic.output, file = file.path(outpath, paste0("bivar.dic.output_", core.no, ".RData")))
    
    return(NULL)
  }
  
  # Parallel Processing
  library(parallel)
  if(!file.exists("runs")) dir.create(file.path(extension, "runs"))
  outpath <- file.path(extension, "runs")
  
  library(timeDate); library(abind)
  begin <- Sys.time()
  parallel::mclapply(1:num.cores, para.sim, mc.cores = num.cores,
                     mc.silent = FALSE)       # Parallelization Step
  end <- Sys.time() - begin
  
  # Combine Results
  params.output_ALL <- NULL
  brier.output_ALL <- NULL
  eff.parms.output_ALL <- NULL
  deviance.output_ALL <- NULL
  dic.output_ALL <- NULL
  
  for (i in 1:num.cores) {
    load(file = file.path(outpath, paste0('bivar.params.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.brier.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.eff.parms.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.deviance.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.dic.output_', i, ".RData")))
    params.output_ALL <- abind(params.output_ALL, params.output, along = 4)
    brier.output_ALL <- rbind(brier.output_ALL, brier.output)
    eff.parms.output_ALL <- rbind(eff.parms.output_ALL, eff.parms.output)
    deviance.output_ALL <- rbind(deviance.output_ALL, deviance.output)
    dic.output_ALL <- rbind(dic.output_ALL, dic.output)
    
  }
  
  list(`Predictions` = params.output_ALL,
       `Brier` = brier.output_ALL,
       `Deviance` = deviance.output_ALL,
       `effNumParms` = eff.parms.output_ALL,
       `DIC` = dic.output_ALL
  )
}


##########################################################


# Storage
outpath2 <- here("data", "storage")  # Simulation Storage Folder
if(!file.exists(outpath2)) dir.create(file.path(extension, "storage"))

# Pararellization

# (sim.begin <- Sys.time())
# simulation <- meanfunc.sim(num.cores = 8, seeding = 1)
# (simulation$time <- Sys.time() - sim.begin)    # Total Simulation Time
# save(simulation, file = file.path(outpath2, paste0("supp_simulation_weakconf.RData")))