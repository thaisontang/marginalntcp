# BOOTSTRAP RESAMPLING

# SKIN DVHs & GU TOXICITY
# Estimated Run Time: 26.0 hours (1000 resamples using 8 cores)

# Packages
installer(c(c("tidyverse", "compiler"),
            c("gam", "akima", "fields", "abind",
              "monoreg", "pheatmap")))

# Bootstrap Resampling of Functional Monotonic Regression
extension <- here("simulated_illustration", "data")

# Uncomment for Troubleshooting Purposes
# num.cores = 8
# seeding = 1
# r <- 1; core.no <- 1


meanfunc.boot <- function(num.cores = 8, seeding = 1) {
  
  # Parallel Processing Set-Up
  total <- 1000
  R <- total %/% num.cores                      # Resampling per Core
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
    
    
    # PACKAGES
    installer(c(c("tidyverse", "compiler"),
                c("scam", "gam", "akima", "fields", "abind", "monoreg",
                  "pheatmap")))
    enableJIT(3)
    
    # Grid of Doses and Volumes
    dose.range <- c(0, 70) # Range of Dose Distribution
    dvals <- seq(dose.range[1], dose.range[2], by = 1)  # 71 (dose indices)
    gvals <- seq(0, 1, by = 0.02)   # 51 (volume indices)
    num.d <- length(dvals) # Number of Dose Indices
    num.g <- length(gvals) # Number of Volume Indices
    
    # Threshold
    q_skin <- 0.2
    d.star_skin <- 20
    d.star_skin <- dvals[findInterval(d.star_skin, dvals)]
    
    # Sample Size
    n <- nrow(wide.data)
    
    
    # EMPTY VECTORS
    param.names <- c("Linear", "Polynomial", "Additive", "Bivariable")
    
    # Bootstrap Pointwise Causal NTCP Estimates
    params.output <- array(NA, dim = c(num.d, num.g, length(param.names), R))
    
    # Bootstrap Stochastic Causal NTCP Estimates
    params2.output <- matrix(NA, nrow = R, ncol = length(param.names))
    
    # Bootstrap Empirical Means of the Outcome
    params3.output <- matrix(NA, nrow = R, ncol = 1)
    
    # Bootstrap In-Sample Fit Statistics
    deviance.output <- matrix(NA, nrow = R, ncol = length(param.names))
    colnames(deviance.output) <- param.names
    brier.output <- deviance.output 
    eff.parms.output <- dic.output <- deviance.output[,-5]
    
    from <- seq(1, total.sim, by = R)[core.no]
    to <- seq(R, total.sim, by = R)[core.no]
    
    for(r in 1:R) {
      
      iter <- from + (r - 1)
      set.seed(iter + seeding)
      
      boot.id <- sample(1:n, n, replace = T)
      boot.wide <- wide.data[boot.id,]
      boot.long <- long.data[unlist(lapply(boot.id, function(id) seq((id - 1)*71 + 1, id*71, 1))),]
      boot.d.star.data_skin <- boot.long[which(boot.long$Dose == d.star_skin),]
      
      # Outcome
      y.out <- boot.long$`GI Toxicity`
      
      # Covariates Design Matrix (without Intercept)
      cov.labs <- c("Age65_numeric")    # can include additional covariates
      cov.mat <- data.matrix(boot.long[,cov.labs])
      colnames(cov.mat) <- cov.labs
      
      
      ## FITTED MODELS ##
      
      # Initialization
      nobs <- nrow(boot.long)
      dose.ind <- (boot.long$Dose - min(dvals))/(diff(range(dvals)))  # scaled to [0, 1]
      vol.ind <- boot.long$Volume_Skin
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
      colnames(mod.preds.train) <- colnames(mod.preds.test) <- param.names 
      
      # Model Fitting under different Functional Forms
      
      
      # COVARIATE-ADJUSTED MODELS
      
      # MODEL 1: LINEAR (LOGISTIC)
      sim.mod1 <- monofunc(covs = covs.noadditive, spec.use = spec.minadapt, parm.use = parm.nomono,
                           iter = iter, seeding = seeding)
      
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
      sim.mod2 <- monofunc(covs = covs.noadditive.poly, spec.use = spec.minadapt, parm.use = parm.nomono,
                           iter = iter, seeding = seeding)
      
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
      sim.mod3 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.twoadditive,
                           iter = iter, seeding = seeding)
      
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
      sim.mod4 <- monofunc(covs = covs.additive, spec.use = spec.minadapt, parm.use = parm.bivariate,
                           iter = iter, seeding = seeding)
      
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
      
      
      
      
      
      # Empty Vectors to Collect Estimates
      fit.risk <- array(NA, dim = c(num.d, num.g, mod.size))
      fit.srisk <- lapply(1:mod.size, function(x) { NA })
      
      # P(X = x)
      emp.dist <- boot.long[,cov.labs] %>% group_by_all() %>% summarize(prop = n()/nrow(boot.long))
      
      # f(V_d | X) / f(G_d | X)
      # The code below is for the special case of a binary covariate
      
      # Volumes relevant to Stochastic Intervention; F(G_d = q | X)
      obs.g.0 <- boot.d.star.data_skin$Volume_Skin[cov.dstar == 0 & boot.d.star.data_skin$Volume_Skin <= q_skin]
      obs.g.1 <- boot.d.star.data_skin$Volume_Skin[cov.dstar == 1 & boot.d.star.data_skin$Volume_Skin <= q_skin]
      
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
          `g` = boot.long$Volume_Skin[which(boot.long[,"Dose"] == d.star_skin)],
          `risk`= mod.preds.train[,z][which(boot.long[,"Dose"] == d.star_skin)]
        ) %>% cbind(boot.long[which(boot.long[,"Dose"] == d.star_skin), cov.labs])
        
        # P(Y = 1 | d^*, g, x) relevant to Stochastic Intervention
        # The code below is for the special case of a binary covariate
        est.risks.0 <- cond.means$risk[cov.dstar == 0 & cond.means$g <= q_skin]
        est.risks.1 <- cond.means$risk[cov.dstar == 1 & cond.means$g <= q_skin]
        
        # Estimated
        est.stochast.0 <- sum(est.risks.0*est.isws.0)
        est.stochast.1 <- sum(est.risks.1*est.isws.1)
        
        fit.srisk[[z]] <- mean(c(est.stochast.0, est.stochast.1)[boot.d.star.data_skin[,cov.labs] %>% pull + 1])
      }
      obs.risk <- prop.table(table(boot.wide$`GI Toxicity`))[2]
      
      
      # Pointwise Causal NTCP Estimates
      params.output[,,,r] <- array(fit.risk, dim = c(num.d, num.g, length(param.names), 1))
      
      # Stochastic Causal NTCP Estimates
      nan.recode <- function(x) {na_if(x, "NaN")}
      params2.output[r,] <- nan.recode(unlist(fit.srisk))
      
      # Empirical Mean Estimates
      params3.output[r,] <- nan.recode(unname(obs.risk))
      
      # In-Sample Statistics
      eff.parms.output[r,] <- nan.recode(eff.parms)
      deviance.output[r,] <- nan.recode(deviance)
      brier.output[r,] <- nan.recode(brier)
      dic.output[r,] <- nan.recode(dic)
      
      # cat(r,"\n")
    }
    
    save(params.output, file = file.path(outpath, paste0("boot_bivar.params.output_", core.no, ".RData")))
    save(params2.output, file = file.path(outpath, paste0("boot_bivar.params2.output_", core.no, ".RData")))
    save(params3.output, file = file.path(outpath, paste0("boot_bivar.params3.output_", core.no, ".RData")))
    save(brier.output, file = file.path(outpath, paste0("boot_bivar.brier.output_", core.no, ".RData")))
    save(deviance.output, file = file.path(outpath, paste0("boot_bivar.deviance.output_", core.no, ".RData")))
    save(eff.parms.output, file = file.path(outpath, paste0("boot_bivar.eff.parms.output_", core.no, ".RData")))
    save(dic.output, file = file.path(outpath, paste0("boot_bivar.dic.output_", core.no, ".RData")))
    
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
  params.output2_ALL <- NULL
  params.output3_ALL <- NULL
  brier.output_ALL <- NULL
  eff.parms.output_ALL <- NULL
  deviance.output_ALL <- NULL
  dic.output_ALL <- NULL
  
  for (i in 1:num.cores) {
    load(file = file.path(outpath, paste0('boot_bivar.params.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.params2.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.params3.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.brier.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.eff.parms.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.deviance.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.dic.output_', i, ".RData")))
    
    params.output_ALL <- abind(params.output_ALL, params.output, along = 5)
    params.output2_ALL <- rbind(params.output2_ALL, params2.output)
    params.output3_ALL <- rbind(params.output3_ALL, params3.output)
    brier.output_ALL <- rbind(brier.output_ALL, brier.output)
    eff.parms.output_ALL <- rbind(eff.parms.output_ALL, eff.parms.output)
    deviance.output_ALL <- rbind(deviance.output_ALL, deviance.output)
    dic.output_ALL <- rbind(dic.output_ALL, dic.output)
  }
  
  list(`Predictions` = params.output_ALL,
       `Stochast Predictions` = params.output2_ALL,
       `Obs Means` = params.output3_ALL,
       `Brier` = brier.output_ALL,
       `Deviance` = deviance.output_ALL,
       `effNumParms` = eff.parms.output_ALL,
       `DIC` = dic.output_ALL
  )
}


##########################################################


# Storage
outpath2 <- here("simulated_illustration", "data", "storage")  # Resampling Storage Folder
if(!file.exists(outpath2)) dir.create(file.path(extension, "storage"))

# Pararellization

# Pointwise & Stochastic Causal Risks
# (sim.begin <- Sys.time())
# bootresamp <- meanfunc.boot(num.cores = 8, seeding = 1)
# (bootresamp$time <- Sys.time() - sim.begin)    # Total Bootstrap Resampling Time
# save(bootresamp, file = file.path(outpath2, paste0("ntcp_bootresamp_age65_skin.RData")))
