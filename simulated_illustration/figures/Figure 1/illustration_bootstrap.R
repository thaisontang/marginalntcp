library(here)
res.dir <- here("results", "Manuscript", "Data Application", "Exploratory",
                "Bootstrap Resampling", "Bladder")

source(here("analysis code", "functions.R"))

# Packages
installer(c(c("tidyverse", "compiler"),
            c("gam", "akima", "fields", "abind",
              "monoreg", "pheatmap")))

# Bootstrap Resampling of Functional Monotonic Regression
wd <- here("data")
extension <- wd

# ~10 hours for 104 replicates
num.cores = 8
seeding = 1
r <- 1; core.no <- 1

# Loading Required Data
source(here("analysis code", "Manuscript Code", "Data Application",
            "Exploratory", "2022_03_31 - Data Load.R"))

# Run for Current Working Dataset (Skin DVH)
source(here("analysis code", "Manuscript Code", "Data Application",
            "Exploratory", "Marginal Model - Skin DVH.R"))


meanfunc.boot <- function(num.cores = 8, seeding = 1) {
  
  # Parallel Processing Set-Up
  total <- 500
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
    
    
    
    # Threshold
    q <- 0.2
    d.star <- 50
    
    # Reference Dose and DVH Values
    # dose.range <- range(unlist((ListSubset(dvh, sub.vec = which(dvh.pat$ROI %in% oars.id)))$Dose))
    dose.range <- c(0, 70)
    dvals <- seq(dose.range[1], dose.range[2], by = 1)  # 71
    gvals <- seq(0, 1, by = 0.02)   # 51
    num.d <- length(dvals); num.g <- length(gvals)
    d.star <- dvals[findInterval(d.star, dvals)]
    n <- nrow(wide.data)
    
    # Empirical PMF
    epmf <- function(orig.vec) {
      uniq.vec <- sort(unique(orig.vec))
      return(data_frame(
        `x` = uniq.vec,
        `pmf` = (ecdf(orig.vec)(uniq.vec) - lag(ecdf(orig.vec)(uniq.vec), default = 0))
      ))
    }
    
    # EMPTY VECTORS
    param.names <- c("Logistic", "Polynomial Logistic",
                     "Additive (Doubly)", "Bivariable")
    params.output <- array(NA, dim = c(num.d, num.g, 4 + 1, length(param.names), R))
    params2.output <- matrix(NA, nrow = R, ncol = length(param.names))
    params3.output <- matrix(NA, nrow = R, ncol = 1)
    deviance.output <- matrix(NA, nrow = R, ncol = length(param.names))
    colnames(deviance.output) <- param.names
    brier.output <- deviance.output 
    eff.parms.output <- dic.output <- deviance.output[,-5]
    
    stochast.save <- lapply(vector("list", R), function(x) list(
      `CondMeans00_Logistic` = NULL,
      `CondMeans00_PolyLog` = NULL,
      `CondMeans00_Additive` = NULL,
      `CondMeans00_Bivariable` = NULL,
      `CondMeans01_Logistic` = NULL,
      `CondMeans01_PolyLog` = NULL,
      `CondMeans01_Additive` = NULL,
      `CondMeans01_Bivariable` = NULL,
      `CondMeans10_Logistic` = NULL,
      `CondMeans10_PolyLog` = NULL,
      `CondMeans10_Additive` = NULL,
      `CondMeans10_Bivariable` = NULL,
      `CondMeans11_Logistic` = NULL,
      `CondMeans11_PolyLog` = NULL,
      `CondMeans11_Additive` = NULL,
      `CondMeans11_Bivariable` = NULL,
      
      `ISWs00_Logistic` = NULL,
      `ISWs00_PolyLog` = NULL,
      `ISWs00_Additive` = NULL,
      `ISWs00_Bivariable` = NULL,
      `ISWs01_Logistic` = NULL,
      `ISWs01_PolyLog` = NULL,
      `ISWs01_Additive` = NULL,
      `ISWs01_Bivariable` = NULL,
      `ISWs10_Logistic` = NULL,
      `ISWs10_PolyLog` = NULL,
      `ISWs10_Additive` = NULL,
      `ISWs10_Bivariable` = NULL,
      `ISWs11_Logistic` = NULL,
      `ISWs11_PolyLog` = NULL,
      `ISWs11_Additive` = NULL,
      `ISWs11_Bivariable` = NULL,
      
      `Age65` = NULL,
      `Stage` = NULL
    ))
    data.save <- vector("list", R)
    
    from <- seq(1, total.sim, by = R)[core.no]
    to <- seq(R, total.sim, by = R)[core.no]
    
    for(r in 1:R) {
      
      iter <- from + (r - 1)
      set.seed(iter + seeding)
      
      boot.id <- sample(1:n, n, replace = T)
      boot.wide <- wide.data[boot.id,]
      boot.long <- boot.wide %>% left_join(skin.dvh, by = "studynum") %>%
        pivot_longer(cols = any_of(colnames(skin.dvh)[-1]), names_to = "Dose", values_to = "Volume") %>%
        mutate(Dose = as.numeric(Dose),
               y = as.numeric(anyderm) - 1)
      boot.d.star <- boot.long[which(boot.long$Dose == d.star),]
      boot.wide.age65 <- boot.wide$age65
      boot.wide.stage <- boot.wide$tstage
      boot.long.age65 <- boot.long$age65
      boot.long.stage <- boot.long$tstage
      boot.d.star.age65 <- boot.d.star$age65
      boot.d.star.stage <- boot.d.star$tstage
      
      
      
      ## FITTED MODELS ##
      
      # Initialization
      nobs <- nrow(boot.long)
      x1 <- (boot.long$Dose - min(dvals))/(diff(range(dvals)))
      x2 <- boot.long$Volume
      grid.x1 <- (dvals - min(dvals))/diff(range(dvals))
      grid.x2 <- gvals
      grid.dat <- expand.grid(`Dose` = grid.x1, `Volume` = grid.x2,
                              `Stage` = c(0, 1), `Age65` = c(0, 1))
      
      # Training and Testing
      train <- data.matrix(cbind(rep(1, nobs), x1, x2, x1^2, x2^2, x1^3, x2^3, as.numeric(boot.long.age65) - 1, as.numeric(boot.long.stage) - 1))
      test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-(3:4)], grid.dat[,-(3:4)]^2, grid.dat[,-(3:4)]^3, grid.dat[,3:4]))
      alldat <- rbind(train, test)
      colnames(train) <- colnames(test) <- colnames(alldat) <- c("Intercept",
                                                                 "Dose", "Volume", "Dose-Squared", "Volume-Squared", "Dose-Cubed",
                                                                 "Volume-Cubed", "Age65", "Stage")
      
      # Covariates
      covs.additive <- list(`train` = train, `test` = test, y = boot.long$y, nobs = nobs, add.id = 2:3, covs.id = c(1, 8:9))
      covs.noadditive <- list(`train` = train, `test` = test, y = boot.long$y, nobs = nobs, add.id = 2:3, covs.id = c(1:3, 8:9))
      covs.noadditive.poly <- list(`train` = train, `test` = test, y = boot.long$y, nobs = nobs, add.id = 2:3, covs.id = c(1:5, 8:9))
      
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
      
      # Summarize in a List
      mod.list <- list(
        sim.mod1, sim.mod2, sim.mod3, sim.mod4
      )
      mod.preds.train <- lapply(mod.list, function(model) {colMeans(model$risk[,1:nobs])})
      mod.preds.test <- lapply(mod.list, function(model) {colMeans(model$risk[,-(1:nobs)])})
      names(mod.list) <- names(mod.preds.train) <- names(mod.preds.test) <- c(
        "LR", "PLR", "Additive", "Bivariable"
      )
      
      
      # Empirical Estimators
      p.stage.emp <- prop.table(table(boot.wide.age65, boot.wide.stage))
      p.stage.emp.age65_0.stage_0 <- p.stage.emp[1,1]
      p.stage.emp.age65_0.stage_1 <- p.stage.emp[1,2]
      p.stage.emp.age65_1.stage_0 <- p.stage.emp[2,1]
      p.stage.emp.age65_1.stage_1 <- p.stage.emp[2,2]
      
      # ePMF Functions
      fx.00 <- epmf(boot.d.star$Volume[
        which((as.numeric(boot.d.star.age65) - 1 == 0) &
                (as.numeric(boot.d.star.stage) - 1 == 0))])
      fx.01 <- epmf(boot.d.star$Volume[
        which((as.numeric(boot.d.star.age65) - 1 == 0) &
                (as.numeric(boot.d.star.stage) - 1 == 1))])
      fx.10 <- epmf(boot.d.star$Volume[
        which(as.numeric(boot.d.star.age65) - 1 == 1 &
                as.numeric(boot.d.star.stage) - 1 == 0)])
      fx.11 <- epmf(boot.d.star$Volume[
        which(as.numeric(boot.d.star.age65) - 1 == 1 &
                as.numeric(boot.d.star.stage) - 1 == 1)])
      
      # Empirical CDF at g = q
      Fq.00 <- mean(boot.d.star$Volume[
        which((as.numeric(boot.d.star.age65) - 1 == 0) &
                (as.numeric(boot.d.star.stage) - 1 == 0))
      ] <= q)
      Fq.01 <- mean(boot.d.star$Volume[
        which((as.numeric(boot.d.star.age65) - 1 == 0) &
                (as.numeric(boot.d.star.stage) - 1 == 1))
      ] <= q)
      Fq.10 <- mean(boot.d.star$Volume[
        which(as.numeric(boot.d.star.age65) - 1 == 1 &
                as.numeric(boot.d.star.stage) - 1 == 0)
      ] <= q)
      Fq.11 <- mean(boot.d.star$Volume[
        which(as.numeric(boot.d.star.age65) - 1 == 1 &
                as.numeric(boot.d.star.stage) - 1 == 1)
      ] <= q)
      
      # Empirical DVHs
      obs.g.00 <- boot.d.star$Volume[
        as.numeric(boot.d.star.age65) - 1 == 0 &
          as.numeric(boot.d.star.stage) - 1 == 0 &
          boot.d.star$Volume <= q]
      obs.g.01 <- boot.d.star$Volume[
        as.numeric(boot.d.star.age65) - 1 == 0 &
          as.numeric(boot.d.star.stage) - 1 == 1 &
          boot.d.star$Volume <= q]
      obs.g.10 <- boot.d.star$Volume[
        as.numeric(boot.d.star.age65) - 1 == 1 &
          as.numeric(boot.d.star.stage) - 1 == 0 &
          boot.d.star$Volume <= q]
      obs.g.11 <- boot.d.star$Volume[
        as.numeric(boot.d.star.age65) - 1 == 1 &
          as.numeric(boot.d.star.stage) - 1 == 1 & 
          boot.d.star$Volume <= q]
      
      # f(g | x)
      est.isws.00 <- 1/length(obs.g.00)
      est.isws.01 <- 1/length(obs.g.01)
      est.isws.10 <- 1/length(obs.g.10)
      est.isws.11 <- 1/length(obs.g.11)
      
      # Empty Vectors
      fit.risk <- array(NA, dim = c(num.d, num.g, 4 + 1, length(mod.list)))
      fit.srisk <- lapply(1:length(mod.list), function(x) { NA })
      
      # Predicted Risk by DVH
      for(z in 1:length(mod.list)) {
        pred <- mod.preds.test[[z]]
        test_eval <- as_tibble(cbind(test, pred)) %>%
          dplyr::select(Dose, Volume, Stage, Age65, pred) %>%
          group_by(Dose, Volume) %>%
          mutate(`cpred` = (Age65 == 0)*(Stage == 0)*p.stage.emp[1,1]*pred +
                   (Age65 == 0)*(Stage == 1)*p.stage.emp[1,2]*pred +
                   (Age65 == 1)*(Stage == 0)*p.stage.emp[2,1]*pred +
                   (Age65 == 1)*(Stage == 1)*p.stage.emp[2,2]*pred)
        test_eval_00 <- test_eval %>% dplyr::filter(Age65 == 0, Stage == 0)
        test_eval_01 <- test_eval %>% dplyr::filter(Age65 == 0, Stage == 1)
        test_eval_10 <- test_eval %>% dplyr::filter(Age65 == 1, Stage == 0)
        test_eval_11 <- test_eval %>% dplyr::filter(Age65 == 1, Stage == 1)
        test_eval <- test_eval %>%
          summarise(cpred = sum(cpred)) %>% ungroup
        pred_mat <- list(
          `00` = matrix(test_eval_00$pred, length(dvals), length(gvals)),
          `01` = matrix(test_eval_01$pred, length(dvals), length(gvals)),
          `10` = matrix(test_eval_10$pred, length(dvals), length(gvals)),
          `11` = matrix(test_eval_11$pred, length(dvals), length(gvals))
        )
        cpred_mat <- matrix(test_eval$cpred, length(dvals), length(gvals), byrow = T)
        
        fit.risk[,,1,z] <- pred_mat$`00`
        fit.risk[,,2,z] <- pred_mat$`01`
        fit.risk[,,3,z] <- pred_mat$`10`
        fit.risk[,,4,z] <- pred_mat$`11`
        fit.risk[,,5,z] <- cpred_mat
        
        # Stochastic Risk (Conditional on X)
        cond.means <- data.frame(
          `g` = boot.long$Volume[which(boot.long[,"Dose"] == d.star)],
          `risk`= mod.preds.train[[z]][which(boot.long[,"Dose"] == d.star)],
          `x1` = boot.long.age65[which(boot.long[,"Dose"] == d.star)],
          `x2` = boot.long.stage[which(boot.long[,"Dose"] == d.star)]
        )
        
        # P(Y = 1 | d^*, g, x)
        est.risks.00 <- cond.means$risk[
          as.numeric(boot.d.star.age65) - 1 == 0 &
            as.numeric(boot.d.star.stage) - 1 == 0 &
            cond.means$g <= q]
        est.risks.01 <- cond.means$risk[
          as.numeric(boot.d.star.age65) - 1 == 0 &
            as.numeric(boot.d.star.stage) - 1 == 1 &
            cond.means$g <= q]
        est.risks.10 <- cond.means$risk[
          as.numeric(boot.d.star.age65) - 1 == 1 &
            as.numeric(boot.d.star.stage) - 1 == 0 &
            cond.means$g <= q]
        est.risks.11 <- cond.means$risk[
          as.numeric(boot.d.star.age65) - 1 == 1 &
            as.numeric(boot.d.star.stage) - 1 == 1 &
            cond.means$g <= q]
        
        # Estimated
        est.stochast.00 <- sum(est.risks.00*est.isws.00)
        est.stochast.01 <- sum(est.risks.01*est.isws.01)
        est.stochast.10 <- sum(est.risks.10*est.isws.10)
        est.stochast.11 <- sum(est.risks.11*est.isws.11)
        
        # Save Conditional Values
        stochast.save[[r]][[z]] = est.risks.00       # Conditional Means 
        stochast.save[[r]][[4 + z]] = est.risks.01      
        stochast.save[[r]][[8 + z]] = est.risks.10
        stochast.save[[r]][[12 + z]] = est.risks.11
        stochast.save[[r]][[16 + z]] = est.isws.00   # ISW
        stochast.save[[r]][[20 + z]] = est.isws.01   
        stochast.save[[r]][[24 + z]] = est.isws.10   
        stochast.save[[r]][[28 + z]] = est.isws.11   
        stochast.save[[r]][[29]] = boot.d.star.age65  # Covariates
        stochast.save[[r]][[30]] = boot.d.star.stage
        
        fit.srisk[[z]] <- sum(c(est.stochast.00,
                                est.stochast.01, est.stochast.10, est.stochast.11)*as.numeric(t(p.stage.emp)))
      }
      obs.risk <- prop.table(table(boot.wide$anyderm))[2]
      
      
      # Combine True Information
      mod.preds <- fit.risk
      # identical(mod.preds[,,1,5], unname(true.risk$`x = 0`))
      
      # Array Manipulation Toy
      # arr <- array(1:18, dim = c(3, 3, 2)); arr[,,1]; arr[,,2]
      
      # Fit Statistics
      p.fit <- mod.preds.train
      brier.fun <- function(p) {mean((p - boot.long$y)^2, na.rm = T)}
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
      
      # Save
      params.output[,,,,r] <- array(mod.preds, dim = c(num.d, num.g, 4 + 1, length(param.names), 1))
      params2.output[r,] <- nan.recode(unlist(fit.srisk))
      params3.output[r,] <- nan.recode(unname(obs.risk))
      eff.parms.output[r,] <- nan.recode(eff.parms)
      deviance.output[r,] <- nan.recode(deviance)
      brier.output[r,] <- nan.recode(brier)
      dic.output[r,] <- nan.recode(dic)
      
      # Data Save
      data.save[[r]] <- boot.long
      
      # cat(r,"\n")
    }
    
    save(params.output, file = file.path(outpath, paste0("boot_bivar.params.output_", core.no, ".RData")))
    save(params2.output, file = file.path(outpath, paste0("boot_bivar.params2.output_", core.no, ".RData")))
    save(params3.output, file = file.path(outpath, paste0("boot_bivar.params3.output_", core.no, ".RData")))
    save(brier.output, file = file.path(outpath, paste0("boot_bivar.brier.output_", core.no, ".RData")))
    save(deviance.output, file = file.path(outpath, paste0("boot_bivar.deviance.output_", core.no, ".RData")))
    save(eff.parms.output, file = file.path(outpath, paste0("boot_bivar.eff.parms.output_", core.no, ".RData")))
    save(dic.output, file = file.path(outpath, paste0("boot_bivar.dic.output_", core.no, ".RData")))
    
    save(data.save, file = file.path(outpath, paste0("boot_bivar.datalist.output_", core.no, ".RData")))
    save(stochast.save, file = file.path(outpath, paste0("boot_bivar.stochastests.output_", core.no, ".RData")))
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
  datasave.output_ALL <- vector("list", 0)
  stochast.save.output_ALL <- vector("list", 0)
  
  for (i in 1:num.cores) {
    load(file = file.path(outpath, paste0('boot_bivar.params.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.params2.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.params3.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.brier.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.eff.parms.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.deviance.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.dic.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.datalist.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('boot_bivar.stochastests.output_', i, ".RData")))
    
    params.output_ALL <- abind(params.output_ALL, params.output, along = 5)
    params.output2_ALL <- rbind(params.output2_ALL, params2.output)
    params.output3_ALL <- rbind(params.output3_ALL, params3.output)
    brier.output_ALL <- rbind(brier.output_ALL, brier.output)
    eff.parms.output_ALL <- rbind(eff.parms.output_ALL, eff.parms.output)
    deviance.output_ALL <- rbind(deviance.output_ALL, deviance.output)
    dic.output_ALL <- rbind(dic.output_ALL, dic.output)
    datasave.output_ALL <- c(datasave.output_ALL, data.save)
    stochast.save.output_ALL <- c(stochast.save.output_ALL, stochast.save)
  }
  
  list(`Predictions` = params.output_ALL,
       `Stochast Predictions` = params.output2_ALL,
       `Obs Means` = params.output3_ALL,
       `Brier` = brier.output_ALL,
       `Deviance` = deviance.output_ALL,
       `effNumParms` = eff.parms.output_ALL,
       `DIC` = dic.output_ALL,
       `Data` = datasave.output_ALL,
       `Conditional Stochastic Predictions` = stochast.save.output_ALL
  )
}


##########################################################


# Storage
outpath2 <- here("data", "storage")  # Resampling Storage Folder
if(!file.exists(outpath2)) dir.create(file.path(extension, "storage"))

# Pararellization

# Pointwise & Stochastic Causal Risks
# Start Time @ 2022-09-14 02:01:57 EDT - 1000 Bootstrap Resamples (15.3 hours)
(sim.begin <- Sys.time())
bootresamp <- meanfunc.boot(num.cores = 8, seeding = 1943)
(bootresamp$time <- Sys.time() - sim.begin)    # Total Bootstrap Resampling Time
# save(bootresamp, file = file.path(outpath2, paste0("pmh_bootresamp_causal_skin_age65_TWOCONFS_correctstochast.RData")))
