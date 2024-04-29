library(here)
res.dir <- here("main simulation", "results") # Results Directory
source(here("functions", "functions.R"))

# Packages
installer(c(c("tidyverse", "compiler"),
            c("gam", "akima", "fields", "abind",
              "monoreg", "pheatmap")))

# Simulation of Functional Monotonic Regression
wd <- here("data")
extension <- wd

# Generate Grid of True Risks
source(here("main simulation", "true_risk.R"))

meanfunc.sim <- function(num.cores = 8, seeding = 1) {
  
  # Parallel Processing Set-Up
  total <- 104
  R <- total %/% num.cores   # Simulations per Core
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
    bmax <- 2      
    amin <- 35     
    amax <- 45  
    bmin_x0 <- bmin
    bmax_x0 <- bmax
    beta.parms.x0 <- c(4, 6)/3
    beta.parms.x1 <- rev(beta.parms.x0)
    amin_upper <- amin
    amax_lower <- amax
    p.grade <- 0.4
    n <- 500 # n <- 100 # Change depending on scenario
    q <- 0.8
    d.star <- 40
    
    # Reference Dose and DVH Values
    num.d <- 26
    num.g <- 26
    resolution <- 0.001
    dvals <- seq(30, 50, length.out = num.d)
    gvals <- seq(0, 1, length.out = num.g)
    d.star <- dvals[findInterval(d.star, dvals)]
    g.seq <- seq(0, 1, by = 0.005)
    q.seq <- g.seq[g.seq <= q]
    
    # Functional Mean Model
    alpha <- -18;
    beta <- 0.45
    beta2 <- 0.5
    
    ## DGM FUNCTIONS ##
    dg.calc <- function(mu, sd, d) {
      dnorm(d, mean = mu, sd = sd)
    }
    g.calc <- function(mu, sd, d) {
      1 - pnorm(d, mean = mu, sd = sd)
    }
    fun.mod <- function(mu, x) { expit(alpha + beta*mu + beta2*x) }
    nan.recode <- function(x) {na_if(x, "NaN")}
    
    # Limits of Mu
    mufind <- function(sigma, dstar, g, muroot) {
      dvhfun <- function(mu) {1 - pnorm(dstar, mu, sigma) - g}
      if(class(suppressWarnings(try(uniroot(dvhfun, lower = muroot[1], upper = muroot[2])$root, silent = T))) != "try-error") {
        return(uniroot(dvhfun, lower = muroot[1], upper = muroot[2])$root)
      } else {
        return(NA)
      }
    }
    
    # Limits of Sigma
    sigmafind <- function(mu, dstar, g, sigroot) {
      dvhfun <- function(sigma) {1 - pnorm(dstar, mu, sigma) - g}
      if(class(suppressWarnings(try(uniroot(dvhfun, lower = sigroot[1] - 0.01, upper = sigroot[2] + 0.01)$root, silent = T))) != "try-error") {
        ans <- uniroot(dvhfun, lower = sigroot[1] - 0.01, upper = sigroot[2] + 0.01)$root
        if(ans < sigroot[1]) {ans <- sigroot[1]}
        if(ans > sigroot[2]) {ans <- sigroot[2]}
        return(ans)
      } else {
        return(NA)
      }
    }
    
    # Finding Valid Values of Mu
    sig.lims <- function(g, dstar, x, prnt = F) {
      muroot <- c(amin, amax_lower)                    # x = 0
      sigroot <- c(bmin_x0, bmax_x0)
      if(x == 1) {
        muroot <- c(amin_upper, amax)
        sigroot <- c(bmin, bmax)
      }       # x = 1
      
      mu.sols <- sort(sapply(sigroot, function(sigval) {mufind(sigma = sigval, dstar = dstar, g = g, muroot = muroot)}), na.last = T)
      bound.sols <- ifelse(!is.na(sapply(muroot, function(muval) {sigmafind(mu = muval, dstar = dstar, g = g, sigroot = sigroot)})), muroot, NA)
      if(prnt) {contour(x = muvals, y = sigvals, z = dgmat(dstar))}
      if(sum(!is.na(c(mu.sols, bound.sols))) == 0) {
        return(c(NA, NA))
      } else {
        return(range(c(mu.sols, bound.sols), na.rm = T))
      }
    }
    
    # Kernel Functions for Mu (Calculating Risk)
    mu.expit.kernel <- function(mu, dstar, g, x) {
      beta.parms <- beta.parms.x0	
      if(x == 1) {beta.parms <- beta.parms.x1}	
      alpha.val <- beta.parms[1]; beta.val <- beta.parms[2]
      sigroot <- c(bmin, bmax)
      sig.val <- sigmafind(mu, dstar = dstar, g = g, sigroot = sigroot)
      if(is.na(sig.val)) {return(NA)}
      norm.fun <- function(d, mu) {
        dnorm((d - mu)/(sig.val))*(1 - ((d - mu)/(sig.val))^2)
      }
      J <- (sig.val)^2/abs(integrate(Vectorize(norm.fun), 0, dstar, mu = mu)$value)
      # if(J > 10^15) {return(NA)}
      fun.mod(mu = mu, x = x)*((mu - amin)^(alpha.val - 1))*((amax - mu)^(beta.val - 1))*J
    }
    mu.kernel <- function(mu, dstar, g, x) {
      beta.parms <- beta.parms.x0	
      if(x == 1) {beta.parms <- beta.parms.x1}	
      alpha.val <- beta.parms[1]; beta.val <- beta.parms[2]
      sigroot <- c(bmin, bmax)
      sig.val <- sigmafind(mu, dstar = dstar, g = g, sigroot = sigroot)
      if(is.na(sig.val)) {return(NA)}
      norm.fun <- function(d, mu) {
        dnorm((d - mu)/(sig.val))*(1 - ((d - mu)/(sig.val))^2)
      }
      J <- (sig.val)^2/abs(integrate(Vectorize(norm.fun), 0, dstar, mu = mu)$value)
      # if(J > 10^15) {return(NA)}
      ((mu - amin)^(alpha.val - 1))*((amax - mu)^(beta.val - 1))*J
    }
    
    # mu.lims <- sapply(g.grid, function(g) {
    #   sig.lims(g = g, dstar = dstar)
    # })
    # g.ll <- splinefun(g.grid, mu.lims[1,])
    # g.ul <- splinefun(g.grid, mu.lims[2,])
    
    # Kernel Functions for g (Calculating Stochastic Risk)
    g.dens.make <- function(g, dstar, x) {
      if(g %in% c(0, 0.5, 1)) {return(NA)}
      
      mu.range <- sig.lims(g = g, dstar = dstar, x = x)
      
      if(sum(is.na(mu.range)) > 0) {return(NA)}
      mu.min <- mu.range[1]
      mu.max <- mu.range[2]
      
      integrate(Vectorize(mu.kernel), lower = mu.min, upper = mu.max,
                dstar = dstar, g = g, x = x)$value/4/(60^(2 + 4/3 - 1))/(gamma(2)*gamma(4/3)/gamma(2 + 4/3))
    }
    g.expit.kernel.make <- function(g, dstar, x) {
      if(g %in% c(0, 0.5, 1)) {return(NA)}
      mu.range <- sig.lims(g, dstar = dstar, x = x)
      if(sum(is.na(mu.range)) > 0) {return(NA)}
      mu.min <- mu.range[1]
      mu.max <- mu.range[2]
      
      integrate(Vectorize(mu.expit.kernel), lower = mu.min, upper = mu.max,
                dstar = dstar, g = g, x = x)$value/4/(60^(2 + 4/3 - 1))/(gamma(2)*gamma(4/3)/gamma(2 + 4/3))
    }
    
    # Posthoc Adjustment of Density Functions of g
    posthoc.adjust <- function(g.func.vals, g.s = g.seq) {
      g.func.vals[which(is.na(g.func.vals))] <- 0
      critvals <- which(g.s %in% c(0, 0.5, 1))
      if(0 %in% g.s) {
        g.func.vals[which(g.s == 0)] <- (spline(g.s[-critvals], g.func.vals[-critvals], xout = g.s, method = "natural")$y)[which(g.s == 0)]
      }
      
      # g = 0.5
      if(0.5 %in% g.s) {
        midvals <- which((g.s > 0.48) & (g.s < 0.52))
        g.func.vals[midvals] <- (spline(g.s[-sort(c(midvals, critvals))], g.func.vals[-sort(c(midvals, critvals))], xout = g.s, method = "natural")$y)[midvals]
      }
      
      # g = 1
      if(1 %in% g.s) {
        g.func.vals[which(g.s == 1)] <- (spline(g.s[-critvals], g.func.vals[-critvals], xout = g.s, method = "natural")$y)[which(g.s == 1)]
      }
      
      # Rescaling g = 0 and g = 1
      tf <- splinefun(g.seq[-which(g.seq %in% c(0, 0.5, 1))],
                      g.func.vals[-which(g.seq %in% c(0, 0.5, 1))],
                      method = "natural")
      tf.fin <- function(y) {tf(y)/integrate(Vectorize(tf), 0, 1)$value}
      
      # integrate(tf.fin, 0, 1)$value
      return(list(
        `fun` = tf.fin,
        `SF` = integrate(Vectorize(tf), 0, 1)$value
      ))
    }
    
    posthoc.scale <- function(g.func.vals, scale, g.s = g.seq, x) {
      g.func.vals[which(is.na(g.func.vals))] <- 0
      critvals <- which(g.s %in% c(0, 0.5, 1))
      dens <- dens.0
      if(x == 1) dens <- dens.1
      if(0 %in% g.s) {
        g.func.vals[which(g.s == 0)] <- (spline(g.s[-critvals], g.func.vals[-critvals], xout = g.s, method = "natural")$y)[which(g.s == 0)]
      }
      
      # g = 0.5
      if(0.5 %in% g.s) {
        midvals <- which((g.s > 0.48) & (g.s < 0.52))
        g.func.vals[midvals] <- (spline(g.s[-sort(c(midvals, critvals))], g.func.vals[-sort(c(midvals, critvals))], xout = g.s, method = "natural")$y)[midvals]
      }
      
      # g = 1
      if(1 %in% g.s) {
        g.func.vals[which(g.s == 1)] <- (spline(g.s[-critvals], g.func.vals[-critvals], xout = g.s, method = "natural")$y)[which(g.s == 1)]
      }
      
      # Rescaling g = 0 and g = 1
      tf <- splinefun(g.seq[-which(g.seq %in% c(0, 0.5, 1))],
                      g.func.vals[-which(g.seq %in% c(0, 0.5, 1))],
                      method = "natural")
      tf.fin <- function(y) {tf(y)/scale}
      return(tf.fin)
    }
    
    # True Risk Function
    risk.fun <- function(dstar, g, x) {
      parms <- c(amin, amax_lower)
      if(x == 1) {parms <- c(amin_upper, amax)}
      
      mu.range <- sig.lims(g, dstar = dstar, x = x)
      if(sum(is.na(mu.range)) > 0) {return(NA)}
      
      mu.min <- mu.range[1]
      mu.max <- mu.range[2]
      
      # Set as NA if Equation Cannot Be Solved
      if(is.na(mu.kernel(mu.min, dstar = dstar, g = g, x = x)) | is.na(mu.kernel(mu.max, dstar = dstar, g = g, x = x)) |
         is.infinite(mu.kernel(mu.min, dstar = dstar, g = g, x = x)) | is.infinite(mu.kernel(mu.max, dstar = dstar, g = g, x = x)) |
         sum(c(mu.kernel(mu.min, dstar = dstar, g = g, x = x), mu.kernel(mu.max, dstar = dstar, g = g, x = x)) > 10000) > 1 |
         is.na(mu.expit.kernel(mu.min, dstar = dstar, g = g, x = x)) | is.na(mu.expit.kernel(mu.max, dstar = dstar, g = g, x = x)) |
         is.infinite(mu.expit.kernel(mu.min, dstar = dstar, g = g, x = x)) | is.infinite(mu.expit.kernel(mu.max, dstar = dstar, g = g, x = x)) |
         sum(c(mu.expit.kernel(mu.min, dstar = dstar, g = g, x = x), mu.expit.kernel(mu.max, dstar = dstar, g = g, x = x)) > 10000) > 1
      ) {return(NA)}
      if(exists("a")) {rm(a)}
      if(exists("b")) {rm(b)}
      a <- tryCatch(
        integrate(Vectorize(mu.expit.kernel), mu.min, mu.max, dstar = dstar, g = g, x = x,
                       stop.on.error = F)$value,
        error = function(e) NA
      )
      b <- tryCatch(
        integrate(Vectorize(mu.kernel), mu.min, mu.max, dstar = dstar, g = g, x = x,
                  stop.on.error = F)$value,
        error = function(e) NA
      )
      # if(!exists("a") | !exists("b")) {return(NA)}
      a/b
    }
    
    # True Pointwise Causal Risk
    crisk_point.fun <- function(dstar, g) {
      sum(sapply(0:1, function(x) {nan.recode(risk.fun(g = g, dstar = dstar, x = x))})*c(1 - p.grade, p.grade))
    }
    
    # True Stochastic Causal Risk
    scrisk.fun <- function(dstar, x) {
      parms <- c(amin, amax_lower)
      if(x == 1) {parms <- c(amin_upper, amax)}
      
      # Kernel for g
      dens.var <- posthoc.adjust(sapply(g.seq, function(g) {g.dens.make(g = g, dstar = dstar, x = x)}))
      g.kernel <- dens.var$fun
      
      # Expit Kernel for g
      g.expit.kernel <- posthoc.scale(sapply(g.seq, function(g) {g.expit.kernel.make(g = g, dstar = dstar, x = x)}),
                                      scale = dens.var$SF, x = x)
      
      out <- integrate(Vectorize(g.expit.kernel), 0, q)$value/integrate(Vectorize(g.kernel), 0, q)$value
      if(integrate(Vectorize(g.kernel), 0, q)$value == 0) {return(0)}
      return(out)
    }
    
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
    epmf <- function(orig.vec) {
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
                     "Additive (Doubly)", "Bivariable",
                     "True")
    params.output <- array(NA, dim = c(num.d, num.g, 2 + 1, length(param.names), R))
    params2.output <- matrix(NA, nrow = R, ncol = length(param.names))
    deviance.output <- matrix(NA, nrow = R, ncol = length(param.names))
    colnames(deviance.output) <- param.names
    brier.output <- deviance.output 
    eff.parms.output <- dic.output <- deviance.output[,-5]
  
    
    from <- seq(1, total.sim, by = R)[core.no]
    to <- seq(R, total.sim, by = R)[core.no]
    
    for(r in 1:R) {
      
      iter <- from + (r - 1)
      set.seed(iter + seeding)
      
      x <- rbinom(n, 1, p.grade)
      wide.data <- data.frame(
        id = 1:n,
        x = x,
        mus = amin + (amax - amin)*(rbeta(n, beta.parms.x1[1], beta.parms.x1[2])*(x == 1) +	
                                      rbeta(n, beta.parms.x0[1], beta.parms.x0[2])*(x == 0)),      # Simulated Means,
        sds = runif(n, bmin, bmax)*(x == 1) + runif(n, bmin_x0, bmax_x0)*(x == 0)                # Simulated Standard Deviations
      )
      wide.data$y <- rbinom(n, 1, prob = fun.mod(mu = wide.data$mus, x = wide.data$x))
      
      long.data <- data.frame(
        id = rep(1:length(wide.data$mus), each = length(dvals)),
        x = rep(wide.data$x, each = length(dvals)),
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
      d.star.data <- long.data[which(long.data$dose == d.star),]
      
      
      ## FITTED MODELS ##
      
      # Tumor Stage (Empirical)
      p.stage.emp <- as.numeric(prop.table(table(wide.data$x))); names(p.stage.emp) <- c("x = 0", "x = 1")
      
      # Initialization
      nobs <- nrow(long.data)
      x1 <- (long.data$dose - min(dvals))/(diff(range(dvals)))
      x2 <- long.data$dvh
      grid.x1 <- (dvals - min(dvals))/diff(range(dvals))
      grid.x2 <- gvals
      grid.dat <- expand.grid(`Dose` = grid.x1, `Volume` = grid.x2, `Stage` = c(0, 1))
      
      # Training and Testing
      train <- data.matrix(cbind(rep(1, nobs), x1, x2, x1^2, x2^2, x1^3, x2^3, long.data$x))
      test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat[,-3], grid.dat[,-3]^2, grid.dat[,-3]^3, grid.dat[,3]))
      alldat <- rbind(train, test)
      colnames(train) <- colnames(test) <- colnames(alldat) <- c("Intercept",
                                                                 "Dose", "Volume", "Dose-Squared", "Volume-Squared", "Dose-Cubed",
                                                                 "Volume-Cubed", "Stage")
      
      # Covariates
      covs.additive <- list(`train` = train, `test` = test, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1, 8))
      covs.noadditive <- list(`train` = train, `test` = test, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1:3, 8))
      covs.noadditive.poly <- list(`train` = train, `test` = test, y = long.data$y, nobs = nobs, add.id = 2:3, covs.id = c(1:5, 8))
      
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
      # str(mod.preds.test)
      
      
      # Empirical Estimators
      
      # Assignment Model (Beta Inflated Regression)
      # exp.mod <- gamlss(dvh ~ x, data = d.star.data, family = BEINF)
      # predict(exp.mod, what = "mu", type = "response")
      
      # ePMF Functions
      fx.0 <- epmf(d.star.data$dvh[which(d.star.data$x == 0)])
      fx.1 <- epmf(d.star.data$dvh[which(d.star.data$x == 1)])
      
      # Empirical CDF at g = q
      Fq.0 <- mean(d.star.data$dvh[which(d.star.data$x == 0)] <= q)
      Fq.1 <- mean(d.star.data$dvh[which(d.star.data$x == 1)] <= q)
      
      # Empirical DVHs
      obs.g.0 <- d.star.data$dvh[d.star.data$x == 0 & d.star.data$dvh <= q]
      obs.g.1 <- d.star.data$dvh[d.star.data$x == 1 & d.star.data$dvh <= q]
      
      # f(g | x)
      est.isws.0 <- 1/length(obs.g.0)
      est.isws.1 <- 1/length(obs.g.1)
      
      # Empty Vectors
      true.crisk <- matrix(NA, nrow = num.d, ncol = num.g)
      true.risk <- lapply(1:2, function(x) {
        nullmat <- matrix(NA, nrow = num.d, ncol = num.g)
        rownames(nullmat) <- dvals
        colnames(nullmat) <- gvals
        return(nullmat)
      }); names(true.risk) <- c("x = 0", "x = 1")
      true.crisk <- matrix(NA, num.d, num.g)
      true.scrisk <- NA
      true.srisk <- rep(NA, 2); names(true.srisk) <- c("x = 0", "x = 1")
      
      fit.risk <- array(NA, dim = c(num.d, num.g, 2 + 1, length(mod.list)))
      fit.srisk <- lapply(1:length(mod.list), function(x) { NA })
      
      # Predicted Risk by DVH
      for(z in 1:length(mod.list)) {
        pred <- mod.preds.test[[z]]
        test_eval <- as_tibble(cbind(test, pred)) %>% dplyr::select(Dose, Volume, Stage, pred) %>%
          group_by(Dose, Volume) %>%
          mutate(`cpred` = (Stage == 0)*p.stage.emp[1]*pred + (Stage == 1)*p.stage.emp[2]*pred)
        test_eval_0 <- test_eval %>% dplyr::filter(Stage == 0)
        test_eval_1 <- test_eval %>% dplyr::filter(Stage == 1)
        test_eval <- test_eval %>%
          summarise(cpred = sum(cpred)) %>% ungroup
        pred_mat <- list(
          `x = 0` = matrix(test_eval_0$pred, num.d, num.g),
          `x = 1` = matrix(test_eval_1$pred, num.d, num.g)
        )
        cpred_mat <- matrix(test_eval$cpred, num.d, num.g, byrow = T)
        
        fit.risk[,,1,z] <- pred_mat$`x = 0`
        fit.risk[,,2,z] <- pred_mat$`x = 1`
        fit.risk[,,3,z] <- cpred_mat
        
        # Stochastic Risk (Conditional on X)
        cond.means <- data.frame(
          `g` = long.data$dvh[which(long.data[,"dose"] == d.star)],
          `risk`= mod.preds.train[[z]][which(long.data[,"dose"] == d.star)],
          `x` = long.data$x[which(long.data[,"dose"] == d.star)]
        )
        
        # P(Y = 1 | d^*, g, x)
        est.risks.0 <- cond.means$risk[d.star.data$x == 0 & cond.means$g <= q]
        est.risks.1 <- cond.means$risk[d.star.data$x == 1 & cond.means$g <= q]
        
        # Estimated
        est.stochast.0 <- sum(est.risks.0*est.isws.0)
        est.stochast.1 <- sum(est.risks.1*est.isws.1)
        
        fit.srisk[[z]] <- mean(c(est.stochast.0, est.stochast.1)[d.star.data$x + 1])
      }
      
      # Combine
      true.risk_arr <- array(NA, dim = c(num.d, num.g, 2 + 1))
      true.risk_arr[,,1] <- true.risk$`x = 0`
      true.risk_arr[,,2] <- true.risk$`x = 1`
      true.risk_arr[,,3] <- true.crisk
      
      # Combine True Information
      mod.preds <- abind(fit.risk, array(true.risk_arr, dim = c(num.d, num.g, 3)), along = 4)
      # identical(mod.preds[,,1,5], unname(true.risk$`x = 0`))
      
      # True Probabilities Based on Training Data
      subdat.true.risk <- rep(NA, nrow(long.data))
      for(i in 1:nrow(long.data)) {
        row <- long.data[i,]
        cat("Row:", i, "\n")
        subdat.true.risk[i] <- risk.fun(g = row$dvh, dstar = row$dose, x = row$x)
      }
      subdat.true.risk <- nan.recode(subdat.true.risk)
      
      # Array Manipulation Toy
      # arr <- array(1:18, dim = c(3, 3, 2)); arr[,,1]; arr[,,2]
      
      # Fit Statistics
      p.fit <- mod.preds.train
      brier.fun <- function(p) {mean((p - long.data$y)^2, na.rm = T)}
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
      params.output[,,,,r] <- array(mod.preds, dim = c(num.d, num.g, 3, length(param.names), 1))
      params2.output[r,] <- nan.recode(c(unlist(fit.srisk), true.scrisk))
      eff.parms.output[r,] <- nan.recode(eff.parms)
      deviance.output[r,] <- nan.recode(c(deviance, 2*(-1*sum(log(subdat.true.risk), na.rm = T))))
      brier.output[r,] <- nan.recode(c(brier, brier.fun(subdat.true.risk)))
      dic.output[r,] <- nan.recode(dic)
      
      # cat(r,"\n")
    }
    
    save(params.output, file = file.path(outpath, paste0("bivar.params.output_", core.no, ".RData")))
    save(params2.output, file = file.path(outpath, paste0("bivar.params2.output_", core.no, ".RData")))
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
  params.output2_ALL <- NULL
  brier.output_ALL <- NULL
  eff.parms.output_ALL <- NULL
  deviance.output_ALL <- NULL
  dic.output_ALL <- NULL
  
  for (i in 1:num.cores) {
    load(file = file.path(outpath, paste0('bivar.params.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.params2.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.brier.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.eff.parms.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.deviance.output_', i, ".RData")))
    load(file = file.path(outpath, paste0('bivar.dic.output_', i, ".RData")))
    params.output_ALL <- abind(params.output_ALL, params.output, along = 5)
    params.output2_ALL <- rbind(params.output2_ALL, params2.output)
    brier.output_ALL <- rbind(brier.output_ALL, brier.output)
    eff.parms.output_ALL <- rbind(eff.parms.output_ALL, eff.parms.output)
    deviance.output_ALL <- rbind(deviance.output_ALL, deviance.output)
    dic.output_ALL <- rbind(dic.output_ALL, dic.output)
    
  }
  
  list(`Predictions` = params.output_ALL,
       `Stochast Predictions` = params.output2_ALL,
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

# Pointwise & Stochastic Causal Risks
# Logistic, Polynomial Logistic, Additive (Doubly), Bivariable (200 MCMC Iterations, n = 100)
(sim.begin <- Sys.time())
simulation <- meanfunc.sim(num.cores = 8, seeding = 57)
(simulation$time <- Sys.time() - sim.begin)    # Total Simulation Time
save(simulation, "simulation_n=500")