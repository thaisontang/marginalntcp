### UNPACKING SIMULATION FILE ###
wd <- here("data")
setwd(wd)
outpath2 <- here("data", "storage")
res.dir <- here("supplementary simulation", "results")

# Unpack Files of Interest
file.names <- c("supp_simulation_strongconf.RData")
simulation.list <- lapply(1:length(file.names), function(id) {
  load(file = file.path(outpath2, file.names[[id]]))
  simulation
}); names(simulation.list) <- file.names
a.obj <- simulation.list$`supp_simulation_strongconf.RData`$Predictions

# Separate as Misspecified + Correctly Specified Models
sim_bivar <- simulation
true.rep.arr <- array(NA, dim = c(grid.num.d,grid.num.g, 1, nrow(sim_bivar$Brier)))
true.rep.arr[,,,1:nrow(sim_bivar$Brier)] <- true.risk
correct.idx <- which(grepl("Adjusted", colnames(sim_bivar$Brier), fixed = T))
misspec.idx <- which(grepl("Unadjusted", colnames(sim_bivar$Brier), fixed = T))
sim_bivar_correct <- list(
  `Predictions` = abind(sim_bivar$Predictions[,,correct.idx,], true.rep.arr, along = 3),
  `Brier` = sim_bivar$Brier[,c(correct.idx)],
  `Deviance` = sim_bivar$Deviance[,c(correct.idx)],
  `effNumParms` = sim_bivar$effNumParms[,correct.idx],
  `DIC` = sim_bivar$DIC[,correct.idx]
)
sim_bivar_misspec <- list(
  `Predictions` = abind(sim_bivar$Predictions[,,misspec.idx,], true.rep.arr, along = 3),
  `Brier` = sim_bivar$Brier[,c(misspec.idx)],
  `Deviance` = sim_bivar$Deviance[,c(misspec.idx)],
  `effNumParms` = sim_bivar$effNumParms[,misspec.idx],
  `DIC` = sim_bivar$DIC[,misspec.idx]
)

source(here("analysis code", "functions.R"))
installer(c("tidyverse", "plot3D", "colorspace", "rasterpdf"))

### PROCESSING FOR VISUALIZATION ###
sci.exponent <- function(x) {
  floor(log(x)/log(10))
}

latexprint <- function(y, col.include = F, row.include = F, dig = NULL, as.char = T, latex = T){
  library(xtable)
  options(scipen=999)   # Remove Scientific Notation
  # Change to character, keeping formatting
  if(as.char != F) x = cbind(y[,1], sapply(1:ncol(y[,-1]),
                                           function(z) formatC(y[,z+1], format = "f", digits = 10))) #
  #if(as.char != F) x = cbind(sapply(1:ncol(y), function(z) formatC(y[,z], format = "f")))#, digits = dig)))
  rownames(x) <- rownames(y); colnames(x) <- colnames(y)
  if(col.include != F) x = rbind(colnames(x), x)
  if(row.include != F) x = cbind(rownames(x), x)
  if(row.include == T & col.include == T) x[1,1] = ""
  options(scipen=3)     
  # write(print(xtable(x, digits = dig), comment = F, tabular.environment = "longtable",
  #             floating = F, include.rownames = F, include.colnames = F),
  # #       file = paste0(extension,"Code/Output/RandomPrint.txt"))
  if(latex == F) return(x)
  if(latex == T) return(print(xtable(x, digits = dig), comment = F, tabular.environment = "longtable",
                              floating = F, include.rownames = F, include.colnames = F))
}

# Plotting Functions
persp.plot <- function(zvec, xvec, yvec, ztrue, title) {
  
  # True Matrix Preprocessing
  ind.omit1 <- which(is.na(ztrue), arr.ind = T)
  ind.omit2 <- which(is.na(zvec), arr.ind = T)
  ex.vec1 <- which(is.na(ztrue))
  ex.vec2 <- which(is.na(zvec))
  
  xvec.sub1 <- xvec[-ex.vec1]
  yvec.sub1 <- yvec[-ex.vec1]
  xvec.sub2 <- xvec[-ex.vec2]
  yvec.sub2 <- yvec[-ex.vec2]
  zvec.sub <- as.numeric(zvec)[-ex.vec2]
  ztrue.sub <- as.numeric(ztrue)[-ex.vec1]
  if(length(ex.vec1) == 0) {
    xvec.sub1 <- xvec
    yvec.sub1 <- yvec
    ztrue.sub <- as.numeric(ztrue)
  }
  if(length(ex.vec2) == 0) {
    xvec.sub2 <- xvec
    yvec.sub2 <- yvec
    zvec.sub <- as.numeric(zvec)
  }
  
  zmat <- interp::interp(x = xvec.sub2, y = yvec.sub2, z = zvec.sub,
                         ny = nrow(zvec), nx = nrow(zvec))
  zmat.true <- interp::interp(x = xvec.sub1, y = yvec.sub1, z = ztrue.sub,
                              ny = nrow(zvec), nx = nrow(zvec))
  
  upper.flag <- (zvec > ztrue) + 1
  upper.flag <- upper.flag[-nrow(upper.flag),]  # Remove Last Row
  # pheatmap(upper.flag, cluster_rows = F, cluster_cols = F)
  upper.flag <- as.numeric(upper.flag)
  # pheatmap(matrix(upper.flag, 50 - 1, 50), cluster_rows = F, cluster_cols = F)
  
  pdf(file = here(res.dir, title), width = 6, height = 5, paper = "special")
  op <- par(mar=c(2,2,0,0), oma=c(0,0,0,0), mgp=c(2.5,1,0), cex=0.75)
  plot3D::persp3D(zmat$x, zmat$y, zmat$z, col = alpha("orange", 0.3), expand = 0.5,
                  zlim=c(0,1), xlim=c(0,1), ylim=c(0,1), lit = F,
                  contour = list(col = "orange", side = c("zmin")),
                  ticktype='detailed', theta=-45, phi=25, ltheta=25,
                  xlab='Dose (Scaled)', ylab='Volume', zlab='Risk')
  plot3D::persp3D(zmat$x, zmat$y, zmat.true$z, col = alpha("blue", 0.15), expand = 0.5,
                  contour = list(col = "blue", side = c("zmin")),        
                  zlim=c(0,1), xlim=c(0,1), ylim=c(0,1), lit = F, add = T)
  par(op)
  dev.off()
}

mat_pairs <- function(mat, xvals, yvals) {
  grid.vals <- expand.grid(xvals, yvals)
  data.fit2 <- data.frame(xvals = grid.vals$Var1, yvals = grid.vals$Var2, `zval` = as.vector(mat))
  data.fit2[complete.cases(data.fit2),]
}

# MC Results Simulation
replicates <- simulation$Predictions
no.dec = T
latex = F


############### VISUALIZATIONS ##################

#### CORRECT SCENARIO ####

# Mean Curve
sim_B <- nrow(sim_bivar_correct$Brier)
grid.num.d <- 12
grid.num.g <- 12

min.range <- 30
max.range <- 50

dvals.grid <- seq(min.range, max.range, length.out = grid.num.d + 2)[-c(1, grid.num.d + 2)]
gvals.grid <- seq(0, 1, length.out = grid.num.g + 2)[-c(1, grid.num.g + 2)]
param.names <- c("Logistic", "Polynomial Logistic",
                 "Additive (Doubly)", "Bivariable")
stats.labels <- c("Total Absolute Bias", "Total MCSD",
                  "Total rMSE", "MCE")
# Empty Storage Vectors

summ.results_test <- matrix(NA, nrow = length(param.names),
                            ncol = length(stats.labels))
colnames(summ.results_test) <- stats.labels; rownames(summ.results_test) <- param.names

# Testing Set
n <- 500
nobs <- grid.num.d*n
grid.x1 <- (dvals.grid - min.range)/(max.range - min.range)
grid.x2 <- gvals.grid
grid.dat <- expand.grid(data.frame(`Dose` = grid.x1, `Volume` = grid.x2))
test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat))
colnames(test) <- c("Intercept", "Dose", "Volume")

listLabel <- c("PointwiseRisk")
grid.num.pars <- length(param.names)
mean.mat <- lapply(1, function(x) {
  nullmat <- matrix(NA, nrow = grid.num.d, ncol = grid.num.g)
  rownames(nullmat) <- dvals.grid
  colnames(nullmat) <- gvals.grid
  return(nullmat)
}); names(mean.mat) <- listLabel
true.mean.mat <- true.ul.mat <- true.ll.mat <- mean.mat
true.sd.mat <- true.rmse.mat <- true.mce.mat <- mean.mat
ul.mat <- ll.mat <- mean.mat
sd.mat <- mean.mat

for(d in 1:grid.num.d) {
  for(g in 1:grid.num.g) {
    for(k in 1) {
      if(k == 1) {
        true.mean.mat[[k]][d,g] <- mean(sim_bivar_correct$Predictions[d,g,5,])
      }
      true.sd.mat[[k]][d,g] <- 0
      true.ul.mat[[k]][d,g] <- 0
      true.ll.mat[[k]][d,g] <- 0
      # true.rmse.mat[[k]][d,g] <- true.sd.mat[[k]][d,g]
      # true.mce.mat[[k]][d,g] <- true.sd.mat[[k]][d,g]/sqrt(sim_B)
    }
  }
}
true.rmse.mat <- true.sd.mat
true.mce.mat <- lapply(1, function(l) {
  true.sd.mat[[l]]/sqrt(sim_B)
}); names(true.mce.mat) <- listLabel
summ.results_test[length(param.names),] <- c(
  0,                                 # True Bias
  mean(true.sd.mat$PointwiseRisk, na.rm = T),      # True MCSD
  mean(true.rmse.mat$PointwiseRisk, na.rm = T),    # True rMSE
  mean(true.mce.mat$PointwiseRisk, na.rm = T)     # True MCE
)

for(j in 1:grid.num.pars) {
  
  for(d in 1:grid.num.d) {
    for(g in 1:grid.num.g) {
      for(k in 1) {
        if(sum(!is.na(sim_bivar_correct$Predictions[d,g,j,])) > 0) {
          mean.mat[[k]][d,g] <- mean(sim_bivar_correct$Predictions[d,g,j,], na.rm = T)
          # mat.bias[d,g] <- mean(sim_bivar_correct$Predictions[d,g,j,], na.rm = T) - true.mean.mat[[k]][d,g]
          sd.mat[[k]][d,g] <- sd(sim_bivar_correct$Predictions[d,g,j,], na.rm = T)
          ul.mat[[k]][d,g] <- quantile(sim_bivar_correct$Predictions[d,g,j,], probs = 0.975, na.rm = T)
          ll.mat[[k]][d,g] <- quantile(sim_bivar_correct$Predictions[d,g,j,], probs = 0.025, na.rm = T)
          # rmse.mat[[k]][d,g] <- sqrt(bias.mat[[k]][d,g]^2 + sd.mat[[k]][d,g]^2)
          # mce.mat[[k]][d,g] <- sd.mat[[k]][d,g]/sqrt(sim_B)
        } else {
          mean.mat[[k]][d,g] <- NA
          # mat.bias[d,g] <- NA
          ul.mat[[k]][d,g] <- NA
          ll.mat[[k]][d,g] <- NA
          # rmse.mat[[k]][d,g] <- NA
          # mce.mat[[k]][d,g] <- NA
        }
      }
    }
  }
  mat.bias <- lapply(1, function(l) {
    mean.mat[[l]] - true.mean.mat[[l]]
  }); names(mat.bias) <- listLabel
  rmse.mat <- lapply(1, function(l) {
    sqrt(mat.bias[[l]]^2 + sd.mat[[l]]^2)
  }); names(rmse.mat) <- listLabel
  mce.mat <- lapply(1, function(l) {
    sd.mat[[l]]/sqrt(sim_B)
  }); names(mce.mat) <- listLabel
  
  # 2D Contour Plots
  
  # Mean Estimates (Causal Risk)
  pdf(here(res.dir, paste0("Mean vs. True (Causal) - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 6), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  contour(x = dvals.grid, y = gvals.grid, z = true.mean.mat$PointwiseRisk, col = "blue", xlim = range(dvals.grid), lwd = 1.5,
          xlab = "Dose (Gy)", ylab = "Volume")
  contour(x = dvals.grid, y = gvals.grid, z = mean.mat$PointwiseRisk, nlevels = 10, col = "orange", add = T, lwd = 1.5)
  abline(h = c(0, 1), v = range(dvals.grid), col = "gray")
  legend("topright", lty = c(1,1), c("Truth", "Fit"), col = c("blue", "orange"),
         inset = c(-0.13, 0.4), xpd = T)
  par(op)
  rasterpdf::dev.off()
  
  # 3D Contour Plots
  
  # Mean vs. True
  added_lab <- param.names[j]
  persp.plot(xvec = test[,"Dose"],
             yvec = test[,"Volume"],
             zvec = mean.mat$PointwiseRisk,
             ztrue = true.mean.mat$PointwiseRisk,
             title = paste0("Mean vs. True (Causal) - 3D Perspective (Iterations = ", nrow(sim_bivar_correct$Brier), ") - ", added_lab,".pdf"))
  
  # Bias - 2D Contour Plots
  pdf(here(res.dir, paste0("Bias - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  # par(xpd = T)
  resolution <- 0.001
  bias.mat <- mat_pairs(mat.bias$PointwiseRisk, xvals = dvals.grid, yvals = gvals.grid)
  leg.range <- max(abs(c(floor(min(bias.mat$zval)*10)/10,
                         ceiling(max(bias.mat$zval)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals.grid, y=gvals.grid, z=mat.bias$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorspace::diverging_hsv(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(-leg.range, leg.range), xlim = c(min(dvals.grid), max(dvals.grid)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals.grid, y = gvals.grid, z = mat.bias$PointwiseRisk, add = T, labcex = 1.2)
                 })
  mtext("Bias", side = 4, line = 2.5, las = 0)
  
  par(op)
  dev.off()
  
  # Monte Carlo Standard Deviation
  pdf(here(res.dir, paste0("MCSD - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  resolution <- 0.001
  leg.range <- max(abs(c(floor(min(sd.mat$PointwiseRisk, na.rm = T)*10)/10,
                         ceiling(max(sd.mat$PointwiseRisk, na.rm = T)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals.grid, y=gvals.grid, z=sd.mat$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorRampPalette(c("white", "red"))(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(0, leg.range), xlim = c(min(dvals.grid), max(dvals.grid)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals.grid, y = gvals.grid, z = sd.mat$PointwiseRisk, add = T, labcex = 1.2)
                 })
  mtext("MCSD", side = 4, line = 2.5, las = 0)
  
  par(op)
  dev.off()
  
  # Root Mean-Squared Error
  pdf(here(res.dir, paste0("MSE - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  resolution <- 0.001
  leg.range <- max(abs(c(floor(min(rmse.mat$PointwiseRisk, na.rm = T)*10)/10,
                         ceiling(max(rmse.mat$PointwiseRisk, na.rm = T)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals.grid, y=gvals.grid, z=rmse.mat$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorRampPalette(c("white", "red"))(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(0, leg.range), xlim = c(min(dvals.grid), max(dvals.grid)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals.grid, y = gvals.grid, z = rmse.mat$PointwiseRisk, add = T, labcex = 1.2)
                 })
  mtext("rMSE", side = 4, line = 2.5, las = 0)
  par(op)
  dev.off()
  
  # Absolute Bias (One-Number Summary)
  set.idx <- which(!is.na(abs(mat.bias$PointwiseRisk)), arr.ind = T)
  summ.results_test[j,"Total Absolute Bias"] <- mean(abs(mat.bias$PointwiseRisk), na.rm = T)
  
  # Total MCSD (One-Number Summary)
  summ.results_test[j,"Total MCSD"] <- mean(sd.mat$PointwiseRisk[set.idx], na.rm = T)
  
  # Total rMSE (One-Number Summary)
  summ.results_test[j,"Total rMSE"] <- mean(rmse.mat$PointwiseRisk, na.rm = T)
  
  # Total MCE (One-Number Summary)
  summ.results_test[j,"MCE"] <- mean(mce.mat$PointwiseRisk[set.idx], na.rm = T) # ? - Double Check
}



############## Summary Statistics ###############

# Training Set Statistics
# Brier Score, Deviance, Eff. Num of Parameters, DIC
stats.labels_train <- c("Brier", "Deviance", "Eff. Num", "DIC")
summ.results_train <- matrix(NA, nrow = length(param.names),
                             ncol = length(stats.labels_train))
colnames(summ.results_train) <- stats.labels_train
rownames(summ.results_train) <- param.names

summ.results_train[,"Brier"] <- round(colMeans(sim_bivar_correct$Brier, na.rm = T)*100, 2)
summ.results_train[,"Deviance"] <- c(round(colMeans(sim_bivar_correct$Deviance, na.rm = T), 2))
summ.results_train[,"Eff. Num"] <- c(round(colMeans(sim_bivar_correct$effNumParms, na.rm = T), 2))
summ.results_train[,"DIC"] <- c(round(colMeans(sim_bivar_correct$DIC, na.rm = T), 2))


(summ.results_all <- cbind(summ.results_test, summ.results_train))
library(xtable)
xtable(summ.results_all, digits = 3)

# mce.rounder <- function(ests, mce, alpha = qnorm(0.975)) {
#   cis <- matrix(rep(as.numeric(ests), 2), nrow = 2, byrow = T) + alpha*matrix(rep(c(-1, 1), each = ncol(ests)), nrow = 2, byrow = T)*mce
#   rownames(cis) <- rep(rownames(ests), 2)
#   # eq = T
#   # counter <- 0
#   # while(eq) {
#   #   round(cis, counter)
#   #   if(all.equal()
#   # }
#   # return(round(2))
#   cis
# }

# t(sapply(1:nrow(noncurve.summary), function(x) {
#   mce.rounder(noncurve.summary[x,-5], noncurve.summary[x, 5], alpha = 1)
# }))

# noncurve.summary


# Rename
library(stringi)
rename.files <- T
dir.spec <- here(res.dir, "To Present", "2023_12_02", "Correct")
if(rename.files) {
  old.names <- dir.files <- list.files(dir.spec)
  dir.files <- stri_replace_all_fixed(dir.files, "Bias - 2D Contour (104 Replicates) - ", "bias_2d_contour_", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "MCSD - 2D Contour (104 Replicates) - ", "mcsd_2d_contour_", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "Mean vs. True (Causal) - 2D Contour (104 Replicates) - ", "mean_vs_true_2d_", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "Mean vs. True (Causal) - 3D Perspective (Iterations = 104) - ", "mean_vs_true_3d_", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "MSE - 2D Contour (104 Replicates) - ", "mse_2d_contour_", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "Additive (Doubly)", "additive", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "Polynomial Logistic", "polylog", vectorize_all=FALSE)
  dir.files <- stri_replace_all_fixed(dir.files, "Logistic", "logistic", vectorize_all=FALSE)
  new.names <- dir.files <- stri_replace_all_fixed(dir.files, "Bivariable", "bivariable", vectorize_all=FALSE)
  cbind(old.names, new.names) %>% View()
  file.rename(here(dir.spec, old.names), here(dir.spec, new.names))
}



##### MISSPECIFIED SCENARIO #####

# Mean Curve
sim_B <- nrow(sim_bivar_misspec$Brier)
grid.num.d <- 12
grid.num.g <- 12
min.range <- 30
max.range <- 50
dvals.grid <- seq(min.range, max.range, length.out = grid.num.d + 2)[-c(1, grid.num.d + 2)]
gvals.grid <- seq(0, 1, length.out = grid.num.g + 2)[-c(1, grid.num.g + 2)]
param.names <- c("Logistic", "Polynomial Logistic",
                 "Additive (Doubly)", "Bivariable")
stats.labels <- c("Total Absolute Bias", "Total MCSD",
                  "Total rMSE", "MCE")
# Empty Storage Vectors

summ.results_test <- matrix(NA, nrow = length(param.names),
                            ncol = length(stats.labels))
colnames(summ.results_test) <- stats.labels; rownames(summ.results_test) <- param.names

# Testing Set
n <- 500
nobs <- grid.num.d*n
grid.x1 <- (dvals.grid - min.range)/(max.range - min.range)
grid.x2 <- gvals.grid
grid.dat <- expand.grid(data.frame(`Dose` = grid.x1, `Volume` = grid.x2))
test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat))
colnames(test) <- c("Intercept", "Dose", "Volume")

listLabel <- c("PointwiseRisk")
grid.num.pars <- length(param.names)
mean.mat <- lapply(1, function(x) {
  nullmat <- matrix(NA, nrow = grid.num.d, ncol = grid.num.g)
  rownames(nullmat) <- dvals.grid
  colnames(nullmat) <- gvals.grid
  return(nullmat)
}); names(mean.mat) <- listLabel
true.mean.mat <- true.ul.mat <- true.ll.mat <- mean.mat
true.sd.mat <- true.rmse.mat <- true.mce.mat <- mean.mat
ul.mat <- ll.mat <- mean.mat
sd.mat <- mean.mat

for(d in 1:grid.num.d) {
  for(g in 1:grid.num.g) {
    for(k in 1) {
      if(k == 1) {
        true.mean.mat[[k]][d,g] <- mean(sim_bivar_misspec$Predictions[d,g,5,])
      }
      true.sd.mat[[k]][d,g] <- 0
      true.ul.mat[[k]][d,g] <- 0
      true.ll.mat[[k]][d,g] <- 0
      # true.rmse.mat[[k]][d,g] <- true.sd.mat[[k]][d,g]
      # true.mce.mat[[k]][d,g] <- true.sd.mat[[k]][d,g]/sqrt(sim_B)
    }
  }
}
true.rmse.mat <- true.sd.mat
true.mce.mat <- lapply(1, function(l) {
  true.sd.mat[[l]]/sqrt(sim_B)
}); names(true.mce.mat) <- listLabel
summ.results_test[length(param.names),] <- c(
  0,                                 # True Bias
  mean(true.sd.mat$PointwiseRisk, na.rm = T),      # True MCSD
  mean(true.rmse.mat$PointwiseRisk, na.rm = T),    # True rMSE
  mean(true.mce.mat$PointwiseRisk, na.rm = T)     # True MCE
)

for(j in 1:grid.num.pars) {
  
  for(d in 1:grid.num.d) {
    for(g in 1:grid.num.g) {
      for(k in 1) {
        if(sum(!is.na(sim_bivar_misspec$Predictions[d,g,j,])) > 0) {
          mean.mat[[k]][d,g] <- mean(sim_bivar_misspec$Predictions[d,g,j,], na.rm = T)
          # mat.bias[d,g] <- mean(sim_bivar_misspec$Predictions[d,g,j,], na.rm = T) - true.mean.mat[[k]][d,g]
          sd.mat[[k]][d,g] <- sd(sim_bivar_misspec$Predictions[d,g,j,], na.rm = T)
          ul.mat[[k]][d,g] <- quantile(sim_bivar_misspec$Predictions[d,g,j,], probs = 0.975, na.rm = T)
          ll.mat[[k]][d,g] <- quantile(sim_bivar_misspec$Predictions[d,g,j,], probs = 0.025, na.rm = T)
          # rmse.mat[[k]][d,g] <- sqrt(bias.mat[[k]][d,g]^2 + sd.mat[[k]][d,g]^2)
          # mce.mat[[k]][d,g] <- sd.mat[[k]][d,g]/sqrt(sim_B)
        } else {
          mean.mat[[k]][d,g] <- NA
          # mat.bias[d,g] <- NA
          ul.mat[[k]][d,g] <- NA
          ll.mat[[k]][d,g] <- NA
          # rmse.mat[[k]][d,g] <- NA
          # mce.mat[[k]][d,g] <- NA
        }
      }
    }
  }
  mat.bias <- lapply(1, function(l) {
    mean.mat[[l]] - true.mean.mat[[l]]
  }); names(mat.bias) <- listLabel
  rmse.mat <- lapply(1, function(l) {
    sqrt(mat.bias[[l]]^2 + sd.mat[[l]]^2)
  }); names(rmse.mat) <- listLabel
  mce.mat <- lapply(1, function(l) {
    sd.mat[[l]]/sqrt(sim_B)
  }); names(mce.mat) <- listLabel
  
  # 2D Contour Plots
  
  # Mean Estimates (Causal Risk)
  pdf(here(res.dir, paste0("Mean vs. True (Causal) - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 6), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  contour(x = dvals.grid, y = gvals.grid, z = true.mean.mat$PointwiseRisk, col = "blue", xlim = range(dvals.grid), lwd = 1.5,
          xlab = "Dose (Gy)", ylab = "Volume")
  contour(x = dvals.grid, y = gvals.grid, z = mean.mat$PointwiseRisk, nlevels = 10, col = "orange", add = T, lwd = 1.5)
  abline(h = c(0, 1), v = range(dvals.grid), col = "gray")
  legend("topright", lty = c(1,1), c("Truth", "Fit"), col = c("blue", "orange"),
         inset = c(-0.13, 0.4), xpd = T)
  par(op)
  rasterpdf::dev.off()
  
  # 3D Contour Plots
  
  # Mean vs. True
  added_lab <- param.names[j]
  persp.plot(xvec = test[,"Dose"],
             yvec = test[,"Volume"],
             zvec = mean.mat$PointwiseRisk,
             ztrue = true.mean.mat$PointwiseRisk,
             title = paste0("Mean vs. True (Causal) - 3D Perspective (Iterations = ", nrow(sim_bivar_misspec$Brier), ") - ", added_lab,".pdf"))
  
  # Bias - 2D Contour Plots
  pdf(here(res.dir, paste0("Bias - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  # par(xpd = T)
  resolution <- 0.001
  bias.mat <- mat_pairs(mat.bias$PointwiseRisk, xvals = dvals.grid, yvals = gvals.grid)
  leg.range <- max(abs(c(floor(min(bias.mat$zval)*10)/10,
                         ceiling(max(bias.mat$zval)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals.grid, y=gvals.grid, z=mat.bias$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorspace::diverging_hsv(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(-leg.range, leg.range), xlim = c(min(dvals.grid), max(dvals.grid)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals.grid, y = gvals.grid, z = mat.bias$PointwiseRisk, add = T, labcex = 1.2)
                 })
  mtext("Bias", side = 4, line = 2.5, las = 0)
  
  par(op)
  dev.off()
  
  # Monte Carlo Standard Deviation
  pdf(here(res.dir, paste0("MCSD - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  resolution <- 0.001
  leg.range <- max(abs(c(floor(min(sd.mat$PointwiseRisk, na.rm = T)*10)/10,
                         ceiling(max(sd.mat$PointwiseRisk, na.rm = T)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals.grid, y=gvals.grid, z=sd.mat$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorRampPalette(c("white", "red"))(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(0, leg.range), xlim = c(min(dvals.grid), max(dvals.grid)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals.grid, y = gvals.grid, z = sd.mat$PointwiseRisk, add = T, labcex = 1.2)
                 })
  mtext("MCSD", side = 4, line = 2.5, las = 0)
  
  par(op)
  dev.off()
  
  # Root Mean-Squared Error
  pdf(here(res.dir, paste0("MSE - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  resolution <- 0.001
  leg.range <- max(abs(c(floor(min(rmse.mat$PointwiseRisk, na.rm = T)*10)/10,
                         ceiling(max(rmse.mat$PointwiseRisk, na.rm = T)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals.grid, y=gvals.grid, z=rmse.mat$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorRampPalette(c("white", "red"))(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(0, leg.range), xlim = c(min(dvals.grid), max(dvals.grid)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals.grid, y = gvals.grid, z = rmse.mat$PointwiseRisk, add = T, labcex = 1.2)
                 })
  mtext("rMSE", side = 4, line = 2.5, las = 0)
  par(op)
  dev.off()
  
  # Absolute Bias (One-Number Summary)
  set.idx <- which(!is.na(abs(mat.bias$PointwiseRisk)), arr.ind = T)
  summ.results_test[j,"Total Absolute Bias"] <- mean(abs(mat.bias$PointwiseRisk), na.rm = T)
  
  # Total MCSD (One-Number Summary)
  summ.results_test[j,"Total MCSD"] <- mean(sd.mat$PointwiseRisk[set.idx], na.rm = T)
  
  # Total rMSE (One-Number Summary)
  summ.results_test[j,"Total rMSE"] <- mean(rmse.mat$PointwiseRisk, na.rm = T)
  
  # Total MCE (One-Number Summary)
  summ.results_test[j,"MCE"] <- mean(mce.mat$PointwiseRisk[set.idx], na.rm = T) # ? - Double Check
}



############## Summary Statistics ###############

# Training Set Statistics
# Brier Score, Deviance, Eff. Num of Parameters, DIC
stats.labels_train <- c("Brier", "Deviance", "Eff. Num", "DIC")
summ.results_train <- matrix(NA, nrow = length(param.names),
                             ncol = length(stats.labels_train))
colnames(summ.results_train) <- stats.labels_train
rownames(summ.results_train) <- param.names

summ.results_train[,"Brier"] <- round(colMeans(sim_bivar_misspec$Brier, na.rm = T)*100, 2)
summ.results_train[,"Deviance"] <- c(round(colMeans(sim_bivar_misspec$Deviance, na.rm = T), 2))
summ.results_train[,"Eff. Num"] <- c(round(colMeans(sim_bivar_misspec$effNumParms, na.rm = T), 2))
summ.results_train[,"DIC"] <- c(round(colMeans(sim_bivar_misspec$DIC, na.rm = T), 2))


(summ.results_all <- cbind(summ.results_test, summ.results_train))
library(xtable)
xtable(summ.results_all, digits = 3)

# mce.rounder <- function(ests, mce, alpha = qnorm(0.975)) {
#   cis <- matrix(rep(as.numeric(ests), 2), nrow = 2, byrow = T) + alpha*matrix(rep(c(-1, 1), each = ncol(ests)), nrow = 2, byrow = T)*mce
#   rownames(cis) <- rep(rownames(ests), 2)
#   # eq = T
#   # counter <- 0
#   # while(eq) {
#   #   round(cis, counter)
#   #   if(all.equal()
#   # }
#   # return(round(2))
#   cis
# }

# t(sapply(1:nrow(noncurve.summary), function(x) {
#   mce.rounder(noncurve.summary[x,-5], noncurve.summary[x, 5], alpha = 1)
# }))

# noncurve.summary