### UNPACKING SIMULATION FILE ###
wd <- here("data")
setwd(wd)
outpath2 <- here("data", "storage")
res.dir <- here("main simulation", "results")

# Unpack Files of Interest
# file.names <- c("simulation_n=500_01-52.RData",
#                 "simulation_n=500_52-104")
file.names <- "simulation_n=500.RData"
simulation.list <- lapply(1:length(file.names), function(id) {
  load(file = file.path(outpath2, file.names[[id]]))
  simulation
}); names(simulation.list) <- file.names
a.obj <- simulation.list[[1]]$Predictions

# Unpack and Combine
library(abind)
simulation <- vector("list", length(simulation.list[[1]]));
names(simulation) <- names(simulation.list[[1]])
for(i in 1:length(file.names)) {
  simulation$Predictions <- abind(simulation$Predictions,
             simulation.list[[i]]$Predictions, along = 5)
  simulation$`Stochast Predictions` <- rbind(simulation$`Stochast Predictions`,
             simulation.list[[i]]$`Stochast Predictions`)
  simulation$Brier <- rbind(simulation$Brier,
                            simulation.list[[i]]$Brier)
  simulation$Deviance <- rbind(simulation$Deviance,
                               simulation.list[[i]]$Deviance)
  simulation$effNumParms <- rbind(simulation.list[[i]]$effNumParms,
                                  simulation.list[[i]]$effNumParms)
  simulation$DIC <- rbind(simulation$DIC,
                          simulation.list[[i]]$DIC)
  simulation$time <- rbind(simulation$time,
                           simulation.list[[i]]$time)
}

sim_bivar <- simulation

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

# MC Results Simulation
replicates <- simulation$Predictions
no.dec = T
latex = F


############### VISUALIZATIONS ##################

# Mean Curve
sim_B <- nrow(sim_bivar$Brier)
num.d <- 26
num.g <- 26
dvals <- seq(30, 50, length.out = num.d)
gvals <- seq(0, 1, length.out = num.g)
param.names <- c("Logistic", "Polynomial Logistic",
                 "Additive (Doubly)", "Bivariable",
                 "True")
stats.labels <- c("Total Absolute Bias", "Total MCSD",
                  "Total rMSE", "MCE")
# Empty Storage Vectors

summ.results_test <- matrix(NA, nrow = length(param.names),
                            ncol = length(stats.labels))
colnames(summ.results_test) <- stats.labels; rownames(summ.results_test) <- param.names

# Testing Set
# nobs <- nrow(long.data)
nobs <- num.d*n
grid.x1 <- (dvals - min(dvals))/(max(dvals) - min(dvals))
grid.x2 <- gvals
grid.dat <- expand.grid(data.frame(`Dose` = grid.x1, `Volume` = grid.x2))
test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat))
colnames(test) <- c("Intercept", "Dose", "Volume")

listLabel <- c("x = 0", "x = 1", "PointwiseRisk")
num.pars <- length(param.names)
mean.mat <- lapply(1:3, function(x) {
  nullmat <- matrix(NA, nrow = num.d, ncol = num.g)
  rownames(nullmat) <- dvals
  colnames(nullmat) <- gvals
  return(nullmat)
}); names(mean.mat) <- listLabel
true.mean.mat <- true.ul.mat <- true.ll.mat <- mean.mat
true.sd.mat <- true.rmse.mat <- true.mce.mat <- mean.mat
ul.mat <- ll.mat <- mean.mat
sd.mat <- mean.mat

for(d in 1:num.d) {
  for(g in 1:num.g) {
    for(k in 1:3) {
      if(k == 1) {
        true.mean.mat[[k]][d,g] <- true.risk$`x = 0`[d,g]
      }
      if(k == 2) {
        true.mean.mat[[k]][d,g] <- true.risk$`x = 1`[d,g]
      }
      if(k == 3) {
        true.mean.mat[[k]][d,g] <- true.crisk[d,g]
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
true.mce.mat <- lapply(1:3, function(l) {
  true.sd.mat[[l]]/sqrt(sim_B)
}); names(true.mce.mat) <- listLabel
summ.results_test[length(param.names),] <- c(
  0,                                 # True Bias
  mean(true.sd.mat$PointwiseRisk, na.rm = T),      # True MCSD
  mean(true.rmse.mat$PointwiseRisk, na.rm = T),    # True rMSE
  mean(true.mce.mat$PointwiseRisk, na.rm = T)     # True MCE
)

for(j in 1:(num.pars - 1)) {
  
  for(d in 1:num.d) {
    for(g in 1:num.g) {
      for(k in 1:3) {
        if(sum(!is.na(sim_bivar$Predictions[d,g,k,j,])) > 0) {
          mean.mat[[k]][d,g] <- mean(sim_bivar$Predictions[d,g,k,j,], na.rm = T)
          # mat.bias[d,g] <- mean(sim_bivar$Predictions[d,g,k,j,], na.rm = T) - true.mean.mat[[k]][d,g]
          sd.mat[[k]][d,g] <- sd(sim_bivar$Predictions[d,g,k,j,], na.rm = T)
          ul.mat[[k]][d,g] <- quantile(sim_bivar$Predictions[d,g,k,j,], probs = 0.975, na.rm = T)
          ll.mat[[k]][d,g] <- quantile(sim_bivar$Predictions[d,g,k,j,], probs = 0.025, na.rm = T)
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
  mat.bias <- lapply(1:3, function(l) {
    mean.mat[[l]] - true.mean.mat[[l]]
  }); names(mat.bias) <- listLabel
  rmse.mat <- lapply(1:3, function(l) {
    sqrt(mat.bias[[l]]^2 + sd.mat[[l]]^2)
  }); names(rmse.mat) <- listLabel
  mce.mat <- lapply(1:3, function(l) {
    sd.mat[[l]]/sqrt(sim_B)
  }); names(mce.mat) <- listLabel
  
  # 2D Contour Plots
  
  # Mean Estimates (Causal Risk)
  pdf(here(res.dir, paste0("Mean vs. True (Causal) - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 6), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  contour(x = dvals, y = gvals, z = true.mean.mat$PointwiseRisk, col = "blue", xlim = range(dvals), lwd = 1.5,
          xlab = "Dose (Gy)", ylab = "Volume")
  contour(x = dvals, y = gvals, z = mean.mat$PointwiseRisk, nlevels = 10, col = "orange", add = T, lwd = 1.5)
  abline(h = c(0, 1), v = range(dvals), col = "gray")
  legend("topright", lty = c(1,1), c("Truth", "Fit"), col = c("blue", "orange"),
         inset = c(-0.13, 0.4), xpd = T)
  par(op)
  rasterpdf::dev.off()
  
  # Mean Estimates (by Strata)
  pdf(here(res.dir, paste0("Mean by Tumor Stage - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 8), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  contour(x = dvals, y = gvals, z = mean.mat$`x = 1`, col = "blue", xlim = range(dvals), lwd = 1.5,
          xlab = "Dose (Gy)", ylab = "Volume")
  contour(x = dvals, y = gvals, z = mean.mat$`x = 0`, nlevels = 10, col = "orange", add = T, lwd = 1.5)
  abline(h = c(0, 1), v = range(dvals), col = "gray")
  legend("topright", lty = c(1,1), c(paste0("3+ (N = ", table(wide.data$x)[2],")"),
                                     paste0("1-2 (N = ", table(wide.data$x)[1],")")), col = c("orange", "blue"),
         inset = c(-0.19, 0.4), xpd = T,
         title = "Tumor Stage")
  rasterpdf::dev.off()
  
  # 3D Contour Plots
  
  # Mean vs. True
  added_lab <- param.names[j]
  persp.plot(xvec = test[,"Dose"],
             yvec = test[,"Volume"],
             zvec = mean.mat$PointwiseRisk,
             ztrue = true.mean.mat$PointwiseRisk,
             title = paste0("Mean vs. True (Causal) - 3D Perspective (Iterations = ", nrow(sim_bivar$Brier), ") - ", added_lab,".pdf"))
  
  # Mean By Tumor Stage
  added_lab <- param.names[j]
  persp.plot(xvec = test[,"Dose"],
             yvec = test[,"Volume"],
             zvec = mean.mat$`x = 0`,
             ztrue = mean.mat$`x = 1`,
             title = paste0("Mean by Tumor Stage - 3D Perspective (Iterations = ", nrow(sim_bivar$Brier), ") - ", added_lab,".pdf"))
  
  
  
  # Bias - 2D Contour Plots
  pdf(here(res.dir, paste0("Bias - 2D Contour (104 Replicates) - ", param.names[j],".pdf")),
      width = 12, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0 + 3), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  
  # par(xpd = T)
  resolution <- 0.001
  bias.mat <- mat_pairs(mat.bias$PointwiseRisk, xvals = dvals, yvals = gvals)
  leg.range <- max(abs(c(floor(min(bias.mat$zval)*10)/10,
                         ceiling(max(bias.mat$zval)*10)/10)))
  leg.range <- 0.5
  filler.contour(x=dvals, y=gvals, z=mat.bias$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorspace::diverging_hsv(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(-leg.range, leg.range), xlim = c(min(dvals), max(dvals)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals, y = gvals, z = mat.bias$PointwiseRisk, add = T, labcex = 1.2)
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
  filler.contour(x=dvals, y=gvals, z=sd.mat$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorRampPalette(c("white", "red"))(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(0, leg.range), xlim = c(min(dvals), max(dvals)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals, y = gvals, z = sd.mat$PointwiseRisk, add = T, labcex = 1.2)
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
  filler.contour(x=dvals, y=gvals, z=rmse.mat$PointwiseRisk, nlevels = 64,
                 color.palette = function(n) colorRampPalette(c("white", "red"))(n), lwd = 1.5,
                 xlab = "Dose (Gy)", ylab = "Volume",
                 zlim = c(0, leg.range), xlim = c(min(dvals), max(dvals)), ylim = c(-0.01, 1.01),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(x = dvals, y = gvals, z = rmse.mat$PointwiseRisk, add = T, labcex = 1.2)
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

# Dotplots of Stochastic Causal Risk

stochest.dat <- as_tibble(cbind(
  1:nrow(sim_bivar$`Stochast Predictions`),
  sim_bivar$`Stochast Predictions`
))
colnames(stochest.dat) <- c("ID", param.names)
stochest.dat <- stochest.dat %>% select(-True) %>%
  pivot_longer(cols = Logistic:Bivariable,
               names_to = "Model",
               values_to = "Estimate")
stochest.summ <- stochest.dat %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Estimate),
    LL = quantile(Estimate, probs = 0.025),
    UL = quantile(Estimate, probs = 0.975)
  ) %>% ungroup
stochest.true <- true.scrisk
stochest.plot <- stochest.summ %>%
  ggplot(aes(x = factor(Model, levels = param.names[-5]), y = Mean)) +
  geom_hline(yintercept = stochest.true, col = "red", size = 0.25) +
  geom_hline(yintercept = c(0, 1), linetype = "solid", color = "gray", size = 0.25) +
  geom_point(size = 2) + theme_bw() +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  labs(x = "Model", y = "Stochastic Causal Risk") +
  geom_errorbar( aes(ymin = `LL`, ymax = `UL`),
                 width = 0.05) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 18))
  )
ggsave(plot = stochest.plot, filename = here(res.dir, "Sampling Distribution - Stochastic Causal Risk.pdf"),
       device = "pdf", width = 6, height = 4)


############## Summary Statistics ###############

# Training Set Statistics
# Brier Score, Deviance, Eff. Num of Parameters, DIC
stats.labels_train <- c("Brier", "Deviance", "Eff. Num", "DIC")
summ.results_train <- matrix(NA, nrow = length(param.names),
                             ncol = length(stats.labels_train))
colnames(summ.results_train) <- stats.labels_train
rownames(summ.results_train) <- param.names

summ.results_train[,"Brier"] <- round(colMeans(sim_bivar$Brier, na.rm = T)*100, 2)
summ.results_train[,"Deviance"] <- round(colMeans(sim_bivar$Deviance, na.rm = T), 2)
summ.results_train[,"Eff. Num"] <- c(round(colMeans(sim_bivar$effNumParms, na.rm = T), 2), NA)
summ.results_train[,"DIC"] <- c(round(colMeans(sim_bivar$DIC, na.rm = T), 2), NA)


(summ.results_all <- cbind(summ.results_test, summ.results_train))
library(xtable)
xtable(summ.results_all, digits = 3)

# Combine n = 100 and n = 500
# Unpack Files of Interest
file.names <- c("stochastic_causal_refine_n=100.RData",
                "stochastic_causal_refine_n=500.RData")
simulation.list <- lapply(1:length(file.names), function(id) {
  load(file = file.path(outpath2, file.names[[id]]))
  simulation
}); names(simulation.list) <- file.names
a.obj <- simulation.list[[1]]$Predictions

# Unpack and Combine
library(abind)
simulation <- vector("list", length(simulation.list[[1]]));
names(simulation) <- names(simulation.list[[1]])
for(i in 1:length(file.names)) {
  simulation$Predictions <- abind(simulation$Predictions,
                                  simulation.list[[i]]$Predictions, along = 5)
  simulation$`Stochast Predictions` <- rbind(simulation$`Stochast Predictions`,
                                             simulation.list[[i]]$`Stochast Predictions`)
  simulation$Brier <- rbind(simulation$Brier,
                            simulation.list[[i]]$Brier)
  simulation$Deviance <- rbind(simulation$Deviance,
                               simulation.list[[i]]$Deviance)
  simulation$effNumParms <- rbind(simulation.list[[i]]$effNumParms,
                                  simulation.list[[i]]$effNumParms)
  simulation$DIC <- rbind(simulation$DIC,
                          simulation.list[[i]]$DIC)
  simulation$time <- rbind(simulation$time,
                           simulation.list[[i]]$time)
}
sim500 <- simulation

load(file = file.path(outpath2, "stochastic_causal_refine_n=100.RData"))
sim100 <- simulation


# Combined Dotplots of Stochastic Causal Risk

stochest.dat <- as_tibble(cbind(
  1:nrow(sim100$`Stochast Predictions`),
  sim100$`Stochast Predictions`
)) %>% bind_rows(as_tibble(cbind(
  1:nrow(sim500$`Stochast Predictions`),
  sim500$`Stochast Predictions`
))) %>% mutate(
  n = factor(rep(c(100, 500), each = nrow(sim100$`Stochast Predictions`)))
)
colnames(stochest.dat) <- c("ID", param.names, "n")
stochest.dat <- stochest.dat %>% select(-True) %>%
  pivot_longer(cols = Logistic:Bivariable,
               names_to = "Model",
               values_to = "Estimate")
stochest.summ <- stochest.dat %>%
  group_by(Model, n) %>%
  summarise(
    Mean = mean(Estimate),
    LL = quantile(Estimate, probs = 0.025),
    UL = quantile(Estimate, probs = 0.975)
  ) %>% ungroup
stochest.true <- true.scrisk
dodge <- position_dodge(width = 0.5)
parm.labs <- c("Logistic", "Polynomial", "Additive", "Bivariable")
stochest.plot <- stochest.summ %>%
  ggplot(aes(x = factor(Model, levels = param.names[-5], labels = parm.labs), y = Mean, color = n)) +
  geom_hline(yintercept = stochest.true, col = "red", size = 0.2) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "black", size = 0.2) +
  geom_point(size = 1.8, position = dodge) + theme_bw() +
  scale_y_continuous(limits = c(0.2, 0.6), breaks = seq(0, 1.0, by = 0.1)) +
  scale_color_manual(values = c("orange", "blue")) +
  labs(x = "Model", y = "Stochastic Causal Risk", color = "Sample Size") +
  geom_errorbar( aes(ymin = `LL`, ymax = `UL`, color = n), size = 0.2,
                 width = 0.16, position = dodge) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 8, margin = margin(t = 11)),
        axis.title.y = element_text(size = 8, margin = margin(r = 14)),
        axis.text = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, -2),
        legend.box.margin = margin(0, 0, 0, -2)
  )
ggsave(plot = stochest.plot, filename = here(res.dir, "Sampling Distribution - Stochastic Causal Risk (n = 100, 500).pdf"),
       device = "pdf", width = 4.5, height = 1.8)




