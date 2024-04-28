### UNPACKING SIMULATION FILE ###
res.dir <- here("simulated_illustration", "figures")
installer(c("tidyverse", "plot3D", "colorspace", "rasterpdf"))


# Load Bootstrap Results
outpath2 <- here("simulated_illustration", "data", "storage")  # Resampling Storage Folder
load(file = file.path(outpath2, paste0("ntcp_bootresamp_age65_bladder.RData")))

# Combine Bootstrap pointwise.preds across 8 Cores
pointwise.preds <- R.utils::wrap(bootresamp$Predictions, map = list(1, 2, 3, NA))

# VISUALIZATION FUNCTIONS

# Restrict Contour Region (based on Observed Data)
signif.floor <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- floor(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}
signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}

dvh.bound <- long.data %>% group_by(Dose) %>%
  summarise(`minVol` = signif.floor(min(Volume_Bladder), 2),
            `maxVol` = signif.ceiling(max(Volume_Bladder), 2)) %>%
  mutate(dose.idx = 1:num.d,
         minVol.idx = findInterval(minVol, gvals),
         maxVol.idx = ifelse(findInterval(maxVol, gvals) != 51, findInterval(maxVol, gvals) + 1, findInterval(maxVol, gvals))) %>%
  ungroup %>% dplyr::select(dose.idx:maxVol.idx) %>%
  rowwise %>% transmute(dose.idx, Vol.idx = list(seq(minVol.idx, maxVol.idx, 1))) %>%
  unnest(Vol.idx) %>% data.matrix


poly.seq <- long.data %>% group_by(Dose) %>%
  summarise(`minVol` = signif.floor(min(Volume_Bladder), 2),
            `maxVol` = signif.ceiling(max(Volume_Bladder), 2) + 0.02)
poly.seq <- rbind(cbind(poly.seq$Dose, poly.seq$minVol),
                  cbind(rev(poly.seq$Dose), rev(poly.seq$maxVol)))

bound.mat <- function(mat) {
  nullmat <- matrix(NA, nrow(mat), ncol(mat))
  nullmat[dvh.bound] <- mat[dvh.bound]
  return(nullmat)
}

# Visualization Test Set
dose.ind <- (dvals - min(dvals))/(max(dvals) - min(dvals))
vol.ind <- gvals
grid.dat <- expand.grid(`Dose` = dose.ind, `Volume` = vol.ind)
test <- data.matrix(cbind(rep(1, nrow(grid.dat)), grid.dat))
colnames(test) <- c("Intercept", "Dose", "Volume")

param.names <- c("Linear", "Polynomial", "Additive", "Bivariable")
num.pars <- length(param.names)

for(j in 1:num.pars) {
  
  # Mean Pointwise NTCP Estimates
  mean.mat <- apply(pointwise.preds[,,j,], c(1, 2), function(x) {mean(x, na.rm = T)})
  
  # 2D Contour Plots
  
  # Mean Estimates (Causal Risk)
  pdf(here(res.dir, paste0("Figure2", letters[j],".pdf")),
      width = 10, height = 6)
  op <- par(las = 1, mar = c(4.5, 5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
  contour(x = dvals, y = gvals, z = bound.mat(mean.mat), nlevels = 10, col = "black",
          xlim = range(dvals), lwd = 1.5, xlab = "Dose (Gy)", ylab = "Volume")
  abline(h = c(0, 1), v = range(dvals), col = "gray")
  lines(poly.seq, lty = "dashed")
  # legend("topright", lty = c(1), c("Fit"), col = c("blue"),
  #        inset = c(-0.13, 0.4), xpd = T)
  par(op)
  rasterpdf::dev.off()

}


# Dotplots of Stochastic Causal Risk

stochest.dat <- as_tibble(cbind(
  1:nrow(bootresamp$`Stochast Predictions`),
  bootresamp$`Stochast Predictions`/bootresamp$`Obs Means`[,rep(1, num.pars)]
))
colnames(stochest.dat) <- c("ID", param.names)
stochest.dat <- stochest.dat %>%
  pivot_longer(cols = Linear:Bivariable,
               names_to = "Model",
               values_to = "Estimate")
stochest.summ <- stochest.dat %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Estimate),
    LL = quantile(Estimate, probs = 0.025),
    UL = quantile(Estimate, probs = 0.975)
  ) %>% ungroup 

# Training Set Statistics
# Brier Score, Deviance, Eff. Num of Parameters, DIC
stats.labels_train <- c("Brier", "Deviance", "Eff. Num", "DIC")
summ.results_train <- matrix(NA, nrow = length(param.names),
                             ncol = length(stats.labels_train))
colnames(summ.results_train) <- stats.labels_train
rownames(summ.results_train) <- param.names

summ.results_train[,"Brier"] <- round(colMeans(bootresamp$Brier, na.rm = T)*100, 2)
summ.results_train[,"Deviance"] <- round(colMeans(bootresamp$Deviance, na.rm = T), 2)
summ.results_train[,"Eff. Num"] <- round(colMeans(bootresamp$effNumParms, na.rm = T), 2)
summ.results_train[,"DIC"] <- round(colMeans(bootresamp$DIC, na.rm = T), 2)

# Table 1 (Bladder)
cbind(stochest.summ[c(3:4, 1:2),-1], summ.results_train)