library(here)
source(here("analysis code", "functions.R"))
outpath2 <- here("data", "storage")
# Packages
installer(c(c("tidyverse", "compiler"),
            c("gam", "akima", "fields", "abind",
              "monoreg", "pheatmap")))

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
nexp <- 1000

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

d.star <- dvals.grid[findInterval(d.star, dvals.grid)]
g.seq <- seq(0, 1, by = 0.005)

# Functional Mean Model
out.parms <- c(-21, 0.5, 1, 3)

## DGM FUNCTIONS ##
dg.calc <- function(mu, sd, d) {
  dnorm(d, mean = mu, sd = sd)
}
g.calc <- function(mu, sd, d) {
  1 - pnorm(d, mean = mu, sd = sd)
}
fun.mod <- function(mu, id) { expit(as.numeric(c(1, mu, x1 = true.data$x1[id], true.data$x2[id]) %*% out.parms)) }
nan.recode <- function(x) {na_if(x, "NaN")}
expit <- function(q) {1/(1 + exp(-q))}
logit <- function(p) {log(p/(1 - p))}

# Matrix to Long-Format Data
mat_pairs <- function(mat, xvals, yvals) {
  grid.vals <- expand.grid(xvals, yvals)
  data.fit2 <- data.frame(xvals = grid.vals$Var1, yvals = grid.vals$Var2, `zval` = as.vector(mat))
  data.fit2[complete.cases(data.fit2),]
}

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
sig.lims <- function(g, dstar, id, prnt = F) {
  muroot <- c(amin, amax)
  sigroot <- c(bmin, bmax)
  
  mu.sols <- sort(sapply(sigroot, function(sigval) {mufind(sigma = sigval, dstar = dstar, g = g, muroot = muroot)}), na.last = T)
  bound.sols <- ifelse(!is.na(sapply(muroot, function(muval) {sigmafind(mu = muval, dstar = dstar, g = g, sigroot = sigroot)})), muroot, NA)
  if(sum(!is.na(c(mu.sols, bound.sols))) == 0) {
    return(c(NA, NA))
  } else {
    return(range(c(mu.sols, bound.sols), na.rm = T))
  }
}

# Kernel Functions for Mu (Calculating Risk)
mu.expit.kernel <- function(mu, dstar, g, id) {
  # g = 0.2; mu = mean(g.ll(g), g.ul(g)); dstar = d.star; id = 1
  alpha.val <- true.data$alpha.breg[id]
  beta.val <- true.data$beta.breg[id]
  sigroot <- c(bmin, bmax)
  sig.val <- sigmafind(mu, dstar = dstar, g = g, sigroot = sigroot)
  if(is.na(sig.val)) {return(NA)}
  norm.fun <- function(d, mu) {
    dnorm((d - mu)/(sig.val))*(1 - ((d - mu)/(sig.val))^2)
  }
  J <- (sig.val)^2/abs(integrate(Vectorize(norm.fun), 0, dstar, mu = mu)$value)
  # if(J > 10^15) {return(NA)}
  fun.mod(mu = mu, id = id)*((mu - amin)^(alpha.val - 1))*((amax - mu)^(beta.val - 1))*J
}
mu.kernel <- function(mu, dstar, g, id) {
  # g = 0.2; mu = mean(g.ll(g), g.ul(g)); dstar = d.star; id = 1
  alpha.val <- true.data$alpha.breg[id]
  beta.val <- true.data$beta.breg[id]
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

# g.grid <- seq(0.01, 0.99, by = 0.01)
# mu.lims <- sapply(g.grid, function(g) {
#   sig.lims(g = g, dstar = d.star, id = 1)
# })
# g.ll <- splinefun(g.grid, mu.lims[1,])
# g.ul <- splinefun(g.grid, mu.lims[2,])
# plot(gvals, g.ll(gvals), type = "l", lwd = 2, ylim = c(30, 50),
#      xlab = "g", ylab = bquote(mu))
# lines(gvals, g.ul(gvals), type = "l", lwd = 2)


# True Risk Function
risk.fun <- function(dstar, g, id = id) {
  parms <- c(amin, amax)
  
  mu.range <- sig.lims(g, dstar = dstar, id = id)
  if(sum(is.na(mu.range)) > 0) {return(NA)}
  
  mu.min <- mu.range[1]
  mu.max <- mu.range[2]
  
  # Set as NA if Equation Cannot Be Solved
  if(is.na(mu.kernel(mu.min, dstar = dstar, g = g, id = id)) | is.na(mu.kernel(mu.max, dstar = dstar, g = g, id = id)) |
     is.infinite(mu.kernel(mu.min, dstar = dstar, g = g, id = id)) | is.infinite(mu.kernel(mu.max, dstar = dstar, g = g, id = id)) |
     sum(c(mu.kernel(mu.min, dstar = dstar, g = g, id = id), mu.kernel(mu.max, dstar = dstar, g = g, id = id)) > 10000) > 1 |
     is.na(mu.expit.kernel(mu.min, dstar = dstar, g = g, id = id)) | is.na(mu.expit.kernel(mu.max, dstar = dstar, g = g, id = id)) |
     is.infinite(mu.expit.kernel(mu.min, dstar = dstar, g = g, id = id)) | is.infinite(mu.expit.kernel(mu.max, dstar = dstar, g = g, id = id)) |
     sum(c(mu.expit.kernel(mu.min, dstar = dstar, g = g, id = id), mu.expit.kernel(mu.max, dstar = dstar, g = g, id = id)) > 10000) > 1
  ) {return(NA)}
  if(exists("a")) {rm(a)}
  if(exists("b")) {rm(b)}
  a <- tryCatch(
    integrate(Vectorize(mu.expit.kernel), mu.min, mu.max, dstar = dstar, g = g, id = id,
              stop.on.error = F)$value,
    error = function(e) NA
  )
  b <- tryCatch(
    integrate(Vectorize(mu.kernel), mu.min, mu.max, dstar = dstar, g = g, id = id,
              stop.on.error = F)$value,
    error = function(e) NA
  )
  # if(!exists("a") | !exists("b")) {return(NA)}
  a/b
}

# True Dataset
true.data <- tibble(
  id = 1:n,
  x1 = rbinom(n, 1, p.grade),
  x2 = rnorm(n),
  eta.breg = as.numeric(cbind(1, x1, x2) %*% beta.parms),
  alpha.breg = up.scale*expit(eta.breg),
  beta.breg = up.scale*(1 - expit(eta.breg))
)

# GENERATING TRUE VALUES

# True Pointwise Causal Risk
true.covs.mat <- expand.grid(dvals.grid, gvals.grid, 1:nrow(true.data))
colnames(true.covs.mat) <- c("Dose", "Volume", "ID")

# Discontinued (Takes a Long Time)
# true.risk <- sapply(1:grid.num.d, function(g) {
#   sapply(1:grid.num.g, function(d) {
#     mean(sapply(1:nrow(true.data), function(id) {
#       nan.recode(risk.fun(g = gvals.grid[g], dstar = dvals.grid[d], id = id))
#     }), na.rm = T)
#   })
# })
# rownames(true.risk) <- dvals.grid; colnames(true.risk) <- gvals.grid

# 1 Hour (Can Parallelize)
true.covs.mat$TruePred <- sapply(1:nrow(true.covs.mat), function(id) {
  nan.recode(risk.fun(g = true.covs.mat$Volume[id], dstar = true.covs.mat$Dose[id], id = true.covs.mat$ID[id]))
})
true.covs.mat <- as_tibble(true.covs.mat) %>%
  group_by(Dose, Volume) %>% summarize(
    TruePred = mean(TruePred, na.rm = T)
  ) %>% ungroup
true.risk <- matrix(true.covs.mat$TruePred, nrow = grid.num.d, ncol = grid.num.g, byrow = T)
rownames(true.risk) <- dvals.grid; colnames(true.risk) <- gvals.grid
# save(true.risk, file = file.path(outpath2, paste0("simulation_supp_strongconf_truerisk.RData")))