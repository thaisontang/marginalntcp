## DGM FUNCTIONS ##
dg.calc <- function(mu, sd, d) {
  dnorm(d, mean = mu, sd = sd)
}
g.calc <- function(mu, sd, d) {
  1 - pnorm(d, mean = mu, sd = sd)
}
fun.mod <- function(mu, x) { expit(alpha + beta*mu + beta2*x) }
nan.recode <- function(x) {na_if(x, "NaN")}

expit <- function(q) {1/(1 + exp(-q))}
logit <- function(p) {log(p/(1 - p))}

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

g.grid <- seq(0.01, 0.99, by = 0.01)
mu.lims <- sapply(g.grid, function(g) {
  sig.lims(g = g, dstar = dstar)
})
g.ll <- splinefun(g.grid, mu.lims[1,])
g.ul <- splinefun(g.grid, mu.lims[2,])

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


# Initial Parameter Values
bmin <- 1      
bmax <- 2      # bmax <- 10
amin <- 35     # amin <- 5
amax <- 45     # amax <- 76.5
bmin_x0 <- bmin
bmax_x0 <- bmax
beta.parms.x0 <- c(4, 6)/3
beta.parms.x1 <- rev(beta.parms.x0)
amin_upper <- amin
amax_lower <- amax
p.grade <- 0.4
n <- 500  # n <- 100
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

true.risk <- lapply(1:2, function(x) {
  nullmat <- matrix(NA, nrow = num.d, ncol = num.g)
  rownames(nullmat) <- dvals
  colnames(nullmat) <- gvals
  return(nullmat)
}); names(true.risk) <- c("x = 0", "x = 1")


# True Quantities
true.risk <- lapply(1:2, function(x) {matrix(NA, nrow = num.d, ncol = num.g)})
for(d in 1:num.d) {
  for(g in 1:num.g) {
    for(x in 1:2) {
      xval <- c(0, 1)[x]
      
      # Print
      cat("d:", d, "/", num.d, "     ", "g:", g, "/", num.g, "     x:", x, "/ 2\n")
      
      # True Risk (Conditional)
      true.risk[[x]][d,g] <- nan.recode(risk.fun(g = gvals[g], dstar = dvals[d], x = xval))
    }
  }
}
names(true.risk) <- c("x = 0", "x = 1")

# True Densities
dens.0 <- posthoc.adjust(sapply(g.seq, function(g) {g.dens.make(g = g, dstar = d.star, x = 0)}), g.s = g.seq)$fun
dens.1 <- posthoc.adjust(sapply(g.seq, function(g) {g.dens.make(g = g, dstar = d.star, x = 1)}), g.s = g.seq)$fun

true.crisk <- true.risk[[1]]*(1 - p.grade) + true.risk[[2]]*p.grade
rownames(true.crisk) <- dvals; colnames(true.crisk) <- gvals
true.srisk <- c(scrisk.fun(dstar = d.star, x = 0), scrisk.fun(dstar = d.star, x = 1)); names(true.srisk) <- c("x = 0", "x = 1")
true.scrisk <- unname(true.srisk[1]*(1 - p.grade) + true.srisk[2]*p.grade)