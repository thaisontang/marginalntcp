# BASE FUNCTIONS

# Package Installation
installer <- function(vector) {
  for(i in vector) {
    if(!require(i, character.only = T)) install.packages(i)
    library(i, character.only = T)
  }
}

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


# Quick Rugplots
rugplot <- function(x, from = NULL, to = NULL, label = "") {
  if(is.null(from)) from <- min(x, na.rm = T)
  if(is.null(to)) to <- max(x, na.rm = T)
  cat("Missing Obs: ", sum(is.na(x)))
  plot(density(x, from = from, to = to, na.rm = T), lwd = 2, main = "", col = "navy", xlab = label)
  rug(x, col = "red")
}

# BASE PACKAGES
installer(c("tidyverse", "here"))


# filled.contour customize
filler.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                         length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n, 
                                                                                                                 "YlOrRd", rev = TRUE), col = color.palette(length(levels) - 
                                                                                                                                                              1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
                            xaxs = "i", yaxs = "i", las = 1, axes = TRUE, 
                            frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", 
                                 "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, 
              asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

# Calculate DVH from Dose Densities
DVHmapper <- function(obj, doses) {
  Dose <- doses
  Volume <- matrix(NA, nrow = length(obj[[1]]), ncol = length(doses))
  
  for(y in 1:length(obj[[1]])) {
    ecdf.lint <- approxfun(obj$Dose[[y]], obj$CDF[[y]], yleft = 0, yright = 1)
    if(y %% 10 == 0) print(y)
    Volume[y,] <- ecdf.lint(Dose)
  }
  colnames(Volume) <- Dose
  return(Volume)
}

# Subsetting a List by Row Index
ListSubset <- function(list.obj, sub.vec = NULL){
  if(is.null(sub.vec)) sub.vec <- 1:length(list.obj[[1]])
  new.list <- lapply(1:length(list.obj), function(x)list.obj[[x]][sub.vec])
  names(new.list) <- names(list.obj)
  return(new.list)
}



monofunc <- function(covs,
                     spec.use = spec.short,         # Choice of Number of Iterations/Burn-In
                     parm.use = parm.nomono         # Choice of Model Parametrization
) {
  (begin <- Sys.time())
  results <- monosurv(niter=spec.use$niter, burnin=spec.use$burnin, adapt=spec.use$adapt,
                      refresh = 10, thin = 5,
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

