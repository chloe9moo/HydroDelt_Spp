# functions for 01_FindSppThreshold.R

get_dens_ratio <- function(gf_model, pred_var) { 
  #this is modified from gradientForest::split.density.plot.method2 script!
  #pull out necessary values
  i <- pred_var
  imp <- importance(gf_model)[i]
  resA <- gf_model$res[gf_model$res$var == i, ]
  splits <- resA$split
  w <- pmax(resA$improve.norm, 0)
  X <- na.omit(gf_model$X[,i])
  rX <- range(X)
  dX <- diff(rX)
  
  #prep needed functions
  normalize.density <- function(d,integral=1,integrate=T) {
    # scale y values so that density integrates to integral
    Id <- if(integrate) integrate.density(d) else 1; 
    d$y <- d$y/Id*integral; 
    d
  }
  
  integrate.density <- function(d) {
    integrate(approxfun(d,rule=2),lower=min(d$x),upper=max(d$x))$value
  }
  
  scale.density <- function(d,scale=1/mean(d$y)) {
    d$y <- d$y*scale
    d
  }
  
  is.binned <- function(obj) {
    compact <- obj$call$compact
    if (is.null(compact)) 
      FALSE
    else
      eval(compact)
  }
  
  normalize.histogram <- function(ci, integral=1, bin=F, nbin=101) {
    # scale y values so that histogram integrates to integral
    # optionally aggregate the y's into binned x ranges 
    if (bin) {
      brks <- seq(min(ci$x),max(ci$x),len=nbin)
      xx <- cut(ci$x,breaks=brks,inc=T)
      yy <- tapply(ci$y,xx,sum)
      yy[is.na(yy)] <- 0
      ci <- list(x=0.5*(brks[-1]+brks[-nbin]),y=yy)
    }
    dx <- min(diff(ci$x)); 
    Id <- sum(ci$y*dx); 
    ci$y <- ci$y/Id*integral; 
    ci
  }
  
  # importance density I(x)
  suppressWarnings({
    dImp <- density(splits, weight = w/sum(w), from = rX[1], to = rX[2])
  })
  if((dX/dImp$bw) > 50)  dImp <- density(splits, weight = w/sum(w), from = rX[1], to = rX[2], bw=dX/50)
  dImpNorm <- normalize.density(dImp, imp, integrate=T)
  
  # data density d(x)
  dObs <- density(X, from = rX[1], to = rX[2])
  if((dX/dObs$bw) > 50)  dObs <- density(X, from = rX[1], to = rX[2], bw=dX/50)
  lambda <- 0.9
  dObs$y <- lambda*dObs$y + (1-lambda)/diff(range(dObs$x))
  dObsNorm <- normalize.density(dObs,imp,integrate=T)
  
  # raw importances I
  ci <- cumimp(gf.mod, i, standardize=F)   
  ci$y <- diff(c(0, ci$y))
  ci <- normalize.histogram(ci, imp, bin=F | !is.binned(gf.mod), nbin=101)
  
  # standardized density f(x) = I(x)/d(x)
  dStd <- dImp
  dStd$y <- dImp$y/dObs$y
  dStdNorm <- try(normalize.density(dStd,imp,integrate=T)) # this sometimes does not converge
  if (class(dStdNorm) == "try-error") dStdNorm <- normalize.histogram(dStd,imp)
  
  # get ratio value
  ratio_val <- mean(dStdNorm$y)/mean(dStd$y)
  return(ratio_val)
}

# this is how the splits density plot is made FYI, ratio was taken from abline:
# plot(ci, type='h', col="grey60", xlim=range(splits), lwd=1,
#      #ylim = c(0, max(dImpNorm$y, dObsNorm$y, dStdNorm$y, ci$y)), lend=2,
#      ylim = c(0, max(dImpNorm$y, dObsNorm$y, dStdNorm$y)*1.1), lend=2,
#      xlab = i, ylab = "")
# lines(dImpNorm, col = "black", lwd = 2)
# lines(dObsNorm, col = "red", lwd = 2)
# lines(dStdNorm, col = "blue", lwd = 2)
# abline(h = mean(dStdNorm$y)/mean(dStd$y), lty = 2, col = "blue")
