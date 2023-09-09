# functions for 01_FindSppThreshold.R
##contains:
##1. function to calculate ratio of density splits / density observations
##2. function to output stacked ggplots for looking at different types of thresholds against cum. imp. curves, and presence / absence data (requires ggpubr!)
if(!"ggpubr" %in% (.packages())) { cat("warning: package 'ggpubr' is required to run 'check_threshold' plotting function and is not loaded.\n") }
if(!"gradientForest" %in% (.packages())) { cat("warning: package 'gradientForest' is required to run these functions and is not loaded.\n") }
if(!"tidyverse" %in% (.packages())) { cat("warning: package 'tidyverse' is required to run these functions and is not loaded.\n") }

#this is modified from gradientForest::split.density.plot.method2 script!
get_dens_ratio <- function(gf_model, pred_var) { 

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

check_threshold <- function(thresh_df, pa_df, spp_name, pred_var, plot_it = TRUE, save_location = FALSE) { #make function to make it easier
  
  if(nrow(thresh_df %>% filter(species == spp_name)) == 0) { stop(paste("Species name not found.")) }
  if(nrow(thresh_df %>% filter(env_var == pred_var)) == 0) { stop(paste("Variable name not found.")) }
  spp.cu <- thresh_df %>% filter(species == spp_name & env_var == pred_var)
  
  #get thresholds
  if( sum(spp.cu$slope, na.rm = T) != 0 ) {
    t.x <- spp.cu %>% select(contains("t_")) %>% select(-contains("_prop")) %>% distinct() %>% select(!where(is.na)) %>%
      pivot_longer(everything(), values_to = "thresh", names_to = "type")
  }

  #get raw values
  t.spp <- pa_df %>% 
    select(one_of(spp_name, pred_var), site_x, site_y)
  
  #plot theme for consistency
  thresh_theme <- list(
    {if(exists("t.x")) geom_vline(data = t.x, aes(xintercept = thresh, color = type), linewidth = 1.5, alpha = 0.5)},
    theme_linedraw(),
    scale_x_continuous(name = pred_var,
                       limits = c(floor(min(t.spp[, pred_var], na.rm = T)), ceiling(max(t.spp[, pred_var], na.rm = T))),
                       n.breaks = abs( ceiling(max(t.spp[, pred_var], na.rm = T)) - floor(min(t.spp[, pred_var], na.rm = T)) + 1 )
    ),
    theme(axis.title = element_text(size = 9), legend.title = element_text(size = 9)),
    labs(color = "Threshold Type:")
  )
  
  #boxplot for presence/absence distro around threshold
  p_box <- ggplot() +
    geom_boxplot(data = t.spp, aes(.data[[pred_var]], .data[[spp_name]]), outlier.shape = NA) +
    geom_jitter(data = t.spp, aes(.data[[pred_var]], .data[[spp_name]]), width = 0, height = 0.3, size = 1, alpha = 0.5) +
    ylab("pres. / abs.") +
    thresh_theme
  # p_box
  
  #cumulative importance plot
  p_cimp <- ggplot() +
    geom_path(data = spp.cu, aes(X, Y)) +
    xlab(pred_var) + ylab("cumul. imp.") +
    thresh_theme
  # p_cimp
  
  #slope plot
  p_slope <- ggplot() +
    geom_path(data = spp.cu , aes(X, slope)) +
    xlab(pred_var) + ylab("cumul. imp. slope") +
    thresh_theme
  # p_slope
  
  #stack
  p <- ggpubr::ggarrange(p_cimp, p_slope, p_box, ncol = 1, nrow = 3, common.legend = T, legend = "bottom", align = "v")
  p_all <- ggpubr::annotate_figure(p, top = text_grob(spp_name, face = "bold"))
  # p_all
  
  if(save_location != FALSE) {
    # ggsave(paste0(save_location, "/ThreshPlot_", spp_name, ".png"), plot = p_all)
    png(filename = paste0(save_location, "/ThreshPlot_", spp_name, "_", pred_var, ".png"), width = 5, height = 6, units = "in", res = 200)
    print(p_all)
    dev.off()
  }
  
  if(plot_it) { return(p_all) }
}

# check_threshold(thresh_df = CU_all, pa_df = full.in, spp_name = "Sialis", pred_var = "site_x",
#                 save_location = paste0(PATH, "/11_Thresholds/Plots"))

# #map it ----
# library(maps); library(sf)
# 
# ozark <- st_as_sf(map(database = "state", region = c('missouri', 'arkansas', 'oklahoma'), plot = FALSE, fill = TRUE))
# occ.dat <- t.spp %>% st_as_sf(., coords = c("site_x", "site_y"), crs = st_crs(4326), remove = F)
# 
# ggplot() +
#   geom_sf(data = ozark) +
#   geom_sf(data = occ.dat, aes(fill = .data[[spp]]), shape = 21) +
#   # geom_vline(xintercept = t.x, color = "red") +
#   theme_minimal()