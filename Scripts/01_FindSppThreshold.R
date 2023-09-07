## AUTOMATED SPP THRESHOLD DETERMINATION
# by C.E.Moore V.
# 
# OBJ: Find individual species' threshold based on [criteria] 
# >> determine proportion of species in assemblage with a threshold

library(tidyverse); library(gradientForest)
library(ggpubr)

# set variables + load ----
PATH <- getwd()
source(paste0(PATH, "/Scripts/01a_FindThreshFunc.R"))

num.vars <- 5 #how many variables to pull out

file.name <- paste0(PATH, "/10_GFOutput/gf.bugs.23.GW_full_cat.rds") #save name separate for later

gf.mod <- readRDS(file.name)

spp_imp <- as.data.frame(importance(gf.mod, type = "Species")) %>% rownames_to_column("species") #pull out species performance
colnames(spp_imp)[2] <- "rel_err"

full.in <- cbind(gf.mod$Y, gf.mod$X) #pull out full predictor + response data and combine (for later)

iVars <- names(importance(gf.mod))[1:num.vars]


# find thresholds, all spp x imp. vars. ----
sub_val <- 5 #set how far out to crop for second threshold method
CU <- list()

for (i in 1:length(iVars)) { #for every important variable ...
  CU[[i]] <- cumimp(gf.mod, iVars[[i]], "Species") #get cumulative importance for all species
  
  CU[[i]] <- lapply(seq_along(CU[[i]]), function(x) { #get various threshold values for all species

    spp.cu <- as.data.frame(CU[[i]][[x]])
    colnames(spp.cu) <- c("X", "Y")

    spp.cu <- spp.cu %>%
      mutate(species = names(CU[[i]])[[x]], #get species name for later
             env_var = iVars[[i]], #get variable name for later
             slope = (Y - lag(Y)) / (X - lag(X))) #get slope of cumulative importance curve
    
    if (sum(spp.cu$slope, na.rm = T) == 0) { #deal with spp, def. no threshold
      spp.cu <- spp.cu %>%
        mutate(t_max_start = NA, t_max_end = NA, t_5p = NA, t_m2_end = NA, t_m2_start = NA, t_m2_5p = NA)
      return(spp.cu)
    }
    
    #threshold = largest change in cumulative importance, start and end of the change
    idx <- which.max(spp.cu$slope) #index of largest slope
    
    spp.cu$t_max_start <- spp.cu[idx-1,]$X
    spp.cu$t_max_end <- spp.cu[idx,]$X
    
    #threshold = find point a bit before largest change, to account for slight increases before the jump
    if (sub_val < idx) { 
      tmp <- spp.cu %>% rownames_to_column("index") %>% slice((idx-sub_val):(idx+sub_val)) 
    } else {
      tmp <- spp.cu %>% rownames_to_column("index") %>% slice(0:(idx+sub_val)) }
    
    tmp <- tmp %>%
      mutate(norm_y = (Y - min(Y)) / (max(Y) - min(Y)),
             find_peak = abs(norm_y - 0.05))
    idx <- tmp[which.min(tmp$find_peak),]$index
    
    spp.cu$t_5p <- spp.cu[idx,]$X
    
    #other potential thresholds that are close in magnitude to largest change?
    tmp.t <- max(spp.cu$slope, na.rm = T) * 0.9 #within 90% of largest value
    spp.cu <- spp.cu %>%
      mutate(n_peaks = sum(spp.cu$slope > tmp.t, na.rm = T))
    
    if (unique(spp.cu$n_peaks) > 1) {
      spp.cu$t_m2_end <- (spp.cu %>% filter(slope > tmp.t & slope != max(slope, na.rm = T)) %>% arrange(desc(slope)))[1,]$X
      
      idx <- which(spp.cu$t_m2_end == spp.cu$X)
      spp.cu$t_m2_start <- spp.cu[(idx-1),]$X
      
      if (sub_val < idx) { 
        tmp <- spp.cu %>% rownames_to_column("index") %>% slice((idx-sub_val):(idx+sub_val)) 
      } else {
        tmp <- spp.cu %>% rownames_to_column("index") %>% slice(0:(idx+sub_val)) }
      
      tmp <- tmp %>%
        mutate(norm_y = (Y - min(Y)) / (max(Y) - min(Y)),
               find_peak = abs(norm_y - 0.05))
      idx <- tmp[which.min(tmp$find_peak),]$index
      
      spp.cu$t_m2_5p <- spp.cu[idx,]$X
      
    } else {
      spp.cu$t_m2_end <- NA
      spp.cu$t_m2_start <- NA
      spp.cu$t_m2_5p <- NA
    }
    
    return(spp.cu)
    })
  
  CU[[i]] <- bind_rows(CU[[i]]) #combine
  
  CU[[i]] <- CU[[i]] %>%
    left_join(., spp_imp, by = "species") #add error just in case it's useful later
  
  cat(paste0(iVars[[i]], " thresholded..\n"))
  if(i == length(iVars)) cat("\nDone!")
}

# quick check threshold creation // see below for nicer plotting
# plot(spp.cu$X, spp.cu$Y, type = "l")
# abline(v = unique(spp.cu$t_max_start), col = "red")
# abline(v = unique(spp.cu$t_max_end), col = "blue")
# abline(v = unique(spp.cu$t_5p), col = "green")
# abline(v = unique(spp.cu$t_m2_end), col = "purple")
# abline(v = unique(spp.cu$t_m2_start), col = "orange")
# abline(v = unique(spp.cu$t_m2_5p), col = "pink")
# dev.off()

CU_all <- bind_rows(CU)

write_csv(CU_all, paste0(PATH, "/11_Thresholds/full_thresh_", 
                         gsub(".rds", "", gsub(paste0(PATH, "/10_GFOutput/"), "", file.name)), #matching name to GF model names
                         ".csv"))

# proportion of spp with threshold ----
## Method 1. Spp above ratio dens split importance:dens obs ----
get_dens_ratio(gf.mod, "TminCat")

# indiv. species plots to check things ----
check_threshold <- function(spp_name, pred_var) { #make function to make it easier
  
  if(nrow(CU_all %>% filter(species == spp_name)) == 0) { stop(paste("Species name not found.")) }
  if(nrow(CU_all %>% filter(env_var == pred_var)) == 0) { stop(paste("Variable name not found.")) }
  spp.cu <- CU_all %>% filter(species == spp_name & env_var == pred_var)
  
  #get thresholds
  t.x <- spp.cu %>% select(contains("t_")) %>% distinct() %>% select(!where(is.na)) %>%
    pivot_longer(everything(), values_to = "thresh", names_to = "type")
  
  #get raw values
  t.spp <- full.in %>% 
    select(one_of(spp_name, pred_var), site_x, site_y)
  
  #plot theme for consistency
  thresh_theme <- list(
    theme_linedraw(),
    scale_x_continuous(name = pred_var, 
                       limits = c(floor(min(t.spp[, pred_var], na.rm = T)), ceiling(max(t.spp[, pred_var], na.rm = T))),
                       n.breaks = abs( ceiling(max(t.spp[, pred_var], na.rm = T)) - floor(min(t.spp[, pred_var], na.rm = T)) + 1 )
    )
  )
  
  #boxplot for presence/absence distro around threshold
  p_box <- ggplot() +
    geom_boxplot(data = t.spp, aes(.data[[pred_var]], .data[[spp_name]])) +
    geom_jitter(data = t.spp, aes(.data[[pred_var]], .data[[spp_name]]), width = 0.2) +
    geom_vline(data = t.x, aes(xintercept = thresh, color = type), linewidth = 2, alpha = 0.6) +
    ylab("presence / absence") +
    thresh_theme
  # p_box
  
  #cumulative importance plot
  p_cimp <- ggplot() +
    geom_path(data = spp.cu, aes(X, Y)) +
    geom_vline(data = t.x, aes(xintercept = thresh, color = type), linewidth = 2, alpha = 0.6) +
    xlab(pred_var) + ylab("cumulative importance") +
    thresh_theme
  # p_cimp
  
  #slope plot
  p_slope <- ggplot() +
    geom_path(data = spp.cu , aes(X, slope)) +
    geom_vline(data = t.x, aes(xintercept = thresh, color = type), linewidth = 2, alpha = 0.6) +
    xlab(pred_var) + ylab("cumulative importance slope") +
    thresh_theme
  # p_slope
  
  #stack
  p <- ggpubr::ggarrange(p_cimp, p_slope, p_box, ncol = 1, nrow = 3, common.legend = T, legend = "right", align = "v")
  p_all <- ggpubr::annotate_figure(p, top = text_grob(spp_name, face = "bold"))
  # p_all
  
  return(p_all)
}

check_threshold(spp_name = "Maccaffertium", pred_var = "TminCat")

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


