## AUTOMATED SPP THRESHOLD DETERMINATION
# by C.E.Moore V.
# 
# OBJ: Find individual species' threshold based on change in cumulative importance 

library(tidyverse); library(gradientForest)
# library(parallel) #<<multimodel runs

# set variables----
PATH <- getwd()

#which model(s) to analyze? (what is the name of the file)
file.name <- paste0(PATH, "/10_GFOutput/gf.bugs.23.GW_full_cat.rds")  #single model run
# file.list <- list.files(paste0(PATH, "/10_GFOutput"), pattern = "gf.", full.names = T)
#how many variables (ordered by importance) to threshold?
num.vars <- 5
#what is the allowable slope cut off for consolidating thresholds - ie, percent of targeted slope?
prx.val <- 0.60 #<< could be part of the threshold if its at least 60% of the highest slope value
#what is the allowable slope cut off for alternative thresholds besides highest slope?
sl2.val <- 0.80 #<< other peaks are highlighted if slope is greater then at least 80% of highest slope value

# load necessary files // start multi model run ----
# n.cores <- detectCores()/2
# mclapply(file.list, mc.cores = n.cores, function(file.name) { #start parallel run, ~3 min for 18 models on 8 cores, max ~14 GB RAM
  
gf.mod <- readRDS(file.name) #gf model of interest

spp_imp <- as.data.frame(importance(gf.mod, type = "Species")) %>% rownames_to_column("species") #pull out species performance
colnames(spp_imp)[2] <- "rel_err"

iVars <- names(importance(gf.mod))[1:num.vars] #pull out variables of interest


# find thresholds, all spp x imp. vars. ----
CU <- list()
# i <- 4; x <- 104 #code testing ONLY!

for (i in 1:length(iVars)) { #for every important variable ...
  
  CU[[i]] <- cumimp(gf.mod, iVars[[i]], "Species") #get cumulative importance for all species
  
  CU[[i]] <- lapply(seq_along(CU[[i]]), function(x) { #get various threshold values for all species

    spp.cu <- as.data.frame(CU[[i]][[x]])
    colnames(spp.cu) <- c("X", "Y")

    spp.cu <- spp.cu %>%
      mutate(species = names(CU[[i]])[[x]], #get species name for later
             env_var = iVars[[i]], #get variable name for later
             slope = (Y - lead(Y)) / (X - lead(X))) #get slope of cumulative importance curve
    
    if (sum(spp.cu$slope, na.rm = T) == 0) { #deal with spp with flat cumulative importance curve
      spp.cu <- spp.cu %>%
        mutate(t_max_start = NA, t_max_end = NA, t_max_mn = NA, n_peaks = 0, t_m2_start = NA, t_m2_end = NA, t_m2_mn = NA)
      return(spp.cu)
    }
    
    #find slopes that are close together to consolidate
    idx <- which.max(spp.cu$slope) #index of largest slope
    n <- (max(spp.cu$X, na.rm = T) - min(spp.cu$X, na.rm = T)) * 0.05 #10% window (of total env gradient)
    
    tmp <- spp.cu %>%
      filter((spp.cu[idx,]$X - n) < X & X < (spp.cu[idx,]$X + n))
    
    tmp.t <- max(tmp$slope, na.rm = T) * prx.val
    
    #get index of start and end of change
    t.idx <- which(spp.cu$X %in% (tmp %>% filter(slope > tmp.t & slope != max(slope, na.rm = T)))$X)
    t.idx <- c(t.idx, idx) #combine all potential start / end points
    t.idx <- c(start = min(t.idx, na.rm = T), end = max(t.idx, na.rm = T)+1)
    
    #add to df
    spp.cu$t_max_start <- spp.cu[t.idx[["start"]],]$X
    spp.cu$t_max_end <- spp.cu[t.idx[["end"]],]$X
    spp.cu$t_max_mn <- ( unique(spp.cu$t_max_start) + unique(spp.cu$t_max_end) ) / 2
    
    #get total cumulative importance change and threshold proportion to total
    d <- max(spp.cu$Y, na.rm = T) - min(spp.cu$Y, na.rm = T)
    
    spp.cu$t_max_prop <- ( 
      spp.cu[which(spp.cu$X == unique(spp.cu$t_max_end)),]$Y - spp.cu[which(spp.cu$X == unique(spp.cu$t_max_start)),]$Y
    ) / d * 100
    
    #find other potential thresholds throughout rest of curve
    tmp <- spp.cu
    tmp[t.idx[["start"]]:t.idx[["end"]], ] <- NA
    
    tmp.t <- max(spp.cu$slope, na.rm = T) * sl2.val #find other peaks
    spp.cu <- spp.cu %>% mutate(n_peaks = sum(tmp$slope > tmp.t, na.rm = T)+1) #count other alternatives
    
    if (unique(spp.cu$n_peaks) > 1) { #if more than 1 threshold is found...
      
      #threshold second time
      idx <- which.max(tmp$slope) #index of largest slope
      
      tmp2 <- tmp %>%
        filter((tmp[idx,]$X - n) < X & X < (tmp[idx,]$X + n))
      
      tmp.t <- max(tmp2$slope, na.rm = T) * prx.val
      
      t.idx <- which(tmp$X %in% (tmp2 %>% filter(slope > tmp.t & slope != max(slope, na.rm = T)))$X)
      t.idx <- c(t.idx, idx) #combine all potential start / end points
      t.idx <- c(start = min(t.idx, na.rm = T), end = max(t.idx, na.rm = T)+1)
      
      #add to df
      spp.cu$t_m2_end <- tmp[t.idx[["start"]],]$X
      spp.cu$t_m2_start <- tmp[t.idx[["end"]],]$X
      spp.cu$t_m2_mn <- ( unique(spp.cu$t_m2_start) + unique(spp.cu$t_m2_end) ) / 2
      
      #get total cumulative importance change and threshold proportion to total
      spp.cu$t_m2_prop <- ( spp.cu[t.idx[["end"]],]$Y - spp.cu[t.idx[["start"]],]$Y ) / d * 100
      
    } else {
      spp.cu$t_m2_end <- NA; spp.cu$t_m2_start <- NA; spp.cu$t_m2_mn <- NA; spp.cu$t_m2_prop <- NA
    }

    # #threshold = find point a bit before largest change, to account for slight increases before the jump
    # sub_val <- 5 #set how far out to crop for second threshold method
    # if (sub_val < idx) { 
    #   tmp <- spp.cu %>% rownames_to_column("index") %>% slice((idx-sub_val):(idx+sub_val)) 
    # } else {
    #   tmp <- spp.cu %>% rownames_to_column("index") %>% slice(0:(idx+sub_val)) }
    # 
    # tmp <- tmp %>%
    #   mutate(norm_y = (Y - min(Y)) / (max(Y) - min(Y)),
    #          find_peak = abs(norm_y - 0.05))
    # idx <- tmp[which.min(tmp$find_peak),]$index
    # 
    # spp.cu$t_5p <- spp.cu[idx,]$X
    # 
    # spp.cu$t_5p_prop <- ( spp.cu[spp.cu$X == unique(spp.cu$t_max_end), ]$Y - spp.cu[idx,]$Y ) / d * 100
    
    return(spp.cu)
    })
  
  CU[[i]] <- bind_rows(CU[[i]]) #combine
  
  CU[[i]] <- CU[[i]] %>%
    left_join(., spp_imp, by = "species") #add error just in case it's useful later
  
  cat(paste0(iVars[[i]], " thresholds calculated..\n"))
  if(i == length(iVars)) cat("\nDone!")
}

# quick check threshold creation // see below for nicer plotting
# plot(spp.cu$X, spp.cu$Y, type = "l")
# points(tmp$X, tmp$Y, col = "blue", lwd = 2)
# abline(v = unique(spp.cu$t_max_start), col = "red")
# abline(v = unique(spp.cu$t_max_end), col = "blue")
# abline(v = unique(spp.cu$t_max_mn), col = "green")
# abline(v = unique(spp.cu$t_m2_end), col = "purple")
# abline(v = unique(spp.cu$t_m2_start), col = "orange")
# abline(v = unique(spp.cu$t_m2_mn), col = "pink")
# dev.off()

CU_all <- bind_rows(CU)

write_csv(CU_all, paste0(PATH, "/11_Thresholds/full_thresh_", 
                         gsub(".rds", "", gsub(paste0(PATH, "/10_GFOutput/"), "", file.name)), #matching name to GF model names
                         ".csv"))

# }) #end parallel

# indiv. species plots to check things ----
# library(ggpubr)
# source(paste0(PATH, "/Scripts/XX_find_thresh_func.R"))
# ##need to adjust to plot larger env gradients (e.g., WS Area)
# CU_all <- read_csv(paste0(PATH, "/11_Thresholds/full_thresh_gf.bugs.23.GW_HIT_cat.csv"))
# gf.mod <- readRDS(file.list[[2]])
# full.in <- cbind(gf.mod$Y, gf.mod$X) #pull out full predictor + response data and combine (for later)
# spp_imp <- as.data.frame(importance(gf.mod, type = "Species")) %>% rownames_to_column("species")
# 
# check_threshold(thresh_df = CU_all, pa_df = full.in, spp_name = "Orthocladiinae", pred_var = "fh10")

#save top spp plots
# top_spp <- as.list((spp_imp %>% arrange(desc(rel_err)))[1:5, "species"]) 
# lapply(top_spp, function(x) {
#   for (var in iVars) {
#     check_threshold(thresh_df = CU_all, pa_df = full.in, spp_name = x, pred_var = var, plot_it = F,
#                     save_location = paste0(PATH, "/11_Thresholds/Plots"))
#   }
# })
