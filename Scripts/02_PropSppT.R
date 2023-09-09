## CALC. PROP. SPECIES W/ THRESHOLD
# by C.E.Moore V.
# 
# OBJ: Calculate proportion of species with a threshold

library(tidyverse); library(gradientForest)
source(paste0(PATH, "/Scripts/01a_FindThreshFunc.R"))





# proportion of spp with threshold -
# CU_all <- read_csv(paste0(PATH, "/11_Thresholds/full_thresh_", gsub(".rds", "", gsub(paste0(PATH, "/10_GFOutput/"), "", file.name)), ".csv"))
# 
# ## Method 1. Spp above ratio dens split importance:dens obs ----
# #from Chen et al. 2023
# #make table of ratios
# r.df <- data.frame(env_var = iVars, ratio = NA)
# for (var in iVars) {
#   r.df[r.df$env_var == var, ]$ratio <- get_dens_ratio(gf.mod, var)
# }
# 
# #find species with threshold density ratio > 1
# # tmp <- CU_all %>% filter(env_var == "site_x")
# # # CU_all %>% group_by(env_var)
# # tmp %>%
# #   group_by(species)