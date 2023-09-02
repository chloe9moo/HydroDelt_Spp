## AUTOMATED SPP THRESHOLD DETERMINATION
# by C.E.Moore V.
# 
# OBJ: Find individual species threshold based on [criteria] >> determine proportion of species in assemblage with a threshold

library(tidyverse); library(gradientForest)

# set variables ----
PATH <- getwd()
num.vars <- 5 #how many variables to pull out
#!!! probably set a way to call a specific variable later

# load in gf model, pull out info ----
gf.mod <- readRDS(paste0(PATH, "/10_GFOutput/gf.bugs.23.GW_full_cat.rds"))

imp.vars <- names(importance(gf.mod))[1:num.vars]

CU <- list()
for (varX in imp.vars) {
  # print(varX)
  CU[[varX]] <- cumimp(gf.mod, varX, "Species")
}

#!!! look up how others are defining thresholds, specifically chen et al., since they compared thresholds



# plotting ----
imp.sp <- sapply(CU$site_x, function(cu) max(cu$y))
imp.sp <- order(-imp.sp)[1:min(length(CU$site_x), length(imp.sp))]
spp1 <- names(CU$site_x)[imp.sp][5]
spp1.xy <- CU$site_x[spp1] %>% as.data.frame()
colnames(spp1.xy) <- c("x","y")

ggplot(spp1.xy) +
  geom_point(aes(x, y)) +
  geom_path(aes(x, y)) +
  labs(x = "site_x", y = "cumulative importance") +
  theme_classic()






