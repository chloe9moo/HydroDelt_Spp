#GRADIENT FORESTS 

#mostly updates

library(tidyverse); library(gradientForest)

PATH <- getwd()

#read in data ----
# fish <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/ARMOOK_Fishes_bysite23a_StreamCat_Flow_Full_15km.csv")) %>% select(-`...1`)
# bug <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/BenthicInsect_23_StreamCat_Flow_Full_15km.csv"))
fish <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/fish_bio_env_updated_20231110.csv"))
bug <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/bug_bio_env_updated_20231110.csv"))

all <- list(fish, bug) #combine for common manipulations

all <- lapply(all, function(df){ #gets rid of duplicate columns, keeps 1
  df <- df[!duplicated(t(df))]
  return(df)
})

spp_list <- read_csv(paste0(PATH, "/XX_archive/species_list_original.csv"))
# rm(fish, bug)

#get variable list ----
vars <- list.files(paste0(PATH, "/10_GFOutput/"), full.names = TRUE)
vars <- lapply(vars, function(x){
  tmp <- readRDS(x)
  tmp <- names(tmp$X)
  return(tmp)
})

names(vars) <- gsub(".rds", "", list.files(paste0(PATH, "/10_GFOutput/")))

vars.all <- data.frame(var_type = names(vars), var = I(vars), stringsAsFactors = FALSE) %>%
  unnest(var)
vars.all <- vars.all %>%
  mutate(var_type = gsub("gf\\.(bugs|fish)\\.23\\.(GW|GWlump|Int|RO)[._](full|HIT|LULC)(_cat)?", "\\3", var_type)) %>%
  distinct() %>%
  filter(var_type != "full") %>%
  filter(!(var_type == "LULC" & grepl("[Cc]at", var))) #going to do watershed vars only

##add in new variables ----
#mcmanamay hydro alt vars
hydro.alt <- names(bug)[grepl("pn", names(bug))]
hydro.alt <- c(hydro.alt[grepl(paste0("pn", toupper(vars.all[vars.all$var_type == "HIT", ]$var), collapse = "|"), hydro.alt)], "pnHA_rank", "pnSeasonal")

vars.all <- bind_rows(
  vars.all,
  data.frame(var_type = rep("H_ALT", length(hydro.alt)), var = hydro.alt)
)

#predicted stream temp
strm.temp <- data.frame(var_type = rep("LULC", length(names(bug)[grepl("temp", names(bug))])), var = names(bug)[grepl("temp", names(bug))])
strm.temp <- strm.temp %>%
  filter(grepl("8.5", var)) %>% #since it's only historical, doesn't really matter because the two models are the same
  filter(grepl("cv|ssn", var)) %>% #only doing seasonal avgs
  filter(grepl("mn|cv", var))

vars.all <- bind_rows(vars.all, strm.temp)

rm(hydro.alt, strm.temp, vars)
  
#get ref type + attach to datasets (also attached updated env vars) ----
hit <- read_csv(paste0(PATH, "/02_EnvDat/Updated_HIT_4.23.a.csv"))
hit <- hit[!duplicated(t(hit))]
hit <- hit %>% mutate(type = ifelse(GAGESII_HDI.HYDRO_DISTURB_INDX > 13, "non-ref", "ref"))

all <- lapply(all, function(x){
  x <- x %>%
    select(-any_of(colnames(hit))) %>% #replace updated vars
    left_join(., hit, by = c("STAID" = "STAID_0"))
  return(x)
})

length(unique(c(all[[1]]$STAID, all[[2]]$STAID))) #number of gages used in data

#get gage types within each flow type
n <- bind_rows(all[[1]] %>% select(STAID, type, Flow_type), all[[2]] %>% select(STAID, type, Flow_type)) %>% distinct()
n2 <- n %>% select(STAID, type) %>% distinct()
table(n$Flow_type, n$type)

#run GF models ----
#make function for repeated runs
gf_sub <- function(full.data = all[[1]], #biological data + env data dataframe
                   species.list = spp_list$species, #species of interest (in case less than the full dataframe)
                   env.variable.list = vars.all, #variables of interest; previously pulled out HIT, LULC, and full run
                   flow.type = "Int", #flow type to run, set as NA if all
                   gage.type = "ref", #gage type to run, set as NA if all
                   gf.classification = TRUE, #classification or regression? likely don't change ever for our purposes
                   #gradientForest options:
                   ntree = 999,
                   transform = NULL,
                   compact = TRUE,
                   trace = TRUE,  
                   nbin=201,
                   corr.threshold = 0.5) {
  
  ##remove sites with < 5 spp ----
  full.data <- full.data %>% 
    mutate(loc_tot = rowSums(select(., any_of(species.list)))) %>%
    filter(loc_tot > 5) %>%
    select(-loc_tot)
  
  ##subset by flow + reference ----
  #flow
  if(!any(grepl(flow.type, c("Int", "RO", "GW"), ignore.case = TRUE) | is.na(flow.type))) { stop("flow type not recognized.") }
  if(!is.na(flow.type)) { full.data <- full.data[grepl(flow.type, full.data$Flow_type, ignore.case = TRUE), ] }
  
  #reference
  if(!any(grepl(gage.type, c("ref", "non-ref"), ignore.case = TRUE) | is.na(gage.type))) { stop("gage type not recognized.") }
  if(!is.na(gage.type)) { full.data <- full.data[grepl(paste0("^", gage.type, "$"), full.data$type, ignore.case = TRUE), ] }
  
  ##pull out env vars ----
  ##- remove NA env var rows
  full.data <- full.data %>% filter(!if_any(all_of(env.variable.list), is.na))
  env_col <- full.data[, grepl(paste(paste0("^", env.variable.list, "$"), collapse = "|"), colnames(full.data))]
  if(any(!env.variable.list %in% colnames(env_col))) { 
    warning(paste("The following environmental variables were not in the dataframe:", paste(env.variable.list[!env.variable.list %in% colnames(env_col)], collapse = ", "))) 
  }
  ## remove constant env vars
  col_rem <- sapply(env_col, function(col) all(col == col[1]))
  env_col <- env_col[, !col_rem]
  
  ##pull out species ----
  spp_col <- full.data[, grepl(paste(species.list, collapse = "|"), colnames(full.data))]
  ##- remove species w/ < 10 records
  spp_col <- spp_col[, colSums(spp_col) >= 10]
  # ##- remove species w/ < 5% of max collection records?? ##not sure they did this in the end...
  # spp_col <- spp_col[, colSums(spp_col) >= max(colSums(spp_col))*0.05]
  
  ##set model parameters ----
  n.sites <- nrow(spp_col)
  n.spp <- ncol(spp_col)
  l <- floor(log2(n.sites * 0.368/2))
  
  ##run as classification or regression? ----
  if(gf.classification) {
    spp_col <- spp_col %>% mutate(across(where(is.numeric), factor))
  }
  
  ##run model ----
  message("Running GF on ", n.spp, " species and ", n.sites, " sites\nUsing ", flow.type, " flow and ", gage.type, " gages..\n")
  set.seed(31)
  
  gf.out <- gradientForest(cbind(env_col, spp_col),
                           predictor.vars = colnames(env_col), response.vars = colnames(spp_col),
                           ntree = ntree, transform = transform, compact = compact, trace = trace,  
                           nbin = nbin, maxLevel = l, corr.threshold = corr.threshold)
  
  return(gf.out)
}

#models to run: 
#fish + bug (2)
# all[[1]], all[[2]]
#all flows, int, ro, gw
flow.opt <- c("Int", "RO", "GW")
#all gages, ref, non-ref
gage.opt <- c(NA, "ref", "non-ref")
#all vars, HIT, LULC, H_ALT
var.opt <- list(vars.all$var, 
                vars.all[vars.all$var_type == "HIT", ]$var, 
                vars.all[vars.all$var_type == "LULC", ]$var, 
                vars.all[vars.all$var_type == "H_ALT", ]$var)

#run models
dir.create(paste0(PATH, "/10_GFOutput/", gsub("-", "_", Sys.Date()))) #date directory for future comparisons

for(i in 1:length(all)) {
  for(f in flow.opt) {
    for(g in gage.opt) {
      for(v in 1:length(var.opt)) {
        
        if(v == 1) { var <- "all" }
        if(v == 2) { var <- "HIT" }
        if(v == 3) { var <- "LULC" }
        if(v == 4) { var <- "HyALT" }
        if(i == 1) { bio <- "fish"} else { bio <- "bugs" }
        if(is.na(g)) { ga <- "all" } else { ga <- g }
        
        cat("\nRunning", bio, ":", f, "flow,", ga, "gage,", var, "variable type...\n")
        file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_", f, "_", ga, "gage_", var, ".rds")
        if(file.exists(file.name)) { cat("Model exists, next.\n"); next }
        
        gf.out <- NULL

        gf.out <- gf_sub(full.data = all[[i]], #biological data + env data dataframe
                         species.list = spp_list$species, #species of interest (in case less than the full dataframe)
                         env.variable.list = var.opt[[v]], #variables of interest; previously pulled out HIT, LULC, and full run
                         flow.type = f, #flow type to run, set as NA if all
                         gage.type = g #gage type to run, set as NA if all
        )

        saveRDS(gf.out, file = file.name)

        Sys.sleep(1)

      }
    }
  }
}





# most_important.gf.GW_lump <- names(importance(gf.out))[0:25]
# most_important.gf.GW_lump
# 
# #GF model plots
# plot(gf.out , plot.type = "O", xlab="")
# plot(gf.fish.23.Int.full, plot.type = "O", xlab="")
# 
# ##note: because of the way the function is set up, need to run the following line before running split density plot
# gf.out$call$compact <- TRUE
# gf.out$call$nbin <- nbin
# plot(gf.out, plot.type = "S", imp.vars = most_important.gf.GW_lump,
#      leg.posn = "topright", cex.legend = 0.9, cex.axis = 0.9,
#      cex.lab = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
# 
# plot(gf.out, plot.type = "C", imp.vars = most_important.gf.GW_lump,
#      show.overall = F, legend = T, leg.posn = "topleft",
#      leg.nspecies = 5, cex.lab = 0.9, cex.legend = 0.9, ylim=c(0,0.4),
#      cex.axis = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5,
#                                                             0.5, 0), mar = c(2.5, 1, 0.5, 0.5), omi = c(0,0.3, 0, 0)))
# 
# plot(gf.out, plot.type = "C", imp.vars = most_important.gf.GW_lump,
#      show.species = F, common.scale = T, cex.axis = 1,
#      cex.lab = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5,
#                                                            0.5, 0), mar = c(2.5, 1, 1, 0.5), omi = c(0,0.3, 0, 0)))
# 
# par(oma=c(7,3,1,1),mar=c(2,2,2,2),mfrow=c(2,1))
# plot(gf.out, plot.type = "P", show.names = T, horizontal = F,
#      cex.axis = 1, cex.labels = 0.6, line = 2.5, ylim=c(0,1))
# 
# 
# plot(gf.out, plot.type = "Split.Density")
