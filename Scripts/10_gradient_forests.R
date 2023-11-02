#gradient forests
library(tidyverse); library(gradientForest)

PATH <- getwd()

#read in data ----
fish <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/ARMOOK_Fishes_bysite23a_StreamCat_Flow_Full_15km.csv")) %>% select(-`...1`)
bug <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/BenthicInsect_23_StreamCat_Flow_Full_15km.csv"))
all <- list(fish, bug) #combine for common manipulations

all <- lapply(all, function(df){ #gets rid of duplicate columns, keeps 1
  df <- df[!duplicated(t(df))]
  return(df)
})

spp_list <- read_csv(paste0(PATH, "/species_list_original.csv"))
# rm(fish, bug)

#get variable list ----
vars <- list.files(paste0(PATH, "/10_GFOutput/"), full.names = TRUE)
vars <- lapply(vars, function(x){
  tmp <- readRDS(x)
  tmp <- names(tmp$X)
  return(tmp)
})

names(vars) <- gsub(".rds", "", list.files(paste0(PATH, "/10_GFOutput/")))
vars.all <- unique(unlist(vars))

#get ref type + attach to datasets (also attached updated env vars)
hit <- read_csv(paste0(PATH, "/02_EnvDat/Updated_HIT_4.23.a.csv"))
hit <- hit[!duplicated(t(hit))]
hit <- hit %>% mutate(type = ifelse(GAGESII_HDI.HYDRO_DISTURB_INDX > 13, "non-ref", "ref"))

all <- lapply(all, function(x){
  x <- x %>%
    select(-any_of(colnames(hit))) %>% #replace updated vars
    left_join(., hit, by = c("STAID" = "STAID...2"))
  return(x)
})

length(unique(c(all[[1]]$STAID, all[[2]]$STAID))) #number of gages used in data

#get gage types within each flow type
n <- bind_rows(all[[1]] %>% select(STAID, type, Flow_type), all[[2]] %>% select(STAID, type, Flow_type)) %>% distinct()
n2 <- n %>% select(STAID, type) %>% distinct()
table(n$Flow_type, n$type)

#run GF models ----
flow <- c("Int", "RO", "GW")

gf_sub <- function(full.data = all[[1]],
                   species.list = spp_list$species, #previously pulled out HIT, LULC, and full run, also did Cat level runs for each
                   env.variable.list = vars.all,
                   flow.type = "Int",
                   gage.type = "ref",
                   gf.classification = TRUE,
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













most_important.gf.GW_lump <- names(importance(gf.out))[0:25]
most_important.gf.GW_lump

#GF model plots
plot(gf.out , plot.type = "O", xlab="")

##note: because of the way the function is set up, need to run the following line before running split density plot
gf.out$call$compact <- TRUE
gf.out$call$nbin <- nbin
plot(gf.out, plot.type = "S", imp.vars = most_important.gf.GW_lump,
     leg.posn = "topright", cex.legend = 0.9, cex.axis = 0.9,
     cex.lab = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))

plot(gf.out, plot.type = "C", imp.vars = most_important.gf.GW_lump,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.9, cex.legend = 0.9, ylim=c(0,0.4),
     cex.axis = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                            0.5, 0), mar = c(2.5, 1, 0.5, 0.5), omi = c(0,0.3, 0, 0)))

plot(gf.out, plot.type = "C", imp.vars = most_important.gf.GW_lump,
     show.species = F, common.scale = T, cex.axis = 1,
     cex.lab = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                           0.5, 0), mar = c(2.5, 1, 1, 0.5), omi = c(0,0.3, 0, 0)))

par(oma=c(7,3,1,1),mar=c(2,2,2,2),mfrow=c(2,1))
plot(gf.out, plot.type = "P", show.names = T, horizontal = F,
     cex.axis = 1, cex.labels = 0.6, line = 2.5, ylim=c(0,1))


plot(gf.out, plot.type = "Split.Density")
