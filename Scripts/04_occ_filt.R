# OCCURRENCE FILTERING + TAXONOMY FINALIZING

library(tidyverse); library(sf)

PATH <- getwd()

#load in originals, with updated names, combine new datasets ----
occ.dat.orig <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_updatedspp_nofilt_updated20240621.csv"),  col_types = cols(lat = col_number(),
                                                                                                                        long = col_number(),
                                                                                                                        taxa_count = col_number(),
                                                                                                                        .default = "c"))

#remove old names ----
filt.dat <- occ.dat.orig[, !grepl("_orig", names(occ.dat.orig))] %>% distinct()

#remove names flagged for removal ----
#remove hybrids
filt.dat <- filt.dat[is.na(filt.dat$flag_hybrid), ] 
filt.dat <- filt.dat[, names(filt.dat) != "flag_hybrid"]
#remove unidentified etc.
tmp <- apply(filt.dat, 1, function(row) any(grepl("remove", row, ignore.case = TRUE))) #find rows where any column has 'remove' in it
filt.dat <- filt.dat[!tmp, ]

#remove non-target taxa ----
##phylum
filt.dat <- filt.dat[filt.dat$phylum %in% c("Arthropoda", "Chordata"), ]
##class
filt.dat <- filt.dat[!filt.dat$class %in% c("Mammalia", "Amphibia", "Maxillopoda", "Aves", "Branchiopoda", "Collembola", "Copepoda", "Ostracoda", "Arachnida", "Malacostraca"), ]
##family
filt.dat <- filt.dat[!filt.dat$family %in% c("Chelydridae"),]

#fix some coordinates
filt.dat <- filt.dat %>% mutate(long = ifelse(long > 0, long*-1, long))

#fix bio_type for grouping later
filt.dat <- filt.dat %>%
  mutate(bio_type = case_when(phylum == "Chordata" ~ "fish",
                              class == "Insecta" ~ "bug",
                              T ~ "neither"))

#remove sites outside target region ----
##ecoregion shapefile made in XX_data_gathering.R
eco <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "ecoreg_l3_interior_highlands_crop") %>%
  st_transform(5070) #change crs for buffer

eco.buff <- st_buffer(eco, dist = 16000) #16 km buffer for wiggle room for later filters

##convert bio dat to spatial
sf.dat <- st_as_sf(filt.dat, coords = c("long", "lat"), remove = FALSE, crs = 4269) %>%
  st_transform(5070)

##remove points outside int. high region
sf.dat <- st_filter(sf.dat, eco.buff)

##check map
# ggplot() +
#   geom_sf(data = sf.dat) +
#   geom_sf(data = eco, fill = NA, color = "gray", lwd = 3) +
#   theme_minimal()

filt.dat <- sf.dat %>% st_drop_geometry()

##save up to regional filter
write_csv(filt.dat, paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_", gsub("-", "", Sys.Date()), ".csv"))

rm(eco, eco.buff, occ.dat.orig, sf.dat)

#find + snap nearest nhd, flow, and gages ----
##NOTE: ON APRIL 11, 2024, REMOVED SOME NHD STREAMS IN QGIS THAT DIDN'T HAVE STREAMCAT DATA (N=30)
##ALSO EDITED 2 FLOWLINES BY REMOVING A PART OF THE MULTILINE (FOR JOINING), FOR ONE FLOWLINE, ITS WEIRDLY CROPPED BUT IT'S BEYOND THE RANGE OF THE STUDY AREA
PATH <- getwd()

filt.dat <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_20240626.csv"), col_types = cols(lat = col_number(),
                                                                                                       long = col_number(),
                                                                                                       taxa_count = col_number(),
                                                                                                       .default = "c"))

#if redoing the snapping steps:
# filt.dat <- filt.dat %>% select(-c(contains("dist"), contains("flw"), COMID, contains("gage_no"), contains("_new")))

#make spatial w/ projected crs
sf.albers <- st_as_sf(filt.dat, coords = c("long", "lat"), remove = FALSE, crs = 4269)
sf.albers <- st_transform(sf.albers, crs = 5070)

##snap to NHD stream ----
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% #nhd streams prev. clipped to region (good to use clipped because of MEM usage)
  st_transform(5070) %>%
  filter(!WBAreaType %in% c("LakePond", "Reservoir"))

#assign COMID to each site + calc distance to nearest line (meters)
nr.line <- st_nearest_feature(sf.albers, nhd, check_crs = TRUE) #find nearest nhd strm
nhd_near <- nhd[nr.line,]

dist <- as.vector(st_distance(sf.albers, nhd_near, by_element = TRUE)) #get distance to nearest nhd strm

df <- data.frame(COMID = nhd_near$COMID, dist2strm_m_nhd = dist)

#add to occurrence df
filt.dat <- bind_cols(filt.dat, df)

#sites at the same stream (same NHD code?)
# sites <- filt.dat %>% select(site_id, date, long, lat, COMID) %>% distinct() %>% arrange(COMID)
# View(sites[duplicated(sites$COMID)|duplicated(sites$COMID, fromLast = TRUE),])

###check distances between 'unique' sites ----
# sites <- filt.dat %>% 
#   select(bio_type, long, lat) %>%
#   distinct() %>%
#   #it's faster to convert here rather than getting distinct sites from sf.albers object
#   st_as_sf(., coords = c("long", "lat"), remove = FALSE, crs = 4269) %>%
#   st_transform(5070)
# 
# site.dist <- st_distance(sites)
# # t <- site.dist[1:100,1:100]
# t <- site.dist
# t <- as.data.frame(t)
# t$site_1 <- seq_len(nrow(t))
# #pivot longer for filtering etc.
# t <- pivot_longer(t, cols = -site_1, names_to = "site_2", values_to = "site_dist_m")
# #remove same site distances
# t <- t[t$site_1 != t$site_2, ] 
# #remove large distances that aren't worrisome (less than 5 km)
# t$site_dist_m <- round(as.numeric(t$site_dist_m), digits = 3)
# t <- t[t$site_dist_m < 1000, ]
# 
# #look at all closely matched sites
# tmp <- unique(t$site_1, t$site_2)
# tmp <- unique(tmp)
# 
# sites.xy <- sites[tmp, ]
# sites.xy$close_site <- "close_site"
# 
# tmp <- left_join(filt.dat, sites.xy %>% st_drop_geometry()) %>% 
#   filter(!is.na(close_site)) %>%
#   mutate(taxa_name = case_when(bio_type == "fish" ~ species,
#                                bio_type == "bug" & !is.na(genus) ~ genus,
#                                bio_type == "bug" & is.na(genus) & !is.na(tribe) ~ tribe,
#                                bio_type == "bug" & is.na(genus) & is.na(tribe) ~ subfamily,
#                                T ~ NA)) %>%
#   select(-c(phylum, class, order, family, genus, species, superfamily, tribe, subfamily)) %>%
#   filter(!is.na(taxa_name))
# 
# write_csv(tmp, paste0(PATH, "/01_BioDat/site_check_for_overlap.csv"))

rm(nhd_near, df, dist, nr.line, df, x, tmp, sites, sites.xy, t, tmp, zero.dist, zero.sites)

##snap to classified flow regimes ----
flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% select(-PopupInfo) %>%
  st_transform(5070)
flw <- flw %>% filter(Flow_type != "BigR")
  
nr.line <- st_nearest_feature(sf.albers, flw, check_crs = TRUE) #find nearest flw line
flw_near <- flw[nr.line,]

dist <- as.vector(st_distance(sf.albers, flw_near, by_element = TRUE)) #get distance to nearest flw line

df <- data.frame(flw_name = flw_near$Name, 
                 flw_type = flw_near$Flow_type, 
                 dist2strm_m_flw = dist)

#add to occurrence df
filt.dat <- bind_cols(filt.dat, df)

rm(flw_near, flw, df)

##snap to gage ----
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))
g.sf <- st_as_sf(g.info, coords = c("dec_long_va", "dec_lat_va"), remove = FALSE, crs = 4269) %>%
  st_transform(5070)

snap_gage <- function(x, yr_to_clip, df_gage_obj) { #yr_to_clip is so I can attach either gages with 15 yrs of data OR with 10 years of data
  
  res_list <- vector("list", length(yr_to_clip))
  
  for(i in seq_along(yr_to_clip)) {
    yr_min <- yr_to_clip[i]
    
    #filter gages
    g.sub.sf <- df_gage_obj %>% filter(hit_ttl_yrs >= yr_min)
    
    #find nearest gage within filtered gage
    nr.gage <- st_nearest_feature(x, g.sub.sf, check_crs = TRUE) 
    gage_near <- g.sub.sf[nr.gage,]
    
    #get distance
    dist <- as.vector(st_distance(x, gage_near, by_element = TRUE))
    
    #put into a dataframe for return
    res_list[[i]] <- data.frame(gage_no = gage_near$site_no, dist2gage_m = dist)
    names(res_list[[i]]) <- paste0(names(res_list[[i]]), "_", yr_min, "yr")
  }
  
  df <- do.call(cbind, res_list)
  
  return(df)
  
}

near_gage <- lapply(sf.albers, function(x) snap_gage(x, c(10, 15), g.sf))

#add to occurrence df
occ.list[[1]] <- bind_cols(occ.list[[1]], near_gage[[1]])
occ.list[[2]] <- bind_cols(occ.list[[2]], near_gage[[2]])

occ.list[[1]] <- occ.list[[1]] %>% mutate(across(contains("dist2"), ~ round(.x, digits = 2)))
occ.list[[2]] <- occ.list[[2]] %>% mutate(across(contains("dist2"), ~ round(.x, digits = 2)))

# t <- occ.list[[1]] %>% filter(dist2gage_m_15yr <= 15000) %>% select(COMID, gage_no_15yr) %>% distinct()

#save occurrence datasets with new info
write_csv(occ.list[[1]], paste0(PATH, "/01_BioDat/occ_bug_inthigh_long.csv"))
write_csv(occ.list[[2]], paste0(PATH, "/01_BioDat/occ_fish_inthigh_long.csv"))

rm(i, sf.albers, near_gage, g.sf)


##get midpoint coordinate for COMIDs ----
#reduce NHD datset
comids <- lapply(occ.list, function(x) x$COMID)
comids <- unique( c(comids[[1]], comids[[2]]) )
nhd.sub <- nhd[nhd$COMID %in% comids, ]

#find midpoint along all lines
nhd.linestring <- st_cast(nhd.sub, "LINESTRING")

nhd.midpoint <- st_line_sample(nhd.linestring, sample = 0.5) %>% 
  st_transform(4269) #get NAD83 x, y

nhd.midpoint <- data.frame(COMID = nhd.linestring$COMID,
                           long_new = st_coordinates(nhd.midpoint)[, "X"],
                           lat_new = st_coordinates(nhd.midpoint)[, "Y"])

nhd.midpoint <- nhd.midpoint[!duplicated(nhd.midpoint$COMID),] #note, this will just take first duplicate, which in this case is fine but other cases a more nuanced approach might be needed

#add new lat long to occurrence dataframe
occ.list <- lapply(occ.list, function(x) { left_join(x, nhd.midpoint) })

rm(nhd, nhd.midpoint, nhd.linestring)

# #plot check
# ggplot() +
#   geom_sf(data = nhd.sub[2,]) +
#   geom_sf(data = st_as_sf(nhd.midpoint[2,], coords = c("long_new", "lat_new"), crs = 4269), color = "red") +
#   theme_minimal()

#save
write_csv(occ.list[[1]], paste0(PATH, "/01_BioDat/occ_bug_inthigh_long.csv"))
write_csv(occ.list[[2]], paste0(PATH, "/01_BioDat/occ_fish_inthigh_long.csv"))

####get distance between old and new coordinates for later (JIC) ----
old.pt <- lapply(occ.list, st_as_sf, coords = c("long", "lat"), remove = FALSE, crs = 4269)
new.pt <- lapply(occ.list, st_as_sf, coords = c("long_new", "lat_new"), remove = FALSE, crs = 4269)

dist <- vector("list", length(old.pt))
for(i in seq_along(old.pt)) {
  
  dist[[i]] <- st_distance(old.pt[[i]], new.pt[[i]], by_element = TRUE) #get distance between coordinates
  
  dist[[i]] <- bind_cols(occ.list[[i]], dist[[i]])
  
  colnames(dist[[i]])[ncol(dist[[i]])] <- "dist_new_old_coord"
  
  dist[[i]] <- dist[[i]] %>%
    select(COMID, matches("lat|long"), dist_new_old_coord) %>%
    distinct()
  
}

lapply(dist, function(x) max(x$dist_new_old_coord)) #get furthest distance b/w old/new coords

write_csv(dist[[1]], paste0(PATH, "/01_BioDat/occ_bug_dist_updated_coords.csv"))
write_csv(dist[[2]], paste0(PATH, "/01_BioDat/occ_fish_dist_updated_coords.csv"))

rm(new.pt, old.pt, dist)

##remake spatial w/ projected crs + new consolidated coordinates
sf.albers <- lapply(occ.list, st_as_sf, coords = c("long_new", "lat_new"), remove = FALSE, crs = 4269)
sf.albers <- lapply(sf.albers, st_transform, crs = 5070)


#get species list (all tax) ----
spp <- filt.dat %>%
  select(phylum, class, order, family, subfamily, tribe, genus, species) %>%
  distinct()
write_csv(spp, paste0(PATH, "/01_BioDat/species_list_complete_20240205.csv"))