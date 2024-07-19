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

#snap to NHD stream ----
##NOTE: ON APRIL 11, 2024, REMOVED SOME NHD STREAMS IN QGIS THAT DIDN'T HAVE STREAMCAT DATA (N=30)
##ALSO EDITED 2 FLOWLINES BY REMOVING A PART OF THE MULTILINE (FOR JOINING), FOR ONE FLOWLINE, ITS WEIRDLY CROPPED BUT IT'S BEYOND THE RANGE OF THE STUDY AREA
PATH <- getwd()

filt.dat <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_20240716.csv"), col_types = cols(lat = col_number(),
                                                                                                       long = col_number(),
                                                                                                       taxa_count = col_number(),
                                                                                                       .default = "c"))

#if redoing the snapping steps:
# filt.dat <- filt.dat %>% select(-c(contains("dist"), contains("flw"), COMID, contains("gage_no"), contains("_new")))

#make spatial w/ projected crs
sf.albers <- st_as_sf(filt.dat, coords = c("long", "lat"), remove = FALSE, crs = 4269)
sf.albers <- st_transform(sf.albers, crs = 5070)

#nhd load in
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% #nhd streams prev. clipped to region (good to use clipped because of MEM usage)
  st_transform(5070) %>%
  filter(!WBAreaType %in% c("LakePond", "Reservoir"))

#assign COMID to each site + calc distance to nearest line (meters)
nr.line <- st_nearest_feature(sf.albers, nhd, check_crs = TRUE) #find nearest nhd strm
nhd_near <- nhd[nr.line,]

dist <- as.vector(st_distance(sf.albers, nhd_near, by_element = TRUE)) #get distance to nearest nhd strm

df <- data.frame(COMID = nhd_near$COMID, dist2strm_m_nhd = round(dist, digits = 3))

#add to occurrence df
filt.dat <- bind_cols(filt.dat, df)

#save jic
write_csv(filt.dat, paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_", gsub("-", "", Sys.Date()), ".csv"))

#sites at the same stream (same NHD code?)
# sites <- filt.dat %>% select(site_id, date, long, lat, COMID) %>% distinct() %>% arrange(COMID)
# View(sites[duplicated(sites$COMID)|duplicated(sites$COMID, fromLast = TRUE),])

rm(nhd, nr.line, nhd_near, dist, df)

#create unique site ids for snapping to flow regimes + gages ----
PATH <- getwd()
filt.dat <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_20240716.csv"), col_types = cols(lat = col_number(),
                                                                                                       long = col_number(),
                                                                                                       taxa_count = col_number(),
                                                                                                       dist2strm_m_nhd = col_number(),
                                                                                                       .default = "c")) %>%
  mutate(dist2strm_m_nhd = round(dist2strm_m_nhd, 3))

##get distances between 'unique' sites ----
sites <- filt.dat %>%
  select(bio_type, long, lat) %>%
  distinct() %>%
  #it's faster to convert here rather than getting distinct sites from sf.albers object
  st_as_sf(., coords = c("long", "lat"), remove = FALSE, crs = 4269) %>%
  st_transform(5070)

site.dist <- st_distance(sites)
# t <- site.dist[1:100,1:100]
t <- site.dist
t <- as.data.frame(t)
t$site_1 <- seq_len(nrow(t))
#pivot longer for filtering etc.
t <- pivot_longer(t, cols = -site_1, names_to = "site_2", values_to = "site_dist_m")
#remove same site distances
t <- t[t$site_1 != t$site_2, ]
#remove large distances that aren't worrisome (less than 5 km)
t$site_dist_m <- round(as.numeric(t$site_dist_m), digits = 3)
t <- t[t$site_dist_m < 5000, ]

#look at all closely matched sites
tmp <- unique(t$site_1, t$site_2)
tmp <- unique(tmp)

sites.xy <- sites[tmp, ]
sites.xy$close_site <- "close_site"

tmp <- left_join(filt.dat, sites.xy %>% st_drop_geometry()) %>%
  mutate(taxa_name = case_when(bio_type == "fish" ~ species,
                               bio_type == "bug" & !is.na(genus) ~ genus,
                               bio_type == "bug" & is.na(genus) & !is.na(tribe) ~ tribe,
                               bio_type == "bug" & is.na(genus) & is.na(tribe) ~ subfamily,
                               T ~ NA)) %>%
  select(-c(phylum, class, order, family, genus, species, superfamily, tribe, subfamily)) %>%
  filter(!is.na(taxa_name))

write_csv(tmp, paste0(PATH, "/01_BioDat/site_check_for_overlap.csv"))

#checking duplicates / overlap
# tmp <- tmp %>%
#   filter(!is.na(close_site)) %>%
#   mutate(taxa_count = ifelse(taxa_count > 0, 1, 0)) %>%
#   select(-contains("date")) %>%
#   distinct() %>%
#   arrange(COMID, taxa_name)
# 
# df <- tmp[duplicated(tmp[c("COMID", "taxa_name")]) | duplicated(tmp[c("COMID", "taxa_name")], fromLast = TRUE), ]

## make site ids ----
#if two sites at the same COMID are within 5 km of each other, combine.
#all others are unique sites.
sites <- sites %>% mutate(site_ind = row_number())
site_df <- data.frame(site_ind = as.numeric(unique( c(t$site_1, t$site_2) )),
                      site_id_new = NA) %>%
  left_join(., st_drop_geometry(sites))

site_match <- t %>%
  #to get rid of duplicates (e.g., switched site1 and site2)
  mutate(site_2 = as.integer(site_2),
         site_a = pmin(site_1, site_2),
         site_b = pmax(site_1, site_2)) %>%
  select(-site_1, -site_2) %>%
  relocate(site_a, site_b) %>%
  distinct()

site_match_out <- data.frame()

#track site_id creation
track_site_bug <- NA
track_site_fish <- NA

# seq_len(nrow(site_df))

for(i in seq_len(nrow(site_df))) {
  #get close site index
  ind_site1 <- site_df[i, ]$site_ind
  #get close site location info
  site1 <- sites[ind_site1, ]
  comid1 <- filt.dat[filt.dat$bio_type == site1$bio_type & filt.dat$long == site1$long & filt.dat$lat == site1$lat, ]$COMID
  comid1 <- unique(comid1)
  #stop check
  if(length(comid1) > 1) { message("too many comids for target site."); break }
  
  #find all pairs for target ind
  search.df <- bind_rows(site_match[site_match$site_a == ind_site1, ],
                         site_match[site_match$site_b == ind_site1, ])
  
  #paired sites
  sites2 <- site_df[site_df$site_ind %in% search.df$site_b, ]
  
  #check bio type match (skip combining if fish + bug combo)
  check_df <- pmap(sites2, function(site_ind, site_id_new, bio_type, long, lat) {
    bio_check <- bio_type == site1$bio_type
    
    #check for COMID match in close sites
    comid2 <- filt.dat[filt.dat$bio_type == bio_type & filt.dat$long == long & filt.dat$lat == lat, ]$COMID
    comid2 <- unique(comid2)
    # return(comid2)
    
    #stop check
    if(length(comid2) > 1) { message("more than 1 comid returned for match site."); break }

    comid_check <- comid2 == comid1

    return(data.frame(site_a = ind_site1, site_b = site_ind, bio_match = bio_check, COMID_match = comid_check))
  })
  
  check_df <- bind_rows(check_df)
  check_df <- check_df[check_df$site_a != check_df$site_b, ]
  
  #create new site ID for combining
  ##first check for previous id assignment
  if(!is.na(site_df[i, ]$site_id_new)) {
    
    site_name <- site_df[i, ]$site_id_new
    
  } else {
    
    site_name <- paste0(site1$bio_type, "_site_NUMBER")
    
    if(all(is.na(get(paste0("track_site_", site1$bio_type))))) {
      site_name <- sub("NUMBER", "1", site_name)
      #add new name to track vector
      assign(paste0("track_site_", site1$bio_type), site_name) 
    } else {
      #get next site number from vector
      n <- length(get(paste0("track_site_", site1$bio_type))) + 1
      site_name <- sub("NUMBER", n, site_name)
      #add new name to site tracking vector
      assign(paste0("track_site_", site1$bio_type), c(get(paste0("track_site_", site1$bio_type)), site_name))
    }
    
  }
  
  #check for bio + comid match, assign name
  comb_sites <- check_df[check_df$bio_match & check_df$COMID_match, ]
  #if matches of bio and COMID
  if(nrow(comb_sites) != 0) {
    #add new site name to target site
    site_df[i, ]$site_id_new <- site_name
    #add new site name to matched sites
    site_df[site_df$site_ind %in% comb_sites$site_b, ]$site_id_new <- site_name
    #add new site name to matched pairs tracking df
    check_df$site_id_new <- site_name
    check_df <- check_df %>%
      mutate(site_id_new = case_when(bio_match & COMID_match ~ site_name, T ~ NA))
  }

  ##bind to output for tracking + saving
  site_match_out <- bind_rows(site_match_out, check_df)
  
  rm(check_df, comb_sites, search.df, site1, sites2, comid1, ind_site1)
  
  if(i == round(nrow(site_df)*0.125, digits = 0)) { message("12.5% complete at ", Sys.time()) }
  if(i == round(nrow(site_df)*0.25, digits = 0)) { message("25% complete at ", Sys.time()) }
  if(i == round(nrow(site_df)*0.375, digits = 0)) { message("37.5% complete at ", Sys.time()) }
  if(i == round(nrow(site_df)*0.5, digits = 0)) { message("50% complete at ", Sys.time()) }
  if(i == round(nrow(site_df)*0.625, digits = 0)) { message("62.5% complete at ", Sys.time()) }
  if(i == round(nrow(site_df)*0.75, digits = 0)) { message("75% complete at ", Sys.time()) }
  if(i == round(nrow(site_df)*0.875, digits = 0)) { message("87.5% complete at ", Sys.time()) }
  if(i == nrow(site_df)) { message("done!") }
  
}

#save jic
site_match_out <- left_join(site_match_out, site_match)
write_csv(site_match_out, paste0(PATH, "/01_BioDat/site_check_new_name_pairs.csv"))
write_csv(site_df, paste0(PATH, "/01_BioDat/site_check_new_name_siteindex.csv"))

##update site coordinates for now combined sites ----
site_df <- site_df[!is.na(site_df$site_id_new),]

sites_new <- sites %>% 
  mutate(site_ind = row_number()) %>%
  filter(site_ind %in% site_df$site_ind)

sites_new <- left_join(sites_new, site_df)

#get midpoint of grouped sites
sites_new <- sites_new %>%
  group_by(site_id_new) %>%
  mutate(geometry = st_centroid(st_union(geometry))) %>%
  ungroup() %>%
  st_transform(4269) %>%
  mutate(long_new = st_coordinates(.)[, "X"],
         lat_new = st_coordinates(.)[, "Y"]) %>%
  select(-site_ind)

#save jic
write_csv(sites_new, paste0(PATH, "/01_BioDat/site_check_new_coords.csv"))

##add site ids and coordinates to full dataset ----
tmp <- filt.dat %>% 
  left_join(., st_drop_geometry(sites_new))

tmp.na.bug <- tmp %>%
  filter(is.na(site_id_new) & bio_type == "bug") %>%
  group_by(long, lat, COMID) %>%
  mutate(site_id_new = cur_group_id()) %>%
  ungroup() %>%
  arrange(site_id_new) %>%
  mutate(site_id_new = paste0("bug_site_", site_id_new + length(track_site_bug)))
tmp.na.fish <- tmp %>%
  filter(is.na(site_id_new) & bio_type == "fish") %>%
  group_by(long, lat, COMID) %>%
  mutate(site_id_new = cur_group_id()) %>%
  ungroup() %>%
  arrange(site_id_new) %>%
  mutate(site_id_new = paste0("fish_site_", site_id_new + length(track_site_fish)))

#add updated dat back in
filt.dat <- bind_rows(tmp.na.bug,
                      tmp.na.fish,
                      tmp[!is.na(tmp$site_id_new),])

filt.dat <- filt.dat %>%
  mutate(long_new = case_when(is.na(long_new) ~ long, T ~ long_new),
         lat_new = case_when(is.na(lat_new) ~ lat, T ~ lat_new)) %>%
  distinct()

#save JIC
write_csv(filt.dat, paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_updatedsites_", gsub("-", "", Sys.Date()), ".csv"))

##date range of new sites (for later use) ----
#for each COMID, find range of sample dates (if available) to check dates
date_range <- filt.dat %>%
  select(date, date_min, date_max, site_id_new) %>%
  mutate(across(contains("date"), ~ as.Date(.x)),
         min_date = pmin(date, date_min, date_max, na.rm = TRUE),
         max_date = pmax(date, date_min, date_max, na.rm = TRUE)) %>%
  select(-c(date_min, date_max)) %>%
  group_by(site_id_new) %>%
  summarise(min_date = ifelse(all(is.na(min_date)), NA, min(min_date, na.rm = TRUE)),
            max_date = ifelse(all(is.na(max_date)), NA, max(max_date, na.rm = TRUE)),
            unique_dates = n_distinct(date),
            approx_nspp = n()) %>%
  mutate(across(c(min_date, max_date), ~ as.Date(.x))) %>%
  # filter(approx_nspp >= 5) %>%
  mutate(yr_range = as.numeric((max_date - min_date)/365),
         yr_range = round(yr_range, 3))

#save jic
write_csv(date_range, paste0(PATH, "/01_BioDat/site_check_date_range.csv"))

#clean up workspace
rm(sf.albers, site_df, sites, sites_new, sites.xy, t, tmp, tmp.na.bug, tmp.na.fish, i, n, site.dist, track_site_bug, track_site_fish)

#snap to classified flow regimes ----
filt.dat <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_updatedsites_20240717.csv"), col_types = cols(lat = col_number(),
                                                                                                                    long = col_number(),
                                                                                                                    lat_new = col_number(),
                                                                                                                    long_new = col_number(),
                                                                                                                    taxa_count = col_number(),
                                                                                                                    dist2strm_m_nhd = col_number(),
                                                                                                                    .default = "c"))

#make spatial w/ projected crs
sf.albers <- st_as_sf(filt.dat, coords = c("long_new", "lat_new"), remove = FALSE, crs = 4269)
sf.albers <- st_transform(sf.albers, crs = 5070)

flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% select(-PopupInfo) %>%
  st_transform(5070)
flw <- flw %>% filter(Flow_type != "BigR")

nr.line <- st_nearest_feature(sf.albers, flw, check_crs = TRUE) #find nearest flw line
flw_near <- flw[nr.line,]

dist <- as.vector(st_distance(sf.albers, flw_near, by_element = TRUE)) #get distance to nearest flw line

df <- data.frame(flw_name = flw_near$Name, 
                 flw_type = flw_near$Flow_type, 
                 dist2strm_m_flw = round(dist, digits = 3))

#add to occurrence df
filt.dat <- bind_cols(filt.dat, df)

rm(flw_near, flw, df)

#snap to gage ----
#gage load in
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))
g.info <- g.info[g.info$final_filter == "keep", ]
#make spatial
g.sf <- st_as_sf(g.info, coords = c("dec_long_va", "dec_lat_va"), remove = FALSE, crs = 4269) %>%
  st_transform(5070)

#find nearest gage within filtered gage set of matching flow type
##split by flw type
filt.flw <- split(filt.dat, filt.dat$flw_type)
g.sf <- split(g.sf, g.sf$flw_type)

#function to apply to each regime
snap_to_gage <- function(occ.dat, gage.dat, by_flw = TRUE) {
  
  #make spatial w/ projected crs
  sf.albers <- st_as_sf(occ.dat, coords = c("long_new", "lat_new"), remove = FALSE, crs = 4269)
  sf.albers <- st_transform(sf.albers, crs = 5070)
  
  #find nearest gage
  nr.gage <- st_nearest_feature(sf.albers, gage.dat, check_crs = TRUE) 
  gage_near <- gage.dat[nr.gage,]
  
  #get distance
  dist <- as.vector(st_distance(sf.albers, gage_near, by_element = TRUE))
  
  #put into a dataframe for binding
  if(by_flw) {
    df <- data.frame(flw_gage_no = gage_near$site_no, dist2_flw_gage_m = round(dist, 3))
  } else {
    df <- data.frame(gage_no = gage_near$site_no, dist2_gage_m = round(dist, 3))
  }

  #bind to occ df
  out <- bind_cols(occ.dat, df)
  
  return(out)
}

RO <- snap_to_gage(filt.flw[["RO"]], g.sf[["RO"]])
GW <- snap_to_gage(filt.flw[["GW"]], g.sf[["GW"]])
Int <- snap_to_gage(filt.flw[["Int"]], g.sf[["Int"]])

filt.dat <- bind_rows(RO, GW, Int)

#find nearest gage disregarding gage flow type
filt.dat <- snap_to_gage(filt.dat, bind_rows(g.sf), by_flw = FALSE)

#save jic
write_csv(filt.dat, paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_allinfo_", gsub("-", "", Sys.Date()), ".csv"))

#save site info only
sites <- filt.dat %>%
  select(bio_type, site_id_new, long, lat, long_new, lat_new, COMID, contains("dist"), flw_type, contains("gage")) %>%
  filter(dist2_flw_gage_m <= 15000 & #gage distance within 15 km
           dist2strm_m_nhd <= 1000 & #less than 1 km from a NHD stream
           dist2strm_m_flw <= 1000) %>%
  distinct()
write_csv(sites, paste0(PATH, "/01_BioDat/sites_alltax_inthigh_allinfo_", gsub("-", "", Sys.Date()), ".csv"))

#final occurrence filter + widen data ----
final.dat <- filt.dat %>%
  filter(dist2_flw_gage_m <= 15000 & #gage distance within 15 km
         dist2strm_m_nhd <= 1000 & #less than 1 km from a NHD stream
         dist2strm_m_flw <= 1000) #less than 1 km from a classified flow regime stream

#set final names, remove ambiguous names
final.dat <- final.dat %>%
  mutate(taxa_name = case_when(bio_type == "fish" ~ species,
                               tribe == "Aciliini" ~ tribe,
                               tribe == "Hydrophilini" ~ tribe,
                               tribe == "Pseudochironomini" ~ tribe,
                               bio_type == "bug" & !is.na(genus) ~ genus,
                               bio_type == "bug" & is.na(genus) & !is.na(tribe) ~ tribe,
                               bio_type == "bug" & is.na(genus) & is.na(tribe) & !is.na(subfamily) ~ subfamily,
                               T ~ NA)) %>%
  filter(!is.na(taxa_name))

#set final lat, long, site id names
final.dat <- final.dat %>%
  select(-site_id, -lat, -long) %>%
  rename(site_id = site_id_new,
         long = long_new,
         lat = lat_new) %>%
  distinct()

#date range of final filter
#not filtering by date right now (17JUL2024)
# test <- final.dat %>%
#   filter(date > "1960-01-01" | is.na(date))
# date_range <- test %>%
#   select(date, date_min, date_max, site_id) %>%
#   mutate(across(contains("date"), ~ as.Date(.x)),
#          min_date = pmin(date, date_min, date_max, na.rm = TRUE),
#          max_date = pmax(date, date_min, date_max, na.rm = TRUE)) %>%
#   select(-c(date_min, date_max)) %>%
#   group_by(site_id) %>%
#   summarise(min_date = ifelse(all(is.na(min_date)), NA, min(min_date, na.rm = TRUE)),
#             max_date = ifelse(all(is.na(max_date)), NA, max(max_date, na.rm = TRUE)),
#             unique_dates = n_distinct(date),
#             approx_nspp = n()) %>%
#   mutate(across(c(min_date, max_date), ~ as.Date(.x))) %>%
#   # filter(approx_nspp >= 5) %>%
#   mutate(yr_range = as.numeric((max_date - min_date)/365),
#          yr_range = round(yr_range, 3))

#remove 'not in region' fishes
final.dat <- final.dat %>%
  filter(!taxa_name %in% c("Dicentrarchus punctatus", "Erimystax dissimilis", "Etheostoma duryi", "Lepomis punctatus", "Macrhybopsis aestivalis",
                           "Notropis chrosomus", "Notropis rubellus", "Phoxinus phoxinus"))

##save final dataset (long format)
write_csv(final.dat, paste0(PATH, "/01_BioDat/occ_alltax_finalfilter_long_", gsub("-", "", Sys.Date()), ".csv"))

##wide format data ----
#prep function
widen_final_occ <- function(occ_dat, bio_target) {
  
  taxa.dat <- occ_dat %>%
    filter(bio_type == bio_target) %>%
    select(-c(contains("date"), source, bio_type, phylum, class, order, family, superfamily, subfamily, tribe, genus, species, contains("dist"), gage_no, flw_name)) %>%
    mutate(taxa_name = sub(" ", "_", taxa_name),
           taxa_count = 1) %>%
    distinct() %>%
    arrange(taxa_name) %>%
    pivot_wider(names_from = taxa_name, values_from = taxa_count) %>%
    mutate(across(-c(COMID, site_id, long, lat, flw_type, flw_gage_no), ~ ifelse(is.na(.x), 0, .x)))
  
  ##remove sites with fewer than 5 species
  species.list <- names(taxa.dat)[!names(taxa.dat) %in% c("COMID", "site_id", "long", "lat", "flw_type", "flw_gage_no")]
  taxa.dat <- taxa.dat %>% 
    mutate(loc_tot = rowSums(select(., any_of(species.list)))) 
  
  message(nrow(taxa.dat[taxa.dat$loc_tot < 5, ]), " sites with fewer than 5 taxa removed...")
  
  taxa.dat <- taxa.dat %>%
    filter(loc_tot >= 5) %>%
    select(-loc_tot)
  
  ##remove species that occur at every site (due to gradient forest warning)
  spp_col <- taxa.dat[, grepl(paste(species.list, collapse = "|"), colnames(taxa.dat))]
  spp_col <- spp_col[, colSums(spp_col) != nrow(spp_col)]
  removed_spp <- species.list[!species.list %in% names(spp_col)]
  message(bio_target, " taxa removed: ", removed_spp)
  
  #add back to site info
  taxa.dat <- bind_cols(taxa.dat[, !grepl(paste(species.list, collapse = "|"), colnames(taxa.dat))], spp_col)
  
  return(taxa.dat)
  
}
###FISH ----

fish <- widen_final_occ(final.dat, bio_target = "fish")

write_csv(fish, paste0(PATH, "/01_BioDat/occ_fish_filtered_wide_", gsub("-", "", Sys.Date()), ".csv"))

###BUGS ----

bug <- widen_final_occ(final.dat, bio_target = "bug")

write_csv(bug, paste0(PATH, "/01_BioDat/occ_bug_filtered_wide_", gsub("-", "", Sys.Date()), ".csv"))

##check distro of flw types
table(bug$flw_type)
table(fish$flw_type)

#get species list (all tax) ----
spp <- final.dat %>%
  select(bio_type, phylum, class, order, family, subfamily, tribe, genus, species, taxa_name) %>%
  distinct()
write_csv(spp, paste0(PATH, "/01_BioDat/species_list_complete_20240717.csv"))
