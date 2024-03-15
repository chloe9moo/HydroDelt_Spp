# OCCURRENCE FILTERING + TAXONOMY FINALIZING

library(tidyverse); library(sf)

PATH <- getwd()

#load in originals, with updated names, combine new datasets ----
occ.dat.orig <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/occ_alltax_updatedspp_nofilt_20231128.csv"),  col_types = cols(lat = col_number(),
                                                                                                                                 long = col_number(),
                                                                                                                                 taxa_count = col_number(),
                                                                                                                                 .default = "c")) %>%
  select(-bio_site_no)

wqp <- read_csv(paste0(PATH, "/01_BioDat/occ_wqp_alltax_updatedspp_nofilt_20240129.csv"), col_types = cols(lat = col_number(),
                                                                                                           long = col_number(),
                                                                                                           taxa_count = col_number(),
                                                                                                           .default = "c")) %>% 
  select(-ResultCommentText) %>%
  mutate(phylum = case_when(order == "Podocopa" ~ "Arthropoda", TRUE ~ phylum))

old_dat <- lapply(list.files(path = paste0(PATH, "/01_BioDat/archive_occ_dat/"), pattern = "corrected_updatedsppnames", full.names = TRUE), 
                  read_csv, col_types = cols(lat = col_number(), long = col_number(), .default = "c"))
old_dat <- bind_rows(old_dat)
old_dat <- old_dat %>%
  mutate(kingdom = "Animalia", phylum = case_when(!is.na(species) ~ "Chordata", T ~ phylum)) #these should all be fish


occ.dat.orig <- bind_rows(occ.dat.orig, wqp, old_dat) %>%
  distinct()
  
write_csv(occ.dat.orig, paste0(PATH, "/01_BioDat/occ_alltax_updatedspp_nofilt_updated20240205.csv"))
rm(wqp, old_dat)

occ.dat.orig <- occ.dat.orig %>%
  #fix some names
  mutate(#taxa_name_check = ifelse(is.na(taxa_name_orig), 1, 0), #use this to organize for checking names
         taxa_name_orig = case_when(is.na(taxa_name_orig) ~ genus_orig,
                                    T ~ taxa_name_orig),
         tribe = case_when(genus == "Acilius" ~ "Aciliini",
                           taxa_name_orig == "Chironomini" ~ "Chironomini",
                           taxa_name_orig == "Eriopterini" ~ "Eriopterini",
                           taxa_name_orig == "Hydroporini" ~ "Hydroporini",
                           taxa_name_orig == "Limoniini" ~ "Limoniini",
                           taxa_name_orig == "Palpomyiini" ~ "Palpomyiini",
                           taxa_name_orig == "Pentaneurini" ~ "Pentaneurini",
                           taxa_name_orig == "Stratiomyini" ~ "Stratiomyini",
                           taxa_name_orig == "Tanytarsini" ~ "Tanytarsini",
                           taxa_name_orig == "Pseudochironomini" ~ "Pseudochironomini",
                           genus == "Hydroporus" ~ "Hydroporini", 
                           genus == "Limonia" ~ "Limoniini",
                           genus == "Chironomus" ~ "Chironomini",
                           genus == "Djalmabatista" ~ "Procladiini",
                           TRUE ~ tribe),
         genus = case_when(genus == "Mastor" ~ "Amblyscirtes",
                           taxa_name_orig == "Anisoptera" ~ "Anisoptera",
                           genus == "Gomphurus" ~ "Gomphus",
                           taxa_name_orig == "Nematocera" ~ NA,
                           taxa_name_orig == "Neolophylax" ~ "Neophylax",
                           taxa_name_orig == "Phycentropus" ~ "Phylocentropus",
                           taxa_name_orig == "Smittia" ~ "Smittia",
                           grepl("ammocoete", taxa_name_orig, ignore.case = TRUE) ~ NA,
                           TRUE ~ genus),
         family = case_when(taxa_name_orig == "Anisoptera" ~ "Libellulidae",
                            tribe == "Chironomini" ~ "Chironomidae",
                            taxa_name_orig == "Dropsychid" ~ "Hydropsychidae",
                            taxa_name_orig == "Eriopterini" ~ "Tipulidae",
                            taxa_name_orig == "Forcipomyiinae" ~ "Ceratopogonidae",
                            taxa_name_orig == "Hydroporini" ~ "Dytiscidae",
                            taxa_name_orig == "Limoniini" ~ "Tipulidae",
                            taxa_name_orig == "Nematocera" ~ NA,
                            taxa_name_orig == "Neolophylax" ~ "Uenoidae",
                            taxa_name_orig == "Orthocladiinae" ~ "Chironomidae",
                            taxa_name_orig == "Palpomyiini" ~ "Ceratopogonidae",
                            taxa_name_orig == "Pentaneurini" ~ "Chironomidae",
                            taxa_name_orig == "Phycentropus" ~ "Dipseudopsidae",
                            taxa_name_orig == "Smittia" ~ "Chironomidae",
                            taxa_name_orig == "Stratiomyini" ~ "Stratiomyidae",
                            taxa_name_orig == "Tanypodinae" ~ "Chironomidae",
                            taxa_name_orig == "Tanytarsini" ~ "Chironomidae",
                            genus == "Pentagenia" ~ "Palingeniidae",
                            grepl("ammocoete", taxa_name_orig, ignore.case = TRUE) ~ "Petromyzontidae",
                            grepl("ammocoete", species_orig, ignore.case = TRUE) ~ "Petromyzontidae",
                            TRUE ~ family),
         order = case_when(order == "Acari" ~ NA,
                           genus == "Podura" ~ "Poduromorpha", 
                           TRUE ~ order),
         subfamily = case_when(taxa_name_orig == "Ceratopogoninae" ~ "Ceratopogoninae",
                               taxa_name_orig == "Chironominae" ~ "Chironominae",
                               taxa_name_orig == "Forcipomyiinae" ~ "Forcipomyiinae",
                               taxa_name_orig == "Gerrinae" ~ "Gerrinae",
                               taxa_name_orig == "Hemerodromiinae" ~ "Hemerodromiinae",
                               taxa_name_orig == "Hydroporini" ~ "Hydroporinae",
                               taxa_name_orig == "Limoniinae" ~ "Limoniinae",
                               taxa_name_orig == "Orthocladiinae" ~ "Orthocladiinae",
                               taxa_name_orig == "Psychodinae" ~ "Psychodinae",
                               taxa_name_orig == "Tanypodinae" ~ "Tanypodinae",
                               taxa_name_orig == "Tanytarsini" ~ "Chironominae",
                               taxa_name_orig == "Diamesinae" ~ "Diamesinae",
                               genus == "Hydroporus" ~ "Hydroporinae",
                               genus == "Chironomus" ~ "Chironominae",
                               tribe %in% c("Chironomini", "Tanytarsini") ~ "Chironominae",
                               genus == "Djalmabatista" ~ "Tanypodinae",
                               TRUE ~ subfamily),
         suborder = case_when(taxa_name_orig == "Nematocera" ~ "Nematocera",
                              TRUE ~ suborder)) 

spp <- occ.dat.orig %>%
  select(kingdom, phylum, class, subclass, superorder, order, suborder, infraorder,
         superfamily, family, subfamily, tribe, genus, subgenus, species, subspecies,
         contains("_orig")) %>%
  distinct()
write_csv(spp, paste0(PATH, "/01_BioDat/species_list_complete_20240205.csv"))

occ.dat <- occ.dat.orig %>%
  select(lat, long, contains("date"), taxa_count,
         kingdom, phylum, 
         class, subclass,
         superorder, order, suborder, infraorder,
         superfamily, family, subfamily, 
         tribe, genus, subgenus, species, subspecies, source,
         -contains("orig")) %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) #get rid of zero count ids

# table(occ.dat$phylum)

#get fish dat ----
fish <- occ.dat %>%
  filter(phylum == "Chordata") %>% #verts
  filter(!class %in% c("Aves", "Amphibia", "Mammalia") & order != "Testudines") %>% #remove non-fish verts
  select(-contains("sub"), -contains("super"), -contains("infra"), -c(kingdom, phylum, class, tribe)) %>%
  filter(!is.na(species)) #will def. be going to species, so removing ambiguous names here
  # rowid_to_column("key")

#double check names with original names!
#double check sources with original runs
tmp1 <- occ.dat.orig %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) %>% #get rid of zero count ids
  filter(phylum == "Chordata") %>% #verts
  filter(!class %in% c("Aves", "Amphibia", "Mammalia") & order != "Testudines") %>%
  select(matches("family|genus|species|_orig"), -matches("super|infra"), -subgenus, -subfamily, -subspecies, -order_orig) %>%
  distinct()
#   select(source)
# table(tmp$source)

#get species list with record count
tmp <- fish %>% group_by(order, family, genus, species) %>%
  summarise(num_records_ttl = n())

write_csv(tmp, paste0(PATH, "/01_BioDat/fish_species_list.csv"))
write_csv(fish, paste0(PATH, "/01_BioDat/occ_fish_all_long_20240205.csv"))

#get invert dat ----
t <- read_csv(paste0(PATH, "/01_BioDat/tribe_info.csv")) %>% rename(tribe_new = tribe) #for updating tribes

bug <- occ.dat %>%
  filter(class == "Insecta") %>%
  select(-c(kingdom, phylum, class, subclass, superorder, suborder, infraorder, superfamily, subgenus, species, subspecies)) %>%
  #update tribes in case that is the lowest tax level chosen
  left_join(., t) %>%
  mutate(tribe_new = case_when(is.na(tribe_new) ~ tribe,
                               T ~ tribe_new)) %>%
  select(-tribe) %>% rename(tribe = tribe_new) %>%
  filter(!is.na(family)) #remove ambiguous names (min. req. for now is family)
  

#double check names with original names!
#double check sources with original runs
tmp1 <- occ.dat.orig %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) %>% #get rid of zero count ids
  filter(class == "Insecta") %>%
  select(matches("order|family|subfamily|tribe|genus|_orig|_check"), -matches("super|infra"), -suborder, -subgenus, -species_orig) %>%
  distinct()
  # select(source)
# table(tmp$source)

#get species list
tmp <- bug %>% group_by(order, family, subfamily, tribe, genus) %>%
  summarise(num_records_ttl = n())

write_csv(tmp, paste0(PATH, "/01_BioDat/bug_species_list.csv"))
write_csv(bug, paste0(PATH, "/01_BioDat/occ_bug_all_long_20240205.csv"))

#region filter ----
eco <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "ecoreg_l3_interior_highlands_crop")

sf.occ.dat <- lapply(list(bug, fish), function(x) { #convert bio dat to spatial dat
  x <- st_as_sf(x, coords = c("long", "lat"), remove = FALSE, crs = 4269)
})

eco.buff <- st_buffer(eco, dist = 16000) #16 km buffer for wiggle room

sf.occ.dat <- lapply(sf.occ.dat, function(x) { #remove occurrence points falling outside int. high. region
  x <- st_filter(x, eco.buff) 
})

# ggplot() +
#   geom_sf(data = sf.occ.dat[[2]]) +
#   geom_sf(data = eco, fill = NA, color = "gray", lwd = 3) +
#   theme_minimal()

for(i in 1:length(sf.occ.dat)) { #save regional filter dataset
  x <- sf.occ.dat[[i]] %>% st_drop_geometry()

  if(i == 1) { taxa <- "bug"} else { taxa <- "fish" }
  
  write_csv(x, paste0(PATH, "/01_BioDat/occ_", taxa, "_inthigh_long_", gsub("-", "", Sys.Date()), ".csv"))
}
rm(x, taxa, i)

rm(list = ls())

#fix coordinate error ----
PATH <- getwd()

file.list <- list.files(paste0(PATH, "/01_BioDat/"), pattern = "inthigh_long", full.names = TRUE)
occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                         long = col_number(),
                                                         taxa_count = col_number(),
                                                         .default = "c"))

#make unique sites spatial w/ projected crs
# occ.adj <- lapply(occ.list, function(x) {
#   sites <- x %>%
#     select(long, lat) %>% distinct() %>%
#     st_as_sf(., coords = c("long", "lat"), crs = 4269, remove = FALSE) %>%
#     st_transform(5070)
#   
#   #get buffer (size = amount of uncertainty in coordinates) + group by intersecting buffers
#   buff <- sites %>%
#     st_buffer(., dist = 5) %>% #5 meters, the original crs has 2m accuracy, but went up for bit more wiggle room 
#     st_union()
#   buff <- st_cast(buff, "POLYGON") #break down into ind. polygons
#   
#   #get centroid of unioned buffers + coordinates for using later
#   b.centroid <- st_as_sf(st_centroid(buff)) %>%
#     st_transform(4269) %>%
#     mutate(group_id = row_number(),
#            long_new = st_coordinates(.)[, "X"],
#            lat_new = st_coordinates(.)[, "Y"]) 
#   
#   s.intersect <- as.data.frame(st_intersects(sites, buff)) #find which sites fall within each unioned buffer polygon
#   
#   names(s.intersect) <- c("row_num", "group_id")
#   
#   #add new coordinates to occ df
#   sites <- do.call(cbind, list(sites, s.intersect))
#   
#   sites <- sites %>% left_join(., st_drop_geometry(b.centroid)) %>% select(-c(row_num, group_id)) %>% st_drop_geometry()
#   
#   df <- x %>% left_join(., sites)
#   
#   return(df)
# })
# 
# #checking it works
# # y <- unique(which.max(table(s.intersect$group_id)))
# # t <- filter(s.intersect, group_id == y)
# # ggplot() +
# #   geom_sf(data = buff[y]) +
# #   geom_sf(data = sites[sites$row_num %in% t$row_num, ], fill = "black", shape = 21, alpha = 0.6, size = 4) +
# #   geom_sf(data = b.centroid[y,], color = "red", shape = 4, size = 3, stroke = 1.5) +
# #   theme_minimal()
# # ggsave(paste0(PATH, "/99_figures/ex_site_uncertainty.png"), width = 5, height = 5, bg = "white")
# 
# #save updates
# for(i in seq_along(occ.adj)) { #save regional filter dataset
#   if(i == 1) { taxa <- "bug"} else { taxa <- "fish" }
#   
#   write_csv(occ.adj[[i]], paste0(PATH, "/01_BioDat/occ_", taxa, "_inthigh_long_", gsub("-", "", Sys.Date()), ".csv"))
# }
# 
# rm(list = ls())

#find + snap nearest nhd, flow, and gages ----
PATH <- getwd()
file.list <- list.files(paste0(PATH, "/01_BioDat/"), pattern = "inthigh_long", full.names = TRUE)
occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                         long = col_number(),
                                                         lat_new = col_number(),
                                                         long_new = col_number(),
                                                         taxa_count = col_number(),
                                                         .default = "c"))
#if redoing the snapping steps:
occ.list <- lapply(occ.list, function(x) {
  x %>%
    select(-c(contains("dist"), contains("flw"), COMID, contains("gage_no"), contains("_new")))
})

#make spatial w/ projected crs
sf.albers <- lapply(occ.list, st_as_sf, coords = c("long", "lat"), remove = FALSE, crs = 4269)
sf.albers <- lapply(sf.albers, st_transform, crs = 5070)

##snap to NHD stream ----
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% #nhd streams prev. clipped to region (good to use clipped because of MEM usage)
  st_transform(5070)

#assign COMID to each site + calc distance to nearest line (meters)
near_strms <- lapply(sf.albers, function(x) {
  
  nr.line <- st_nearest_feature(x, nhd, check_crs = TRUE) #find nearest nhd strm
  nhd_near <- nhd[nr.line,]
  
  dist <- as.vector(st_distance(x, nhd_near, by_element = TRUE)) #get distance to nearest nhd strm
  
  df <- data.frame(COMID = nhd_near$COMID, dist2strm_m_nhd = dist)
  
  return(df)
  
})

#add to occurrence df
occ.list[[1]] <- bind_cols(occ.list[[1]], near_strms[[1]])
occ.list[[2]] <- bind_cols(occ.list[[2]], near_strms[[2]])

rm(near_strms)

###get central coordinate for grouped COMIDs ----

nhd.center <- lapply(occ.list, function(x) {
  n.sub <- nhd %>%
    filter(COMID %in% x$COMID)
  
  n.center <- st_centroid(n.sub) %>%
    st_transform(4269) %>%
    select(COMID) %>%
    mutate(long_new = st_coordinates(.)[, "X"],
           lat_new = st_coordinates(.)[, "Y"]) %>%
    st_drop_geometry()
  
  return(n.center)
})

occ.list[[1]] <- left_join(occ.list[[1]], nhd.center[[1]])
occ.list[[2]] <- left_join(occ.list[[2]], nhd.center[[2]])

rm(nhd, nhd.center)

# #plot check
# ggplot() +
#   geom_sf(data = n.sub[2,]) +
#   geom_sf(data = n.center[2,], color = "red") +
#   theme_minimal()

##remake spatial w/ projected crs + new consolidated coordinates
sf.albers <- lapply(occ.list, st_as_sf, coords = c("long_new", "lat_new"), remove = FALSE, crs = 4269)
sf.albers <- lapply(sf.albers, st_transform, crs = 5070)

##snap to classified flow regime ----
flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% select(-PopupInfo) %>%
  st_transform(5070)
flw <- flw %>% filter(Flow_type != "BigR")

near_flw <- lapply(sf.albers, function(x) {
  
  nr.line <- st_nearest_feature(x, flw, check_crs = TRUE) #find nearest flw line
  flw_near <- flw[nr.line,]
  
  dist <- as.vector(st_distance(x, flw_near, by_element = TRUE)) #get distance to nearest flw line
  
  df <- data.frame(flw_name = flw_near$Name, 
                   flw_type = flw_near$Flow_type, 
                   dist2strm_m_flw = dist)
  
  return(df)
  
})

#add to occurrence df
occ.list[[1]] <- bind_cols(occ.list[[1]], near_flw[[1]])
occ.list[[2]] <- bind_cols(occ.list[[2]], near_flw[[2]])

rm(near_flw, flw)

##snap to gage ----
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))
g.sf <- st_as_sf(g.info, coords = c("dec_long_va", "dec_lat_va"), remove = FALSE, crs = 4269) %>%
  st_transform(5070)

snap_gage <- function(x, yr_to_clip, df_gage_obj) {
  
  res_list <- list()
  
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
for(i in seq_along(occ.list)) { 
  x <- occ.list[[i]]
  
  if(i == 1) { taxa <- "bug"} else { taxa <- "fish" }
  
  write_csv(x, paste0(PATH, "/01_BioDat/occ_", taxa, "_inthigh_long_", gsub("-", "", Sys.Date()), ".csv"))
}
rm(x, taxa, i, sf.albers, near_gage, g.sf)


#some name fixes I noticed later ----
bug <- read_csv(paste0(PATH, "/01_BioDat/occ_bug_inthigh_long_20240208.csv"), col_types = cols(lat = col_number(),
                                                                                               long = col_number(),
                                                                                               lat_new = col_number(),
                                                                                               long_new = col_number(),
                                                                                               taxa_count = col_number(),
                                                                                               .default = "c"))
bug <- bug %>%
  mutate(subfamily = case_when(genus == "Thienemannimyia" ~ "Tanypodinae",
                               genus == "Orthocladius" ~ "Orthocladiinae",
                               tribe == "Hydroporini" ~ "Hydroporinae",
                               tribe %in% c("Tanytarsini", "Chironomini") ~ "Chironominae",
                               tribe == "Pentaneurini" ~ "Tanypodinae",
                               tribe == "Tabanini" ~ "Tabaninae",
                               tribe == "Palpomyiini" ~ "Ceratopogoninae",
                               tribe == "Stratiomyini" ~ "Stratiomyinae",
                               T ~ subfamily))

write_csv(bug, paste0(PATH, "/01_BioDat/occ_bug_inthigh_long_20240208.csv"))

fish <- read_csv(paste0(PATH, "/01_BioDat/occ_fish_inthigh_long_20240208.csv"), col_types = cols(lat = col_number(),
                                                                                                 long = col_number(),
                                                                                                 lat_new = col_number(),
                                                                                                 long_new = col_number(),
                                                                                                 taxa_count = col_number(),
                                                                                                 .default = "c"))

fish <- fish %>%
  mutate(species = ifelse(species == "Lepomis punctatus", "Lepomis miniatus", species))

write_csv(fish, paste0(PATH, "/01_BioDat/occ_fish_inthigh_long_20240208.csv"))

rm(list = ls())

# filter + widen for analysis ----
PATH <- getwd()

file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "inthigh_long", full.names = TRUE)
occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                         long = col_number(),
                                                         lat_new = col_number(),
                                                         long_new = col_number(),
                                                         taxa_count = col_number(),
                                                         .default = "c"))

occ.adj <- lapply(occ.list, function(x) { 
  x %>% 
    mutate(across(contains("dist"), as.numeric)) %>%
    filter(dist2gage_m_15yr <= 15000 & dist2strm_m_flw <= 5000) %>% #make sure within 15km of a gage + 5 km classified flow line
    filter(!gage_no_15yr %in% c("06906800", "07020550", "07332500")) %>% #these gages are way off (>15k) from any predicted flow regimes
    select(matches(c("subfamily", "tribe", "genus", "species")), lat_new, long_new, COMID, flw_type, gage_no_15yr, dist2gage_m_15yr, dist2strm_m_flw) %>%
    rename(lat = lat_new, long = long_new) %>%
    distinct() %>%
    group_by(lat, long, COMID) %>%
    mutate(site_id = paste0("s_", cur_group_id())) %>%
    ungroup()
}) 


fish <- occ.adj[[2]] %>%
  select(-genus) %>%
  distinct() %>%
  arrange(species) %>%
  mutate(species = sub(" ", "_", species),
         p = 1) %>%
  pivot_wider(names_from = species, values_from = p) %>%
  mutate(across(-c(site_id, COMID, flw_type, gage_no_15yr), ~ ifelse(is.na(.x), 0, .x)))

write_csv(fish, paste0(PATH, "/01_BioDat/occ_fish_15k_wide_", gsub("-", "", Sys.Date()), ".csv"))

bug_names <- occ.adj[[1]] %>%
  select(subfamily, tribe, genus) %>%
  filter(!(is.na(subfamily) & is.na(tribe) & is.na(genus))) %>%
  group_by(subfamily, tribe, genus) %>%
  summarise(ttl_records = n(), .groups = "drop") #%>%
  # distinct()

bug <- occ.adj[[1]] %>%
  #not really sure if I want to do it this way (to remove ambiguous names or group together)
  # mutate(taxa = case_when(tribe == "Aciliini" ~ tribe,
  #                         tribe == "Chironomini" ~ genus,
  #                         tribe == "Eriopterini" ~ genus,
  #                         tribe == "Hydrophilini" ~ tribe,
  #                         tribe == "Hydroporini" ~ genus,
  #                         tribe == "Limoniini" ~ genus,
  #                         tribe == "Nymphulini" ~ genus,
  #                         tribe == "Palpomyiini" ~ genus,
  #                         subfamily == "Tanypodinae" ~ genus,
  #                         tribe == "Pseudochironomini" ~ tribe,
  #                         tribe == "Stratiomyini" ~ tribe,
  #                         tribe == "Tabanini" ~ genus,
  #                         tribe == "Tanytarsini" ~ genus,
  #                         is.an(subfamily) & is.na(genus) ~ genus,
  #                         ))
  mutate(taxa = case_when(tribe == "Aciliini" ~ tribe,
                          tribe == "Hydrophilini" ~ tribe,
                          tribe == "Pseudochironomini" ~ tribe,
                          !is.na(genus) ~ genus,
                          !is.na(tribe) ~ tribe,
                          !is.na(subfamily) ~ subfamily)) %>%
  filter(!is.na(taxa)) %>%
  select(-c(subfamily, tribe, genus)) %>%
  distinct() %>%
  arrange(taxa) %>%
  mutate(p = 1) %>%
  pivot_wider(names_from = taxa, values_from = p) %>%
  mutate(across(-c(site_id, COMID, flw_type, gage_no_15yr), ~ifelse(is.na(.x), 0, .x)))

write_csv(bug, paste0(PATH, "/01_BioDat/occ_bug_15k_wide_", gsub("-", "", Sys.Date()), ".csv"))

#summarize datasets ----
# occ <- bind_rows(occ.list)
# occ <- occ.list[[2]] %>%
#   filter(dist2gage_m_15yr <= 15000) %>%
#   group_by(order, family, genus, species) %>%
#   summarise(record_ttl_15kbuff = n())
# write_csv(occ, paste0(PATH, "/01_BioDat/fish_species_list_15kcrop.csv"))
# #
# occ <- occ.list[[1]] %>%
#   filter(dist2gage_m_15yr <= 15000) %>%
#   group_by(order, family, subfamily, tribe, genus) %>%
#   summarise(record_ttl_15kbuff = n())
# write_csv(occ, paste0(PATH, "/01_BioDat/bug_species_list_15kcrop.csv"))
