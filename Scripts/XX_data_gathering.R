##species list, region shapefile, nhd retrieval, sites and gage points

library(tidyverse); library(sf); library(maps)

PATH <- getwd()

gf.list <- list.files(paste0(PATH, "/10_GFOutput"), ".rds", full.names = T)

# species list ----
spp.list <- lapply(gf.list, function(gf.name) {
  
  gf.mod <- readRDS(gf.name)
  spp.df <- data.frame(species = colnames(gf.mod$Y), type = NA)
  
  spp.df$type <- ifelse(grepl("fish", gf.name), "fish", "bug")
  
  return(spp.df)
  
})

spp.list <- bind_rows(spp.list) %>%
  distinct() %>%
  arrange(type, species)

write_csv(spp.list, paste0(PATH, "/species_list.csv"))

#ecoregion shapefile ----
#ecoregions = arkansas valley, boston mountains, ouachita mountains, ozark highlands, south central plains
eco.full <- st_read(paste0(PATH, "/02_SpatialData/us_eco_l3"), "us_eco_l3")
eco.ih <- nhd.full %>%
  filter(US_L3NAME %in% c("Arkansas Valley", "Boston Mountains", "Ouachita Mountains", "Ozark Highlands", "South Central Plains"))

st_write(eco.ih, paste0(PATH, "/02_SpatialData/ecoreg_l3_interior_highlands.shp"))
rm(eco.full)

st1 <- st_as_sf(map(database = "state", region = c('missouri', 'arkansas', 'oklahoma'), plot = FALSE, fill = TRUE)) %>%
  st_transform(., st_crs(4269)) %>%
  st_union()
eco.ih <- eco.ih %>%
  st_transform(., st_crs(4269)) %>%
  st_intersection(., st1)

st_write(eco.ih, paste0(PATH, "/02_SpatialData/ecoreg_l3_interior_highlands_crop.shp"))

#HUC10 watersheds ----
eco.ih <- read_sf(paste0(PATH, "/02_SpatialData/ecoreg_l3_interior_highlands_crop.shp"))

#go in each folder, pull out shapefile, combine
wbd.dir <- dir(path = paste0(PATH, "/02_SpatialData"), pattern = "WBD_", full.names = T)

huc.full <- lapply(wbd.dir, function(dir.name) {
  h10 <- read_sf(paste0(dir.name, "/Shape"), query = "select tnmid,areasqkm,states,huc10,hutype from WBDHU10")
  return(h10)
})

huc.full <- do.call(rbind, huc.full)

huc.ih <- st_filter(huc.full, eco.ih)

write_sf(huc.ih, paste0(PATH, "/02_SpatialData/huc10_interior_highlands_crop.shp"))
rm(huc.full)

#nhd in highlands ----

# DID IN QGIS

# sites as points ----
spp.l <- read_csv(paste0(PATH, "/species_list.csv"))
fish.dat <- read_csv(paste0(PATH, "/01_Data/ARMOOK_Fishes_bysite23a_StreamCat_Flow_Full_15km.csv")) %>% select(-c('...1', etheostoma.zonale, notropis.boops)) %>%
  st_as_sf(., coords = c("site_x", "site_y"), crs = 26915) #UTM zone 15
bug.dat <- read_csv(paste0(PATH, "/01_Data/BenthicInsect_23_StreamCat_Flow_Full_15km.csv")) %>%
  st_as_sf(., coords = c("site_x", "site_y"), crs = 4269)

dat.list <- list(fish.dat, bug.dat)
  
dat.list <- lapply(dat.list, function(x) {
  x <- x %>%
    st_drop_geometry() %>%
    select(STAID, contains(spp.l$species)) %>%
    mutate(STAID = paste0("0", as.character(STAID))) %>%
    group_by(STAID) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~if_else(.x > 0, 1, 0))) %>%
    mutate(rich = rowSums(select(., where(is.numeric)) %>% as.matrix())) %>%
    select(STAID, rich)
  
  return(x)
})
                  
g.rich <- full_join(dat.list[[1]], dat.list[[2]], by = "STAID") %>%
  rename(rich_fish = rich.x, rich_bug = rich.y)

g.info <- dataRetrieval::whatNWISdata(siteNumbers = unique(g.rich$STAID)) #get gage locations
g.info <- g.info %>%
  select(site_no, dec_long_va, dec_lat_va) %>%
  distinct()

g.rich <- left_join(g.rich, g.info, by = c("STAID"="site_no"))
write_csv(g.rich, paste0(PATH, "/01_Data/USGS_gage_locations.csv"))

#study extent map ----
states <- st_as_sf(map(database = "state", region = c('missouri', 'arkansas', 'oklahoma'), plot = FALSE, fill = TRUE))
eco.ih <- st_read(paste0(PATH, "/02_SpatialData/ecoreg_l3_interior_highlands_crop.shp"))
huc.ih <- st_read(paste0(PATH, "/02_SpatialData/huc10_interior_highlands_crop.shp"))
g.rich <- read_csv(paste0(PATH, "/01_Data/USGS_gage_locations.csv")) %>% st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269)

ggplot() +
  geom_sf(data = states, fill = "lightgray", color = "darkgray", lwd = 1.2) +
  # geom_sf(data = huc.full, fill = NA) +
  geom_sf(data = huc.ih, fill = NA) +
  geom_sf(data = eco.ih, aes(linetype = US_L3NAME), fill = NA, lwd = 1.5) +
  # geom_sf(data = fish.dat, fill = "blue", alpha = 0.6, shape = 21, size = 5) +
  # geom_sf(data = bug.dat, fill = "red", alpha = 0.6, shape = 21, size = 5) +
  geom_sf(data = g.rich, aes(fill = rich_fish), shape = 21, size = 5) +
  scale_fill_viridis_c() +
  theme_minimal()
