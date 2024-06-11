## MISC. PLOTS ##

library(tidyverse); library(sf)

PATH <- getwd()

source(paste0(PATH, "/Scripts/XX_colors.R"))

# study area ----
#eco region of interest
eco <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/ecoreg_l3_interior_highlands_crop.shp")) #%>%
#   st_union()

#armook 
hlnd <- st_as_sf(maps::map("state", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)

#USGS gages
g.info <- read_csv(paste0(PATH, "/02_EnvDat/all_hit_usgs_gage_info.csv")) %>%
  select(site_no, dec_long_va, dec_lat_va) %>%
  distinct() %>%
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269)
g.info <- g.info %>%
  st_filter(., eco)

g.buff <- g.info %>% st_buffer(., dist = 15000) %>% st_union()

#get flow types
#assign flow_type to each site + calc distance to nearest line (meters)
flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% select(-PopupInfo) %>%
  st_transform(4269)
flw <- flw %>% filter(Flow_type != "BigR")
nr.line <- st_nearest_feature(g.info, flw, check_crs = TRUE) #find nearest strm not river
flw_near <- flw[nr.line,]
g.info$flw_name <- flw_near$Name #add strm names
g.info$flw_type <- flw_near$Flow_type #add flow types
# dist <- as.vector(st_distance(alb.sites, flw_near, by_element = TRUE)) #get distance from site to strm assigned
# temp.sites$dist2strm_m_flw <- dist

rm(flw, flw_near, nr.line)

#bio dat
occ.dat <- lapply(list.files(paste0(PATH, "/01_BioDat"), pattern = "^occ_(bug|fish)", full.names = T), function(x) {
  df <- read_csv(x) %>%
    select(lat, long, phylum) %>%
    distinct()
})
occ.dat <- bind_rows(occ.dat) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4269) %>%
  mutate(taxa = ifelse(phylum == "Arthropoda", "aquatic insect", "fish")) %>%
  st_filter(., eco)

occ.buff <- occ.dat %>%
  st_filter(., st_make_valid(g.buff))

#plot
ggplot() +
  geom_sf(data = hlnd) +
  geom_sf(data = eco, fill = NA) +
  geom_sf(data = occ.dat, aes(fill = taxa), shape = 21, size = 2, alpha = 0.5) +
  # geom_sf(data = occ.buff, aes(fill = taxa), shape = 21, size = 2, alpha = 0.5) +
  geom_sf(data = g.info, aes(shape = flw_type, color = flw_type), size = 3) +
  # geom_sf(data = g.buff, fill = NA) +
  scale_fill_manual(values = tax.pal, name = "Survey Site") +
  scale_color_manual(values = flow.pal, name = "Flow Regime") +
  scale_shape_manual(values = c(16, 17, 18), name = "Flow Regime") +
  coord_sf(xlim = c(-98, -89)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        legend.text = element_text(size = 12))

ggsave(paste0(PATH, "/99_figures/study_map_allocc.png"), width = 7, height = 6, bg = "white")

# species threshold comp plots ----
thresh_df <- lapply(list.files(paste0(PATH, "/11_Thresholds"), pattern = "full_thresh", full.names = TRUE), read_csv)

#get type of model and add to dataframe
for(i in 1:length(thresh_df)) {
  
  thresh_df[[i]]$name <- sub(".csv", "", 
                             sub("full_thresh_gf.", "", 
                                 list.files(paste0(PATH, "/11_Thresholds"), pattern = "full_thresh")[[i]]))
  
  thresh_df[[i]] <- thresh_df[[i]] %>%
    mutate(taxa = case_when(grepl("bug", name) ~ "bug",
                            grepl("fish", name) ~ "fish"),
           flow_type = case_when(grepl("GW", name) ~ "GW",
                                 grepl("Int", name) ~ "Int",
                                 grepl("RO", name) ~ "RO"),
           var_type = case_when(grepl("HIT", name) ~ "HIT",
                                grepl("LULC", name) ~ "LULC",
                                TRUE ~ "full"))
  
}

thresh_df <- bind_rows(thresh_df)

#boxplot comparing flow type thresholds (fish only)
fish <- thresh_df %>% 
  filter(taxa == "fish")

#get thresholds:
tmp <- fish %>%
  filter(var_type == "HIT") %>%
  select(species, env_var, flow_type, t_max_mn, t_m2_mn, t_max_prop, t_m2_prop) %>%
  distinct() %>%
  mutate(prop = case_when(t_m2_prop > t_max_prop & !is.na(t_m2_prop) ~ t_m2_prop, #find larger prop threshold
                          TRUE ~ t_max_prop),
         thresh = case_when(prop == t_max_prop ~ t_max_mn,
                            TRUE ~ t_m2_mn)) %>% 
  filter(prop >= 25) %>% #only a threshold if the change is greater than 25% of the overall cum. imp. curve
  select(-contains("t_")) %>%
  #remove variables only with one flow type showing up (for visualization purposes)
  group_by(env_var) %>%
  filter(length(unique(intersect(flow_type, unique(fish$flow_type)))) >= 2) %>%
  ungroup()

ggplot() +
  geom_boxplot(data = tmp, aes(x = flow_type, y = thresh, fill = flow_type)) +
  geom_jitter(data = tmp, aes(x = flow_type, y = thresh), alpha = 0.4, width = 0.3) +
  facet_wrap(~ env_var, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = flow.pal) +
  theme_bw() + xlab("") + ylab("species threshold value")
ggsave(paste0(PATH, "/99_figures/spp_thresh_boxplot_flowcomp.png"), bg = "white", width = 12, height = 5)

rm(list = ls())

#site diversity ----
library(tidyverse); library(sf)
PATH <- getwd()
source(paste0(PATH, "/Scripts/XX_colors.R"))

file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                         long = col_number(),
                                                         COMID = col_character(),
                                                         gage_no_15yr = col_character(),
                                                         dist2gage_m_15yr = col_number(),
                                                         dist2strm_m_flw = col_number()))
# occ.list <- list(occ.list[[2]])

div_list <- read_csv(paste0(PATH, "/98_result_tables/site_div_alltax_alldiv_raw.csv"))

#add COMID and flw type to diversity result
xy_flw <- vector("list", 2)
for(i in seq_along(occ.list)) {
  if(i == 1) { taxa <- "bug" } else { taxa <- "fish" }
  # taxa <- "fish"
  xy_flw[[i]] <- occ.list[[i]] %>% select(site_id, lat, long, COMID, flw_type) %>% mutate(taxa = taxa)
}

xy_flw <- bind_rows(xy_flw)

div <- left_join(div_list, xy_flw)

##boxplot of flow type x diversity estimate ----
tmp <- div %>% 
  pivot_longer(cols = matches("^(n_|f_)"), names_to = "div_type", values_to = "value") %>%
  mutate(div_type = case_when(div_type == "n_sp" ~ "richness",
                              div_type == "n_fsp" ~ "n. functionally unique spp",
                              div_type == "f_eve" ~ "functional evenness",
                              div_type == "f_disp" ~ "functional dispersion",
                              T ~ NA)) %>%
  filter(!is.na(div_type))

ggplot() +
  geom_boxplot(data = tmp %>% filter(taxa == "fish"), aes(x = flw_type, y = value, fill = flw_type), na.rm = T) +
  # geom_jitter(data = tmp, aes(x = flw_type, y = value), alpha = 0.2, width = 0.3) +
  # facet_grid(div_type ~ taxa, scales = "free") +
  facet_wrap(~ div_type, scales = "free_y") +
  scale_fill_manual(values = flow.pal) +
  theme_bw() +
  theme(axis.title = element_blank())
ggsave(paste0(PATH, "/99_figures/site_div_fish_boxplot_comp.png"), width = 7, height = 5)

##density plot flow type x diversity ----
ggplot() +
  geom_density(data = div, aes(f_disp, fill = flw_type), alpha = 0.6) +
  facet_wrap(~ taxa) +
  scale_fill_manual(values = flow.pal, name = "flow type") +
  xlab("functional dispersion") +
  theme_bw()
ggsave(paste0(PATH, "/99_figures/site_div_alltax_density_fdispcomp.png"), width = 8, height = 3)

## funct ~ rich scatterplot ----

ggplot(data = tmp) +
  geom_density(aes(value, fill = div_type)) +
  facet_wrap(~ div_type, scales = "free") 

ggplot(data = filter(div, n_sp > 1), aes(x = n_sp, y = f_disp, color = flw_type)) +
  geom_point(aes(shape = flw_type), size = 3, alpha = 0.1) +
  geom_smooth(method = "glm", linewidth = 1.5) +
  facet_wrap(~ taxa, scales = "free_x") +
  scale_color_manual(values = flow.pal) +
  coord_cartesian(ylim = c(0, 0.26)) +
  xlab("number of species per site") + ylab("functional dispersion") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(fill = NA, size = 3, alpha = 1)))

ggsave(paste0(PATH, "/99_figures/site_div_alltax_rich-v-disp_comp.png"), width = 8, height = 5)

##map of site diversity ----
#nhd load in
# nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>%
#   mutate(COMID = as.character(COMID))
# nhd <- left_join(nhd, div)
div.sf <- st_as_sf(div, coords = c("long", "lat"), crs = 4269) 

#armook 
hlnd <- st_as_sf(maps::map("state", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)

p1 <- ggplot() +
  geom_sf(data = hlnd, fill = "lightgrey", color = "black") +
  geom_sf(data = div.sf %>% filter(taxa == "fish") %>% arrange(n_sp), 
          aes(fill = n_sp), shape = 21, size = 2, alpha = 0.7) +
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = c(-98, -89)) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank())
p2 <- ggplot() +
  geom_sf(data = hlnd, fill = "lightgrey", color = "black") +
  geom_sf(data = div.sf %>% filter(taxa == "fish") %>% arrange(f_disp), 
          aes(fill = f_disp), shape = 21, size = 2, alpha = 0.7) +
  scale_fill_viridis_c(direction = -1) +
  coord_sf(xlim = c(-98, -89)) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank())
ggpubr::ggarrange(p1, p2, nrow = 1)

ggsave(paste0(PATH, "/99_figures/site_div_fish_map_comparison.png"), width = 10, height = 5)

#site diversity relationship ----
div.sf <- st_as_sf(div %>% filter(taxa == "fish"), coords = c("long", "lat"), crs = 4269) 

#break down relationship into cats
div.sf <- div.sf %>%
  mutate(across(.cols = c(n_sp, f_disp),
                ~ case_when(.x <= quantile(.x, c(0.25)) ~ "low",
                            .x > quantile(.x, c(0.25)) & .x <= quantile(.x, c(0.50)) ~ "low-med",
                            .x > quantile(.x, c(0.50)) & .x <= quantile(.x, c(0.75)) ~ "med-high",
                            .x > quantile(.x, c(0.75)) ~ "high")),
         across(.cols = c(n_sp, f_disp), ~ factor(.x, levels = c("low", "low-med", "med-high", "high"), ordered = TRUE)),
         rich_fdisp = interaction(n_sp, f_disp, sep = ":"))

#make legend
l <- ggplot() +
  geom_tile(data = div.sf, aes(n_sp, f_disp, fill = rich_fdisp)) +
  scale_fill_viridis_d(option = "viridis") +
  scale_x_discrete(expand = c(0,0), name = "richness") +
  scale_y_discrete(expand = c(0,0), name = "functional diversity") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
l

#plot
p1 <- ggplot() +
  geom_sf(data = hlnd, fill = "grey", color = "black", linewidth = 0.6) +
  geom_sf(data = div.sf, aes(color = rich_fdisp), shape = 19, size = 5, alpha = 0.8) +
  scale_color_viridis_d() +
  coord_sf(xlim = c(-98, -89)) +
  theme(axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = "none")
p1

# ggsave(paste0(PATH, "/99_figures/site_div_fish_interaction_legend.png"), plot = l, width = 5, height = 4)
ggsave(paste0(PATH, "/99_figures/site_div_fish_interaction_map.png"), plot = p1, width = 10, height = 9)
