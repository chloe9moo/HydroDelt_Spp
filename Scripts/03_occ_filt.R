# OCCURRENCE FILTERING + TAXONOMY FINALIZING

library(tidyverse); library(sf)

PATH <- getwd()

occ.dat.orig <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_updatedspp_nofilt_20231128.csv"),  col_types = cols(date_min = col_character(), date_max = col_character(),
                                                                                                            subclass = col_character(), superorder = col_character(), 
                                                                                                            suborder = col_character(), infraorder = col_character(), 
                                                                                                            superfamily = col_character(), subgenus = col_character(), 
                                                                                                            species = col_character(), subspecies = col_character())) %>%
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
                            TRUE ~ family),
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
                               TRUE ~ subfamily),
         suborder = case_when(taxa_name_orig == "Nematocera" ~ "Nematocera",
                              TRUE ~ suborder)) 

occ.dat <- occ.dat.orig %>%
  select(bio_site_no, lat, long, contains("date"), taxa_count,
         kingdom, phylum, 
         class, subclass,
         superorder, order, suborder, infraorder,
         superfamily, family, subfamily, 
         tribe, genus, subgenus, species, subspecies,
         -contains("orig")) %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) #get rid of zero count ids

# table(occ.dat$phylum)

#get fish dat ----
fish <- occ.dat %>%
  filter(phylum == "Chordata") %>% #verts
  filter(!class %in% c("Aves", "Amphibia") & order != "Testudines") %>% #remove non-fish verts
  select(-contains("sub"), -contains("super"), -contains("infra"), -c(kingdom, phylum, class, tribe)) %>%
  filter(!is.na(species))
  # rowid_to_column("key")

#double check names with original names!
#double check sources with original runs
tmp <- occ.dat.orig %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) %>% #get rid of zero count ids
  filter(phylum == "Chordata") %>% #verts
  filter(!class %in% c("Aves", "Amphibia") & order != "Testudines") %>%
  select(matches("family|genus|species|_orig"), -matches("super|infra"), -subgenus, -subfamily, -subspecies, -order_orig) %>%
  distinct()
#   select(source)
# table(tmp$source)

#get species list
tmp <- fish %>% select(order, family, genus, species) %>%
  distinct()

write_csv(tmp, paste0(PATH, "/01_BioDat/fish_species_list.csv"))
write_csv(fish, paste0(PATH, "/01_BioDat/occ_fish_all_long_20231128.csv"))

#get invert dat ----
t <- read_csv(paste0(PATH, "/01_BioDat/tribe_info.csv")) %>% rename(tribe_new = tribe) #for updating tribes

bug <- occ.dat %>%
  filter(class == "Insecta") %>%
  select(-c(kingdom, phylum, class, subclass, superorder, suborder, infraorder, superfamily, subgenus, species, subspecies)) %>%
  #update tribes in case that is the lowest tax level chosen
  left_join(., t) %>%
  mutate(tribe_new = case_when(is.na(tribe_new) ~ tribe,
                               T ~ tribe_new)) %>%
  select(-tribe) %>% rename(tribe = tribe_new)

#double check names with original names!
#double check sources with original runs
tmp <- occ.dat.orig %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) %>% #get rid of zero count ids
  filter(class == "Insecta") %>%
  select(matches("order|family|subfamily|tribe|genus|_orig|_check"), -matches("super|infra"), -suborder, -subgenus, -species_orig) %>%
  distinct()
  # select(source)
# table(tmp$source)

#get species list
tmp <- bug %>% select(order, family, subfamily, tribe, genus) %>%
  distinct()

write_csv(tmp, paste0(PATH, "/01_BioDat/bug_species_list.csv"))
write_csv(bug, paste0(PATH, "/01_BioDat/occ_bug_all_long_20231128.csv"))

#region filter ----
eco <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "ecoreg_l3_interior_highlands_crop")

sf.occ.dat <- lapply(list(bug, fish), function(x) { #convert bio dat to spatial dat
  x <- st_as_sf(x, coords = c("long", "lat"), remove = FALSE, crs = 4269)
})

eco.buff <- st_buffer(eco, dist = 10000) #10 km buffer for wiggle room

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
  
  write_csv(x, paste0(PATH, "/01_BioDat/occ_", taxa, "_inthigh_long_", gsub("-", "", Sys.Date())))
}
rm(x, taxa, i)

