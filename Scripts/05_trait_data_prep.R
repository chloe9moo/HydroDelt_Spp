### TRAIT DATA PREP ###

library(tidyverse)
options(readr.show_col_types = FALSE)

PATH <- getwd()

t.cols <- read_csv(paste0(PATH, "/20_Traits/trait_column_selection.csv")) %>%
  filter(!is.na(database_col_name) & database_col_name != "--")

source(paste0(PATH, "/Scripts/XX_trait_functions.R"))

#FISH ----
f.spp <- read_csv(paste0(PATH, "/01_BioDat/fish_species_list.csv"))

f.cols <- t.cols %>%
  separate_rows(., database_col_name, sep = ", ") %>%
  filter(taxa == "fish") %>%
  select(trait_name, trait_type, database_col_name) 

##TOFF ----
f1 <- read_delim(paste0(PATH, "/20_Traits/TOFF_fix_1.csv"), delim = ";", escape_double = FALSE, trim_ws = TRUE) %>%
  select(where(~ any(!is.na(.))), -`,,,,,,,,,,,,,,,,,,,,,,,,`) %>%
  mutate(mea_id = sub("^[^ ]* ", "", mea_id)) %>%
  select(tra_name, gen_name, spe_name, mea_value) %>%
  rename(genus = gen_name, spp_ep = spe_name) %>%
  mutate(species = paste0(genus, " ", spp_ep)) %>%
  filter(tra_name %in% f.cols$database_col_name) %>% #get desired traits
  left_join(., f.cols, by = c("tra_name" = "database_col_name")) %>% #rename traits for consistency 
  select(-tra_name) %>%
  distinct()

#pivot wider (using apply function bc of some format weirdness)
tmp <- apply(f1, 1, function(x) {
  df <- data.frame(genus = x[["genus"]], spp_ep = x[["spp_ep"]], species = x[["species"]], trait = x[["mea_value"]])
  if(x[["trait_type"]] == "continuous") { #convert continous traits to numeric class
    df$trait <- as.numeric(df$trait)
  }
  names(df)[names(df) == "trait"] <- x[["trait_name"]]
  return(df)
})
f1 <- bind_rows(tmp)

f1 <- f1 %>%
  mutate(across(c(where(is.character), -c(genus, spp_ep, species)), 
                ~ case_when(grepl("yes", .x, ignore.case = TRUE) ~ "1",
                            grepl("no", .x, ignore.case = TRUE) ~ "0",
                            T ~ .x)),
         max_tl = max_tl / 10, #convert mm to cm for consistency
         db_orig = "TOFF")

##FISHMORPH ---- for the time being, not using FISHMORPH data

# f2 <- read_delim(paste0(PATH, "/20_Traits/FISHMORPH_Database.csv"), delim = ";", escape_double = FALSE, trim_ws = TRUE) %>%
#   mutate(spp_ep = sub("^[^ ]* ", "", `Genus species`),
#          genus = sub(" .*", "", `Genus species`),
#          db_orig = "FISHMORPH") %>%
#   rename(species = `Genus species`) %>%
#   select(genus, spp_ep, species, f.cols$database_col_name) %>%
#   relocate(genus, spp_ep, species) 
# 
# #adjust for synonyms
# tmp <- f2 %>%
#   filter(species %in% c("Erimyzon oblongus", "Etheostoma euzonum", "Etheostoma spectabile", "Noturus albater")) %>%
#   mutate(species = case_when(species == "Erimyzon oblongus" ~ "Erimyzon claviformis",
#                              species == "Etheostoma euzonum" ~ "Etheostoma erizonum",
#                              species == "Etheostoma spectabile" ~ "Etheostoma fragi",
#                              species == "Noturus albater" ~ "Noturus maydeni"),
#          spp_ep = sub("^[^ ]* ", "", species),
#          genus = sub(" .*", "", species))
# f2 <- bind_rows(f2, tmp)
# 
# f2 <- f2 %>% filter(species %in% f.spp$species)

##FishTraits ----
f3 <- readxl::read_excel(paste0(PATH, "/20_Traits/FishTraits_14.3.xls")) %>%
  rename(genus = GENUS, spp_ep = SPECIES) %>%
  mutate(spp_ep = ifelse(spp_ep == "petenese", "petenense", spp_ep), #fix misspelling
         species = paste0(genus, " ", spp_ep)) %>%
  select(genus, spp_ep, species, 
         matches(paste(f.cols$database_col_name[f.cols$database_col_name != "--" & !is.na(f.cols$database_col_name)], collapse = "|")),
         -OTHERNAMES) %>%
  #handle missing value placeholders
  mutate(across(-c(genus, spp_ep, species), ~if_else(.x %in% c(-1, -555, -999), NA, .x)),
         db_orig = "FishTraits")

#update names + add summary columns
c.key <- setNames(f.cols$trait_name, f.cols$database_col_name)
c.key <- c.key[names(c.key) %in% intersect(names(f3), names(c.key))] #get rid of TOFF variables
c.key <- c.key[!c(duplicated(c.key) | duplicated(c.key, fromLast = TRUE))] #only retain variables that don't need to be summarized/consolidated

f3 <- f3 %>%
  rename_with(~c.key[.], 
              .cols = intersect(names(f3), names(c.key))) %>%
  #spawning season
  mutate(spawn_spring = ifelse(rowSums(!is.na(select(., MAR, APR, MAY))) > 0, rowSums(select(., MAR, APR, MAY), na.rm = TRUE), NA),
         spawn_summer = ifelse(rowSums(!is.na(select(., JUN, JUL, AUG))) > 0, rowSums(select(., JUN, JUL, AUG), na.rm = TRUE), NA),
         spawn_fall = ifelse(rowSums(!is.na(select(., SEP, OCT, NOV))) > 0, rowSums(select(., SEP, OCT, NOV), na.rm = TRUE), NA),
         spawn_winter = ifelse(rowSums(!is.na(select(., DEC, JAN, FEB))) > 0, rowSums(select(., DEC, JAN, FEB), na.rm = TRUE), NA),
         across(.cols = c(spawn_spring, spawn_summer, spawn_fall, spawn_winter), ~ case_when(.x > 0 ~ 1,
                                                                                             is.na(.x) ~ NA,
                                                                                             TRUE ~ 0))) %>%
  select(-c(JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC)) %>%
  #spawn behaviors
  mutate(spawn_g_ns = ifelse(rowSums(!is.na(select(., starts_with("B_2_")))) > 0, rowSums(select(., starts_with("B_2_")), na.rm = TRUE), NA),
         spawn_g_sc = ifelse(rowSums(!is.na(select(., starts_with("B_1_")))) > 0, rowSums(select(., starts_with("B_1_")), na.rm = TRUE), NA),
         spawn_ng_bh = ifelse(rowSums(!is.na(select(., starts_with("A_2_")))) > 0, rowSums(select(., starts_with("A_2_")), na.rm = TRUE), NA),
         spawn_ng_os = ifelse(rowSums(!is.na(select(., starts_with("A_1_")))) > 0, rowSums(select(., starts_with("A_1_")), na.rm = TRUE), NA)) %>%
  select(-c(starts_with("B_"), starts_with("A_"))) %>%
  #carnivore comp
  mutate(troph_carniv = ifelse(rowSums(!is.na(select(., c(FSHCRCRB, BLOOD, EGGS)))) > 0, rowSums(select(., c(FSHCRCRB, BLOOD, EGGS)), na.rm = TRUE), NA),
         troph_carniv = case_when(troph_carniv > 0 ~ 1,
                                  is.na(troph_carniv) ~ NA,
                                  TRUE ~ 0)) %>%
  select(-c(FSHCRCRB, BLOOD, EGGS))

#adjust synonyms
tmp <- f3 %>%
  filter(species == "Etheostoma spectabile") %>%
  mutate(species = "Etheostoma fragi",
         spp_ep = sub("^[^ ]* ", "", species),
         genus = sub(" .*", "", species))
tmp1 <- f3 %>%
  filter(species %in% c("Phoxinus erythrogaster", "Dorosoma petenese", "Erimyzon oblongus",
                        "Etheostoma whipplei", "Etheostoma euzonum", "Etheostoma spectabile",
                        "Lampetra appendix", "Moxostoma macrolepidotum", "Notropis volucellus",
                        "Percina caprodes")) %>%
  mutate(species = case_when(species == "Phoxinus erythrogaster" ~ "Chrosomus erythrogaster",
                             species == "Dorosoma petenese" ~ "Dorosoma petenense",
                             species == "Erimyzon oblongus" ~ "Erimyzon claviformis",
                             species == "Etheostoma whipplei" ~ "Etheostoma artesiae", 
                             species == "Etheostoma euzonum" ~ "Etheostoma erizonum",
                             species == "Etheostoma spectabile" ~ "Etheostoma uniporum",
                             species == "Lampetra appendix" ~ "Lethenteron appendix", 
                             species == "Moxostoma macrolepidotum" ~ "Moxostoma pisolabrum", 
                             species == "Notropis volucellus" ~ "Notropis wickliffi",
                             species == "Percina caprodes" ~ "Percina fulvitaenia"),
         spp_ep = sub("^[^ ]* ", "", species),
         genus = sub(" .*", "", species))
f3 <- bind_rows(f3, tmp, tmp1)

#convert cat vars to character for joining
f3 <- f3 %>%
  mutate(across(.cols = intersect(names(f3), t.cols[t.cols$trait_type == "binary" & t.cols$database_source == "FishTraits", ]$trait_name),
                as.character))

##olden et al. data ----
f2 <- read_csv(paste0(PATH, "/20_Traits/olden_fish_traits.csv")) %>%
  select(genus, species, 
         matches(paste(f.cols$database_col_name[f.cols$database_col_name != "--" & !is.na(f.cols$database_col_name)], collapse = "|"))) %>%
  #fix unit missmatch 
  mutate(MAXBODYL = MAXBODYL / 10) #mm to cm

f.cols2 <- t.cols %>%
  filter(database_source == "OldenTraits") %>%
  select(trait_name, trait_type, database_col_name)  %>%
  mutate(trait_cat = case_when(grepl(",", database_col_name) ~ sub(".+, (.+)", "\\1", database_col_name), #get individual categories from col
                               TRUE ~ sub(".+([0-9])", "\\1", database_col_name)),
         database_col_name = case_when(grepl(",", database_col_name) ~ sub(",[^,]+$", "", database_col_name), #fix original col name column
                                       TRUE ~ sub("[0-9]", "", database_col_name))) %>%
  separate_rows(., database_col_name, sep = ", ")

#update categorical/multi column trait
tmp <- f2 %>%
  select(genus, species, RGUILD) %>%
  pivot_longer(cols = -c(genus, species), names_to = "Trait_group", values_to = "Trait") %>%
  left_join(., f.cols2, by = c("Trait_group" = "database_col_name"), multiple = "all", relationship = "many-to-many") %>% #rename traits for consistency 
  group_by(Trait_group) %>%
  filter(str_detect(Trait, regex(trait_cat, ignore_case = TRUE))) %>% #remove extra rows created by join (by selecting only rows where cat == value)
  mutate(value = 1) %>% #prep for widening dataframe
  group_by(genus, species, trait_name) %>%
  summarise(sum = sum(value)) %>%
  pivot_wider(names_from = trait_name, values_from = sum) %>%
  ungroup() %>%
  mutate(across(-c(genus, species), ~ ifelse(is.na(.x), 0, .x)),
         across(-c(genus, species), as.character))

#update all other names + add guild back in
c.key <- setNames(f.cols2$trait_name, f.cols2$database_col_name)
c.key <- c.key[names(c.key) != "RGUILD"]

f2 <- f2 %>%
  select(-RGUILD) %>%
  rename_with(~c.key[.], 
              .cols = intersect(names(f2), names(c.key))) %>%
  left_join(., tmp) %>%
  mutate(db_orig = "olden")

##combine databases ----
fish <- bind_rows(f1, f2, f3) %>% select(-spp_ep)

write_csv(fish, paste0(PATH, "/20_Traits/trait_dat_raw_fish.csv"))

##summarize 1 row per species ----
fish <- read_csv(paste0(PATH, "/20_Traits/trait_dat_raw_fish.csv"))

fish.sum <- fish %>%
  filter(species %in% f.spp$species) %>% #generous estimate of spp list for now
  mutate(across(.cols = unique(f.cols[f.cols$trait_type == "binary", ]$trait_name),
                as.numeric)) %>%
  group_by(genus, species) %>%
  summarise(across(.cols =  c("fecundity", "longevity", "matur_age", "spawn_length", "egg_size"),
                   ~ median(.x, na.rm = TRUE)),
            max_tl = ifelse(all(is.na(max_tl)), NA, max(max_tl, na.rm = TRUE)),
            across(.cols = unique(f.cols[f.cols$trait_type == "binary", ]$trait_name),
                   ~ ifelse(all(is.na(.)), NA, sum(.x, na.rm = TRUE)))) %>%
  ungroup() %>%
  mutate(across(.cols = unique(f.cols[f.cols$trait_type == "binary", ]$trait_name),
                ~ ifelse(is.na(.x), NA, ifelse(.x > 0 , 1, 0)))) #%>%
  #checking whether can categorize
  # mutate(test = spawn_indiff + spawn_g_ns + spawn_g_sc + spawn_ng_bh + spawn_ng_os + spawn_offshore) %>%
  # select(species, spawn_indiff, spawn_offshore, contains("spawn_g_"), contains("spawn_ng_"), test) 
  #can not, so not doing below anymore
  # mutate(repro_type = case_when(spawn_indiff == 1 ~ "indiff",
  #                               spawn_g_ns == 1 ~ "g_ns",
  #                               spawn_g_sc == 1 ~ "g_sc",
  #                               spawn_ng_bh == 1 ~ "ng_bh",
  #                               spawn_ng_os == 1 ~ "ng_os",
  #                               TRUE ~ NA)) %>%
  # select(-c(spawn_indiff, spawn_g_ns, spawn_g_sc, spawn_ng_bh, spawn_ng_os))

#get categorical variables back, remove duplicates + all NA rows to add back in
tmp <- fish %>%
  select(species, matches(paste(unique(f.cols[f.cols$trait_type == "categorical", ]$trait_name), collapse = "|"))) %>%
  distinct() %>%
  filter(!if_all(.cols = matches(paste(unique(f.cols[f.cols$trait_type == "categorical", ]$trait_name), collapse = "|")),
                 is.na))

fish.sum <- fish.sum %>% left_join(., tmp, by = "species", multiple = "all")

write_csv(fish.sum, paste0(PATH, "/20_Traits/trait_dat_summ_fish.csv"))

rm(f1, f2, f3, tmp, tmp1, f.cols2, c.key)

##pair to more conservative species list ----
fish.sum <- read_csv(paste0(PATH, "/20_Traits/trait_dat_summ_fish.csv"))
fish <- read_csv(paste0(PATH, "/01_BioDat/occ_fish_inthigh_long_20240208.csv")) %>% select(order, family, genus, species) %>% distinct()

fish.traits <- left_join(fish, fish.sum)

##identify missing traits ----
na.check <- fish.traits %>%
  mutate(troph = rowSums(!is.na(select(., starts_with("troph_")))),
         hab.pref = rowSums(!is.na(select(., matches("bed|bould|claysilt|cob|debr|grav|lwd|muck|sand|veg")))),
         hab.flow = rowSums(!is.na(select(., matches("fast|mod|slow")))),
         sp.season = rowSums(!is.na(select(., matches("winter|spring|summer|fall")))),
         sp.type = rowSums(!is.na(select(., matches("_ng_|_g_|indiff|offshore")))),
         across(c(troph, hab.pref, hab.flow, sp.season, sp.type), ~ifelse(.x < 1, NA, .x))) %>%
  select(-contains("hab_"), -contains("troph_"), -matches("winter|spring|summer|fall"), -matches("_ng_|_g_|indiff|offshore")) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(-c(order, family, genus, species), names_to = "trait", values_to = "value") %>%
  filter(is.na(value)) %>%
  arrange(species, trait)

write_csv(na.check, paste0(PATH, "/20_Traits/trait_dat_fish_miss_dat_list_20240109.csv"))

rm(list = ls())

##bring in missing traits from lit search ----
PATH <- getwd()

source(paste0(PATH, "/Scripts/XX_trait_functions.R"))

fish.sum <- read_csv(paste0(PATH, "/20_Traits/trait_dat_summ_fish.csv"))

new.t <- read_csv(paste0(PATH, "/20_Traits/trait_dat_fish_miss_dat_list_updated.csv"))
new.t <- new.t[, 1:3] #get rid of information columns

new.t <- new.t %>%
  separate_longer_delim(., value, ";") %>% 
  mutate(value = gsub(" ", "", value)) %>% 
  filter(!is.na(value)) %>%
  mutate(trait = case_when(grepl("sp.season|sp.type", trait) ~ paste0("spawn_", value),
                           grepl("hab.", trait) ~ paste0("hab_", value),
                           grepl("troph", trait) ~ paste0("troph_", value),
                           T ~ trait)) %>%
  pivot_wider(., names_from = trait, values_from = value) %>%
  mutate(across(c(contains("troph"), contains("hab"), matches("spring|summer|winter|fall"), matches("_ng_|_g_|offshore")), 
                ~ ifelse(!is.na(.x), 1, NA)),
         across(-c(species, spawn_freq, temp_pref, col_pos), as.numeric))

tmp <- bind_rows(fish.sum[, !names(fish.sum) %in% "genus"], new.t)

#compress down to one row per species
tmp <- tmp %>%
  mutate(across(c(spawn_freq, temp_pref, col_pos), ~ str_to_lower(.x))) %>%
  group_by(species) %>%
  summarise(across(.cols = c(fecundity, longevity, matur_age, spawn_length, egg_size, max_tl),
                   ~ median(.x, na.rm = TRUE)),
            across(.cols = c(contains("hab_"), contains("spawn_"), contains("troph_"), -spawn_length, -spawn_freq),
                   ~ ifelse(all(is.na(.x)), NA, sum(.x, na.rm = TRUE))),
            across(.cols = c(spawn_freq, temp_pref, col_pos),
                   ~ find_cat_mode(.x)))

#adjust for binary vars in setting actual NAs for later replacement
#essentially, if all values in a category are NA, we can assume we didn't find that variable
tmp <- tmp %>%
  rowwise() %>%
  mutate(across(matches("fast|mod|slow"), ~ ifelse(all(is.na(c_across(matches("fast|mod|slow")))),
                                                   NA, replace_na(.x, 0))),
         across(c(contains("hab_"), -matches("fast|mod|slow")), ~ ifelse(all(is.na(c_across(c(contains("hab_"), -matches("fast|mod|slow"))))),
                                                                         NA, replace_na(.x, 0))),
         across(matches("winter|spring|summer|fall"), ~ ifelse(all(is.na(c_across(matches("winter|spring|summer|fall")))),
                                                               NA, replace_na(.x, 0))),
         across(contains("troph_"), ~ ifelse(all(is.na(c_across(contains("troph_")))),
                                            NA, replace_na(.x, 0))),
         across(matches("g_ns|g_sc|indiff|ng_bh|ng_os|offshore"), ~ ifelse(all(is.na(c_across(matches("g_ns|g_sc|indiff|ng_bh|ng_os|offshore")))),
                                                                           NA, replace_na(.x, 0))))

#pair to spp names
##read in final occurrence set with multiple hierarchy incl.
f.spp <- read_csv(paste0(paste0(PATH, "/01_BioDat/occ_alltax_finalfilter_long_20240718.csv"))) %>%
  filter(bio_type == "fish") %>% 
  select(order, family, genus, species) %>% 
  distinct()

fish.traits <- left_join(f.spp, tmp)

write_csv(fish.traits, paste0(PATH, "/20_Traits/trait_dat_summ_fish_with_updates.csv"))

#summarize amount of missing data
#by trait
sapply(fish.traits, function(x) round(sum(is.na(x)) / length(x) * 100, digits = 2))
#by spp
apply(fish.traits, 1, function(x) round(sum(is.na(x)) / length(x) * 100, digits = 2))

rm(list = ls()) #clean workspace

##fill in missing traits via RF ----
library(missForest); library(tidyverse)

PATH <- getwd()

fish.traits <- read_csv(paste0(PATH, "/20_Traits/trait_dat_summ_fish_with_updates.csv"))

#get binary traits
bi.vars <- sapply(fish.traits, function(x) {
  if(is.character(x) | is.factor(x)) {
    # tmp <- bind_cols(tmp, x)
    return(FALSE)
  } else {
    if(max(x, na.rm = TRUE) == 1 & min(x, na.rm = TRUE) == 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
})
bi.vars <- names(bi.vars[bi.vars == TRUE])

fish.traits <- fish.traits %>%
  mutate(across(c(order, family, temp_pref, col_pos, spawn_freq, matches(bi.vars)), as.factor))

mutate(higher_lvl = case_when(#NEOPTERA (INFRACLASS)
  order %in% c("Orthoptera", "Plecoptera", "Hemiptera") ~ "Polyneoptera/Paraneoptera",
  ##HOLOMETABOLA (SUPERORDER)
  order %in% c("Coleoptera", "Hymenoptera", "Megaloptera", "Neuroptera") ~ "Neuropteroidea +",
  order %in% c("Diptera", "Lepidoptera", "Trichoptera") ~ "Mecoptera",
  #PALAEOPTERA (INFRACLASS)
  order %in% c("Ephemeroptera", "Odonata") ~ "Palaeoptera (infraclass)",
  T ~ NA)) %>%

tmp <- as.data.frame(fish.traits[,-c(3,4)])

#if xtfrm error appears, need to restart R and try again
fish.imp <- missForest(tmp, #including order+family, as taxonomic relationship vars for improved prediction
                       ntree = 1000, #set to 1000 for full run
                       variablewise = TRUE,
                       verbose = TRUE) #return individual var performance or not

imp.traits <- bind_cols(fish.traits[,c(3,4)], fish.imp$ximp)
imp.traits <- imp.traits %>%
  relocate(order, family, genus, species) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

write_csv(imp.traits, paste0(PATH, "/20_Traits/trait_dat_summ_fish_imputed.csv"))

imp.tr.error <- data.frame(trait = names(fish.imp$ximp), 
                           error_type = names(fish.imp$OOBerror), 
                           error_val = fish.imp$OOBerror)

imp.tr.error <- imp.tr.error %>%
  mutate(adj_error = case_when(error_type == "MSE" ~ sqrt(error_val), #RMSE
                               error_type == "PFC" ~ 1 - error_val), #Proportion true classification
         across(c(error_val, adj_error), ~ round(.x, digits = 3))) 

write_csv(imp.tr.error, paste0(PATH, "/20_Traits/impute_trait_error_rates_fish.csv"))

#plot traits
gg.traits <- fish.traits[,-c(1:4)] %>% 
  select(where(is.numeric)) %>% 
  pivot_longer(cols = everything()) %>%
  filter(!is.na(value))
gg.new <- imp.traits %>% 
  select(where(is.numeric)) %>% 
  pivot_longer(cols = everything()) %>%
  filter(!is.na(value))
gg.tr.mean <- gg.traits %>%
  group_by(name) %>%
  summarise(mean = mean(value)) %>%
  left_join(., 
            imp.tr.error[imp.tr.error$error_type == "MSE",] %>% rename(name = trait)
            ) %>%
  mutate(upper = mean + adj_error,
         lower = mean - adj_error)

ggplot() +
  geom_violin(data = gg.traits, aes(x = 1, y = value, fill = name)) +
  geom_point(data = gg.tr.mean, aes(x = 1, y = mean), size = 3, pch = 15) +
  geom_errorbar(data = gg.tr.mean, aes(x = 1, ymax = upper, ymin = lower), width = 0.25) +
  facet_wrap(~ name, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggplot() +
  geom_density(data = gg.traits, aes(x = value)) +
  geom_density(data = gg.new, aes(x = value)) +
  facet_wrap(~ name, scales = "free") +
  theme_classic() +
  theme(axis.title.x = element_blank())

rm(list = ls())

#fill in missing traits via mean of higher tax. lvls ---
# cats <- names(fish.traits)[sapply(fish.traits, is.character)]
# cats <- cats[!cats %in% c("order", "family", "genus", "species")]
# 
# fish.traits <- fish.traits %>%
#   #fill in across increasingly higher levels
#   replace_na_traits(., "genus", cats) %>%
#   replace_na_traits(., "family", cats) %>%
#   replace_na_traits(., "order", cats) %>%
#   replace_na_traits(., NA, cats)
#   
# #update binary vars so that it's either 0 or 1 (by rounding)
# fish.traits <- fish.traits %>%
#   mutate(across(matches(paste0(t.cols[t.cols$taxa == "fish" & t.cols$trait_type == "binary",]$trait_name, collapse = "|")),
#                 ~ round(.x, digits = 0)))
# 
# write_csv(fish.traits, paste0(PATH, "/20_Traits/trait_dat_summ_fish_na_removed.csv"))

#INVERT ----
i.spp <- read_csv(paste0(PATH, "/01_BioDat/bug_species_list.csv"))

i.cols <- t.cols %>%
  filter(taxa == "invert") %>%
  select(trait_name, trait_type, database_col_name)  %>%
  mutate(trait_cat = case_when(grepl(",", database_col_name) ~ sub(".+, (.+)", "\\1", database_col_name),
                               TRUE ~ sub(".+([0-9])", "\\1", database_col_name)),
         database_col_name = case_when(grepl(",", database_col_name) ~ sub(",[^,]+$", "", database_col_name),
                                       TRUE ~ sub("[0-9]", "", database_col_name))) %>%
  separate_rows(., database_col_name, sep = ", ")

##Freshwater insect CONUS ----
i1 <- read_csv(paste0(PATH, "/20_Traits/Genus_Traits.csv")) %>%
  filter(Trait_group %in% i.cols$database_col_name) %>% #get desired traits
  left_join(., i.cols, by = c("Trait_group" = "database_col_name"), multiple = "all") %>% #rename traits for consistency 
  filter(str_detect(Trait, regex(trait_cat, ignore_case = TRUE))) %>% #remove extra rows created by join
  filter(ifelse(Trait_group == "Rheophily_abbrev", Trait == trait_cat, TRUE)) %>% #because of value overlaps here
  mutate(value = 1) %>% #prep for widening dataframe
  group_by(Genus, trait_name) %>%
  summarise(sum = sum(value)) %>%
  mutate(sum = ifelse(sum > 0, 1, NA)) %>%
  pivot_wider(names_from = trait_name, values_from = sum) %>%
  relocate(Genus, contains("troph"), contains("max_size"), contains("repro"), contains("disp"), contains("emerg"), contains("dev"), 
           contains("therm"), contains("habit"), contains("synch"), contains("resp"), contains("rheo"))

#get higher taxonomic rankings to append
tmp <- read_csv(paste0(PATH, "/20_Traits/Ancillary_Taxonomy.csv")) %>% 
  select(Order, Family, Genus) %>%
  #deal with duplicates for joining
  filter(!(Family == "Pyralidae" & Genus == "Acentria")) %>%
  filter(!(Family == "Glossosomatidae" & Genus == "Agabus")) %>%
  filter(!(Family == "Siphlonuridae" & Genus == "Ameletus")) %>%
  filter(!(Family == "Perlodidae" & Genus == "Doroneuria")) %>%
  filter(!(Family == "Helodidae" & Genus == "Elodes")) %>%
  filter(!(Family == "Cordulegastridae" & Genus == "Epitheca")) %>%
  filter(!(Order == "Plecoptera" & Genus == "Glutops")) %>%
  filter(!(Family == "Hydrophilidae" & Genus == "Helophorus")) %>%
  filter(!(Family == "Hydrophilidae" & Genus == "Hydrochus")) %>%
  filter(!(Family == "Dolicopodidae" & Genus == "Hydrophorus")) %>%
  filter(!(Family == "Heptagenidae" & Genus == "Macdunnoa")) %>%
  filter(!(Family == "Macroveliidae" & Genus == "Mesovelia")) %>%
  filter(!(Family == "Pyralidae" & Genus == "Paraponyx")) %>%
  filter(!(Family == "Ephemeridae" & Genus == "Pentagenia")) %>%
  filter(!(Family %in% c("Pyralidae", "Crambiidae") & Genus == "Petrophila")) %>%
  filter(!(Family == "Coenagrionidae" & Genus == "Plathemis")) %>%
  filter(!(Family == "Heptageniidae" & Genus == "Pseudiron")) %>%
  filter(!(Family == "Simulidae" & Genus == "Simulium")) %>%
  filter(!(Family == "Caenidae" & Genus == "Tricorythodes")) %>%
  distinct()

# filter(spp_names, Genus == i1[13539,]$Genus) #to check duplicates

i1 <- i1 %>% left_join(., tmp)
i1 <- i1 %>%
  relocate(Order, Family, Genus) %>%
  mutate(db_orig = "frshwtr_CONUS")
 
##Poff et al. 2006 traits ----
i2 <- readxl::read_excel(paste0(PATH, "/20_Traits/poffetal_traitmatrix.xls"), sheet = 2) %>%
  mutate(Order = str_to_title(Order)) %>%
  select(Order, Family, Genus, matches(unique(i.cols$database_col_name))) %>%
  pivot_longer(cols = matches(unique(i.cols$database_col_name)), names_to = "database_col_name", values_to = "trait_cat") %>%
  mutate(trait_cat = as.character(trait_cat)) %>%
  left_join(., i.cols) %>%
  mutate(trait_cat = ifelse(trait_cat > 0, 1, NA),
         db_orig = "Poffetal2006") %>% #just use this column to make 1 / 0 for wide format
  select(-c(database_col_name, trait_type)) %>%
  pivot_wider(names_from = trait_name, values_from = trait_cat) %>%
  relocate(Order, Family, Genus, contains("swim"), contains("desic"), contains("long")) %>%
  mutate(Family = case_when(Genus == "Pentagenia" ~ "Palingeniidae", T ~ Family))

##combine datasets + summ ----
bug <- bind_rows(i1, i2)

write_csv(bug, paste0(PATH, "/20_Traits/trait_dat_raw_bug.csv"))

bug.sum <- bug %>%
  mutate(Family = ifelse(Family == "Tipulidae2", "Tipulidae", Family)) %>%
  group_by(Order, Family, Genus) %>%
  summarise(across(.cols = intersect(names(bug), t.cols[t.cols$trait_type == "binary", ]$trait_name),
                   ~ ifelse(all(is.na(.)), NA, sum(.x, na.rm = TRUE)))) %>%
  ungroup() %>%
  mutate(across(.cols = intersect(names(bug), t.cols[t.cols$trait_type == "binary", ]$trait_name),
                ~ ifelse(is.na(.x), NA, ifelse(.x > 0 , 1, 0)))) %>%
  mutate(max_size = case_when(max_size_lar == 1 ~ "large",
                              max_size_med == 1 ~ "med",
                              max_size_small == 1 ~ "small",
                              TRUE ~ NA),
         repro_disp = case_when(repro_disp_low == 1 ~ "low",
                                repro_disp_high == 1 ~ "high",
                                TRUE ~ NA),
         disp_strength = case_when(disp_strong == 1 ~ "strong",
                                   disp_weak == 1 ~ "weak",
                                   TRUE ~ NA),
         gen_num = case_when(dev_univolt == 1 ~ "univolt",
                             dev_multivolt == 1 ~ "multivolt",
                             dev_semivolt == 1 ~ "semivolt",
                             TRUE ~ NA),
         therm_pref = case_when(therm_cold == 1 ~ "cold",
                                therm_cold_cool == 1 ~ "cold-cool",
                                therm_cool_warm == 1 ~ "cool-warm",
                                therm_warm == 1 ~ "warm",
                                therm_hot == 1 ~ "hot",
                              TRUE ~ NA),
         synch = case_when(synch_poorly == 1 ~ "poorly",
                           synch_well == 1 ~ "well",
                           TRUE ~ NA),
         resp_type = case_when(resp_teg == 1 ~ "tegument",
                               resp_gill == 1 ~ "gills",
                               resp_spir == 1 ~ "plast_spir",
                               TRUE ~ NA),
         rheo_type = case_when(rheo_depo == 1 ~ "depo",
                               rheo_depo_eros == 1 ~ "depo_eros",
                               rheo_eros == 1 ~ "eros",
                               TRUE ~ NA),
         swim_abil = case_when(swim_none == 1 ~ "none",
                               swim_strong == 1 ~ "strong",
                               swim_weak == 1 ~ "weak",
                               TRUE ~ NA),
         desic_tol = case_when(desic_absent == 1 ~ "absent",
                               desic_present == 1 ~ "present",
                               TRUE ~ NA),
         life_span = case_when(long_vshort == 1 ~ "vshort",
                               long_short == 1 ~ "short",
                               long_long == 1 ~ "long",
                               TRUE ~ NA)) %>%
  select(-c(synch_well, synch_poorly, resp_teg, resp_gill, resp_spir,
            rheo_depo, rheo_depo_eros, rheo_eros, swim_none, swim_weak,
            swim_strong, desic_absent, desic_present, long_vshort, long_short,
            long_long, dev_semivolt, dev_univolt, dev_multivolt, therm_cold,
            therm_cold_cool, therm_cool_warm, therm_warm, therm_hot, max_size_small,
            max_size_med, max_size_lar, repro_disp_low, repro_disp_high, disp_weak, disp_strong)) %>%
  mutate(across(-c(Order, Family, Genus), as.factor)) %>%
  filter(Genus %in% i.spp$genus | Genus %in% i.spp$tribe | Genus %in% i.spp$subfamily)

write_csv(bug.sum, paste0(PATH, "/20_Traits/trait_dat_summ_bug.csv"))

rm(i1, i2, tmp)

##identify missing traits ----
tmp <- bug.sum %>%
  mutate(troph = rowSums(!is.na(select(., starts_with("troph_")))),
         emerg = rowSums(!is.na(select(., starts_with("emerg_")))),
         habit = rowSums(!is.na(select(., starts_with("habit_")))),
         across(c(troph, emerg, habit), ~ifelse(.x < 1, NA, .x))) %>%
  mutate(taxa_name = ifelse(is.na(Genus), Family, Genus)) %>%
  select(-contains("troph_"), -contains("emerg_"), -contains("habit_"), -c(Order, Family, Genus))

na.check <- i.spp %>%
  mutate(taxa_name = ifelse(is.na(genus), tribe, genus),
         taxa_name = ifelse(is.na(taxa_name), subfamily, taxa_name),
         taxa_name = ifelse(is.na(taxa_name), family, taxa_name)) %>%
  filter(!is.na(taxa_name)) %>%
  left_join(., tmp) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(-c(order, family, subfamily, tribe, genus, taxa_name), names_to = "trait", values_to = "value") %>%
  filter(is.na(value))

write_csv(na.check, paste0(PATH, "/20_Traits/trait_dat_bug_miss_dat_list_", gsub("-", "", Sys.Date()), ".csv"))

rm(list = ls())

##bring in missing traits from lit search ----

#TO BE DONE AT A LATER DATE WHEN SEARCH IS COMPLETE

##pair to species list + fill in missing data ----
PATH <- getwd()

bug.sum <- read_csv(paste0(PATH, "/20_Traits/trait_dat_summ_bug.csv")) %>%
  rename(genus = Genus, family = Family, order = Order) %>%
  mutate(join_track = "joined")
bug <- read_csv(paste0(PATH, "/01_BioDat/occ_bug_inthigh_long.csv"), col_types = cols(.default = "c")) %>%
  select(order, family, subfamily, tribe, genus) %>%
  filter(!(is.na(subfamily) & is.na(tribe) & is.na(genus))) %>%
  mutate(taxa = case_when(tribe == "Aciliini" ~ tribe,
                          tribe == "Hydrophilini" ~ tribe,
                          tribe == "Pseudochironomini" ~ tribe,
                          !is.na(genus) ~ genus,
                          !is.na(tribe) ~ tribe,
                          !is.na(subfamily) ~ subfamily)) %>%
  distinct()

#attempt to join traits to species list
#genus level first
bug.traits <- left_join(bug, bug.sum, by = c("genus"="genus"), na_matches = "never") %>%
  filter(!is.na(join_track)) %>%
  select(-contains(".y")) %>%
  rename_with(~gsub("\\.x", "", .x))

tmp <- anti_join(bug, bug.traits, by = c("genus")) #taxa without traits
tmp.tr <- anti_join(bug.sum, bug.traits, by = c("genus")) #traits not yet joined

tmp.tr <- left_join(tmp, tmp.tr, by = c("taxa"="genus"), na_matches = "never") %>%
  select(-contains(".y")) %>%
  rename_with(~gsub("\\.x", "", .x))

bug.traits <- bind_rows(bug.traits, tmp.tr) %>%
#   filter(!is.na(join_track)) %>%
    select(-join_track)

#essentially, if all values in a category are NA, we can assume we didn't find that variable
bug.traits <- bug.traits %>%
  rowwise() %>%
  mutate(across(contains("troph"), ~ ifelse(all(is.na(c_across(contains("troph")))),
                                                   NA, replace_na(.x, 0))),
         across(contains("emerg"), ~ ifelse(all(is.na(c_across(contains("emerg")))),
                                            NA, replace_na(.x, 0))),
         across(contains("habit"), ~ ifelse(all(is.na(c_across(contains("habit")))),
                                            NA, replace_na(.x, 0)))) %>%
  mutate(habit_crawl = ifelse(habit_crawl == FALSE, 0, NA))

write_csv(bug.traits, paste0(PATH, "/20_Traits/trait_dat_summ_bug_with_updates.csv"))

rm(list = ls())

##fill in missing traits via RF ----
library(missForest); library(tidyverse)

PATH <- getwd()

bug.traits <- read_csv(paste0(PATH, "/20_Traits/trait_dat_summ_bug_with_updates.csv"))

#handle duplicates from collapsing down for taxa level
source(paste0(PATH, "/Scripts/XX_trait_functions.R"))
bug.traits$taxa[duplicated(bug.traits$taxa)]
dups <- bug.traits[duplicated(bug.traits$taxa) | duplicated(bug.traits$taxa, fromLast = TRUE),]
tmp <- dups[, names(dups) %in% c("order", "family", "subfamily", "tribe", "genus", "taxa")] %>% mutate(genus = taxa) %>% distinct()
dups <- compress_traits(data = dups, "taxa", 
                c("max_size", "repro_disp", "disp_strength", "gen_num", "therm_pref", 
                  "synch", "resp_type", "rheo_type", "swim_abil", "desic_tol", "life_span")) %>%
  mutate(across(matches("troph|emerg|habit"), ~ ifelse(.x > 0, 1, 0)))
dups <- left_join(tmp, dups)

#add back in
bug.traits <- bind_rows(bug.traits[!duplicated(bug.traits$taxa) & !duplicated(bug.traits$taxa, fromLast = TRUE),], dups) #should be 10 missing

# #summarize amount of missing data
# sapply(bug.traits, function(x) round(sum(is.na(x)) / length(x) * 100, digits = 2))

#note: randomforest cannot handle categorical predictors with more than 53 categories, so need to break down for incl. taxa relationships in imputation
#breaking down based on relationships in higher taxonomic levels 
#(source: itis.gov, https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-14-52/figures/1)
bug.traits <- bug.traits %>%
  mutate(higher_lvl = case_when(#NEOPTERA (INFRACLASS)
                                order %in% c("Orthoptera", "Plecoptera", "Hemiptera") ~ "Polyneoptera/Paraneoptera",
                                ##HOLOMETABOLA (SUPERORDER)
                                order %in% c("Coleoptera", "Hymenoptera", "Megaloptera", "Neuroptera") ~ "Neuropteroidea +",
                                order %in% c("Diptera", "Lepidoptera", "Trichoptera") ~ "Mecoptera",
                                #PALAEOPTERA (INFRACLASS)
                                order %in% c("Ephemeroptera", "Odonata") ~ "Palaeoptera (infraclass)",
                                T ~ NA)) %>%
  mutate(across(-c(family, genus, taxa), as.factor))
#get counts:
# bug.traits %>%
#   group_by(higher_lvl) %>%
#   summarise(unique_fam = n_distinct(family),
#             unique_taxa = n_distinct(taxa))

bug.list <- split(bug.traits, bug.traits$higher_lvl) #split for meeting randomforest criteria

imputed.results <- lapply(bug.list, function(x) {
  
  x$family <- as.factor(x$family)
  tmp.tr <- as.data.frame(x[, !names(x) %in% c("subfamily", "tribe", "genus", "taxa", "higher_lvl", "index")])
  
  #note: if xtfrm error appears, need to restart R and try again
  imp <- missForest(tmp.tr, #including order+family, as taxonomic relationship var for improved prediction
                    ntree = 1000, #set to 1000 for full run
                    variablewise = TRUE,
                    verbose = TRUE) #return individual var performance or not
  
  imp.traits <- bind_cols(x[, names(x) %in% c("subfamily", "tribe", "genus", "taxa", "higher_lvl", "index")], imp$ximp)
  
  imp.tr.error <- data.frame(trait = names(imp$ximp), 
                             error_type = names(imp$OOBerror), 
                             error_val = imp$OOBerror)
  
  res <- list(traits = imp.traits, error = imp.tr.error)
  
  return(res)
  
})

#bind back together
imp.traits <- data.frame()
imp.tr.error <- data.frame()
for(i in seq_along(imputed.results)) {
  imp.traits <- bind_rows(imp.traits, imputed.results[[i]][["traits"]])
  imp.tr.error <- bind_rows(imp.tr.error, imputed.results[[i]][["error"]])
}

write_csv(imp.traits %>% relocate(order, family, subfamily, tribe, genus, taxa) %>% select(-higher_lvl), paste0(PATH, "/20_Traits/trait_dat_summ_bug_imputed.csv"))


imp.tr.error.sum <- imp.tr.error %>%
  mutate(adj_error = case_when(error_type == "MSE" ~ sqrt(error_val), #RMSE
                               error_type == "PFC" ~ 1 - error_val)) %>% #Proportion true classification
  group_by(trait) %>%
  summarise(mn_error = mean(error_val, na.rm = TRUE),
            mn_prop_true = mean(adj_error, na.rm = TRUE))

write_csv(imp.tr.error.sum, paste0(PATH, "/20_Traits/impute_trait_error_rates_bug_summ.csv"))
write_csv(imp.tr.error, paste0(PATH, "/20_Traits/impute_trait_error_rates_bug_raw.csv"))

##fill in missing traits via mean of higher tax. lvls ----
# cats <- names(bug.traits)[sapply(bug.traits, is.character)]
# 
# #since currently only 1, no 0s, assume missing only if genera has whole row missing
# bug.traits <- bug.traits %>%
#   mutate(across(everything(), as.character)) %>%
#   rowwise() %>%
#   mutate(missing_flag = ifelse(
#     all(is.na(c_across(matches(names(bug.sum)[!names(bug.sum) %in% c("order", "family", "genus")])))),
#     "F", NA)) %>%
#   ungroup() %>%
#   mutate(across(-matches(cats), ~ ifelse(is.na(missing_flag) & is.na(.x), 0, .x)),
#          across(-matches(cats), as.numeric)) %>%
#   select(-missing_flag)
#   
# cats <- cats[!cats %in% c("order", "family", "subfamily", "tribe", "genus", "taxa")]
# 
# bug.traits <- bug.traits %>%
#   #fill in across increasingly higher levels
#   replace_na_traits(., "genus", cats) %>%
#   replace_na_traits(., "tribe", cats) %>%
#   replace_na_traits(., "subfamily", cats) %>%
#   replace_na_traits(., "family", cats) %>%
#   replace_na_traits(., "order", cats) %>%
#   select(-c(order, family, subfamily, tribe, genus)) %>%
#   compress_traits(., "taxa", cats)
# 
# #update binary vars so that it's either 0 or 1 (by rounding)
# bug.traits <- bug.traits %>%
#   mutate(across(where(is.numeric),
#                 ~ round(.x, digits = 0)))
# 
# write_csv(bug.traits, paste0(PATH, "/20_Traits/trait_dat_summ_bug_na_removed.csv"))
# 
rm(list = ls())

#ALL TAXA ----
library(tidyverse)

PATH <- getwd()
source(paste0(PATH, "/Scripts/XX_trait_functions.R"))

#sites
file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                         long = col_number(),
                                                         site_id = col_character(),
                                                         COMID = col_character(),
                                                         gage_no_15yr = col_character(),
                                                         dist2gage_m_15yr = col_number(),
                                                         dist2strm_m_flw = col_number()))
# occ.list <- list(occ.list[[2]]) #for now, without imputed bug data

#wide traits load in
file.list <- list.files(paste0(PATH, "/20_Traits/"), pattern = "trait_dat_summ_(fish|bug)_imputed", full.names = TRUE)
traits <- lapply(file.list, read_csv)

## trait clustering ----
library(cluster); library(vegan)

#this can be used to identify the possible levels within a trait dataset
# lvls <- lapply(t_mod, function(x) {
#   lapply(select(x, where(is.character)), unique) #check potential cat levels
# })
#for ordinal variables
##FISH
oc <- list(
  spawn_freq = c("single", "multiple"), 
  temp_pref = c("cold", "cold/cool", "cool", "cool/warm", "warm") 
  # col_pos = c("Benthic", "Non-benthic") #not ordered
)
##BUG

group_res_list <- comp_trait_groups(trait_dat =  traits[[1]][,-c(1:3)], ord_col = oc, max_k = 100)

###plots to look at groups ----
#dendrogram
library(dendextend)
pltree(group_res_list$hierarchical_cluster_results, cex = 0.6, labels = FALSE)

dend <- group_res_list$hierarchical_cluster_results %>% as.dendrogram %>% hang.dendrogram
par(mar = c(15,2,1,1))
dend %>% color_branches(k=10) %>% color_labels(k=15) %>% plot
dend %>% rect.dendrogram(k=10)

#compare silhouette, r2
clust.df <- group_res_list[["group_comp_results"]]
group_res_list[["hierarchical_cluster_results"]]$height
plot(clust.df$k, clust.df$mn_sil_w, type = "b", pch = 19)
plot(clust.df$k, clust.df$grp_r2, type = "b", pch = 19)

ggplot() +
  geom_point(data = clust.df, aes(x = mn_sil_w, y = grp_r2, fill = k), shape = 21, size = 3) +
  scale_fill_viridis_c() +
  theme_minimal()

#pcoa vizualize
plot_clusters(trait_group_output = group_res_list, num_clust = 20, return_pc = c(1,2),
              ellipse = TRUE, vectors = TRUE, trait_dat = traits[[1]][,-c(1:3)], ord_col = oc)

table(cutree(group_res_list[["hierarchical_cluster_results"]], k = 20))

#assign groups and make site_x_trait matrix
clust <- cutree(group_res_list$hierarchical_cluster_results, k = 20)
clust <- as.data.frame(clust) %>% rownames_to_column("species") %>% mutate(clust = paste0("clust", clust))

clust.wide <- site_x_trait(occ.list[[1]], clust)
write_csv(clust.wide[[1]], paste0(PATH, "/20_Traits/site_x_trait_fish_presence-absence_clust.csv"))
write_csv(clust.wide[[2]], paste0(PATH, "/20_Traits/site_x_trait_fish_abundance_clust.csv"))

## site by trait matrices ----
#FISH
#convert continuous vars
fish <- site_x_trait(occ.list[[2]], traits[[2]][,-c(1:3)])

write_csv(fish[[1]], paste0(PATH, "/20_Traits/site_x_trait_fish_presence-absence_cont2cat.csv"))
write_csv(fish[[2]], paste0(PATH, "/20_Traits/site_x_trait_fish_abundance_cont2cat.csv"))

#don't convert
fish <- site_x_trait(occ.list[[2]], traits[[2]][,-c(1:3)], convert_cont_to_cat = FALSE)

write_csv(fish[[1]], paste0(PATH, "/20_Traits/site_x_trait_fish_presence-absence_cont2mn.csv"))
write_csv(fish[[2]], paste0(PATH, "/20_Traits/site_x_trait_fish_abundance_cont2mn.csv"))

#BUG
#convert continuous vars
bug <- site_x_trait(occ.list[[1]], traits[[1]][,-c(1:5)])

write_csv(bug[[1]], paste0(PATH, "/20_Traits/site_x_trait_bug_presence-absence_cont2cat.csv"))
write_csv(bug[[2]], paste0(PATH, "/20_Traits/site_x_trait_bug_abundance_cont2cat.csv"))

#don't convert
bug <- site_x_trait(occ.list[[1]], traits[[1]][,-c(1:5)], convert_cont_to_cat = FALSE)

write_csv(bug[[1]], paste0(PATH, "/20_Traits/site_x_trait_bug_presence-absence_cont2mn.csv"))
write_csv(bug[[2]], paste0(PATH, "/20_Traits/site_x_trait_bug_abundance_cont2mn.csv"))

# summarize + plot ----
library(vegan); library(ape)
fish.sum <- read_csv(paste0(PATH, "/20_Traits/trait_dat_summ_fish.csv"))

##pcoa ----
#prep dataframe
tmp <- fish.sum %>% 
  select(-c(genus, spp_ep)) %>% 
  mutate(species = sub(" ", "_", species)) %>% 
  column_to_rownames("species") %>%
  mutate(across(matches(paste(t.cols[t.cols$trait_type == "binary" & t.cols$database_source == "FishTraits", ]$trait_name, collapse = "|")),
                as.factor),
         repro_type = as.factor(repro_type)) 
#scale and center numeric data
tmp <- cbind(apply(tmp %>% select(where(is.numeric)), 2, scale, center = TRUE, scale = TRUE),
             tmp %>% select(-where(is.numeric)))

#gower distance
dist <- cluster::daisy(tmp, metric = "gower", stand = TRUE)

pcoa <- pcoa(dist, correction = "lingoes", rn = rownames(tmp))
# pcoa <- cmdscale(dist, eig = T, k = 10)

#check
# ordiplot(pcoa, display = 'sites', type = 'text')

fi_pcoa <- as.data.frame(pcoa$vectors[,1:3]) %>% 
  rownames_to_column("species") %>% 
  rename(PC1 = Axis.1, PC2 = Axis.2, PC3 = Axis.3) %>%
  mutate(species = sub("_", " ", species)) %>%
  left_join(., f.spp) %>%
  relocate(order, family, genus, species) %>%
  group_by(family) %>%
  mutate(fam_n = n()) %>%
  left_join(., tmp %>% rownames_to_column("species") %>% mutate(species = sub("_", " ", species))) %>%
  mutate(repro_type = case_when(grepl("g_ns", repro_type) ~ "guarder, nest spawner",
                                grepl("indiff", repro_type) ~ "bearer, substrate indifferent",
                                grepl("g_sc", repro_type) ~ "guarder, substrate chooser",
                                grepl("ng_bh", repro_type) ~ "nonguarder, brood hider",
                                grepl("ng_os", repro_type) ~ "nonguarder, open substrate",
                                TRUE ~ NA))

#get trait arrows
vec.tr <- envfit(pcoa$vectors, tmp, perm = 1000, na.rm = TRUE)
# vec.tr <- envfit(scores(pcoa), tmp, perm=1000, na.rm = TRUE) cmdscale method

t.scrs <- bind_rows(as.data.frame(scores(vec.tr, display = "vectors")),
                    as.data.frame(scores(vec.tr, display = "factors"))) %>%
  rownames_to_column("trait") %>%
  rename(PC1 = Axis.1, PC2 = Axis.2) %>%
  mutate(cat = case_when(grepl("0", trait) ~ "0",
                         grepl("1", trait) ~ "1",
                         grepl("g_ns", trait) ~ "guarder, nest spawner",
                         grepl("indiff", trait) ~ "bearer, substrate indifferent",
                         grepl("g_sc", trait) ~ "guarder, substrate chooser",
                         grepl("ng_bh", trait) ~ "nonguarder, brood hider",
                         grepl("ng_os", trait) ~ "nonguarder, open substrate",
                         TRUE ~ NA),
         trait_name = gsub("0|1|g_ns|g_sc|indiff|ng_bh|ng_os", "", trait))

#get centroid
centroid <- at_pcoa %>% select(PC1, PC2) %>% mutate(across(everything(), mean)) %>% distinct()

#plot
ggplot() +
  geom_point(data = fi_pcoa, aes(x = PC1, y = PC2, fill = repro_type), size = 3, shape = 21, alpha = 0.7) +
  #family color plot
  # geom_point(data = fi_pcoa %>% filter(fam_n < 15), aes(x = PC1, y = PC2), fill = "gray", size = 1.5, shape = 21, alpha = 0.5) +
  # geom_point(data = fi_pcoa %>% filter(fam_n > 15), aes(x = PC1, y = PC2, fill = family), size = 2, shape = 21, alpha = 0.8) +
  #continuous vectors
  geom_segment(data = t.scrs %>% filter(is.na(cat)),
               aes(x = 0, xend = PC1*2, y = 0, yend = PC2*2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black", linewidth = 0.5) +
  geom_text(data = t.scrs %>% filter(is.na(cat)), 
            aes(x = PC1*2, y = PC2*2, label = trait), size = 4, angle = 0, nudge_y = -0.005, nudge_x = -0.04) +
  #cat vectors
  geom_point(data = t.scrs %>% filter(trait_name == "repro_type"), 
             aes(x = PC1, y = PC2, fill = cat), size = 4, shape = 24) +
  scale_fill_viridis_d(direction = -1, name = "", na.translate = F) +
  theme_bw() 
  # theme(legend.position = "bottom") +
  # guides(fill = guide_legend(nrow = 2))
ggsave(paste0(PATH, "/99_figures/fish_pcoa.png"), width = 8, height = 5)

# pcoa$values$Corr_eig/sum(pcoa$values$Corr_eig)*100


fish <- read_csv(paste0(PATH, "/20_Traits/comb_traits_long_fish.csv"))
f.spp <- read_csv(paste0(PATH, "/01_BioDat/occ_fish_all_long_20231128.csv")) %>% select(order, family, genus, species) %>% #fish list
  distinct()

fish.w <- fish %>% 
  select(-db_orig) %>%
  pivot_wider(names_from = trait, values_from = tra_val, values_fill = NA, values_fn = median)

#which spp with no data?
f.spp %>% filter(!is.na(species)) %>%
  filter(!species %in% fish.w$species) %>% View

#find traits x NA percent
na.summ <- as.data.frame(colSums(is.na(fish.w)) / nrow(f.spp) * 100)
names(na.summ) <- "perc_missing"
na.summ <- na.summ %>%
  rownames_to_column("trait") %>%
  filter(!trait %in% c("genus", "species", "spp_ep"))
na.summ <- bind_rows(na.summ %>% slice_max(perc_missing, n = 10), 
                     na.summ %>% slice_min(perc_missing, n = 10))
##plot
ggplot() +
  geom_col(data = na.summ, aes(x = perc_missing, y = reorder(trait, perc_missing))) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), limits = c(0, 100)) +
  ylab("")
ggsave(paste0(PATH, "/99_figures/traits_pct_missing_4trait.png"), height = 5, width = 5)

#find fish x NA percent
fish.na <- fish.w %>%
  pivot_longer(cols = -c(genus, spp_ep, species), names_to = "trait", values_to = "tra_val") %>%
  group_by(species) %>%
  summarise(perc_missing = sum(is.na(tra_val)) / n() * 100)
fish.na <- bind_rows(fish.na %>% slice_max(perc_missing, n = 10), 
                     fish.na %>% slice_min(perc_missing, n = 10))
##plot
ggplot() +
  geom_col(data = fish.na, aes(x = perc_missing, y = reorder(species, perc_missing))) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), limits = c(0, 100)) +
  ylab("")
ggsave(paste0(PATH, "/99_figures/traits_pct_missing_4fish.png"), height = 5, width = 5)



