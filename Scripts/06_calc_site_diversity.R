## CALC DIVERSITY METRICS AT EACH SITE ##

library(tidyverse); library(FD); library(sf)
options(readr.show_col_types = FALSE)

PATH <- getwd()

file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                         long = col_number(),
                                                         site_id = col_character(),
                                                         COMID = col_character(),
                                                         gage_no_15yr = col_character(),
                                                         dist2gage_m_15yr = col_number(),
                                                         dist2strm_m_flw = col_number())) 
file.list <- list.files(paste0(PATH, "/20_Traits/"), pattern = "na_removed", full.names = TRUE)
traits <- lapply(file.list, read_csv)

#site richness ----
rich_res <- lapply(occ.list, function(x) {
  t_names <- names(x)[!names(x) %in% c("lat", "long", "site_id", "COMID", "flw_type", "gage_no_15yr", "dist2gage_m_15yr", "dist2strm_m_flw")]
  
  tmp <- x %>%
    mutate(rich = rowSums(select(., all_of(t_names)), na.rm = TRUE)) %>%
    select(-matches(t_names))
  
  tmp2 <- x %>%
    pivot_longer(cols = matches(t_names), names_to = "species", values_to = "p") %>%
    filter(p != 0) %>%
    group_by(species) %>%
    summarise(ttl_sites = n())
  
  out_list <- list(tmp, tmp2)
  
  return(out_list)
})

for(i in seq_along(rich_res)) {
  if(i == 1) { taxa <- "bug" } else { taxa <- "fish" }
  
  write_csv(rich_res[[i]][[1]], paste0(PATH, "/98_result_tables/site_div_rich_", taxa, ".csv"))
  write_csv(rich_res[[i]][[2]], paste0(PATH, "/98_result_tables/site_div_species_x_site_", taxa, ".csv"))
}


#trait diversity ----
##prep traits ----
t_mod <- lapply(traits, function(t) {
  t %>%
    select(-c(matches(c("order", "family", "genus", "record_ttl_15kbuff")))) %>%
    { if("species" %in% names(.)) rename(., taxa = species) else . } %>%
    mutate(taxa = gsub(" ", "_", taxa)) %>%
    column_to_rownames("taxa")
})

###get and set ordered traits ----
lvls <- lapply(t_mod, function(x) {
  lapply(select(x, where(is.character)), unique) #check potential cat levels
})

# dput(lvls[[2]])
#have to manually alter order
lvls[[1]] <- list(
                  max_size = c("small", "med", "large"), 
                  repro_disp = c("low", "high"), 
                  disp_strength = c("weak", "strong"), 
                  gen_num = c("univolt", "semivolt", "multivolt"), 
                  therm_pref = c("cold", "cold-cool", "cool-warm", "warm", "hot"), 
                  synch = c("poorly", "well"), 
                  # resp_type = c("tegument", "gills", "plast_spir"), #not ordered
                  rheo_type = c("depo", "depo_eros", "eros"), 
                  swim_abil = c("none", "weak", "strong"), 
                  desic_tol = c("absent", "present"), 
                  life_span = c("vshort", "short", "long")
                  )
lvls[[2]] <- list(
                  spawn_freq = c("Single", "Multiple"), 
                  temp_pref = c("Cold", "Cold/cool", "Cool", "Cool/warm", "Warm") 
                  # col_pos = c("Benthic", "Non-benthic") #not ordered
                  )

for(i in seq_along(lvls)) {
  l <- lvls[[i]]
  
  for (col_name in names(l)) {
    if (col_name %in% names(t_mod[[i]])) {
      t_mod[[i]][[col_name]] <- factor(t_mod[[i]][[col_name]], levels = l[[col_name]], ordered = TRUE)
    }
  }
}


##calc div ----
div_list <- list()
for(i in seq_along(occ.list)) {
  x <- occ.list[[i]]
  t <- t_mod[[i]]
  
  x <- x[, names(x) %in% rownames(t)]
  x <- x[sort(names(x))]
  t <- t[rownames(t) %in% names(x), ]
  t <- t[match(names(x), rownames(t)),]
  fdiv <- dbFD(x = t, a = x, corr = "cailliez", calc.FRic = FALSE)
  
  if(i == 1) { taxa <- "bug" } else { taxa <- "fish" }
  
  #make df
  div_list[[i]] <- data.frame(taxa = rep(taxa, length(fdiv$nbsp)),
                              site_id = occ.list[[i]]$site_id,
                              n_sp = fdiv$nbsp,
                              n_fsp = fdiv$sing.sp,
                              f_eve = fdiv$FEve,
                              # f_div = fdiv$FDiv,
                              f_disp = fdiv$FDis,
                              f_raoq = fdiv$RaoQ)
}
  
div_list <- bind_rows(div_list)
write_csv(div_list, paste0(PATH, "/98_result_tables/site_div_alltax_alldiv_raw.csv"))

#summarize results
div_sum <- div_list %>%
  group_by(taxa) %>%
  summarise(across(-c(site_id), 
                   list(mean = ~ mean(.x, na.rm = T),
                        min = ~ min(.x, na.rm = T),
                        max = ~ max(.x, na.rm = T)),
                   .names = "{.col}_{.fn}"))
write_csv(div_sum, paste0(PATH, "/98_result_tables/site_div_alltax_alldiv_summ.csv"))

