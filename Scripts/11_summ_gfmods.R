### EXAMINE + COMPARE GF MODEL OUTPUTS ###
#make summary tables + figures for comparing variable types, taxa, gage types etc.

library(tidyverse); library(gradientForest)
options(readr.show_col_types = FALSE)

PATH <- getwd()

#select model(s) by run type ----
taxa <- "fish"
bio.type <- "(taxonomic_pa|trait_cont2cat)"
flow.type <- "(GW|RO|Int)"
var.type <- "hydro_fm24"
model.run.date <- "2024_08_13" #folder of model rds's to load in 

#prep functions ----
## to match variable name with category later
var_cat_match <- function(var_name) {
  #find which name matches
  matching_var <- sapply(env.table$variable_col_name, function(x) grepl(x, var_name))
  matching_var <- which(matching_var)
  #find category for that name
  matching_cat <- env.table[matching_var, ]$data_category
  # matching_name <- env.table[matching_var, ]$variable_full_name
  
  df <- data.frame(input_name = var_name, env_table_name = names(matching_var), data_category = matching_cat)
  
  return(df)
}

##function to read in model and set inputs: 
gf_model_load <- function(gf_name) {
  
  n <- sub("gf_", "", sub(".rds", "", gf_name))
  gf.out <- readRDS(paste(PATH, "10_GFOutput", model.run.date, gf_name, sep = "/"))
  
  gf.out$type_name <- n
  
  return(gf.out)                  
}

##function to get model input as columns in result dataframe:
add_info_cols <- function(df, gf_model) {
  df <- df %>%
    mutate(flow = case_when(grepl("GW", gf_model$type_name) ~ "GW",
                            grepl("Int", gf_model$type_name) ~ "Int",
                            grepl("RO", gf_model$type_name) ~ "RO"),
           taxa = case_when(grepl("bug", gf_model$type_name) ~ "bug",
                            grepl("fish", gf_model$type_name) ~ "fish"),
           run_type = case_when(grepl("trait", gf_model$type_name) ~ "trait",
                                grepl("taxonomic", gf_model$type_name) ~ "taxonomic"))
  
  return(df)
}

##function to get nlcd names:
add_nlcd_names <- function(df) {
  df <- df %>%
    mutate(ID = case_when(grepl("nlcd", env_var) ~ str_extract(env_var, "\\d{2}$"),
                   T ~ NA),
           ID = as.numeric(ID)) %>%
    left_join(., FedData::nlcd_colors()[, names(FedData::nlcd_colors()) %in% c("ID", "Class")]) %>%
    mutate(Class = str_to_lower(gsub(" ", "_" , gsub("[,/()]", " ", Class)))) %>%
    rowwise() %>%
    mutate(env_var = case_when(
      grepl("nlcd_pct", env_var) ~ gsub("\\d{2}$", as.character(Class), gsub("_nlcd", "", env_var)),
      TRUE ~ env_var
    )) %>%
    ungroup() %>%
    select(-ID, -Class)
  return(df)
}


#load necessary inputs ----
#for getting variable categories
env.table <- read_csv("02_EnvDat/environmental_variable_info.csv") %>% filter(!is.na(variable_col_name)) 

#get list of models specified to load in
gf.files <- list.files(paste0(PATH, "/10_GFOutput/", model.run.date), pattern = ".rds")

#get colors for plots
source(paste0(PATH, "/Scripts/XX_colors.R"))

#set file save name for all files
save.file.name <- gsub("\\|", "_", gsub("\\)", "", gsub("\\(", "", paste(taxa, bio.type, flow.type, var.type, sep = "_"))))

#get names of models
gf.set <- gf.files[grepl(paste(taxa, bio.type, flow.type, var.type, sep = "_"), gf.files)]

#load models specified
gf.list <- lapply(gf.set, gf_model_load)

#variable importance (weighted R2) ----
#note: same as max. cum imp value
##table ----
variable_importance_table <- function(x) {
  
  if(is.character(x)) {
    message("loading model rds...")
    x <- gf_model_load(x)
  }
  
  imp.vars <- data.frame(importance(x)) %>%
    rownames_to_column("env_var")
  names(imp.vars)[[2]] <- "weighted_r2"
  imp.vars <- arrange(imp.vars, desc(weighted_r2))
  
  imp.vars$type_name <- x$type_name
  
  imp.vars <- add_info_cols(imp.vars, x)
  
  return(imp.vars)
  
}

# vi.table <- variable_importance_table(gf.list[[1]])

vi.table <- lapply(gf.list, variable_importance_table)
vi.table <- bind_rows(vi.table)

write_csv(vi.table, paste0(PATH, "/98_result_tables/variable_importance_R2_", save.file.name, ".csv"))

##plot ----
vi.table <- read_csv(paste0(PATH, "/98_result_tables/variable_importance_R2_", save.file.name, ".csv"))

#get categories for colors
cat.df <- lapply(vi.table$env_var, var_cat_match)
cat.df <- bind_rows(cat.df) %>% select(-env_table_name) %>% distinct()
# any(duplicated(cat.df$input_name))

#prep table for plotting
vi.table <- vi.table %>%
  left_join(., cat.df, by = c("env_var" = "input_name")) %>%
  mutate(env_type = case_when(grepl("spatial", data_category) ~ "spatial",
                              grepl("hydrology", data_category) ~ "hydrology",
                              grepl("lithology", data_category) ~ "lithology",
                              grepl("soil", data_category) ~ "soil",
                              grepl("climate", data_category) ~ "climate",
                              grepl("land cover|land use", data_category) ~ "land cover")) 

if(all(vi.table$env_type == "hydrology")) {
  vi.table <- vi.table %>%
    mutate(env_type = case_when(grepl("stream temperature", data_category) ~ "stream temperature",
                                grepl("frequency", data_category) ~ "frequency",
                                grepl("timing", data_category) ~ "timing",
                                grepl("magnitude", data_category) ~ "magnitude",
                                grepl("duration", data_category) ~ "duration",
                                grepl("rate of change", data_category) ~ "rate of change",
                                grepl("cumulative", data_category) ~ "summarized flow"))
}

if(var.type == "lulc") {
  vi.table <- vi.table %>%
    mutate(env_type = case_when(env_type == "spatial" ~ "spatial",
                                grepl("imp_|AGDRAIN|RDCRS|RDDENS|DAMDENS|_21|_22|_23|_24|_81|_82", env_var) ~ "land use",
                                T ~ "land cover"))
  vi.table <- add_nlcd_names(vi.table)
}

fl.lvl <- vi.table %>%
  group_by(flow) %>%
  arrange(weighted_r2) %>%
  do(tibble(al=levels(reorder(interaction(.$flow, .$env_var, drop=TRUE), .$weighted_r2)))) %>%
  pull(al)

plot.df <- vi.table %>%
  group_by(flow, run_type) %>%
  slice_max(weighted_r2, n = 10) %>%
  mutate(al = factor(interaction(flow, env_var), levels = fl.lvl))

#plot
###bar chart ----
#set color scheme
if(grepl("hydro", var.type)) { pal <- hydro.pal }
if(grepl("baseline|alt", var.type)) { pal <- var.high.pal }
if(grepl("lulc", var.type)) { pal <- lulc.pal }

p1 <- ggplot(data = plot.df, aes(x = weighted_r2, y = al, fill = env_type)) +
  geom_col(color = "black", linewidth = 0.1) +
  facet_grid(vars(flow), vars(run_type), scales = "free") +
  scale_x_continuous(expand = c(0.001,0), limits = c(0, max(vi.table$weighted_r2, na.rm = TRUE)+0.0005)) +
  scale_y_discrete(breaks = fl.lvl, labels = sub("^[^.]+\\.", "", fl.lvl), name = "") +
  scale_fill_manual(values = pal, name = "Variable Type") +
  labs(x = expression(paste(R^2, " weighted importance"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = "bottom",
        panel.spacing.x = unit(1, "lines"))
p1
ggsave(paste0(PATH, "/99_figures/var_imp_barchart_", save.file.name, ".png"), plot = p1, width = 8, height = 8)

###dot plot ----
library(tidytext)

#prep data
plot.df <- vi.table %>%
  select(-type_name) %>%
  pivot_wider(names_from = "flow", values_from = "weighted_r2") %>%
  mutate(mean = rowMeans(select(., GW, Int, RO), na.rm = TRUE),
         env_fac = reorder_within(env_var, mean, run_type)) %>%
  pivot_longer(c(GW, Int, RO, mean), names_to = "flow", values_to = "weighted_r2") %>%
  mutate(facet.group = factor(run_type, levels = c("taxonomic", "trait"))) 

#reduce plot to only top 5 vars in each type
top.vars <- plot.df %>%
  filter(flow != "mean") %>%
  group_by(facet.group, flow) %>%
  slice_max(n = 5, order_by = weighted_r2) %>%
  ungroup() %>%
  select(env_fac) %>% distinct()
top.vars <- as.character(top.vars$env_fac)

p2 <- ggplot(data = plot.df[plot.df$env_fac %in% top.vars, ]) +
  geom_point(aes(x = weighted_r2, y = env_fac, shape = flow, fill = flow), size=4, alpha=0.7) +
  facet_wrap(~facet.group, scales = 'free', nrow = 2) +
  scale_fill_manual(values = c(flow.pal, "mean" = "black"), name = "Flow", breaks = c("mean", "GW", "Int", "RO")) +
  scale_shape_manual('Flow', values = c("GW" = 21, "Int" = 23, "RO" = 24, "mean" = 8), breaks = c("mean", "GW", "Int", "RO")) +
  scale_color_manual(values = var.high.pal) +
  scale_y_reordered() +
  labs(x = expression(paste(R^2, " weighted importance")), y = "") +
  # scale_x_continuous(limits = c(0, 0.02)) +
  # theme_classic() +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.box.margin=margin(0,0,0,0))

p2
ggsave(paste0(PATH, "/99_figures/var_imp_dotplot_", save.file.name, ".png"), plot = p2, width = 8, height = 6)

# df <- data.frame(hydro.pal) %>% rownames_to_column("type") %>% mutate(type = factor(type, levels = c("timing", "rate of change", "magnitude", "frequency", "duration",
#                                                                                                      "stream temperature", "summarized flow")))
# ggplot() +
#   geom_tile(data = df, aes(x = 0, y = type, fill = type), width = 1, height = 1) +
#   scale_fill_manual(values = hydro.pal) +
#   scale_x_continuous(expand = c(0,0)) +
#   coord_fixed(ratio = 1) +
#   theme_minimal() +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 21),
#         legend.position = "none",
#         panel.grid = element_blank())
# ggsave(paste0(PATH, "/99_figures/var_imp_dotplot_tile_legend_", save.file.name, ".png"), width = 6, height = 5)

### box plot ----
plot.df <- vi.table %>%
  mutate(flow = factor(flow, levels = c("RO", "GW", "Int")))
ggplot(data = plot.df) +
  geom_boxplot(aes(x = env_type, y = weighted_r2, fill = flow), outlier.shape = NA, show.legend = F) +
  geom_point(aes(x = env_type, y = weighted_r2, fill = flow), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4, size = 2) +
  facet_wrap(~ run_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = flow.pal) +
  labs(y = expression(paste(R^2, " weighted importance")), x = "") +
  guides(fill = guide_legend("Flow Regime", override.aes = list(shape = 22, alpha = 1, size = 5))) +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 10, color = "black"),
        # axis.text.x = element_text(angle = 20, hjust = 0.9, vjust = 1),
        axis.title = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.box.margin=margin(0,0,0,0))
ggsave(paste0(PATH, "/99_figures/var_imp_boxplot_", save.file.name, ".png"), width = 8, height = 8)

### density plot ----
env_vals <- lapply(gf.list, function(gf.out) {
  
  env <- gf.out$X
  # env <- env[, names(env) %in% top.vars]
  env <- env %>%
    pivot_longer(cols = everything(), names_to = "env_var", values_to = "value") %>%
    add_info_cols(., gf_model = gf.out)
  
  return(env)
  
})

env_vals <- bind_rows(env_vals)

if(var.type == "lulc") {
  env_vals <- add_nlcd_names(env_vals)
}

plot.df <- env_vals[env_vals$env_var %in% c("th1", "phase", "mn_water_temp_spring", "mn_water_temp_fall", "ma4", "ta3"), ]
# plot.df <- plot.df %>% 
#   filter((env_var == "pnTL1" & value < 10) | (env_var == "mh20" & value < 200) | grepl("wtemp", env_var) | (env_var == "ma11" & value > 0))
  # filter((env_var == "ma41" & value < 40) | (env_var == "mh20" & value < 200) | grepl("wtemp|ra5", env_var))

p3 <- ggplot() +
  geom_density(data = plot.df, aes(x = value, fill = flow), alpha = 0.6) +
  facet_wrap(~ env_var, scales = "free") +
  scale_fill_manual(values = flow.pal, name = "Flow") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "") +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 20, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.box.margin=margin(0,0,0,0))
p3

ggsave(paste0(PATH, "/99_figures/env_var_density_", save.file.name, ".png"), plot = p3, width = 15, height = 8)

#cumulative importance ----
##curves ----
if(var.type == "lulc") {
  vi.table <- read_csv(paste0(PATH, "/98_result_tables/variable_importance_R2_", save.file.name, ".csv"))
  plot.df <- vi.table %>%
    select(-type_name) %>%
    pivot_wider(names_from = "flow", values_from = "weighted_r2") %>%
    mutate(mean = rowMeans(select(., GW, Int, RO), na.rm = TRUE),
           env_fac = reorder_within(env_var, mean, run_type)) %>%
    pivot_longer(c(GW, Int, RO, mean), names_to = "flow", values_to = "weighted_r2") %>%
    mutate(facet.group = factor(run_type, levels = c("taxonomic", "trait"))) 
  
  #reduce plot to only top 5 vars in each type
  top.vars <- plot.df %>%
    filter(flow != "mean") %>%
    group_by(facet.group, flow) %>%
    slice_max(n = 5, order_by = weighted_r2) %>%
    ungroup() %>%
    select(env_fac) %>% distinct()
  top.vars <- as.character(top.vars$env_fac)
}

top.vars <- unique(sub("___.*", "", top.vars))
ci.cur <- lapply(gf.list, function(gf.out) { 
  
  #get ci values used for curves
  c_imp <- map(as.character(top.vars), function(v) {
    
    xy <- cumimp(gf.out, predictor = v, "Overall")
    ci <- data.frame(env_var = v, 
                     var_val = xy$x,
                     cum_imp = xy$y)
    
    return(ci)
    
  })
  
  c_imp <- bind_rows(c_imp)
  
  #add info columns for plotting
  c_imp <- add_info_cols(c_imp, gf.out)
  
  return(c_imp)
  
})

ci.cur <- bind_rows(ci.cur)
ci.cur <- add_nlcd_names(ci.cur)

##line plot ----
# plot.df <- ci.cur[ci.cur$env_var %in% tail(levels(top.vars), n = 3), ]
plot.df <- ci.cur[ci.cur$env_var %in% c("th1", "phase", "mn_water_temp_spring", "mn_water_temp_fall", "ma4", "ta3"), ]

p4 <- ggplot(data = plot.df, aes(x = var_val, y = cum_imp)) +
  geom_line(aes(color = flow, linetype = run_type)) +
  scale_color_manual(values = flow.pal, name = "Flow") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "") +
  facet_wrap(~env_var, scales = "free") +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ")"), x = "") +
  theme(
    panel.grid = element_line(color = "lightgray"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    plot.background = element_rect(fill = 'white'),
    axis.text = element_text(size=12),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    legend.box.margin=margin(0,0,0,0)
  )
p4

ggsave(paste0(PATH, "/99_figures/cumulimp_curve_", save.file.name, ".png"), plot = p4, width = 9, height = 6)

#ind. species/trait curves ----
##species performance ----
#NOTE: if regression, R2, if classification, 1-rel. error rate
sp.perf <- lapply(gf.list, function(gf.out) {
  
  sp_perf <- importance(gf.out, type = "Species")
  sp_perf <- as.data.frame(sp_perf) %>%
    rename(perf = sp_perf) %>%
    rownames_to_column("name") %>%
    mutate(rank = row_number()) %>%
    add_info_cols(., gf.out)
  
  return(sp_perf)
  
})

sp.perf <- bind_rows(sp.perf)

plot.df <- sp.perf %>%
  select(-rank) %>%
  pivot_wider(names_from = flow, values_from = perf) %>% 
  mutate(mean = rowMeans(select(., GW, Int, RO), na.rm = TRUE),
         name_frac = reorder_within(name, mean, run_type)) %>%
  pivot_longer(cols = c(GW, Int, RO, mean), names_to = "flow", values_to = "perf") %>%
  mutate(facet.group = factor(run_type, levels = c("taxonomic", "trait"))) %>%
  filter(!is.na(perf))

#reduce plot to only top 5 vars in each type
top.spp <- plot.df %>%
  filter(flow != "mean") %>%
  group_by(facet.group, flow) %>%
  slice_max(n = 5, order_by = perf) %>%
  ungroup() %>%
  select(name_frac) %>% distinct()
top.spp <- as.character(top.spp$name_frac)

ggplot(data = plot.df[plot.df$name_frac %in% top.spp, ]) +
  geom_point(aes(x = perf, y = name_frac, shape = flow, fill = flow), size=3.5, alpha=0.8) +
  facet_wrap(~facet.group, scales = 'free_y', nrow = 2) +
  scale_fill_manual(values = c(flow.pal, "mean" = "black"), name = "Flow", breaks = c("mean", "GW", "Int", "RO")) +
  scale_shape_manual('Flow', values = c("GW" = 21, "Int" = 23, "RO" = 24, "mean" = 8), breaks = c("mean", "GW", "Int", "RO")) +
  # scale_color_manual(values = var.high.pal) +
  scale_y_reordered() +
  labs(x = "Performance of random forests over species/trait", y = "") +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size=12),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.box.margin=margin(0,0,0,0))

ggsave(paste0(PATH, "/99_figures/species_performance_", save.file.name, ".png"), width = 7, height = 7)

##species cumulative importance ----
sp.cur <- lapply(gf.list, function(gf.out) { 
  
  #get sp values used for curves
  sp_imp <- map(as.character(top.vars), function(v) {
    
    xy <- cumimp(gf.out, predictor = v, "Species")
    # xy <- xy[names(xy) %in% top.spp]
    
    xy <- bind_rows(xy, .id = "name")
    
    names(xy) <- c("name", "var_val", "cumul_imp")
    
    xy$env_var <- v
    
    return(xy)
    
  })
  
  sp_imp <- bind_rows(sp_imp)
  
  #add info columns for plotting
  sp_imp <- add_info_cols(sp_imp, gf.out)
  
  return(sp_imp)
  
})

sp.cur <- bind_rows(sp.cur)

write.csv(sp.cur, paste0(PATH, "/98_result_tables/species_cum_imp_", save.file.name, ".csv"))

###plot ----
# plot.df <- sp.cur[sp.cur$env_var %in% tail(levels(top.vars), n = 3) & sp.cur$name %in% head(levels(top.spp), n = 10), ]
plot.df.trait <- sp.cur[sp.cur$env_var %in% c("th1", "phase", "mn_water_temp_spring", "mn_water_temp_fall", "ma4", "ta3") & sp.cur$run_type == "trait", ]
plot.df.tax <- sp.cur[sp.cur$env_var %in% c("th1", "phase", "mn_water_temp_spring", "mn_water_temp_fall", "ma4", "ta3") & sp.cur$run_type == "taxonomic", ]

ggplot(data = plot.df.trait, aes(x = var_val, y = cumul_imp)) +
  geom_line(aes(color = name)) +
  # scale_color_manual(values = flow.pal, name = "Flow Class") +
  # scale_linetype_manual(values = c("imp. variable" = "solid", "not" = "dashed"), name = "") +
  facet_grid(flow ~ env_var, scales = "free") +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ")"), x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", angle = 35, hjust = 0.9, vjust = 1.01),
    panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed")
  )

ggplot(data = plot.df.tax, aes(x = var_val, y = cumul_imp)) +
  geom_line(aes(color = name)) +
  # scale_color_manual(values = flow.pal, name = "Flow Class") +
  # scale_linetype_manual(values = c("imp. variable" = "solid", "not" = "dashed"), name = "") +
  facet_grid(flow ~ env_var, scales = "free") +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ")"), x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", angle = 35, hjust = 0.9, vjust = 1.01),
    panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed")
  )
  
# p4
# ggsave(paste0(PATH, "/99_figures/cumulimp_curve_", save.file.name, ".png"), plot = p4, width = 5, height = 5)

#species/trait threshold comparison ----
source(paste0(PATH, "/Scripts/XX_find_thresh_func.R"))

thresh <- lapply(gf.list, function(gf.out) {
  # gf.out <- gf.list[[2]]
  
  #get R2 vals (and variable names)
  imp <- importance(gf.out, type = "Weighted")
  imp <- data.frame(imp) %>% rownames_to_column("env_var") %>% rename(weighted_R2 = imp)
  
  #get threshold values for each variable (see Chen & Olden 2020)
  thresh <- map(imp$env_var, function(v) {
    # cat(v, "\n")
    # t <- tryCatch(get_var_threshold(gf.out, v), error = function(e) NA)
    t <- get_var_threshold(gf.out, v)
    df <- data.frame(env_var = v, 
                     thresh = t)
    
    return(df)
    
  })
  
  thresh <- bind_rows(thresh)
  
  thresh <- thresh %>%
    left_join(., imp) %>%
    add_info_cols(., gf.out)
  
  return(thresh)
  
})

thresh <- bind_rows(thresh)

write_csv(thresh, file = paste0(PATH, "/98_result_tables/community_threshold_val_", save.file.name, ".csv"))

##comp plot ----
thresh <- read_csv(paste0(PATH, "/98_result_tables/community_threshold_val_", save.file.name, ".csv"))

#get categories for grouping
cat.df <- lapply(thresh$env_var, var_cat_match)
cat.df <- bind_rows(cat.df) %>% select(-env_table_name) %>% distinct()

#get min and max values of env vars for scaling
sc.thresh <- lapply(gf.list, function(gf.out) {
  
  mod_info <- data.frame(x = "A") %>% add_info_cols(., gf.out) %>% select(-x)
  
  sc.t <- map(names(gf.out$X), function(v) {
    #get variable threshold
    t <- filter(thresh, env_var == v & flow == mod_info$flow & taxa == mod_info$taxa & run_type == mod_info$run_type)
    t <- t$thresh
    #get variable values + threshold
    var <- gf.out$X[names(gf.out$X) == v][, 1]
    var <- c(var, t)
    
    #scale + center
    sc.var <- scale(var)
    
    #find scaled threshold
    if(is.na(t)) { new.thresh <- NA } else { new.thresh <- sc.var[which(var == t)] }
    tmp <- data.frame(env_var = v,
                      sc.thresh = new.thresh) %>%
      add_info_cols(., gf.out)
    
    return(tmp)
  })
  
  sc.t <- bind_rows(sc.t)
  
  return(sc.t)
    
})

sc.thresh <- bind_rows(sc.thresh)

#prep table for plotting
thresh <- thresh %>%
  left_join(., sc.thresh) %>%
  left_join(., cat.df, by = c("env_var" = "input_name")) %>%
  mutate(env_type = case_when(grepl("spatial", data_category) ~ "spatial",
                              grepl("hydrology", data_category) ~ "hydrology",
                              grepl("lithology", data_category) ~ "lithology",
                              grepl("soil", data_category) ~ "soil",
                              grepl("climate", data_category) ~ "climate",
                              grepl("land cover|land use", data_category) ~ "land cover"))

if(var.type == "lulc") {
  thresh <- thresh %>%
    mutate(env_type = case_when(env_type == "spatial" ~ "spatial",
                                grepl("imp_|AGDRAIN|RDCRS|RDDENS|DAMDENS|_21|_22|_23|_24|_81|_82", env_var) ~ "land use",
                                T ~ "land cover"))
}

thresh <- add_nlcd_names(thresh)

if(all(thresh$env_type == "hydrology")) {
  thresh <- thresh %>%
    mutate(env_type = case_when(grepl("stream temperature", data_category) ~ "stream temperature",
                                grepl("frequency", data_category) ~ "frequency",
                                grepl("timing", data_category) ~ "timing",
                                grepl("magnitude", data_category) ~ "magnitude",
                                grepl("duration", data_category) ~ "duration",
                                grepl("rate of change", data_category) ~ "rate of change",
                                grepl("cumulative", data_category) ~ "summarized flow"))
}

#set color scheme
if(grepl("hydro", var.type)) { pal <- hydro.pal }
if(grepl("baseline|alt", var.type)) { pal <- var.high.pal }
if(grepl("lulc", var.type)) { pal <- lulc.pal }

#plot
plot.df <- thresh %>%
  select(-weighted_R2, -thresh) %>%
  # mutate(thresh = log(thresh)) %>%
  pivot_wider(names_from = run_type, values_from = sc.thresh) %>%
  mutate(resid = trait - taxonomic) 

if(grepl("hydro", var.type)) {
  plot.df <- plot.df %>%
    mutate(env_type = factor(env_type, levels = c("timing", "rate of change", "magnitude", "frequency", "duration",
                                                  "stream temperature", "summarized flow")))
}

most.diff <- plot.df %>%
  filter(!is.na(resid)) %>%
  arrange(resid) %>%
  slice(c(head(row_number(), 5), tail(row_number(), 5)))

theme_thresh_dot <- list(
  scale_shape_manual('Flow Regime', values = c("GW" = 21, "Int" = 23, "RO" = 24)),
  coord_fixed(ratio = 1),
  # scale_y_continuous(limits = c(-7.6, 7.6)),
  # scale_x_continuous(limits = c(-7.6, 7.6)),
  labs(x = "normalized taxonomic threshold value", y = "normalized trait threshold value"),
  theme(
    panel.grid = element_line(color = "lightgray"),
    axis.line = element_line(color = "black"),
    # panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = 'white'),
    axis.text = element_text(size=10),
    axis.title = element_text(size=10),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    legend.box.margin=margin(0,0,0,0)
    ) 
)

ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = env_type), type = "norm", linewidth = 2, alpha = 0.6) +
  geom_point(data = plot.df, aes(x = taxonomic, y = trait, fill = env_type, shape = flow), size = 3, alpha = 0.5) +
  geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = taxonomic, y = trait, color = env_type), se = F) +
  ggrepel::geom_text_repel(data = most.diff, aes(x = taxonomic, y = trait, label = env_var), min.segment.length = 0.2, force_pull = 0) +
  scale_fill_manual(values = pal, name = "Variable Type") +
  scale_color_manual(values = pal, name = "Variable Type") +
  guides(fill = guide_legend("Variable Type", override.aes = list(shape = 21, linetype = 0, alpha = 1, size = 5))) +
  theme_thresh_dot
ggsave(paste0(PATH, "/99_figures/threshold_comparison_point_byvartype_", save.file.name, ".png"), width = 10, height = 6)  

ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 2, alpha = 0.6) +
  geom_point(data = plot.df, aes(x = taxonomic, y = trait, fill = flow, shape = flow), size = 3, alpha = 0.5) +
  geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = taxonomic, y = trait, color = flow), se = F) +
  ggrepel::geom_text_repel(data = most.diff, aes(x = taxonomic, y = trait, label = env_var), min.segment.length = 0.2, force_pull = 0) +
  scale_fill_manual(values = flow.pal, name = "Flow Regime") +
  scale_color_manual(values = flow.pal, name = "Flow Regime") +
  guides(fill = guide_legend("Flow Regime", override.aes = list(linetype = 0, alpha = 1, size = 5))) +
  theme_thresh_dot
ggsave(paste0(PATH, "/99_figures/threshold_comparison_point_byflow_", save.file.name, ".png"), width = 10, height = 6)  

## residual comparison ----
plot.df <- plot.df %>%
  mutate(resid = trait - taxonomic)
y_lim <- max(abs(floor(min(plot.df$resid, na.rm = T))), ceiling(max(plot.df$resid, na.rm = T)))
sum.df <- plot.df %>%
  group_by(env_type, flow) %>%
  summarise(mn_resid = mean(resid, na.rm = T), sd_resid = sd(resid, na.rm = T))

ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(data = plot.df, aes(x = env_type, y = resid), position = position_jitter(width = 0.2), alpha = 0.4) +
  geom_point(data = sum.df, aes(x = env_type, y = mn_resid, color = flow), size = 5) +
  geom_linerange(data = sum.df, aes(x = env_type, ymin = mn_resid - sd_resid, ymax = mn_resid + sd_resid, color = flow), linewidth = 1) +
  facet_wrap(~ flow) +
  scale_color_manual(values = flow.pal, name = "Flow Regime") +
  # scale_color_manual(values = pal, name = "Variable Type") +
  scale_y_continuous(limits = c(y_lim*-1, y_lim)) +
  labs(x = "variable type", y = "residuals for 1:1 relationship") +
  theme_bw() +
  theme(legend.position = "none")


##flow comp ----
plot.df <- thresh %>%
  select(-weighted_R2, -thresh) %>%
  pivot_wider(names_from = flow, values_from = sc.thresh)

theme_flow_thresh <- list(
  scale_fill_manual(values = pal, name = "Variable Type"),
  scale_color_manual(values = pal, name = "Variable Type"),
  scale_linetype_manual('Run Type', values = c("solid", "longdash")), 
  coord_fixed(ratio = 1),
  guides(fill = guide_legend("Variable Type", override.aes = list(linetype = 0, alpha = 1, size = 5)),
         linetype = guide_legend("Assemblage Type", override.aes = list(linetype = c("solid", "dashed"), color = "black", linewidth = 1))),
  scale_x_continuous(limits = c(-3.5, 3)),
  scale_y_continuous(limits = c(-2.5, 3)),
  theme(
    panel.grid = element_line(color = "lightgray"),
    axis.line = element_line(color = "black"),
    # panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = 'white'),
    axis.text = element_text(size=10),
    axis.title = element_text(size=10),
    legend.key = element_rect(fill = "white"),
    legend.key.width = unit(2, "lines"),
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    legend.box.margin=margin(0,0,0,0), 
    legend.direction = "vertical"
  ) 
)

pA <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 1.5, alpha = 0.6) +
  geom_point(data = plot.df, aes(x = GW, y = Int, fill = env_type), shape = 21, size = 3, alpha = 0.3) +
  geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = GW, y = Int, color = env_type, linetype = run_type), linewidth = 1.5, se = F) +
  labs(x = bquote("normalized" ~ bold("groundwater") ~ "threshold value"), y = bquote("normalized" ~ bold("intermittent") ~ "threshold value")) +
  theme_flow_thresh

pB <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 1.5, alpha = 0.6) +
  geom_point(data = plot.df, aes(x = GW, y = RO, fill = env_type), shape = 21, size = 3, alpha = 0.3) +
  geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = GW, y = RO, color = env_type, linetype = run_type), linewidth = 1.5, se = F) +
  labs(x = bquote("normalized" ~ bold("groundwater") ~ "threshold value"), y = bquote("normalized" ~ bold("runoff") ~ "threshold value")) +
  theme_flow_thresh

pC <- ggplot() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 1.5, alpha = 0.6) +
  geom_point(data = plot.df, aes(x = Int, y = RO, fill = env_type), shape = 21, size = 3, alpha = 0.3) +
  geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = Int, y = RO, color = env_type, linetype = run_type), linewidth = 1.5, se = F) +
  labs(x = bquote("normalized" ~ bold("intermittent") ~ "threshold value"), y = bquote("normalized" ~ bold("runoff") ~ "threshold value")) +
  theme_flow_thresh

ggpubr::ggarrange(pA, pB, pC, nrow = 1, legend = "bottom", common.legend = TRUE)
ggsave(paste0(PATH, "/99_figures/threshold_flowregime_comparison_point_", save.file.name, ".png"), width = 15, height = 6, bg = "white")

##boxplot ----
ggplot() +
  geom_boxplot(data = thresh, aes(x = sc.thresh, y = env_type, fill = run_type)) +
  # geom_jitter(data = thresh, aes(x = sc.thresh, y = env_type, shape = run_type)) +
  # facet_wrap(~ flow) +
  scale_fill_manual(values = type.pal, name = "") +
  labs(y = "", x = "normalized variable threshold value") +
  theme(panel.grid.major = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size=12),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.box.margin=margin(0,0,0,0))
ggsave(paste0(PATH, "/99_figures/threshold_comparison_box_", save.file.name, ".png"), width = 14, height = 7)  

#dot plot
thresh %>%
  # select(-weighted_R2, -thresh) %>%
  # pivot_wider(names_from = "run_type", values_from = "sc.thresh") %>%
  # mutate(env_fac = reorder_within(env_var, trait, flow)) %>%
  # pivot_longer(cols = c(taxonomic, trait), names_to = "run_type", values_to = "sc.thresh") %>%
  ggplot() +
  geom_point(aes(x = sc.thresh, y = env_var, fill = run_type), shape = 21) +
  facet_wrap(~ flow) +
  scale_fill_manual(values = type.pal, name = "") +
  scale_y_reordered() +
  labs(y = "", x = "normalized variable threshold value") +
  theme(panel.grid.major = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size=12),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.box.margin=margin(0,0,0,0))

