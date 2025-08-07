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
if(taxa == "bug") { taxa_label = "aquatic insect" } else { taxa_label = taxa }
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
  # slice_max(weighted_r2, n = 10) %>%
  mutate(al = factor(interaction(flow, env_var), levels = fl.lvl))

#plot
###bar chart ----
#set color scheme
if(grepl("hydro", var.type)) { pal <- hydro.pal }
if(grepl("baseline|alt", var.type)) { pal <- var.high.pal }
if(grepl("lulc", var.type)) { pal <- lulc.pal }

p1 <- ggplot(data = plot.df, aes(y = weighted_r2, x = al, fill = run_type)) +
  geom_col(color = "black", linewidth = 0.1, position = "dodge") +
  # facet_grid(vars(flow), vars(run_type), scales = "free") +
  facet_wrap(~ flow, scales = "free", ncol = 1) +
  # scale_y_continuous(expand = c(0.001,0), limits = c(0, max(vi.table$weighted_r2, na.rm = TRUE)+0.0005)) +
  scale_x_discrete(breaks = fl.lvl, labels = sub("^[^.]+\\.", "", fl.lvl), name = "") +
  scale_fill_manual(values = type.pal, name = "Assemblage Type") +
  labs(y = expression(paste(R^2, " weighted importance"))) +
  theme_bw() +
  ggtitle(label = taxa_label) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "bottom",
        panel.spacing.x = unit(1, "lines"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        strip.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
p1
ggsave(paste0(PATH, "/99_figures/var_imp_barchart_", save.file.name, ".png"), plot = p1, width = 12, height = 15)

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

p2 <- ggplot(data = plot.df) +
  geom_point(aes(y = weighted_r2, x = env_fac, shape = flow, fill = flow), size=4, alpha=0.7) +
  facet_wrap(~facet.group, scales = 'free', nrow = 2) +
  scale_fill_manual(values = c(flow.pal, "mean" = "black"), name = "Flow", breaks = c("mean", "GW", "Int", "RO")) +
  scale_shape_manual('Flow', values = c("GW" = 21, "Int" = 23, "RO" = 24, "mean" = 8), breaks = c("mean", "GW", "Int", "RO")) +
  # scale_color_manual(values = var.high.pal) +
  ggtitle(label = taxa_label) +
  scale_x_reordered() +
  labs(y = expression(paste(R^2, " weighted importance")), x = "") +
  # scale_x_continuous(limits = c(0, 0.02)) +
  # theme_classic() +
  theme(panel.grid.major.x = element_line(color="lightgrey"),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.box.margin=margin(0,0,0,0))

p2
ggsave(paste0(PATH, "/99_figures/var_imp_dotplot_", save.file.name, ".png"), plot = p2, width = 12, height = 8)

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
  mutate(flow = factor(flow, levels = c("RO", "GW", "Int")),
         env_type = gsub(" ", "\n", env_type))
max.r2 <- plot.df %>% group_by(flow, run_type, env_type) %>% slice_max(weighted_r2, n = 1) %>%
  mutate(y_coord = ifelse(weighted_r2 < 0.015, weighted_r2 + 0.001, weighted_r2 - 0.001),
         env_var = gsub("_water_temp_", "_", env_var)) %>%
  mutate(y_coord = ifelse(env_var == "mh13", y_coord + 0.001, y_coord))
ggplot(data = plot.df) +
  geom_boxplot(aes(x = env_type, y = weighted_r2, fill = flow), outlier.shape = NA, show.legend = F) +
  geom_point(aes(x = env_type, y = weighted_r2, fill = flow), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75, seed = 1), alpha = 0.4, size = 2) +
  geom_text(data = max.r2[max.r2$env_type != "stream\ntemperature",],  aes(x = env_type, y = y_coord, label = env_var,  shape = flow),
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75, seed = 1)) +
  geom_text(data = max.r2[max.r2$env_type == "stream\ntemperature",],  aes(x = env_type, y = y_coord, label = env_var,  shape = flow),
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75, seed = 1),
            angle = 330, hjust = 0.9) +
  facet_wrap(~ run_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = flow.pal) +
  labs(y = expression(paste(R^2, " weighted importance")), x = "") +
  guides(fill = guide_legend("Flow Regime", override.aes = list(shape = 22, alpha = 1, size = 5))) +
  ggtitle(taxa_label) +
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        # strip.text = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 10, color = "black"),
        # axis.text.x = element_text(angle = 20, hjust = 0.9, vjust = 1),
        axis.title = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.box.margin=margin(0,0,0,0))
ggsave(paste0(PATH, "/99_figures/var_imp_boxplot_", save.file.name, ".png"), width = 8, height = 8)

### summarize for text ----
vi.summ <- vi.table %>%
  group_by(flow, run_type) %>%
  mutate(rank = dense_rank(desc(weighted_r2))) %>%
  filter(rank <= 5)

vi.table %>%
  group_by(flow) %>%
  summarise(mn_r2 = mean(weighted_r2),
            med_r2 = median(weighted_r2),
            max_r2 = max(weighted_r2))

### density plot ----
# env_vals <- lapply(gf.list, function(gf.out) {
#   
#   env <- gf.out$X
#   # env <- env[, names(env) %in% top.vars]
#   env <- env %>%
#     pivot_longer(cols = everything(), names_to = "env_var", values_to = "value") %>%
#     add_info_cols(., gf_model = gf.out)
#   
#   return(env)
#   
# })
# 
# env_vals <- bind_rows(env_vals)
# 
# if(var.type == "lulc") {
#   env_vals <- add_nlcd_names(env_vals)
# }
# 
# plot.df <- env_vals[env_vals$env_var %in% c("th1", "phase", "mn_water_temp_spring", "mn_water_temp_fall", "ma4", "ta3"), ]
# # plot.df <- plot.df %>% 
# #   filter((env_var == "pnTL1" & value < 10) | (env_var == "mh20" & value < 200) | grepl("wtemp", env_var) | (env_var == "ma11" & value > 0))
#   # filter((env_var == "ma41" & value < 40) | (env_var == "mh20" & value < 200) | grepl("wtemp|ra5", env_var))
# 
# p3 <- ggplot() +
#   geom_density(data = plot.df, aes(x = value, fill = flow), alpha = 0.6) +
#   facet_wrap(~ env_var, scales = "free") +
#   scale_fill_manual(values = flow.pal, name = "Flow") +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "") +
#   theme(panel.grid.major = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA),
#         panel.background = element_rect(fill = "white"),
#         strip.background = element_rect(fill = "white", color = "black"),
#         strip.text = element_text(size = 20, face = "bold"),
#         plot.background = element_rect(fill = 'white'),
#         axis.text = element_text(size=20),
#         axis.title = element_text(size=20),
#         legend.key = element_rect(fill = "white"),
#         legend.title = element_text(size=20),
#         legend.text = element_text(size=20),
#         legend.box.margin=margin(0,0,0,0))
# p3
# 
# ggsave(paste0(PATH, "/99_figures/env_var_density_", save.file.name, ".png"), plot = p3, width = 15, height = 8)

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
ci.cur <- ci.cur %>%
  left_join(select(vi.summ, env_var, flow, run_type, rank))
plot.df <- ci.cur %>% filter(!is.na(rank))

p4 <- ggplot(data = plot.df[plot.df$run_type == "taxonomic",], aes(x = var_val, y = cum_imp)) +
  geom_line(aes(color = flow)) +
  scale_color_manual(values = flow.pal, name = "Flow") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "") +
  scale_y_continuous(n.breaks = 3) +
  scale_x_continuous(n.breaks = 4) +
  facet_wrap(flow~env_var, scales = "free") +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * "), taxonomic composition"), x = "") +
  ggtitle(taxa_label) +
  theme(
    panel.grid = element_line(color = "lightgray"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10),
    plot.background = element_rect(fill = 'white'),
    axis.text = element_text(size=12),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    legend.box.margin=margin(0,0,0,0),
    legend.position = "none"
  )
p4
ggsave(paste0(PATH, "/99_figures/cumulimp_curve_taxonomic_", save.file.name, ".png"), plot = p4, width = 11, height = 7)

p5 <- ggplot(data = plot.df[plot.df$run_type == "trait",], aes(x = var_val, y = cum_imp)) +
  geom_line(aes(color = flow)) +
  scale_color_manual(values = flow.pal, name = "Flow") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "") +
  scale_y_continuous(n.breaks = 3) +
  scale_x_continuous(n.breaks = 4) +
  facet_wrap(flow~env_var, scales = "free") +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * "), trait composition"), x = "") +
  ggtitle(taxa_label) +
  theme(
    panel.grid = element_line(color = "lightgray"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10),
    plot.background = element_rect(fill = 'white'),
    axis.text = element_text(size=12),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    legend.box.margin=margin(0,0,0,0),
    legend.position = "none"
  )
p5
ggsave(paste0(PATH, "/99_figures/cumulimp_curve_trait_", save.file.name, ".png"), plot = p5, width = 11, height = 7)

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

df1 <- plot.df %>%
  filter(run_type == "taxonomic" & flow != "mean") %>%
  mutate(name_flow = paste(name, flow, sep = "___")) %>%
  # compute order within each flow and turn into a factor with those levels
  group_by(flow) %>%
  mutate(name_flow = fct_reorder(name_flow, perf, .desc = FALSE),
         grp = if_else(perf <= 0.4, "low", "high")) %>%
  ungroup()

ggplot(data = df1, aes(x = perf, y = name_flow, fill = flow)) +
  # geom_point() +
  geom_col(color = "black", linewidth = 0.1) +
  # facet_grid(vars(flow), vars(run_type), scales = "free") +
  facet_wrap(grp ~ flow, scales = "free") +
  scale_x_continuous(breaks = c(0, 0.5, 1.0), limits = c(0,1)) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(values = flow.pal, name = "Assemblage Type") +
  labs(x = expression(R^2)) +
  theme_bw() +
  ggtitle(label = taxa_label) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none",
        panel.spacing.x = unit(1, "lines"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        strip.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# ggsave(paste0(PATH, "/99_figures/species_performance_barchart_taxo_", save.file.name, ".png"), width = 14, height = 16)

##too big, make table instead
summ.sp.perf <- sp.perf %>%
  mutate(perf = round(perf, digits = 3),
         name = gsub("_", " ", name)) %>%
  pivot_wider(names_from = flow, values_from = c(perf, rank)) %>%
  arrange(run_type, name) %>%
  select(-taxa) %>%
  relocate(run_type, name, perf_GW, rank_GW, perf_Int, rank_Int, perf_RO, rank_RO)
write_csv(summ.sp.perf, paste0(PATH, "/98_result_tables/taxa_performance_r2_", taxa, "_formatted.csv"))

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
sp.cur.plot <- sp.cur %>%
  left_join(select(vi.summ, env_var, flow, run_type, rank)) %>%
  filter(!is.na(rank) & !is.nan(cumul_imp) & !(name == "Morone_saxatilis" & flow == "GW" & run_type == "taxonomic"))

ggplot(data = sp.cur.plot[sp.cur.plot$run_type == "taxonomic",], aes(x = var_val, y = cumul_imp)) +
  geom_line(aes(color = name)) +
  # scale_color_manual(values = flow.pal, name = "Flow Class") +
  # scale_linetype_manual(values = c("imp. variable" = "solid", "not" = "dashed"), name = "") +
  facet_wrap(flow ~ env_var, scales = "free", nrow = 3) +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ") for taxonomic composition"), x = "") +
  scale_color_discrete(breaks = unique((sp.perf %>% filter(run_type == "taxonomic" & rank <= 5))$name),
                       name = "Best fit species\nfor each flow regime") +
  theme_bw() +
  ggtitle(taxa_label) +
  theme(
    # legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed")
  )
ggsave(paste0(PATH, "/99_figures/cumulimp_curve_spplevel_", taxa, "_taxonomic.png"), width = 13, height = 10)

ggplot(data = sp.cur.plot[sp.cur.plot$run_type == "trait",], aes(x = var_val, y = cumul_imp)) +
  geom_line(aes(color = name)) +
  # scale_color_manual(values = flow.pal, name = "Flow Class") +
  # scale_linetype_manual(values = c("imp. variable" = "solid", "not" = "dashed"), name = "") +
  facet_wrap(flow ~ env_var, scales = "free", nrow = 3) +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ") for trait composition"), x = "") +
  scale_color_discrete(breaks = unique((sp.perf %>% filter(run_type == "trait" & rank <= 5))$name),
                       name = "Best fit modality\nfor each flow regime") +
  theme_bw() +
  ggtitle(taxa_label) +
  theme(
    # legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", angle = 35, hjust = 0.9, vjust = 1.01),
    panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed")
  )
ggsave(paste0(PATH, "/99_figures/cumulimp_curve_spplevel_", taxa, "_trait.png"), width = 13, height = 10)
  
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
    sc.var <- (var - min(var)) / (max(var) - min(var))
    # sc.var <- scale(var)
    
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

write_csv(plot.df, paste0(PATH, "/98_result_tables/normalized_threshold_comparison_residuals_", taxa, ".csv"))

plot.df %>% 
  mutate(diff_zero = abs(resid) - 0) %>% 
  mutate(diff_zero_group = round(diff_zero, 2)) %>% 
  group_by(flow, diff_zero_group, env_type) %>% 
  summarise(num_occ = n()) %>% 
  View()

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
if(TRUE) {
 
  plot.df.r <- plot.df %>% mutate(
    env_type2 = if_else(env_type != "stream temperature" & !env_var %in% c("amplitude", "phase") & !grepl("tau|ar", env_var), TRUE, FALSE),
    index = row_number()) %>%
    mutate(env_type2 = case_when(env_type2 == TRUE & grepl("h", env_var) ~ "high flow",
                                 env_type2 == TRUE & grepl("l", env_var) ~ "low flow",
                                 env_type2 == TRUE & grepl("a", env_var) ~ "avg flow",
                                 env_type == "stream temperature" ~ "stream temp",
                                 T ~ "mag6"))
  most.diff <- plot.df.r %>%
    filter(abs(resid) > 0.45) %>%
    filter(!is.na(resid)) %>%
    group_by(flow, env_type) %>%
    arrange(resid) %>%
    slice(c(head(row_number(), 1), tail(row_number(), 1))) %>%
    distinct() %>%
    mutate(nudge_y = case_when(resid < 0 ~ -0.2, 
                               resid > 0 & resid < 0.8 ~ 0.2,
                               resid >= 0.8 ~ 0))
  sum.df <- plot.df.r %>%
    group_by(env_type, flow) %>%
    summarise(mn_resid = mean(resid, na.rm = T), sd_resid = sd(resid, na.rm = T))
  p1 <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_point(data = plot.df.r[!plot.df.r$index %in% most.diff$index,], aes(x = env_type, y = resid), size = 3, position = position_jitter(width = 0.2), alpha = 0.4) +
    geom_point(data = most.diff, aes(x = env_type, y = resid), size = 3, alpha = 0.4) +
    geom_linerange(data = sum.df, aes(x = env_type, ymin = mn_resid - sd_resid, ymax = mn_resid + sd_resid, color = env_type), linewidth = 1) +
    geom_point(data = sum.df, aes(x = env_type, y = mn_resid, color = env_type), shape = 18, size = 7) +
    ggrepel::geom_text_repel(
      data = most.diff,
      aes(x = env_type, y = resid, label = env_var),
      # max.overlaps = Inf,
      min.segment.length = 0,
      force = 5,                          # still helps spreading
      force_pull = 0,
      nudge_y = most.diff$nudge_y
    ) +
    facet_wrap(~ flow) +
    # scale_color_manual(values = flow.pal, name = "Flow Regime") +
    scale_color_manual(values = hydro.pal, name = "") +
    scale_y_continuous(limits = c(-1, 1), labels = c("higher\ntaxonomic\nthreshold", "", "equivalent\nthresholds", "", "higher\ntrait\nthreshold")) +
    labs(x = "hydrology category", y = "residuals for 1:1 taxonomic:trait threshold relationship") +
    theme_bw() +
    ggtitle(taxa_label) +
    theme(legend.position = "none",
          axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
          axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  
  most.diff <- plot.df.r %>%
    filter(abs(resid) > 0.45) %>%
    filter(!is.na(resid)) %>%
    group_by(flow) %>%
    arrange(resid) %>%
    slice(c(head(row_number(), 1), tail(row_number(), 1))) %>%
    distinct()
  sum.df <- plot.df.r %>%
    group_by(flow) %>%
    summarise(mn_resid = mean(resid, na.rm = T), sd_resid = sd(resid, na.rm = T))
  
  p2 <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_point(data = plot.df.r, aes(x = flow, y = resid), size = 3, position = position_jitter(width = 0.2), alpha = 0.4) +
    # geom_point(data = most.diff, aes(x = flow, y = resid), size = 3, alpha = 0.4) +
    geom_linerange(data = sum.df, aes(x = flow, ymin = mn_resid - sd_resid, ymax = mn_resid + sd_resid, color = flow), linewidth = 1) +
    geom_point(data = sum.df, aes(x = flow, y = mn_resid, color = flow), shape = 18, size = 7) +
    # ggrepel::geom_text_repel(data = most.diff, aes(x = flow, y = resid, label = env_var), min.segment.length = 0.2, force_pull = 0) +
    # facet_wrap(~ flow) +
    scale_color_manual(values = flow.pal, name = "Flow Regime") +
    scale_y_continuous(limits = c(-1, 1), labels = c("higher\ntaxonomic\nthreshold", "", "equivalent\nthresholds", "", "higher\ntrait\nthreshold")) +
    labs(x = "flow regime", y = "residuals for 1:1 taxonomic:trait threshold relationship") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(margin = margin(t = 1.8, unit = "cm")),
          plot.margin = unit(c(1.4, 0.3, 0.1, 0), "cm"))
  
  p3 <- ggpubr::ggarrange(p1, p2, nrow = 1, widths = c(3.5, 1), heights = c(2, 0.9))
  ggsave(paste0(PATH, "/99_figures/threshold_comparison_resid-only_", taxa, ".png"), plot = p3, width = 10, height = 6)
   
}

# ##flow comp ----
# plot.df <- thresh %>%
#   select(-weighted_R2, -thresh) %>%
#   pivot_wider(names_from = flow, values_from = sc.thresh)
# 
# theme_flow_thresh <- list(
#   scale_fill_manual(values = pal, name = "Variable Type"),
#   scale_color_manual(values = pal, name = "Variable Type"),
#   scale_linetype_manual('Run Type', values = c("solid", "longdash")), 
#   coord_fixed(ratio = 1),
#   guides(fill = guide_legend("Variable Type", override.aes = list(linetype = 0, alpha = 1, size = 5)),
#          linetype = guide_legend("Assemblage Type", override.aes = list(linetype = c("solid", "dashed"), color = "black", linewidth = 1))),
#   scale_x_continuous(limits = c(-3.5, 3)),
#   scale_y_continuous(limits = c(-2.5, 3)),
#   theme(
#     panel.grid = element_line(color = "lightgray"),
#     axis.line = element_line(color = "black"),
#     # panel.border = element_rect(color = "black", fill = NA),
#     panel.background = element_rect(fill = "white"),
#     plot.background = element_rect(fill = 'white'),
#     axis.text = element_text(size=10),
#     axis.title = element_text(size=10),
#     legend.key = element_rect(fill = "white"),
#     legend.key.width = unit(2, "lines"),
#     legend.title = element_text(size=10),
#     legend.text = element_text(size=10),
#     legend.box.margin=margin(0,0,0,0), 
#     legend.direction = "vertical"
#   ) 
# )
# 
# pA <- ggplot() +
#   geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
#   # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 1.5, alpha = 0.6) +
#   geom_point(data = plot.df, aes(x = GW, y = Int, fill = env_type), shape = 21, size = 3, alpha = 0.3) +
#   geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = GW, y = Int, color = env_type, linetype = run_type), linewidth = 1.5, se = F) +
#   labs(x = bquote("normalized" ~ bold("groundwater") ~ "threshold value"), y = bquote("normalized" ~ bold("intermittent") ~ "threshold value")) +
#   theme_flow_thresh
# 
# pB <- ggplot() +
#   geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
#   # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 1.5, alpha = 0.6) +
#   geom_point(data = plot.df, aes(x = GW, y = RO, fill = env_type), shape = 21, size = 3, alpha = 0.3) +
#   geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = GW, y = RO, color = env_type, linetype = run_type), linewidth = 1.5, se = F) +
#   labs(x = bquote("normalized" ~ bold("groundwater") ~ "threshold value"), y = bquote("normalized" ~ bold("runoff") ~ "threshold value")) +
#   theme_flow_thresh
# 
# pC <- ggplot() +
#   geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
#   # stat_ellipse(data = plot.df, aes(x = taxonomic, y = trait, color = flow), type = "norm", linewidth = 1.5, alpha = 0.6) +
#   geom_point(data = plot.df, aes(x = Int, y = RO, fill = env_type), shape = 21, size = 3, alpha = 0.3) +
#   geom_smooth(data = plot.df, method = "lm", formula='y~x', aes(x = Int, y = RO, color = env_type, linetype = run_type), linewidth = 1.5, se = F) +
#   labs(x = bquote("normalized" ~ bold("intermittent") ~ "threshold value"), y = bquote("normalized" ~ bold("runoff") ~ "threshold value")) +
#   theme_flow_thresh
# 
# ggpubr::ggarrange(pA, pB, pC, nrow = 1, legend = "bottom", common.legend = TRUE)
# ggsave(paste0(PATH, "/99_figures/threshold_flowregime_comparison_point_", save.file.name, ".png"), width = 15, height = 6, bg = "white")

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

# #dot plot
# thresh %>%
#   # select(-weighted_R2, -thresh) %>%
#   # pivot_wider(names_from = "run_type", values_from = "sc.thresh") %>%
#   # mutate(env_fac = reorder_within(env_var, trait, flow)) %>%
#   # pivot_longer(cols = c(taxonomic, trait), names_to = "run_type", values_to = "sc.thresh") %>%
#   ggplot() +
#   geom_point(aes(x = sc.thresh, y = env_var, fill = run_type), shape = 21) +
#   facet_wrap(~ flow) +
#   scale_fill_manual(values = type.pal, name = "") +
#   scale_y_reordered() +
#   labs(y = "", x = "normalized variable threshold value") +
#   theme(panel.grid.major = element_line(color="lightgrey"),
#         panel.grid.major.x = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA),
#         panel.background = element_rect(fill = "white"),
#         strip.background = element_rect(fill = "white", color = "black"),
#         strip.text = element_text(size = 12, face = "bold"),
#         plot.background = element_rect(fill = 'white'),
#         axis.text = element_text(size=12),
#         legend.position = "bottom",
#         legend.key = element_rect(fill = "white"),
#         legend.title = element_text(size=12),
#         legend.text = element_text(size=12),
#         legend.box.margin=margin(0,0,0,0))

# get direction of change for spp x traits ----
find_thresh_direction <- function(gf.mod, spp, v, 
                                  return.plot = FALSE, 
                                  plot.type = c("point", "box", "density", "bar proportion", "bar freq")) {
  if(FALSE) {
    ##testing
    gf.mod <- gf.list[[5]]
    v = "mn_water_temp_spring"
    # spp = "Luxilus_pilsbryi"
    spp = "troph_invert"
  }
  
  var <- gf.mod$X[, names(gf.mod$X) == v]
  spp.val <- gf.mod$Y[, names(gf.mod$Y) == spp]
  
  df <- data.frame(var_value = var, spp_val = spp.val)
  
  df <- add_info_cols(df, gf.mod)
  
  t <- thresh %>% filter(env_var == v & flow == unique(df$flow) & run_type == unique(df$run_type))
  df <- df %>% mutate(thresh_group = ifelse(var_value > t$thresh, "above", "below"))
  
  ##comparison test
  if(all(spp.val %in% c(1, 0))) {
    
    tbl <- table(df$thresh_group, spp.val)
    ct <- chisq.test(tbl)
    
    res <- data.frame(test_type = "chi sq.", 
                      test_stat = ct$statistic, 
                      p_val = ct$p.value, 
                      direction = names(which.max(tbl[, "1"])))
    
  } else {
    
    wt <- wilcox.test(spp.val ~ df$thresh_group)
    
    summ.df <- df %>%
      group_by(thresh_group) %>%
      summarise(mn_mod = mean(spp_val))
    
    res <- data.frame(test_type = "wilcoxon", 
                      test_stat = wt$statistic, 
                      p_val = wt$p.value, 
                      direction = summ.df$thresh_group[which.max(summ.df$mn_mod)])
    
  }
  
  res$spp <- spp
  res$var <- v
  
  if(return.plot) {
    if(!all(spp.val %in% c(1, 0))) {
      df <- df %>% mutate(transf_spp_val = spp_val * 100)
    }
    
    p <- ggplot(data = df)
    
    if(all(spp.val %in% c(1, 0))) {
      if(plot.type == "bar proportion") {
        p <- p + geom_bar(aes(x = spp_val, fill = thresh_group), color = "black", position = "fill") +
          labs(x = spp, y = "proportion of sites")
      }
      if(plot.type == "bar freq") {
        p <- p + geom_bar(aes(x = spp_val, fill = thresh_group), color = "black", position = "dodge") +
          labs(x = spp, y = "freq. of sites")
      }
      if(plot.type == "box") {
        p <- p + 
          geom_boxplot(aes(x = spp_val, y = var_value), outliers = FALSE) +
          geom_point(aes(x = spp_val, y = var_value), alpha = 0.6, position = position_jitter(width = 0.1)) +
          geom_hline(aes(yintercept = t$thresh), color = "red", linetype = "dashed") +
          labs(x = spp, y = v)
      }
    } else {
      if(plot.type == "point") {
        p <- p +
          geom_point(aes(x = var_value, y = spp_val)) +
          geom_vline(aes(xintercept = t$thresh), color = "red", linetype = "dashed")  +
          labs(x = v, y = spp)
      }
      if(plot.type == "box") {
        p <- p +
          geom_boxplot(aes(x = thresh_group, y = transf_spp_val)) +
          scale_y_continuous(transform = "log1p") +
          labs(x = "threshold distribution", y = paste0("proportion of taxa with ", spp))
      }
      if(plot.type == "density") {
        p <- p +
          geom_density(aes(x = transf_spp_val, fill = thresh_group), alpha = 0.5) +
          scale_x_continuous(transform = "log1p") +
          labs(x = paste0("proportion of taxa with ", spp))
      }
    }
    
    p <- p + theme_bw() +
      ggtitle(paste0("statistic = ", round(res$test_stat, 3), 
                     "; p-value = ", round(res$p_val, 3), 
                     "; direction of change is ", res$direction))
    
    return(p)
    
  } else {
    #if not plotting, return results of stats tests
    res <- add_info_cols(res, gf.mod)
    return(res)
  }
}

# find_thresh_direction(gf.list[[5]], spp = "troph_invert", v = "mn_water_temp_spring",
#                       return.plot = T, plot.type = "density")

t.dir <- lapply(gf.list, function(gf.mod) {
  
  spp.l <- names(gf.mod$Y)
  v.l <- names(gf.mod$X)
  
  res.all <- data.frame()
  for(i in seq_along(spp.l)) {
    
    for(y in seq_along(v.l)) {
      
      r.out <- find_thresh_direction(gf.mod, spp = spp.l[[i]], v = v.l[[y]])
      
      res.all <- do.call(rbind, list(res.all, r.out))
      
    }
  }
  
  return(res.all)
  
})

t.dir <- do.call(rbind, t.dir)

t.dir <- t.dir %>%
  mutate(across(where(is.numeric), ~round(.x, digits = 3)))

write_csv(t.dir, paste0(PATH, "/98_result_tables/threshold_change_direction_", taxa, ".csv"))
##need to account for variable and spp r2 after running this!!! summarize too!!

e.tmp <- vi.table %>% select(env_var, weighted_r2, flow, taxa, run_type, env_type) %>% rename(var = env_var)
t.dir.info <- t.dir %>%
  left_join(e.tmp) %>%
  left_join(rename(sp.perf, spp = name))

t.dir.summ <- t.dir.info %>%
  filter(var %in% max.r2$env_var)

t.dir.summ %>% 
  group_by(flow, run_type, var) %>%
  mutate(ttl_spp_incl = sum(!is.na(perf)),
         ttl_sig_diff = sum(p_val <= 0.05)) %>%
  group_by(flow, run_type, var, ttl_spp_incl, ttl_sig_diff, direction) %>%
  summarise(direction_ttl = sum(p_val <= 0.05)) %>% View()

t.dir.summ %>%
  filter(flow == "GW" 
         & var == "th1" 
         & run_type == "taxonomic" 
         # & p_val <= 0.05
         ) %>% View()
