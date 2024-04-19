#Master project
#Luca Eihlzer

# Select species ----
# Defining the available species
species <- c("betu_pub", "epil_fle", "lari_dec", "rhod_fer", "salix", "shrub", "tree")
mySpecies <- as.character("shrub")
if (length(mySpecies) == 0 || length(mySpecies) > 1){
  mySpecies <- as.character(sample(species, 1))
}
message(paste0(mySpecies, " is the selected species"))

# Ensure that folder "images" exist ----
if (!file.exists("images")) {
  # If it doesn't exist, create it
  dir.create("images")
  print("The 'images' folder has been created.")
} else {
  print("The 'images' folder already exists.")
}

# Load libraries ----
list.of.packages <- c(
  "gridExtra",
  "biomod2", 
  "blockCV", 
  "caret", 
  "dplyr", 
  "earth", 
  "ecospat", 
  "ggimage", 
  "ggplot2", 
  "ggtext", 
  "grImport2",
  "ggspatial", 
  "Hmisc", 
  "lattice", 
  "MASS", 
  "mda",
  "mgcv", 
  "pROC", 
  "randomForest", 
  "readr",
  "remotes", 
  "rpart", 
  "sf", 
  "sp", 
  "svglite",
  "terra", 
  "tidyterra",
  "tidyr", 
  "usdm"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
(.packages())

# Load data ----
setwd("C:/Users/Luca/Documents/_UNIL_/_MP/model")

# npres_df <- readRDS("data/npres_df.rds")
pa_data <- readRDS(paste0("data/",mySpecies,".rds"))
pa_data_gen <- readRDS(paste0("data/",mySpecies,"_gen.rds"))
expl_var <- readRDS("data/expl_var.rds")
expl_var_gen <- readRDS("data/expl_var_gen.rds")

# Scaling (Standardization)
ev_mean <- global(expl_var, "mean", na.rm=TRUE)
ev_sd <- global(expl_var, "sd", na.rm=TRUE)
expl_var <- terra::scale(expl_var)

# Scaling (Standardization)
expl_var_gen <- terra::scale(expl_var_gen, center=ev_mean[,1], scale=ev_sd[,1])

# Predictor selection ----
# Formatting Data
raster_values_tun <- terra::extract(expl_var, pa_data, df = TRUE, ID = FALSE)

biomod_data_tun <- BIOMOD_FormatingData(resp.var = pa_data$occ,
                                        expl.var = raster_values_tun,
                                        resp.xy = sf::st_coordinates(pa_data),
                                        resp.name = "occ",
                                        na.rm = TRUE)

options_tun <- BIOMOD_ModelingOptions(
  RF = list(do.classif = TRUE,
            ntree = 50,
            mtry = "default",
            sampsize = NULL,
            nodesize = 5,
            maxnodes = 15),
  
  GAM = list(algo = "GAM_gam")
)

set.seed(123)
RF_tun <- BIOMOD_Modeling(
  bm.format = biomod_data_tun,
  modeling.id = as.character(format(Sys.time(), "%y%m%d%H%M%S")),
  models = c("RF"),
  bm.options = options_tun,
  CV.strategy = "random",
  CV.nb.rep = 10,
  CV.perc = 0.8,
  metric.eval = c("TSS", "ROC"),
  var.import = 500,
  seed.val = 123)

var_imp <- as.data.frame(get_variables_importance(RF_tun))
ev_mod_name <- unique(var_imp$expl.var)

# Correlation matrix
df_cor <- data.frame(round(cor(raster_values_tun),3))

# Transform the correlation matrix into a distance matrix
cor.clust <- hclust(as.dist(1-abs(df_cor)))
# Plot those distances as a tree
png("images/clust_corr.png")
par(mfrow = c(1, 1), mar=c(1, 4, 2, 1))
plot(cor.clust, main="", ylab=expression(paste("Height as 1-|", rho, "|")))
# Add a red line at the threshold value (distance of 0.3)
abline(h=0.3, lty=2, col="red", lwd=2)
dev.off()

var_imp_mean <-  data.frame(matrix(nrow = length(ev_mod_name), ncol = 1))
rownames(var_imp_mean) <- ev_mod_name
colnames(var_imp_mean) <- c("var.imp.mean")

for (i in 1:length(df_cor)){
  var_imp_mean[i,] <- round(mean(var_imp$var.imp[which(var_imp$expl.var == ev_mod_name[i])]),3)
}

ev_sel <- data.frame(matrix(nrow = length(df_cor), ncol = 1))
colnames(ev_sel) <- c("expl.var")

for (i in 1:length(df_cor)) {
  
  ev_name <- rownames(df_cor)[which(abs(df_cor[,i]) > 0.7)]
  
  if (length(ev_name) == 1) {
    ev_sel[i,] <- ev_name
  } else {
    ev_comp <-  data.frame(matrix(nrow = length(ev_name), ncol = 1))
    colnames(ev_comp) <- c("var.imp.mean")
    rownames(ev_comp) <- ev_name
    ev_comp[,1] <- var_imp_mean[ev_name,1]
    ev_sel[i,] <- rownames(ev_comp)[which.max(ev_comp[,1])]
  }
}

ev_sel <- unique(ev_sel)

# repeat selection in case variables are still correlated
df_cor <- df_cor[ev_sel$expl.var,ev_sel$expl.var,drop=FALSE]
var_imp_mean <- var_imp_mean[ev_sel$expl.var,,drop=FALSE]

ev_sel <- data.frame(matrix(nrow = length(df_cor), ncol = 1))
colnames(ev_sel) <- c("expl.var")

for (i in 1:length(df_cor)) {
  
  ev_name <- rownames(df_cor)[which(abs(df_cor[,i]) > 0.7)]
  
  if (length(ev_name) == 1) {
    ev_sel[i,] <- ev_name
  } else {
    ev_comp <-  data.frame(matrix(nrow = length(ev_name), ncol = 1))
    colnames(ev_comp) <- c("var.imp.mean")
    rownames(ev_comp) <- ev_name
    ev_comp[,1] <- var_imp_mean[ev_name,1]
    ev_sel[i,] <- rownames(ev_comp)[which.max(ev_comp[,1])]
  }
}

ev_sel <- unique(ev_sel)

ev_sel$var.imp.mean <- var_imp_mean[c(ev_sel$expl.var),c("var.imp.mean")]
ev_sel <- ev_sel[order(ev_sel$var.imp.mean, decreasing = TRUE),]
species_npres <- sum(pa_data$occ == 1)
species_nabs <- sum(pa_data$occ == 0)

if (round(species_npres/10,0) <= 2 | round(species_nabs/10,0) <= 2){
  n_ev <- 2
} else {
  n_ev <- 3
}

# ev_sel <- ev_sel[ev_sel$var.imp.mean > 0.01,]

ev_sel <- ev_sel[1:n_ev,]
# ev_sel <- ev_sel[order(ev_sel$expl.var), ]

expl_var_mod <- expl_var[[ev_sel$expl.var]]

# ecospat.cor.plot(expl_var_mod)
# vif(expl_var_mod)

# extract the raster values for the species points as a dataframe
raster_values <- terra::extract(expl_var_mod, pa_data, df = TRUE, ID = FALSE)

head(raster_values)
str(raster_values)

# Spatial cross-validation ----

scv <- cv_spatial(
  pa_data,
  column = "occ",
  r = expl_var,
  k = 5,
  hexagon = FALSE,
  rows_cols = c(10,1),
  selection = "random",
  iteration = 999,
  biomod2 = TRUE,
  seed = 123,
  progress = TRUE,
  report = TRUE,
  plot = TRUE
)

# leave one out

ggsave("images/spat_cv.png")

# Defining the folds for data.split.table
spatial_cv_folds <- as.data.frame(scv$biomod_table)
colnames(spatial_cv_folds) <- c("_allData_RUN1", "_allData_RUN2", "_allData_RUN3", "_allData_RUN4", "_allData_RUN5")

# Model fitting ----

options <- BIOMOD_ModelingOptions(
  GLM = list( test = "none" ),
  
  RF = list(  ntree = 50,
              mtry = n_ev,
              sampsize = NULL,
              maxnodes = 15),
  
  GAM = list( algo = "GAM_gam",
              interaction.level = 1)
)

# Delete occ folder with previous models
unlink("occ", recursive = TRUE)

model_out <- BIOMOD_Modeling(
  bm.format = biomod_data,
  modeling.id = as.character(format(Sys.time(), "%y%m%d%H%M%S")),
  models = c("GLM", "GAM", "RF"),
  bm.options = options,
  CV.strategy = "user.defined",
  CV.user.table = spatial_cv_folds,
  var.import = 200,
  metric.eval = c("ROC", "TSS"),
  seed.val = 123
)

# Predictions ----

model_proj <- BIOMOD_Projection(
  bm.mod = model_out,
  proj.name = mySpecies,
  new.env = expl_var_mod,
  models.chosen = "all",
  metric.binary = "TSS",
  build.clamping.mask = FALSE,
  on_0_1000 = FALSE,
  seed.val = 123
)

model_proj_rast <- rast(paste0("occ/proj_",paste0(mySpecies),"/proj_",paste0(mySpecies),"_occ.tif"))

proj_names <- function(names){
  names <- paste0(names, "_")
  run_num <- substring(names, 16, 16)
  names <- substr(names, 18, nchar(names))
  names <- paste0(names,run_num)
}
names(model_proj_rast) <- proj_names(names(model_proj_rast))

model_proj_rast <- model_proj_rast[[sort(names(model_proj_rast))]]

plot.new()
ggplot() +
  geom_spatraster(data = model_proj_rast) +
  facet_wrap(~lyr, ncol = 5) +
  scale_fill_viridis_c(name = "Probability \nof presence",
                       limits = c(0, 1),
                       breaks = c(0,0.25,0.5,0.75,1)) +
  theme_void()
ggsave(paste0("images/",paste0(mySpecies),"_pred.png"), width = 4, height = 6)


ens_mod <- BIOMOD_EnsembleModeling(bm.mod = model_out,
                                   models.chosen = "all",
                                   em.by = "all",
                                   em.algo = c("EMmean"),
                                   metric.select = c("TSS"),
                                   metric.select.thresh = c(0.4),
                                   metric.eval = c("TSS", "ROC"),
                                   var.import = 50,
                                   EMwmean.decay = "proportional")

ens_mod_proj <- BIOMOD_EnsembleForecasting(bm.em = ens_mod,
                                           bm.proj = model_proj,
                                           models.chosen = "all",
                                           metric.binary = "all",
                                           metric.filter = "all",
                                           on_0_1000 = FALSE)

# load the prediction raster
ens_mod_proj_rast <- terra::as.data.frame(rast(paste0("occ/proj_",paste0(mySpecies),"/proj_",paste0(mySpecies),"_occ_ensemble.tif")), xy=TRUE)
names(ens_mod_proj_rast) <- c("x","y","pa")

# draw a svg north arrow
nor_arr <- '
  <svg viewBox="0 0 54 100 ">
   <line x1="52" y1="100" x2="52" y2="37" style="stroke:black;stroke-width:4" />
   <line x1="52" y1="100" x2="2" y2="37" style="stroke:black;stroke-width:4" />
   <line x1="2" y1="100" x2="2" y2="37" style="stroke:black;stroke-width:4" />
   <line x1="27" y1="100" x2="27" y2="0" style="stroke:black;stroke-width:4" />
   <line x1="42" y1="20" x2="27" y2="0" style="stroke:black;stroke-width:4" />
  </svg>
  '

plot.new()
ggplot() +
  geom_raster(data = ens_mod_proj_rast, aes(x = x, y = y,
                                            fill = pa)) +
  scale_fill_viridis_c(name = "Probability \nof presence",
                       limits = c(0, 1),
                       breaks = c(0,0.25,0.5,0.75,1)) +
  theme_void() +
  theme(legend.justification = c(0, 0.7)) +
  coord_equal() +
  geom_point_svg(
    mapping  = aes(x=2609250, y=1098800),
    svg      = nor_arr,
    size     = 4
  ) +
  annotation_scale(location = "br",
                   style = "ticks",
                   plot_unit = "m",
                   pad_x = unit(0, "mm"),
                   pad_y = unit(3, "mm"),) +
  annotate("text", x=2608700, y=1098600, label= "Glacier")

ggsave(paste0("images/",paste0(mySpecies),"_ens_mod.png"), width = 4, height = 6)

# Generalization ----
# extract the raster values for the species points as a dataframe

expl_var_gen_mod <- expl_var_gen[[ev_sel$expl.var]]

model_proj_gen <- BIOMOD_Projection(
  bm.mod = model_out,
  proj.name = paste0(paste0(mySpecies),"_gen"),
  new.env = expl_var_gen_mod,
  models.chosen = "all",
  build.clamping.mask = FALSE,
  seed.val = 123,
  on_0_1000 = FALSE
)

model_proj_gen_rast <- rast(paste0("occ/proj_",paste0(mySpecies),"_gen/proj_",paste0(mySpecies),"_gen_occ.tif"))

names(model_proj_gen_rast) <- proj_names(names(model_proj_gen_rast))
model_proj_gen_rast <- model_proj_gen_rast[[sort(names(model_proj_gen_rast))]]

plot.new()
ggplot() +
  geom_spatraster(data = model_proj_gen_rast) +
  facet_wrap(~lyr, ncol = 5) +
  scale_fill_viridis_c(name = "Probability \nof presence",
                       limits = c(0, 1),
                       breaks = c(0,0.25,0.5,0.75,1)) +
  theme_void()
ggsave(paste0("images/",paste0(mySpecies),"_pred_gen.png"), width = 4, height = 6)

ens_mod_proj_gen <- BIOMOD_EnsembleForecasting(bm.em = ens_mod,
                                               bm.proj = model_proj_gen,
                                               models.chosen = "all",
                                               metric.binary = "all",
                                               metric.filter = "all",
                                               on_0_1000 = FALSE)

ens_mod_proj_gen_rast <- rast(paste0("occ/proj_",paste0(mySpecies),"_gen","/proj_",paste0(mySpecies),"_gen_occ_ensemble.tif"))

# extract the raster values for the species points as a dataframe
pa_proj_gen <- terra::extract(ens_mod_proj_gen_rast, pa_data_gen, df = TRUE, ID = FALSE)
names(pa_proj_gen) <- c("proj_occ")

ens_mod_proj_gen_df <- terra::as.data.frame(ens_mod_proj_gen_rast, xy=TRUE)
names(ens_mod_proj_gen_df) <- c("x","y","pa")

plot.new()
ggplot() +
  geom_raster(data = ens_mod_proj_gen_df, aes(x = x, y = y,
                                            fill = pa)) +
  scale_fill_viridis_c(name = "Probability \nof presence",
                       limits = c(0, 1),
                       breaks = c(0,0.25,0.5,0.75,1)) +
  theme_void() +
  theme(legend.justification = c(0, 0.5)) +
  coord_equal() +
  geom_point_svg(
    mapping  = aes(x=2615600, y=1103600),
    svg      = nor_arr,
    size     = 4
  ) +
  annotation_scale(location = "br",
                   style = "ticks",
                   plot_unit = "m",
                   pad_x = unit(10, "mm"),
                   pad_y = unit(15, "mm"),) +
  annotate("text", x=2615050, y=1103200, label= "Glacier")

ggsave(paste0("images/",paste0(mySpecies),"_ens_mod_gen.png"), width = 4, height = 6)

# Evaluation ----

# Get evaluation scores
plot.new()
bm_PlotEvalMean(
  bm.out = model_out,
  dataset = "calibration",
  group.by = "algo",
  do.plot = TRUE,
  xlim = c(0,1),
  ylim = c(0,1))
ggsave(paste0("images/",paste0(mySpecies),"_eval_calib.png"), width = 6, height = 5)

bm_PlotEvalMean(
  bm.out = model_out,
  dataset = "validation",
  group.by = "algo",
  do.plot = TRUE,
  xlim = c(0,1),
  ylim = c(0,1))
ggsave(paste0("images/",paste0(mySpecies),"_eval_valid.png"), width = 6, height = 5)

# Response curves

resp_curv <- bm_PlotResponseCurves(
  bm.out = model_out,
  do.plot = FALSE)$tab

resp_curv$pred.name <- sub("occ_allData_", "", resp_curv$pred.name)
resp_curv$run <- as.integer(substr(resp_curv$pred.name, 4, 4))
resp_curv$pred.name <- substr(resp_curv$pred.name, 6, nchar(resp_curv$pred.name))

algo_list <- list("GAM", "GLM", "RF")

for (i in 1:n_ev) {
  for (j in 1:length(algo_list)){
    assign(paste0(algo_list[j],"_mean_",i), resp_curv %>%
             filter(pred.name == paste0(algo_list[j])) %>%
             group_by(id, expl.name) %>%
             summarise(across(where(is.numeric), mean)) %>%
             filter(expl.name == ev_sel$expl.var[i]) %>%
             select(id, expl.name, expl.val, pred.val) %>%
             rename(pred.val.mean = pred.val) %>% 
             mutate(algo = algo_list[j]))
    
    assign(paste0(algo_list[j],"_min_",i), resp_curv %>%
             filter(pred.name == paste0(algo_list[j]) & expl.name == ev_sel$expl.var[i]) %>%
             group_by(id, expl.name) %>%
             summarise(across(where(is.numeric), min)) %>%
             select(id, expl.name, expl.val, pred.val) %>%
             rename(pred.val.min = pred.val) %>% 
             mutate(algo = algo_list[j]))
    
    assign(paste0(algo_list[j],"_max_",i), resp_curv %>%
             filter(pred.name == paste0(algo_list[j]) & expl.name == ev_sel$expl.var[i]) %>%
             group_by(id, expl.name) %>%
             summarise(across(where(is.numeric), max)) %>%
             select(id, expl.name, expl.val, pred.val) %>%
             rename(pred.val.max = pred.val) %>% 
             mutate(algo = algo_list[j]))
    
  }
}

# Unscaling and building response curve data for plot

if (n_ev == 2) {
  resp_curve_1 <- rbind(cbind(GAM_mean_1, GAM_max_1["pred.val.max"], GAM_min_1["pred.val.min"]),
                        cbind(GLM_mean_1, GLM_max_1["pred.val.max"], GLM_min_1["pred.val.min"]),
                        cbind(RF_mean_1, RF_max_1["pred.val.max"], RF_min_1["pred.val.min"]))
  
  ev_1 <- as.character(unique(resp_curve_1$expl.name))
  resp_curve_1$expl.val <- resp_curve_1$expl.val * ev_sd[ev_1,] + ev_mean[ev_1,]
  
  resp_curve_2 <- rbind(cbind(GAM_mean_2, GAM_max_2["pred.val.max"], GAM_min_2["pred.val.min"]),
                        cbind(GLM_mean_2, GLM_max_2["pred.val.max"], GLM_min_2["pred.val.min"]),
                        cbind(RF_mean_2, RF_max_2["pred.val.max"], RF_min_2["pred.val.min"]))
  
  ev_2 <- as.character(unique(resp_curve_2$expl.name))
  resp_curve_2$expl.val <- resp_curve_2$expl.val * ev_sd[ev_2,] + ev_mean[ev_2,]
  
  resp_curve_df <- rbind(resp_curve_1, resp_curve_2)
  
  x_label_df <- data.frame(c(ev_1,ev_2),c(ev_mean[ev_1,],ev_mean[ev_2,]))
  x_label_df <- x_label_df[order(x_label_df[,1]), ]
  
} else if (n_ev == 3) {
  resp_curve_1 <- rbind(cbind(GAM_mean_1, GAM_max_1["pred.val.max"], GAM_min_1["pred.val.min"]),
                        cbind(GLM_mean_1, GLM_max_1["pred.val.max"], GLM_min_1["pred.val.min"]),
                        cbind(RF_mean_1, RF_max_1["pred.val.max"], RF_min_1["pred.val.min"]))
  
  ev_1 <- as.character(unique(resp_curve_1$expl.name))
  resp_curve_1$expl.val <- resp_curve_1$expl.val * ev_sd[ev_1,] + ev_mean[ev_1,]
  
  resp_curve_2 <- rbind(cbind(GAM_mean_2, GAM_max_2["pred.val.max"], GAM_min_2["pred.val.min"]),
                        cbind(GLM_mean_2, GLM_max_2["pred.val.max"], GLM_min_2["pred.val.min"]),
                        cbind(RF_mean_2, RF_max_2["pred.val.max"], RF_min_2["pred.val.min"]))
  
  ev_2 <- as.character(unique(resp_curve_2$expl.name))
  resp_curve_2$expl.val <- resp_curve_2$expl.val * ev_sd[ev_2,] + ev_mean[ev_2,]
  
  resp_curve_3 <- rbind(cbind(GAM_mean_3, GAM_max_3["pred.val.max"], GAM_min_3["pred.val.min"]),
                        cbind(GLM_mean_3, GLM_max_3["pred.val.max"], GLM_min_3["pred.val.min"]),
                        cbind(RF_mean_3, RF_max_3["pred.val.max"], RF_min_3["pred.val.min"]))
  
  ev_3 <- as.character(unique(resp_curve_3$expl.name))
  resp_curve_3$expl.val <- resp_curve_3$expl.val * ev_sd[ev_3,] + ev_mean[ev_3,]
  
  resp_curve_df <- rbind(resp_curve_1, resp_curve_2, resp_curve_3)
  
  x_label_df <- data.frame(c(ev_1,ev_2,ev_3),c(ev_mean[ev_1,],ev_mean[ev_2,],ev_mean[ev_3,]))
  x_label_df <- x_label_df[order(x_label_df[,1]), ]
}

resp_curve_df$algo <- as.character(resp_curve_df$algo)

# undo logarithmic transformation
rows_to_transform <- resp_curve_df$expl.name %in% c("log_vhm", "log_dist_water", "log_rough")
resp_curve_df$expl.val[rows_to_transform] <- exp(resp_curve_df$expl.val[rows_to_transform])
resp_curve_df$expl.name <- as.character(resp_curve_df$expl.name)
resp_curve_df$expl.name[rows_to_transform] <- gsub("log_", "", resp_curve_df$expl.name[rows_to_transform])

rows_to_transform <- x_label_df[,1] %in% c("log_vhm", "log_dist_water", "log_rough")
x_label_df[rows_to_transform,2] <- exp(x_label_df[rows_to_transform,2])

x_label <- x_label_df[,2]

# Variable importance

ev_imp <- get_variables_importance(model_out)
ev_imp <- ev_imp[c("algo","expl.var","var.imp")]

text_ev_imp <- ev_imp %>%
  group_by(algo, expl.var) %>%
  summarise(mean_var_imp = mean(var.imp)) %>%
  mutate(label = paste0("Imp: ", round(mean_var_imp, 2))) %>%
  rename(expl.name = expl.var)

text_ev_imp$y <- rep(c(0.3, 0.27, 0.24), each=n_ev)
text_ev_imp$x <- rep(x_label, times=3)
text_ev_imp <- as.data.frame(text_ev_imp)

rows_to_transform <- text_ev_imp$expl.name %in% c("log_vhm", "log_dist_water", "log_rough")
text_ev_imp$expl.name[rows_to_transform] <- gsub("log_", "", text_ev_imp$expl.name[rows_to_transform])

ggplot(resp_curve_df) +
  geom_ribbon(aes(x = expl.val, ymin = pred.val.min, ymax = pred.val.max, fill = algo), alpha = 0.3) +
  geom_path(aes(x = expl.val, y = pred.val.mean, color = algo)) +
  geom_text(data = text_ev_imp, aes(label = label, x = x, y = y, color = factor(algo)),
            hjust=0, size=3, show.legend=FALSE) +
  coord_cartesian(ylim=c(0,1), clip="off") +
  facet_wrap(~ expl.name, scales = "free_x", ncol = 3) +
  labs(x = "", y = "Probability of presence")

ggsave(paste0("images/",paste0(mySpecies),"_respc.png"), width = 10, height = 5)

# ROC curve for the generalisation

roc <- roc(pa_data_gen$occ, pa_proj_gen[,1])
roc_df <- data.frame(roc$specificities, roc$sensitivities)
names(roc_df) <- c("spec", "sens")

ggplot() +
  geom_path(data = roc_df, mapping = aes(x = spec, y = sens)) +
  scale_x_reverse() +
  labs(x = "Specificity", y = "Sensitivity") +
  annotate("text", x = .5, y = .5, label = paste0('AUC: ', round(roc$auc,2)), size = 5)
ggsave(ggsave(paste0("images/",paste0(mySpecies),"_roc_gen.png"), width = 6, height = 6))

# Results tabs ----

PRE <- round(species_npres/120, 2)

# For the single models 

metrics <- data.frame(get_evaluations(model_out)) %>%
  group_by(algo, metric.eval) %>%
  summarise(across(where(is.numeric), ~ round(mean(.), 2))) %>%
  select(algo, metric.eval, validation) %>%
  arrange(metric.eval, algo)

ev_imp_mean <- ev_imp %>%
  group_by(algo, expl.var) %>%
  summarise(across(where(is.numeric), ~ round(mean(.), 2))) %>%
  arrange(expl.var,algo)

rows_to_transform <- ev_imp_mean$expl.var %in% c("log_vhm", "log_dist_water", "log_rough")
ev_imp_mean$expl.var <- as.character(ev_imp_mean$expl.var)
ev_imp_mean$expl.var[rows_to_transform] <- gsub("log_", "", ev_imp_mean$expl.var[rows_to_transform])

columns <- species
rows <- c("PRE", "AUC.GAM", "AUC.GLM", "AUC.RF", "TSS.GAM", "TSS.GLM", "TSS.RF", "var1", "var1.imp.GAM", "var1.imp.GLM", "var1.imp.RF", "var2", "var2.imp.GAM", "var2.imp.GLM", "var2.imp.RF", "var3", "var.imp.GAM", "var3.imp.GLM", "var3.imp.RF")

if (file.exists("data/results_models.csv")) {
  results_models <- read.csv("data/results_models.csv", stringsAsFactors = FALSE)
  results_models <- results_models[,-1]
  rownames(results_models) <- rows
} else {
  results_models <- data.frame(matrix(nrow = length(rows), ncol = length(columns)))
  rownames(results_models) <- rows
  colnames(results_models) <- columns
}

if(n_ev == 2){
  results_models_mySpecies <- c(PRE, metrics$validation, paste0(ev_imp_mean$expl.var[1]), ev_imp_mean$var.imp[1:3], paste0(ev_imp_mean$expl.var[4]), ev_imp_mean$var.imp[4:6], rep(NA, 4))
} else if (n_ev == 3){
  results_models_mySpecies <- c(PRE, metrics$validation, paste0(ev_imp_mean$expl.var[1]), ev_imp_mean$var.imp[1:3], paste0(ev_imp_mean$expl.var[4]), ev_imp_mean$var.imp[4:6], paste0(ev_imp_mean$expl.var[7]), ev_imp_mean$var.imp[7:9])
}

results_models[,mySpecies] <- results_models_mySpecies
results_models <- as.data.frame(t(results_models))
results_models <- results_models[order(results_models$PRE, decreasing=TRUE),]
results_models <- as.data.frame(t(results_models))
write.csv(results_models, file = "data/results_models.csv")

# For the ensamble modelling¨

metrics <- data.frame(get_evaluations(ens_mod)) %>%
  select(metric.eval, calibration, sensitivity, specificity) %>%
  arrange(metric.eval) %>%
  mutate(calibration = round(calibration, 2))

SE <- round(mean(metrics$sensitivity)/100,2)
SP <- round(mean(metrics$specificity)/100,2)


AUC_gen <- round(as.numeric(gsub("[^0-9.]", "", roc$auc)),2)

conf_mat_gen <- caret::confusionMatrix(data = as.factor(ifelse(pa_proj_gen[,1] >= mean(get_evaluations(ens_mod)$cutoff)/1000, 1, 0)),
                reference = as.factor(pa_data_gen$occ))$table
TA <- conf_mat_gen[1,1]
FA <- conf_mat_gen[1,2]
FP <- conf_mat_gen[2,1]
TP <- conf_mat_gen[2,2]
SE_gen <- TP/(TP+FA)
SP_gen <- TA/(TA+FP)
TSS_gen <- round(SE_gen+SP_gen-1,2)

ev_imp <- get_variables_importance(ens_mod)
ev_imp <- ev_imp[c("algo","expl.var","var.imp")]

ev_imp_mean <- ev_imp %>%
  group_by(expl.var) %>%
  summarise(across(where(is.numeric), ~ round(mean(.), 2))) %>%
  arrange(expl.var) 

rows_to_transform <- ev_imp_mean$expl.var %in% c("log_vhm", "log_dist_water", "log_rough")
ev_imp_mean$expl.var <- as.character(ev_imp_mean$expl.var)
ev_imp_mean$expl.var[rows_to_transform] <- gsub("log_", "", ev_imp_mean$expl.var[rows_to_transform])
ev_imp_mean <-  ev_imp_mean[order(ev_imp_mean$var.imp, decreasing = TRUE),]

rows <- c("PRE", "AUC", "TSS", "SE", "SP", "AUC.gen", "TSS.gen", "SE.gen", "SP.gen","var1", "var1.imp", "var2", "var2.imp", "var3", "var3.imp")

if (file.exists("data/results_ensmod.csv")) {
  results_ensmod <- read.csv("data/results_ensmod.csv", stringsAsFactors = FALSE)
  results_ensmod <- results_ensmod[,-1]
  rownames(results_ensmod) <- rows
} else {
  results_ensmod <- data.frame(matrix(nrow = length(rows), ncol = length(columns)))
  rownames(results_ensmod) <- rows
  colnames(results_ensmod) <- columns
}

if(n_ev == 2){
  results_ensmod_mySpecies <- c(PRE, metrics$calibration, SE, SP, AUC_gen, TSS_gen, round(SE_gen,2), round(SP_gen,2), c(as.vector(t(ev_imp_mean))), rep(NA, 2))
} else if (n_ev == 3){
  results_ensmod_mySpecies <- c(PRE, metrics$calibration, SE, SP, AUC_gen, TSS_gen, round(SE_gen,2), round(SP_gen,2), c(as.vector(t(ev_imp_mean))))
}

results_ensmod[,mySpecies] <- results_ensmod_mySpecies
results_ensmod <- as.data.frame(t(results_ensmod))
results_ensmod <- results_ensmod[order(results_ensmod$PRE, decreasing=TRUE),]
results_ensmod <- as.data.frame(t(results_ensmod))
write.csv(results_ensmod, file = "data/results_ensmod.csv")

# Ending message ----
message("Congratulation, you modeled the distribution for ", paste0(mySpecies))

