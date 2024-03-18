#Master project
#Luca Eihlzer

# Ensure that folder "images" exist ----
if (!file.exists("images")) {
  # If it doesn't exist, create it
  dir.create("images")
  print("The 'images' folder has been created.")
} else {
  print("The 'images' folder already exists.")
}

# Select species ----
# Defining the available species
species <- c("epil_fle", "gram", "herb", "lari_dec", "moss", "rhod_fer", "salix", "shrub", "tree")
mySpecies <- as.character("epil_fle")
if (length(mySpecies) == 0 || length(mySpecies) > 1){
  mySpecies <- as.character(sample(species, 1))
}
message(paste0(mySpecies, " is the selected species"))

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
remotes::install_github("coolbutuseless/ggsvg")
if(!require("ggsvg")){
  remotes::install_github("coolbutuseless/ggsvg")
}
library("ggsvg")
(.packages())

# Load data ----
setwd("C:/Users/Luca/Documents/_UNIL_/_MP/model")

spObs <- readRDS("data/spObs.rds")
pa_data <- readRDS(paste0("data/",mySpecies,".rds"))
pa_data_gen <- readRDS(paste0("data/",mySpecies,"_gen.rds"))
features <- readRDS("data/features.rds")
features_gen <- readRDS("data/features_gen.rds")

# Scaling (Standardization)
feat_mean <- global(features, "mean", na.rm=TRUE)
feat_sd <- global(features, "sd", na.rm=TRUE)
features <- terra::scale(features)

# feat_min <- global(features, "min", na.rm=TRUE)
# feat_max <- global(features, "max", na.rm=TRUE)
# features <- (features - feat_min[,1])/(feat_max[,1]-feat_min[,1])
# 
# 
# 
# # Scaling (Standardization)
features_gen <- terra::scale(features_gen, center=feat_mean[,1], scale=feat_sd[,1])

# Features selection ----
# Formatting Data
raster_values_tun <- terra::extract(features, pa_data, df = TRUE, ID = FALSE)

# vif(raster_values_tun, size = nrow(raster_values_tun))
# 
# raster_values_tun <- subset(raster_values_tun, select = -c(slope, curv))


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

RF_tun <- BIOMOD_Modeling(
  bm.format = biomod_data_tun,
  modeling.id = as.character(format(Sys.time(), "%y%m%d%H%M%S")),
  models = c("RF"),
  bm.options = options_tun,
  CV.strategy = "random",
  CV.nb.rep = 5,
  CV.perc = 0.8,
  var.import = 500,
  seed.val = 123)

var_imp <- as.data.frame(get_variables_importance(RF_tun))
feat_mod_name <- unique(var_imp$expl.var)

# Correlation matrix
df_cor <- data.frame(round(cor(raster_values_tun),3))

# Transform the correlation matrix into a distance matrix
cor.clust <- hclust(as.dist(1-abs(df_cor)))
# Plot those distances as a tree
svg(paste0("images/",paste0(mySpecies),"_clust_corr.svg"))
par(mfrow = c(1, 1), mar=c(1, 4, 2, 1))
plot(cor.clust, main="Cluster of the correlations among variables",
     ylab=expression(paste("Height as 1-|", rho, "|")))
# Add a red line at the threshold value (distance of 0.3)
abline(h=0.3, lty=2, col="red", lwd=2)
dev.off()

var_imp_mean <-  data.frame(matrix(nrow = length(feat_mod_name), ncol = 1))

for (i in 1:length(feat_mod_name)){
  var_imp_mean[i,] <- round(mean(var_imp$var.imp[which(var_imp$expl.var == feat_mod_name[i])]),3)
}
rownames(var_imp_mean) <- feat_mod_name
colnames(var_imp_mean) <- c("var.imp.mean")

feat_sel <- data.frame(matrix(nrow = length(df_cor), ncol = 1))
colnames(feat_sel) <- c("expl.var")

for (i in 1:length(df_cor)) {
  
  feat_name <- rownames(df_cor)[which(df_cor[,i] < 1 & abs(df_cor[,i]) > 0.7)]
  
  if (length(feat_name) == 0) {
    feat_sel[i,] <- colnames(df_cor)[i]
  } else {
    feat_comp <-  data.frame(matrix(nrow = length(feat_name)+1, ncol = 1))
    colnames(feat_comp) <- c("var.imp.mean")
    rownames(feat_comp) <- c(feat_name, colnames(df_cor)[i])
    feat_comp[,1] <- var_imp_mean[rownames(feat_comp),1]
    feat_sel[i,] <- rownames(feat_comp)[which.max(feat_comp[,1])]
  }
}

feat_sel <- unique(feat_sel)
feat_sel$var.imp.mean <- var_imp_mean[c(feat_sel$expl.var),c("var.imp.mean")]
feat_sel <- feat_sel[order(feat_sel$var.imp.mean, decreasing = TRUE),]
species_npres <- sum(pa_data$occ == 1)
species_nabs <- sum(pa_data$occ == 0)

if (round(species_npres/10,0) <= 2 | round(species_nabs/10,0) <= 2){
  n_feat <- 2
} else {
  n_feat <- 3
}

# feat_sel <- feat_sel[feat_sel$var.imp.mean > 0.01,]

feat_sel <- feat_sel[1:n_feat,]

features_mod <- features[[feat_sel$expl.var]]


# ecospat.cor.plot(features_mod)
# vif(features_mod)

# extract the raster values for the species points as a dataframe
raster_values <- terra::extract(features_mod, pa_data, df = TRUE, ID = FALSE)

head(raster_values)
str(raster_values)


# Tuning ----

biomod_data <- BIOMOD_FormatingData(resp.var = pa_data$occ,
                                    expl.var = raster_values,
                                    resp.xy = sf::st_coordinates(pa_data),
                                    resp.name = "occ",
                                    na.rm = TRUE)
caret::trainControl(method="cv", 
             repeats = NA, 
             summaryFunction = "twoClassSummary", 
             classProbs = TRUE, 
             returnData = FALSE,
             number = 5)

bm.tuning <- BIOMOD_Tuning(
  biomod_data,
  bm.options = options_tun,
  models = c("GLM", "GAM", "RF"),
  metric.eval = "ROC",
  ctrl.train = caret::trainControl(method="cv", 
                            repeats = NA, 
                            summaryFunction = caret::twoClassSummary, 
                            classProbs = TRUE, 
                            returnData = FALSE,
                            number = 5),
  ctrl.train.tuneLength = 30,
  GLM.method = "glm", # keep it simple no stepwise AIC or BIC
  GLM.type = c("simple", "quadratic", "polynomial", "s_smoother"),
  GLM.interaction = c(0, 1),
  GAM.method = "gam",
  RF.method = "rf"
)

print(bm.tuning)

# Spatial cross-validation ----

scv <- cv_spatial(
  pa_data,
  column = "occ",
  r = features,
  k = 5,
  hexagon = FALSE,
  # flat_top = FALSE, #Creating hexagonal blocks with topped flat(?)
  # size = 650,
  rows_cols = c(10,1),
  selection = "random",
  iteration = 999,
  biomod2 = TRUE,
  seed = 1,
  progress = TRUE,
  report = TRUE,
  plot = TRUE
)

# leave one out

ggsave("images/spat_cv.svg")

# Defining the folds for data.split.table
spatial_cv_folds <- as.data.frame(scv$biomod_table)
colnames(spatial_cv_folds) <- c("_allData_RUN1", "_allData_RUN2", "_allData_RUN3", "_allData_RUN4", "_allData_RUN5")

# Model fitting ----
# Delete occ folder with previous models
unlink("occ", recursive = TRUE)

model_out <- BIOMOD_Modeling(
  bm.format = biomod_data,
  modeling.id = as.character(format(Sys.time(), "%y%m%d%H%M%S")),
  models = c("GLM", "GAM", "RF"),
  bm.options = bm.tuning$models.options,
  CV.strategy = "user.defined",
  CV.user.table = spatial_cv_folds,
  var.import = 100,
  metric.eval = c("KAPPA", "ROC"),
  seed.val = 123
)

# Cite maxent https://biodiversityinformatics.amnh.org/open_source/maxent/

# Get evaluation scores & variables importance
plot.new()
bm_PlotEvalMean(
  bm.out = model_out,
  dataset = "calibration",
  group.by = "algo",
  do.plot = TRUE,
  xlim = c(0,1),
  ylim = c(0,1))
ggsave(paste0("images/",paste0(mySpecies),"_eval_calib.svg"), width = 6, height = 5)

bm_PlotEvalMean(
  bm.out = model_out,
  dataset = "validation",
  group.by = "algo",
  do.plot = TRUE,
  xlim = c(0,1),
  ylim = c(0,1))
ggsave(paste0("images/",paste0(mySpecies),"_eval_valid.svg"), width = 6, height = 5)

resp_curv <- bm_PlotResponseCurves(
  bm.out = model_out,
  do.plot = FALSE)$tab

resp_curv$pred.name <- sub("occ_allData_", "", resp_curv$pred.name)
resp_curv$run <- as.integer(substr(resp_curv$pred.name, 4, 4))
resp_curv$pred.name <- substr(resp_curv$pred.name, 6, nchar(resp_curv$pred.name))

algo_list <- list("GAM", "GLM", "RF")

for (i in 1:n_feat) {
  for (j in 1:length(algo_list)){
    assign(paste0(algo_list[j],"_mean_",i), resp_curv %>%
      filter(pred.name == paste0(algo_list[j])) %>%
      group_by(id, expl.name) %>%
      summarise(across(where(is.numeric), mean)) %>%
      filter(expl.name == feat_sel$expl.var[i]) %>%
      select(id, expl.name, expl.val, pred.val) %>%
      rename(pred.val.mean = pred.val) %>% 
      mutate(algo = algo_list[j]))
    
    assign(paste0(algo_list[j],"_min_",i), resp_curv %>%
      filter(pred.name == paste0(algo_list[j]) & expl.name == feat_sel$expl.var[i]) %>%
      group_by(id, expl.name) %>%
      summarise(across(where(is.numeric), min)) %>%
        select(id, expl.name, expl.val, pred.val) %>%
        rename(pred.val.min = pred.val) %>% 
        mutate(algo = algo_list[j]))
    
    assign(paste0(algo_list[j],"_max_",i), resp_curv %>%
      filter(pred.name == paste0(algo_list[j]) & expl.name == feat_sel$expl.var[i]) %>%
      group_by(id, expl.name) %>%
      summarise(across(where(is.numeric), max)) %>%
        select(id, expl.name, expl.val, pred.val) %>%
        rename(pred.val.max = pred.val) %>% 
        mutate(algo = algo_list[j]))
    
  }
}

resp_curve_GAM <- rbind(cbind(GAM_mean_1, GAM_max_1["pred.val.max"], GAM_min_1["pred.val.min"]),
                        cbind(GAM_mean_2, GAM_max_2["pred.val.max"], GAM_min_2["pred.val.min"]),
                        cbind(GAM_mean_3, GAM_max_3["pred.val.max"], GAM_min_3["pred.val.min"]))

resp_curve_GLM <- rbind(cbind(GLM_mean_1, GLM_max_1["pred.val.max"], GLM_min_1["pred.val.min"]),
                        cbind(GLM_mean_2, GLM_max_2["pred.val.max"], GLM_min_2["pred.val.min"]),
                        cbind(GLM_mean_3, GLM_max_3["pred.val.max"], GLM_min_3["pred.val.min"]))

resp_curve_RF <- rbind(cbind(RF_mean_1, RF_max_1["pred.val.max"], RF_min_1["pred.val.min"]),
                        cbind(RF_mean_2, RF_max_2["pred.val.max"], RF_min_2["pred.val.min"]),
                        cbind(RF_mean_3, RF_max_3["pred.val.max"], RF_min_3["pred.val.min"]))

resp_curve_df <- rbind(resp_curve_GAM, resp_curve_GLM, resp_curve_RF)

ggplot(resp_curve_df) +
  geom_ribbon(aes(x = expl.val, ymin = pred.val.min, ymax = pred.val.max, fill = algo), alpha = 0.3) +
  geom_line(aes(x = expl.val, y = pred.val.mean, color = algo)) +
  facet_wrap(~ expl.name, scales = "free_x", ncol = 1) +
  labs(x = "Expl.val", y = "Pred.val") +
  theme_minimal()





# caret::unscale()

# resp_curve_1 <- ggplot() +
#   geom_ribbon(aes(x = GAM_mean_1$expl.val, ymin = GAM_min_1$pred.val, ymax = GAM_max_1$pred.val), fill = "#F8766D", alpha = 0.3) +
#   geom_ribbon(aes(x = GLM_mean_1$expl.val, ymin = GLM_min_1$pred.val, ymax = GLM_max_1$pred.val), fill = "#00BA38", alpha = 0.3) +
#   geom_ribbon(aes(x = RF_mean_1$expl.val, ymin = RF_min_1$pred.val, ymax = RF_max_1$pred.val), fill = "#619CFF", alpha = 0.3) +
#   geom_line(aes(x=GAM_mean_1$expl.val,y=GAM_mean_1$pred.val), color = "#F8766D", alpha = 0.7) +
#   geom_line(aes(x=GLM_mean_1$expl.val,y=GLM_mean_1$pred.val), color = "#00BA38", alpha = 0.7) +
#   geom_line(aes(x=RF_mean_1$expl.val,y=RF_mean_1$pred.val), color = "#619CFF", alpha = 0.7) +
#   labs(x = "Varable unit", y = "Probability of presence")
# 
# resp_curve_2 <- ggplot() +
#   geom_ribbon(aes(x = GAM_mean_2$expl.val, ymin = GAM_min_2$pred.val, ymax = GAM_max_2$pred.val), fill = "#F8766D", alpha = 0.3) +
#   geom_ribbon(aes(x = GLM_mean_2$expl.val, ymin = GLM_min_2$pred.val, ymax = GLM_max_2$pred.val), fill = "#00BA38", alpha = 0.3) +
#   geom_ribbon(aes(x = RF_mean_2$expl.val, ymin = RF_min_2$pred.val, ymax = RF_max_2$pred.val), fill = "#619CFF", alpha = 0.3) +
#   geom_line(aes(x=GAM_mean_2$expl.val,y=GAM_mean_2$pred.val), color = "#F8766D", alpha = 0.7) +
#   geom_line(aes(x=GLM_mean_2$expl.val,y=GLM_mean_2$pred.val), color = "#00BA38", alpha = 0.7) +
#   geom_line(aes(x=RF_mean_2$expl.val,y=RF_mean_2$pred.val), color = "#619CFF", alpha = 0.7) +
#   labs(x = "Varable unit", y = "Probability of presence")
# 
# if (n_feat==3){
#   resp_curve_3 <- ggplot() +
#     geom_ribbon(aes(x = GAM_mean_3$expl.val, ymin = GAM_min_3$pred.val, ymax = GAM_max_3$pred.val), fill = "#F8766D", alpha = 0.3) +
#     geom_ribbon(aes(x = GLM_mean_3$expl.val, ymin = GLM_min_3$pred.val, ymax = GLM_max_3$pred.val), fill = "#00BA38", alpha = 0.3) +
#     geom_ribbon(aes(x = RF_mean_3$expl.val, ymin = RF_min_3$pred.val, ymax = RF_max_3$pred.val), fill = "#619CFF", alpha = 0.3) +
#     geom_line(aes(x=GAM_mean_3$expl.val,y=GAM_mean_3$pred.val), color = "#F8766D", alpha = 0.7) +
#     geom_line(aes(x=GLM_mean_3$expl.val,y=GLM_mean_3$pred.val), color = "#00BA38", alpha = 0.7) +
#     geom_line(aes(x=RF_mean_3$expl.val,y=RF_mean_3$pred.val), color = "#619CFF", alpha = 0.7) +
#     labs(x = "Varable unit", y = "Probability of presence")
#   
#   figure <- ggarrange(resp_curve_1, resp_curve_2, resp_curve_3,
#                       labels = c(paste0(feat_sel$expl.var[1]), paste0(feat_sel$expl.var[2])),
#                       ncol = 2, nrow = 1)
#   figure
# } else {
#   figure <- ggarrange(resp_curve_1, resp_curve_2,
#                       labels = c(paste0(feat_sel$expl.var[1]), paste0(feat_sel$expl.var[2])),
#                       ncol = 2, nrow = 1)
#   figure
# }

# Plot utilizzando facet_wrap

print(plot)







bm_PlotResponseCurves(
  bm.out = model_out,
  models.chosen = get_built_models(model_out, algo = "GLM"),
  main = "Response curves for GLM")
ggsave(paste0("images/",paste0(mySpecies),"_respc_GLM.svg"), width = 10, height = 5)


bm_PlotResponseCurves(
  bm.out = model_out,
  models.chosen = get_built_models(model_out, algo = "GAM"),
  main = "Response curves for GAM")
ggsave(paste0("images/",paste0(mySpecies),"_respc_GAM.svg"), width = 10, height = 5)


bm_PlotResponseCurves(
  bm.out = model_out,
  models.chosen = get_built_models(model_out, algo = "RF"),
  main = "Response curves for RF")
ggsave(paste0("images/",paste0(mySpecies),"_respc_RF.svg"), width = 10, height = 5)


# Predictions ----

model_proj <- BIOMOD_Projection(
  bm.mod = model_out,
  proj.name = mySpecies,
  new.env = features_mod,
  models.chosen = "all",
  metric.binary = "ROC",
  build.clamping.mask = FALSE,
  on_0_1000 = FALSE,
  seed.val = 123
)

model_proj_rast <- rast(paste0("occ/proj_",paste0(mySpecies),"/proj_",paste0(mySpecies),"_occ.tif"))
names(model_proj_rast) <- c("GLM_1", "GAM_1", "RF_1", "GLM_2", "GAM_2", "RF_2", "GLM_3", "GAM_3", "RF_3", "GLM_4", "GAM_4", "RF_4", "GLM_5", "GAM_5", "RF_5")
model_proj_rast <- model_proj_rast[[sort(names(model_proj_rast))]]

plot.new()
ggplot() +
  geom_spatraster(data = model_proj_rast) +
  facet_wrap(~lyr, ncol = 5) +
  scale_fill_viridis_c(name = "Probability \nof presence",
                       limits = c(0, 1),
                       breaks = c(0,0.25,0.5,0.75,1)) +
  theme_void()
ggsave(paste0("images/",paste0(mySpecies),"_pred.svg"), width = 4, height = 6)


ens_mod <- BIOMOD_EnsembleModeling(bm.mod = model_out,
                                   models.chosen = "all",
                                   em.by = "all",
                                   em.algo = c("EMmean"),
                                   # metric.select = c("TSS"),
                                   # metric.select.thresh = c(0.7),
                                   metric.eval = c("ROC"),
                                   var.import = 3,
                                   EMci.alpha = 0.05,
                                   EMwmean.decay = "proportional")

get_evaluations(ens_mod)
get_variables_importance(ens_mod)

ens_mod_proj <- BIOMOD_EnsembleForecasting(bm.em = ens_mod, 
                                           bm.proj = model_proj,
                                           models.chosen = "occ_EMmeanByROC_mergedData_mergedRun_mergedAlgo",
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

ggsave(paste0("images/",paste0(mySpecies),"_ens_mod.svg"), width = 4, height = 6)

# Generalization ----

# # extract the raster values for the species points as a dataframe

features_gen_mod <- features_gen[[feat_sel$expl.var]]

model_proj_gen <- BIOMOD_Projection(
  bm.mod = model_out,
  proj.name = paste0(paste0(mySpecies),"_gen"),
  new.env = features_gen_mod,
  models.chosen = "all",
  build.clamping.mask = FALSE,
  seed.val = 123,
  on_0_1000 = FALSE
)


model_proj_gen_rast <- rast(paste0("occ/proj_",paste0(mySpecies),"_gen/proj_",paste0(mySpecies),"_gen_occ.tif"))
names(model_proj_gen_rast) <- c("GLM_1", "GAM_1", "RF_1", "GLM_2", "GAM_2", "RF_2", "GLM_3", "GAM_3", "RF_3", "GLM_4", "GAM_4", "RF_4", "GLM_5", "GAM_5", "RF_5")
model_proj_gen_rast <- model_proj_gen_rast[[sort(names(model_proj_gen_rast))]]

plot.new()
ggplot() +
  geom_spatraster(data = model_proj_gen_rast) +
  facet_wrap(~lyr, ncol = 5) +
  scale_fill_viridis_c(name = "Probability \nof presence",
                       limits = c(0, 1),
                       breaks = c(0,0.25,0.5,0.75,1)) +
  theme_void()
ggsave(paste0("images/",paste0(mySpecies),"_pred_gen.svg"), width = 4, height = 6)

ens_mod_proj_gen <- BIOMOD_EnsembleForecasting(bm.em = ens_mod,
                                               bm.proj = model_proj_gen,
                                               models.chosen = "occ_EMmeanByROC_mergedData_mergedRun_mergedAlgo",
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

ggsave(paste0("images/",paste0(mySpecies),"_ens_mod_gen.svg"), width = 4, height = 6)

roc <- roc(pa_data_gen$occ, pa_proj_gen[,1])
roc_df <- data.frame(roc$specificities, roc$sensitivities)
names(roc_df) <- c("spec", "sens")

ggplot() +
  geom_path(data = roc_df, mapping = aes(x = spec, y = sens)) +
  scale_x_reverse() +
  labs(x = "Specificity", y = "Sensitivity") +
  annotate("text", x = .5, y = .5, label = paste0('AUC: ', round(roc$auc,2)), size = 5)
ggsave(ggsave(paste0("images/",paste0(mySpecies),"_roc_gen.svg"), width = 6, height = 6))


# Figures ----
colnames(spObs)
bp_spObs <- (spObs[ , !names(spObs) %in% c("point" , "DateTimeS", "Elevation",
                                           "DateTime", "POINT_X", "POINT_Y", "POINT_Z")])
presences <- colSums(bp_spObs==1)
absences <- colSums(bp_spObs==0)

bp_names <- c("Trees", "Shrubs", "Graminoids", "Herbs", "Mosses", "Bare soil",
              "Picea abies", "Betula pubescens", "Larix decidua",
              "Populus tremula", "Alnus glutinosa", "Rhododendron ferrugineum",
              "Salix", "Sorbus aucuparia", "Calluna vulgaris", "Juniperus communis",
              "Empetrum nigrum", "Vaccinium myrtillus", "Vaccinium uliginosum",
              "Vaccinium vitis-idaea", "Dryas octopetala", "Epilobium angustifolium",
              "Epilobium fleischeri", "Saxifraga aizoides", "Saxifraga bryoides", 
              "Saxifraga paniculata", "Adenostyles alliariae")
bp_df <- data.frame(matrix(nrow = 27, ncol = 3))
colnames(bp_df) <- c("Species", "Presences", "Absences")

bp_df$Species <- bp_names
bp_df$Presences <- presences
bp_df$Absences <- absences

ggplot(data=bp_df, aes(x=reorder(Species, Presences), y=Presences)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  geom_text(aes(label = Presences), hjust = -0.2, colour = "black", size = 3) +
  xlab("Species") +
  coord_flip()
ggsave("images/species_npres.svg")

# Ending message ----
message("Congratulation, you modeled the distribution for ", paste0(mySpecies))