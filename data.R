# Load libraries ----
list.of.packages <- c(
  "terra", 
  "spatialEco"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
(.packages())

# PREPROCESSING FERPECLE ----
# Load data ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load vegetation dataset (dependent variables)
spObs <- read.csv("data_row/allpoints.csv")
head(spObs) # overview of the data frame

# Load environmental dataset (independent variables)
dem <- rast("data_row/DEM.tif")
slope <- terrain(dem, v = "slope", unit = "radians", neighbors = 8)
aspect <- terrain(dem, v = "aspect", unit = "radians", neighbors = 8)
north_exp <- sin(slope)*cos(aspect) 
east_exp <- sin(slope)*sin(aspect)
rough <- terrain(dem, v = "roughness", unit = "radians", neighbors = 8)
# curv <- rast("data_row/curv.tif")
# glac_age_shp <- vect("data_row/glacier_age_mod.shp") #modified geometry for the 1850 limit because one plot was not included in it
# glac_age <- rasterize(glac_age_shp, dem, field = "Age")
vhm <- project(rast("data_row/VHM_ferp.tif"), y = "epsg:2056")
water_shp <- vect("data_row/rivers_lakes.shp")
water_shp <- project(water_shp, y = "epsg:2056")
water_shp <- as.lines(water_shp)
glacier_shp <- vect("data_row/glacier_2020.shp")
glacier_shp <- project(glacier_shp, y = "epsg:2056")
glacier_shp <- as.lines(glacier_shp)
sol_rad <- rast("data_row/solar_rad.tif")

# Load study area geometry
study_area <- vect("data_row/study_area.shp")

# Data processing for IVs ----

# Resampling
res_new <- 10 # new resolution
res_old <- res(dem)[1]
dem <- aggregate(dem, fact = res_new / res_old, fun = mean)
dem <- crop(dem, study_area)

# Euclidean distance rasters
dist_water <- distance(dem, water_shp)
dist_glac <- distance(dem, glacier_shp)

# Curvature rasters
curv_plan <- curvature(dem, type = c("planform"))
names(curv_plan) <- c("curv_plan")
curv_prof <- curvature(dem, type = c("profile"))
names(curv_prof) <- c("curv_prof")
curv_tot <- curvature(dem, type = c("total"))
names(curv_tot) <- c("curv_tot")

curv_tot[curv_tot == 0] <- min(abs(curv_tot[curv_tot != 0]))/10
curv_tot <- log(abs(curv_tot))
# min_lctot_pos <- min(log(curv_tot[curv_tot > 0]))
# min_lctot_neg <- min(log(abs(curv_tot[curv_tot < 0]))) # the closest cell value to zero is negative
# zdiff_lctot <- abs(min_lctot_pos - min_lctot_neg) # add to positive side
# curv_tot_mod <- terra::ifel(curv_tot > 0, -(log(curv_tot) + abs(min_lctot_pos) + zdiff_lctot),
#                               (log(abs(curv_tot)) + abs(min_lctot_neg))
#                               )

# feat_names <- c("dem", "slope", "north_exp", "east_exp", "rough", "sol_rad", "glac_age", "vhm", "dist_water", "dist_glac")
# feat_names <- c("dem", "slope", "north_exp", "east_exp", "rough", "sol_rad", "vhm", "dist_water", "dist_glac", "curv_plan", "curv_prof", "curv_tot", "log_curv_plan", "log_curv_prof", "log_curv_tot")
feat_names <- c("dem", "slope", "north_exp", "east_exp", "rough", "sol_rad", "vhm", "dist_water", "dist_glac", "curv_tot")

feat_list <- mget(as.character(unlist(feat_names)))
new_feat_list <- list(dem)
for (i in 2:length(feat_list)) {
  new_feat_list[[i]] <- resample(feat_list[[i]], dem, method='bilinear')
}

features <- rast(new_feat_list)
names(features) <- feat_names

# Plot features distribution
# svg("images/features_hist.svg")
par(mfrow = c(3, 5))
for (i in 1:nlyr(features)) {
  hist(features[[i]], main = NULL, xlab = names(features[[i]]))
}

# new features based on the distribution
vhm[vhm == 0] <- min(vhm[vhm !=0])/100
log_vhm <- log(vhm)
log_vhm <- resample(log_vhm, dem, method='bilinear')
names(log_vhm) <- c("log_vhm")
add(features) <- log_vhm

hist(log_vhm)

log_dist_water <- log(features$dist_water)
names(log_dist_water) <- c("log_dist_water")
add(features) <- log_dist_water

hist(log_dist_water)

# hist(rough[rough < 1])
rough[rough <= 0.1] <- 0.1
log_rough <- log(rough)
log_rough <- resample(log_rough, dem, method='bilinear')
names(log_rough) <- c("log_rough")
add(features) <- log_rough

hist(log_rough)
# dev.off()

# Aggregate the features and mask to study area
features <- terra::mask(features, study_area)

# Consistency in crs
project(features, y="epsg:2056")

saveRDS(features, file = paste0("data/features.rds"))

# Plot features (to be improved)
svg("images/features.svg")
plot(features) # visualizza la curv con un logaritmo
dev.off()

# Data processing for DVs ----
spObs[is.na(spObs)] <- 0
# Select the columns to make binary
species_names <- colnames(spObs[,2:28])

# Change the non-zero values to 1 in the selected numeric columns
spObs[species_names] <- ifelse(spObs[species_names] != 0, 1, 0)

saveRDS(spObs, file = "data/spObs.rds")

# Presence-absence data 
for (i in 1:length(species_names)) {
  # Create a new dataframe for each species
  assign("new_df", spObs[c(paste0(species_names[i]), "POINT_X", "POINT_Y", "point")])
  # call it as a new df
  # new_df <- get(paste0(species_names[i]))
  # shuffle
  new_df <- new_df[sample(nrow(new_df), nrow(new_df)), ]
  # Rename the columns of the new dataframes
  colnames(new_df) <- c("occ", "x", "y", "ID")
  
  # Convert to sf
  pa_data <- sf::st_as_sf(new_df, coords = c("x", "y"), crs = "epsg:2056")
  # create object
  saveRDS(pa_data, file = paste0("data/",species_names[i],".rds"))
}

# PREPROCESSING ANNIVIERS ----
# Load data ----
# Load vegetation dataset (dependent variables)
spObs_gen <- read.csv("data_row/allpoints_anni.csv")
head(spObs_gen) # overview of the data frame

# Load environmental dataset (independent variables)
dem_gen <- rast("data_row/DEM_anni.tif")
slope_gen <- terrain(dem_gen , v = "slope", unit = "radians", neighbors = 8)
aspect_gen <- terrain(dem_gen , v = "aspect", unit = "radians", neighbors = 8)
north_exp_gen <- sin(slope_gen)*cos(aspect_gen) 
east_exp_gen <- sin(slope_gen)*sin(aspect_gen)
rough_gen <- terrain(dem_gen , v = "roughness", unit = "radians", neighbors = 8)
# curv_gen <- rast("data_row/curv_anni.tif")
# glac_age_shp_gen <- vect("data_row/glacier_age_anni.shp")
# glac_age_gen <- rasterize(glac_age_shp_gen, dem_gen , field = "Age")
vhm_gen <- project(rast("data_row/VHM_anni.tif"), y = "epsg:2056")
water_shp_gen <- vect("data_row/rivers_lakes_anni.shp")
water_shp_gen <- project(water_shp_gen, y = "epsg:2056")
water_shp_gen <- as.lines(water_shp_gen)
glacier_shp_gen <- vect("data_row/glacier_anni_2020.shp")
glacier_shp_gen <- project(glacier_shp_gen, y = "epsg:2056")
glacier_shp_gen <- as.lines(glacier_shp_gen)
sol_rad_gen <- rast("data_row/solar_rad.tif")

# Load study area geometry
study_area_gen <- vect("data_row/study_area_anni.shp")

# Data processing for IVs ----

# Resampling
res_new <- 10 # new resolution
res_old <- res(dem_gen)[1]
dem_gen <- aggregate(dem_gen, fact = res_new / res_old, fun = mean)
dem_gen <- crop(dem_gen, study_area_gen)

# Euclidean distance rasters
dist_water_gen <- distance(dem_gen, water_shp_gen)
dist_glac_gen <- distance(dem_gen, glacier_shp_gen)

# Curvature rastes
curv_plan_gen <- log(curvature(dem_gen, type = c("planform")))
curv_prof_gen <- log(curvature(dem_gen, type = c("profile")))
curv_tot_gen <- log(curvature(dem_gen, type = c("total")))

# crei una funzione che faccia questa cosa in una riga!
feat_names_gen <- c("dem_gen", "slope_gen", "north_exp_gen", "east_exp_gen", "rough_gen", "sol_rad_gen", "vhm_gen", "dist_water_gen", "dist_glac_gen", "curv_plan_gen", "curv_prof_gen", "curv_tot_gen")
feat_list <- mget(as.character(unlist(feat_names_gen)))
new_feat_list <- list(dem_gen)
for (i in 2:length(feat_list)) {
  new_feat_list[[i]] <- resample(feat_list[[i]], dem_gen, method='bilinear')
}

features_gen <- rast(new_feat_list)
names(features_gen) <- feat_names

# Plot features distribution
# svg("images/features_hist.svg")
par(mfrow = c(3, 5))
for (i in 1:nlyr(features_gen)) {
  hist(features_gen[[i]], main = NULL, xlab = names(features[[i]]))
}

# new features based on the distribution
vhm_gen[vhm_gen == 0] <- min(vhm_gen[vhm_gen !=0])/100 # same as the ferpecle vhm!
log_vhm <- log(vhm_gen)
log_vhm <- resample(log_vhm, dem_gen, method='bilinear')
names(log_vhm) <- c("log_vhm")
add(features_gen) <- log_vhm

hist(log_vhm)

log_dist_water <- log(features_gen$dist_water)
names(log_dist_water) <- c("log_dist_water")
add(features_gen) <- log_dist_water

hist(log_dist_water)

# hist(rough_gen[rough_gen < 1])
rough_gen[rough_gen <= 0.1] <- 0.1
log_rough <- log(rough_gen)
log_rough <- resample(log_rough, dem_gen, method='bilinear')
names(log_rough) <- c("log_rough")
add(features_gen) <- log_rough

hist(log_rough)

# Aggregate the features and mask to study area
features_gen <- terra::mask(features_gen, study_area_gen)

# Consistency in crs
project(features_gen, y="epsg:2056")

saveRDS(features_gen, file = paste0("data/features_gen.rds"))

# Plot features (to be improved)
svg("images/features_gen.svg")
plot(features_gen) # visualizza la curv con un logaritmo
dev.off()

# Data processing for DVs ----
spObs_gen[is.na(spObs_gen)] <- 0
# Select the columns to make binary

# Change the non-zero values to 1 in the selected numeric columns
spObs_gen[species_names] <- ifelse(spObs_gen[species_names] != 0, 1, 0)

# Presence-absence data 
for (i in 1:length(species_names)) {
  # Create a new dataframe for each species
  assign("new_df", spObs_gen[c(paste0(species_names[i]), "POINT_X", "POINT_Y", "point")])
  # call it as a new df
  # new_df <- get(paste0(species_names[i]))
  # shuffle
  new_df <- new_df[sample(nrow(new_df), nrow(new_df)), ]
  # Rename the columns of the new dataframes
  colnames(new_df) <- c("occ", "x", "y", "ID")
  
  # Convert to sf
  pa_data <- sf::st_as_sf(new_df, coords = c("x", "y"), crs = "epsg:2056")
  # create object
  saveRDS(pa_data, file = paste0("data/",species_names[i],"_gen.rds"))
}

