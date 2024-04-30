#Master project
#Luca Eihlzer

# Load libraries ----
list.of.packages <- c(
  "ggplot2", 
  "terra"
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

ev_names <- c("dem", "slope", "north_exp", "east_exp", "rough", "sol_rad", "vhm", "dist_water", "dist_glac")

ev_list <- mget(as.character(unlist(ev_names)))
new_ev_list <- list(dem)
for (i in 2:length(ev_list)) {
  new_ev_list[[i]] <- resample(ev_list[[i]], dem, method='bilinear')
}

expl_var <- rast(new_ev_list)
names(expl_var) <- ev_names

# Plot expl_var distribution
par(mfrow = c(3, 5))
for (i in 1:nlyr(expl_var)) {
  hist(expl_var[[i]], main = NULL, xlab = names(expl_var[[i]]))
}

# logarithmic transformation

vhm[vhm == 0] <- min(vhm[vhm !=0])/100
log_vhm <- log(vhm)
log_vhm <- resample(log_vhm, dem, method='bilinear')
names(log_vhm) <- c("log_vhm")
add(expl_var) <- log_vhm

hist(log_vhm)

log_dist_water <- log(expl_var$dist_water)
names(log_dist_water) <- c("log_dist_water")
add(expl_var) <- log_dist_water

hist(log_dist_water)

# hist(rough[rough < 1])
rough[rough <= 0.1] <- 0.1
log_rough <- log(rough)
log_rough <- resample(log_rough, dem, method='bilinear')
names(log_rough) <- c("log_rough")
add(expl_var) <- log_rough

hist(log_rough)

# Remove layers
expl_var <- subset(expl_var, c("vhm", "dist_water", "rough"), negate=TRUE)

# Aggregate the expl_var and mask to study area
expl_var <- terra::mask(expl_var, study_area)

# Consistency in crs
project(expl_var, y="epsg:2056")

saveRDS(expl_var, file = paste0("data/expl_var.rds"))

# Plot expl_var (to be improved)
png("images/expl_var.png")
plot(expl_var) # visualizza la curv con un logaritmo
dev.off()

png("images/expl_var_hist.png")
par(mfrow = c(3, 4))
for (i in 1:nlyr(expl_var)) {
  hist(expl_var[[i]], main = NULL, xlab = names(expl_var[[i]]))
}
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
  
  # shuffle
  new_df <- new_df[sample(nrow(new_df), nrow(new_df)), ]
  
  # Rename the columns of the new dataframes
  colnames(new_df) <- c("occ", "x", "y", "ID")
  
  # Convert to sf
  pa_data <- sf::st_as_sf(new_df, coords = c("x", "y"), crs = "epsg:2056")
  
  # create object
  saveRDS(pa_data, file = paste0("data/",species_names[i],".rds"))
}

# Figure n.pres ----
colnames(spObs)
npres_spObs <- (spObs[ , !names(spObs) %in% c("point" , "bare_soil", "DateTimeS", "Elevation",
                                           "DateTime", "POINT_X", "POINT_Y", "POINT_Z")])
presences <- colSums(npres_spObs==1)
absences <- colSums(npres_spObs==0)

npres_names <- c("Trees", "Shrubs", "Graminoids", "Herbs", "Mosses",
              "Picea abies", "Betula pubescens", "Larix decidua",
              "Populus tremula", "Alnus glutinosa", "Rhododendron ferrugineum",
              "Salix", "Sorbus aucuparia", "Calluna vulgaris", "Juniperus communis",
              "Empetrum nigrum", "Vaccinium myrtillus", "Vaccinium uliginosum",
              "Vaccinium vitis-idaea", "Dryas octopetala", "Epilobium angustifolium",
              "Epilobium fleischeri", "Saxifraga aizoides", "Saxifraga bryoides", 
              "Saxifraga paniculata", "Adenostyles alliariae")
npres_df <- data.frame(matrix(nrow = 26, ncol = 3))
colnames(npres_df) <- c("Species", "Presences", "Absences")

npres_df$Species <- npres_names
npres_df$Presences <- presences
npres_df$Absences <- absences

ggplot(data=npres_df, aes(x=reorder(Species, Presences), y=Presences)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  geom_text(aes(label = Presences), hjust = -0.2, colour = "black", size = 3) +
  xlab("Species") +
  coord_flip()
ggsave("images/species_npres.png", width = 15, height =10, unit = c("cm"))

npres_df$Species <- colnames(npres_spObs)
saveRDS(npres_df, file = "data/npres_df.rds")

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

# crei una funzione che faccia questa cosa in una riga!
ev_names_gen <- paste0(ev_names, "_gen")
ev_list <- mget(as.character(unlist(ev_names_gen)))
new_ev_list <- list(dem_gen)
for (i in 2:length(ev_list)) {
  new_ev_list[[i]] <- resample(ev_list[[i]], dem_gen, method='bilinear')
}

expl_var_gen <- rast(new_ev_list)
names(expl_var_gen) <- ev_names

# Plot expl_var distribution
par(mfrow = c(3, 5))
for (i in 1:nlyr(expl_var_gen)) {
  hist(expl_var_gen[[i]], main = NULL, xlab = names(expl_var[[i]]))
}

# new expl_var based on the distribution
vhm_gen[vhm_gen == 0] <- min(vhm_gen[vhm_gen !=0])/100 # same as the ferpecle vhm!
log_vhm <- log(vhm_gen)
log_vhm <- resample(log_vhm, dem_gen, method='bilinear')
names(log_vhm) <- c("log_vhm")
add(expl_var_gen) <- log_vhm

hist(log_vhm)

log_dist_water <- log(expl_var_gen$dist_water)
names(log_dist_water) <- c("log_dist_water")
add(expl_var_gen) <- log_dist_water

hist(log_dist_water)

# hist(rough_gen[rough_gen < 1])
rough_gen[rough_gen <= 0.1] <- 0.1
log_rough <- log(rough_gen)
log_rough <- resample(log_rough, dem_gen, method='bilinear')
names(log_rough) <- c("log_rough")
add(expl_var_gen) <- log_rough

hist(log_rough)

# Remove layers
expl_var_gen <- subset(expl_var_gen, c("vhm", "dist_water", "rough"), negate=TRUE)

# Aggregate the expl_var and mask to study area
expl_var_gen <- terra::mask(expl_var_gen, study_area_gen)

# Consistency in crs
project(expl_var_gen, y="epsg:2056")

saveRDS(expl_var_gen, file = paste0("data/expl_var_gen.rds"))

# Plot expl_var (to be improved)
png("images/expl_var_gen.png")
plot(expl_var_gen) # visualizza la curv con un logaritmo
dev.off()

png("images/expl_var_gen_hist.png")
par(mfrow = c(3, 5))
for (i in 1:nlyr(expl_var_gen)) {
  hist(expl_var_gen[[i]], main = NULL, xlab = names(expl_var[[i]]))
}
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

  # shuffle
  new_df <- new_df[sample(nrow(new_df), nrow(new_df)), ]
  # Rename the columns of the new dataframes
  colnames(new_df) <- c("occ", "x", "y", "ID")
  
  # Convert to sf
  pa_data <- sf::st_as_sf(new_df, coords = c("x", "y"), crs = "epsg:2056")
  
  # create object
  saveRDS(pa_data, file = paste0("data/",species_names[i],"_gen.rds"))
}


# Comparing study areas ----

ev_min <- global(expl_var, "min", na.rm=TRUE)
ev_max <- global(expl_var, "max", na.rm=TRUE)
ev_mean <- global(expl_var, "mean", na.rm=TRUE)
ev_sd <- global(expl_var, "sd", na.rm=TRUE)

ev_min_gen <- global(expl_var_gen, "min", na.rm=TRUE)
ev_max_gen <- global(expl_var_gen, "max", na.rm=TRUE)
ev_mean_gen <- global(expl_var_gen, "mean", na.rm=TRUE)
ev_sd_gen <- global(expl_var_gen, "sd", na.rm=TRUE)

sas_df <- cbind(ev_min, ev_min_gen, ev_max, ev_max_gen, ev_mean, ev_mean_gen)

# undo logarithmic transformation
sas_df[c("log_vhm", "log_dist_water", "log_rough"),] <- exp(sas_df[c("log_vhm", "log_dist_water", "log_rough"),])

sas_df <- round(sas_df,2)
colnames(sas_df) <- c("min.FV", "min.AV", "max.FV", "max.AV", "mean.FV", "mean.AV")
rownames(sas_df) <- c("Elevation","Slope", "North exposure", "East exposure", "Solar radiation", "Distance from glacier", "Vegetation heigh", "Distance from water", "Roughness")

write.csv(sas_df, file = "data/sas_comp.csv")
