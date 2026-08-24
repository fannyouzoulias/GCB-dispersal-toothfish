################################################################################
# 02_length_model.R
#
# Spatial model of the total length of mature (stage 4) females during the
# spawning season, fitted with sdmTMB: a Gaussian model on log length with a
# smooth of depth, a linear effect of seabed slope, and a spatial random field.
#
# The island is a hole in the domain, so the field uses a barrier mesh: two
# points on opposite sides of Kerguelen are not neighbours even when they are
# close as the crow flies.
#
# Produces the predicted mean length over the spawning area and the 1000
# posterior draws used downstream to propagate uncertainty into egg numbers.
#
# Requires: 00_setup.R, 01_load_data.R
# Also requires sdmTMBextra (not on CRAN): pbs-assess/sdmTMBextra
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

library(sdmTMBextra)

## Data ------------------------------------------------------------------------
# The table arrives already restricted to mature females over the spawning
# season (weeks 22-30, 2020 onwards) and already carries the seabed slope; see
# 01_load_data.R and the data README for the selection.
data_LT <- spawner_lengths %>%
  mutate(depth_scaled = as.numeric(scale(DEPTH)),
         pente_scaled = as.numeric(scale(pente))) %>%
  as.data.frame()

## Mesh ------------------------------------------------------------------------
# Coordinates are projected to UTM (km) so the range parameter is in kilometres.
data_LT <- sdmTMB::add_utm_columns(data_LT, ll_names = c("LON", "LAT"), units = "km")

# Study domain, with the island as a hole. Projected to UTM zone 42 and scaled
# to kilometres, to match the coordinates add_utm_columns() puts on the data.
study_domain <- st_read(dpath("study_domain.geojson"), quiet = TRUE) %>%
  st_transform(32742) %>%
  mutate(geometry = st_geometry(.) / 1000)
st_crs(study_domain) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"

mesh  <- sdmTMB::make_mesh(data_LT, c("X", "Y"), cutoff = 10)
bspde <- sdmTMBextra::add_barrier_mesh(mesh, study_domain, range_fraction = 0.1,
                                       proj_scaling = 1, plot = FALSE)

## Fit -------------------------------------------------------------------------
fit_LT <- sdmTMB(
  formula = log(LT) ~ s(depth_scaled) + pente_scaled,
  data    = data_LT,
  mesh    = bspde,
  spatial = "on",
  family  = gaussian(link = "identity")
)

summary(fit_LT)
sanity(fit_LT)

## Checks ----------------------------------------------------------------------
# Deviance explained
fit0 <- sdmTMB(log(LT) ~ 1, spatial = "off", data = data_LT,
               family = gaussian(link = "identity"), mesh = bspde)
fit_sp  <- update(fit0, spatial = "on")
fit_cov <- update(fit0, formula. = log(LT) ~ s(depth_scaled) + pente_scaled)

message("deviance explained, full model: ", round(1 - deviance(fit_LT)  / deviance(fit0), 3))
message("deviance explained, spatial field: ", round(1 - deviance(fit_sp)  / deviance(fit0), 3))
message("deviance explained, covariates: ", round(1 - deviance(fit_cov) / deviance(fit0), 3))

data_LT$resids <- residuals(fit_LT, type = "mle-mvn")

p_qq <- ggplot(data_LT, aes(sample = resids)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  stat_qq() + stat_qq_line() +
  labs(title = "Length model: normal Q-Q") + theme_bw()

p_res_map <- ggplot(data_LT, aes(X, Y, colour = resids)) +
  geom_point() + scale_colour_gradient2() + coord_fixed() +
  labs(title = "Length model: spatial residuals") + theme_bw()

save_fig("length_model_diagnostics.png", p_qq | p_res_map, width = 12, height = 5)

## Prediction over the spawning area -------------------------------------------
# The prediction grid is restricted to the spawning polygons
pred_grid <- readr::read_csv(dpath("prediction_grid.csv"), show_col_types = FALSE)

inside <- st_within(st_as_sf(pred_grid, coords = c("LON", "LAT"), crs = 4326),
                    release_area, sparse = FALSE)
pred_grid <- pred_grid[apply(inside, 1, any), ] %>%
  mutate(depth_scaled = (DEPTH - mean(data_LT$DEPTH)) / sd(data_LT$DEPTH),
         pente_scaled = (pente - mean(data_LT$pente)) / sd(data_LT$pente)) %>%
  sdmTMB::add_utm_columns(ll_names = c("LON", "LAT"), units = "km")

# 1000 draws from the joint posterior of the fixed effects and the spatial
# field: one column per draw, one row per grid cell.
pred_LT <- predict(fit_LT, newdata = pred_grid, nsim = 1000)

pred_grid$log_LT_mean <- apply(pred_LT, 1, mean)
pred_grid$log_LT_sd   <- apply(pred_LT, 1, sd)

## Map of predicted length -----------------------------------------------------
p_length <- base_map +
  geom_tile(data = filter(pred_grid, LON < 69),
            aes(x = LON, y = LAT, fill = exp(log_LT_mean))) +
  scale_fill_viridis_c(option = "inferno", name = "TL (cm)    ") +
  iso_layers +
  # the stations the model was fitted to
  geom_point(data = data_LT, aes(x = LON, y = LAT), size = 0.3, colour = "black") +
  geom_sf(data = land_sf, fill = "darkgrey", colour = "darkgrey") +
  coord_sf(xlim = c(62, 70), ylim = c(-51.5, -47.9), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_paper() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

save_fig("predicted_female_length.png", p_length, width = 15, height = 10)

## Save for the egg-production step --------------------------------------------
saveRDS(list(grid = pred_grid, logLT_sims = pred_LT),
        file.path(out_dir, "length_predictions.rds"))
