################################################################################
# 04_abundance_model.R
#
# Spatial model of the number of mature (stage 4) females caught per longline
# station during the spawning season, fitted with sdmTMB: a Tweedie model with
# a smooth of depth, a linear effect of slope, log effort (hooks hauled) as an
# offset-like covariate, and a barrier spatial random field.
#
# Produces the predicted female abundance over the spawning area
#
# Requires: 00_setup.R, 01_load_data.R
# Also requires sdmTMBextra (not on CRAN): pbs-assess/sdmTMBextra
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

library(sdmTMBextra)

## Data ------------------------------------------------------------------------
# One row per station: the raised number of mature females, with the effort that
# produced it and the seabed slope. The selection (stage-4 females, weeks 22-30,
# 2020 onwards) is already applied; see 01_load_data.R and the data README.
data_N <- spawner_counts %>%
  mutate(depth_scaled = as.numeric(scale(DEPTH)),
         pente_scaled = as.numeric(scale(pente))) %>%
  as.data.frame() %>%
  sdmTMB::add_utm_columns(ll_names = c("LON", "LAT"), units = "km")

## Mesh ------------------------------------------------------------------------
# Study domain, with the island as a hole. Projected to UTM zone 42 and scaled
# to kilometres, to match the coordinates add_utm_columns() puts on the data.
study_domain <- st_read(dpath("study_domain.geojson"), quiet = TRUE) %>%
  st_transform(32742) %>%
  mutate(geometry = st_geometry(.) / 1000)
st_crs(study_domain) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"

mesh  <- sdmTMB::make_mesh(data_N, c("X", "Y"), cutoff = 10)
bspde <- sdmTMBextra::add_barrier_mesh(mesh, study_domain, range_fraction = 0.1,
                                       proj_scaling = 1, plot = FALSE)

## Fit -------------------------------------------------------------------------
# Tweedie handles the many zero-catch stations together with a skewed positive
# part, which a Poisson or lognormal model would not.
fit_N <- sdmTMB(
  formula        = nb_corr ~ s(depth_scaled) + pente_scaled + log(HAM_VIR),
  data           = data_N,
  mesh           = bspde,
  spatial        = "on",
  spatiotemporal = "off",
  family         = tweedie(link = "log")
)

summary(fit_N)
sanity(fit_N)

## Checks ----------------------------------------------------------------------
fit0 <- sdmTMB(nb_corr ~ 1, spatial = "off", data = data_N,
               family = tweedie(link = "log"), mesh = bspde)
fit_sp  <- update(fit0, spatial = "on")
fit_cov <- update(fit0,
                  formula. = nb_corr ~ s(depth_scaled) + pente_scaled + log(HAM_VIR))

message("deviance explained, full model: ", round(1 - deviance(fit_N)   / deviance(fit0), 3))
message("deviance explained, spatial field: ", round(1 - deviance(fit_sp)  / deviance(fit0), 3))
message("deviance explained, covariates: ", round(1 - deviance(fit_cov) / deviance(fit0), 3))

## Prediction over the spawning area -------------------------------------------
pred_grid <- readr::read_csv(dpath("prediction_grid.csv"), show_col_types = FALSE)

inside <- st_within(st_as_sf(pred_grid, coords = c("LON", "LAT"), crs = 4326),
                    release_area, sparse = FALSE)
# Effort is held at its mean, so the prediction is an abundance index at
# constant fishing effort rather than a prediction of a particular haul.
pred_grid <- pred_grid[apply(inside, 1, any), ] %>%
  mutate(depth_scaled = (DEPTH - mean(data_N$DEPTH)) / sd(data_N$DEPTH),
         pente_scaled = (pente - mean(data_N$pente)) / sd(data_N$pente),
         HAM_VIR      = mean(data_N$HAM_VIR)) %>%
  sdmTMB::add_utm_columns(ll_names = c("LON", "LAT"), units = "km")

# Draws are on the link scale, so they are exponentiated before averaging.
pred_N <- exp(predict(fit_N, newdata = pred_grid, nsim = 1000))

pred_grid$N_mean <- apply(pred_N, 1, mean)
pred_grid$N_sd   <- apply(pred_N, 1, sd)

## Map of predicted abundance --------------------------------------------------
p_abundance <- base_map +
  geom_tile(data = pred_grid, aes(x = LON, y = LAT, fill = log(N_mean))) +
  scale_fill_viridis_c(option = "inferno",
                       name = "Number of female spawners (log)    ") +
  iso_layers +
  geom_point(data = data_N, aes(x = LON, y = LAT), size = 0.3, colour = "black") +
  geom_sf(data = land_sf, fill = "darkgrey", colour = "darkgrey") +
  coord_sf(xlim = c(62, 70), ylim = c(-51.5, -47.9), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_paper() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

save_fig("predicted_female_abundance.png", p_abundance, width = 15, height = 10)

## Save for the egg-production step --------------------------------------------
saveRDS(list(grid = pred_grid, N_sims = pred_N),
        file.path(out_dir, "abundance_predictions.rds"))
