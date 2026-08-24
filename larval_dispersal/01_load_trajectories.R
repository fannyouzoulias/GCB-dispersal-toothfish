################################################################################
# 01_load_trajectories.R
#
# Builds the particle trajectories from the raw output of the advection model,
# then the release table that weights each particle by the egg production of the
# cell it was released from.
#
# Produces: df_traj    (one row per particle x timestep)
#           df_release (one row per particle: release position and egg weight)
#
# The raw output is about 900 MB over 216 files, so this takes several minutes
# and a few GB of memory. The result is saved next to the data so the later
# scripts can be re-run without redoing it.
#
# Requires: 00_setup.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

raw_dir <- dpath("advection_output")

## 1. Read the raw advection output --------------------------------------------
# One CSV per release week and year, in wide format: one row per particle, one
# column per date and coordinate. Release weeks 23-31 are the spawning weeks
# 22-30 plus the one week the eggs take to rise from about 1500 m to the surface
# mixed layer, during which no horizontal drift is applied.
parts  <- list()
offset <- 0   # keeps particle identifiers unique across years and weeks

for (year in 2000:2023) {
  for (week in 23:31) {

    f <- file.path(raw_dir,
                   paste0("advection_18weeks_from_week", week, "_", year, ".csv"))
    if (!file.exists(f)) { warning("missing file: ", f); next }

    message("reading ", basename(f))
    df <- readr::read_csv(f, show_col_types = FALSE)
    colnames(df)[1] <- "Particle_ID"
    df$Unique_ID <- df$Particle_ID + offset

    parts[[length(parts) + 1]] <- df %>%
      pivot_longer(cols = -c(Particle_ID, Unique_ID),
                   names_to = "Date_Coord", values_to = "Value") %>%
      separate(Date_Coord, into = c("Date", "Coord"), sep = " ", extra = "merge") %>%
      pivot_wider(names_from = "Coord", values_from = "Value") %>%
      mutate(Date = as.Date(Date, format = "%Y%m%d"),
             Week = week,
             Year = year)

    offset <- offset + dplyr::n_distinct(df$Particle_ID)
  }
}

df_traj <- bind_rows(parts) %>%
  dplyr::select(Unique_ID, Date, Longitude, Latitude, Week, Year) %>%
  filter(!is.na(Longitude))

rm(parts)

## 2. Seabed depth under each particle -----------------------------------------
# https://www.gebco.net/. The raster is cropped to the extent of the
# trajectories first, otherwise the lookup over ~10^8 positions is far too slow.
# Coordinates are already WGS84, like the raster.
depth_raster <- raster::raster(dpath("gebco_2024_n-40.0_s-59.0_w55.0_e89.0.tif"))

depth_raster <- raster::crop(
  depth_raster,
  raster::extent(min(df_traj$Longitude, na.rm = TRUE),
                 max(df_traj$Longitude, na.rm = TRUE),
                 min(df_traj$Latitude,  na.rm = TRUE),
                 max(df_traj$Latitude,  na.rm = TRUE))
)

df_traj$Depth <- raster::extract(
  depth_raster,
  as.matrix(dplyr::select(df_traj, Longitude, Latitude))
)

## 3. Beaching -----------------------------------------------------------------
# A particle reaching depth >= 0 has hit the coast. It is frozen at that position
# for the rest of the run rather than deleted, so every particle keeps a final
# position and the denominators stay comparable across years.
df_traj <- df_traj %>%
  arrange(Unique_ID, Date) %>%
  group_by(Unique_ID) %>%
  mutate(
    hit       = match(TRUE, Depth >= 0 & !is.na(Depth)),
    lon_hit   = if_else(is.na(hit), Longitude[1], Longitude[hit]),
    lat_hit   = if_else(is.na(hit), Latitude[1],  Latitude[hit]),
    Longitude = if_else(!is.na(hit) & row_number() >= hit, lon_hit, Longitude),
    Latitude  = if_else(!is.na(hit) & row_number() >= hit, lat_hit, Latitude),
    Depth     = if_else(!is.na(hit) & row_number() >= hit, 0, Depth)
  ) %>%
  ungroup() %>%
  dplyr::select(-hit, -lon_hit, -lat_hit)

## 4. Release table and egg weighting ------------------------------------------
# Each particle is weighted by the modelled egg production of the cell it was
# released from: its release position is matched to the nearest cell of the
# egg-production grid (about 10 x 10 km, coarser than the 0.05-degree release
# grid, so neighbouring particles often share a weight). Week_release is the
# spawning week.
egg_grid <- readr::read_csv(dpath("egg_production_grid.csv"), show_col_types = FALSE)

df_release <- df_traj %>%
  group_by(Unique_ID) %>%
  filter(Date == min(Date, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Week_release = Week - 1)

nn <- FNN::get.knnx(
  data  = as.matrix(egg_grid[, c("LON", "LAT")]),
  query = as.matrix(df_release[, c("Longitude", "Latitude")]),
  k     = 1
)$nn.index

df_release$Eggs_Released <- egg_grid$total_oeufs_mean[nn]

## 5. Save ---------------------------------------------------------------------
# Steps 1-4 are slow and deterministic; saving lets the analysis scripts be
# re-run on their own.
saveRDS(df_traj,    file.path(out_dir, "trajectories.rds"))
saveRDS(df_release, file.path(out_dir, "release_positions.rds"))
