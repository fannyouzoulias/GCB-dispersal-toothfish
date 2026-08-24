################################################################################
# 02_trajectory_maps.R
#
# Figure 5: dispersal pathways over the 18-week pelagic phase, as densities of
# particle passages weighted by the egg production of the release cell.
#   A) all years
#   B) low-retention years  (2010, 2015, 2023)
#   C) high-retention years (2000, 2016, 2022)
# The mean positions of the Polar Front and Subantarctic Front are overlaid.
#
#
# Requires: 00_setup.R, 01_load_trajectories.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

## Weighted trajectories -------------------------------------------------------
# Releases from weeks 24-28 only (spawning peak)
release_weeks <- 24:28

df_w <- df_traj %>%
  mutate(Week_release = Week - 1) %>%
  filter(Week_release %in% release_weeks) %>%
  # Eggs_Released is the modelled egg production
  left_join(dplyr::select(df_release, Unique_ID, Eggs_Released), by = "Unique_ID") %>%
  filter(!is.na(Eggs_Released))

## Density ---------------------------------------------------------------------
# Positions are binned on a regular grid (~5 km).
res_lon <- 0.065
res_lat <- 0.045
x0 <- floor(min(df_w$Longitude, na.rm = TRUE))
y0 <- floor(min(df_w$Latitude,  na.rm = TRUE))

agg_weighted_grid <- function(df_sub) {
  df_sub %>%
    transmute(
      Year, Week_release, Unique_ID, Eggs_Released,
      ix = as.integer((Longitude - x0) / res_lon),
      iy = as.integer((Latitude  - y0) / res_lat)
    ) %>%
    distinct(Year, Week_release, Unique_ID, ix, iy, .keep_all = TRUE) %>%
    group_by(ix, iy) %>%
    summarise(Eggs_Weighted = sum(Eggs_Released, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Longitude = x0 + (ix + 0.5) * res_lon,
      Latitude  = y0 + (iy + 0.5) * res_lat
    )
}

## Panel A: all years ----------------------------------------------------------
# Spawning areas are outlined in white.
grid_all <- agg_weighted_grid(df_w)

p_all <- ggplot(grid_all, aes(x = Longitude, y = Latitude, fill = Eggs_Weighted)) +
  geom_tile(colour = NA) +
  front_layers +
  geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
  geom_sf(data = spawning_zones, colour = "white", fill = NA,
          linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_viridis_c(
    option = "H",
    name   = "Number of particles (weighted)",
    guide  = guide_colorbar(direction = "horizontal", title.position = "top",
                            barwidth = unit(14, "cm"), barheight = unit(0.6, "cm"))
  ) +
  coord_sf(xlim = c(60, 85), ylim = c(-54, -44), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + theme_paper() +
  theme(legend.position = "bottom", legend.box = "horizontal")

save_fig("trajectories_all_years.png", p_all, width = 15, height = 10)

## Panels B and C: contrasted years --------------------------------------------
# Three lowest and three highest retention years (see 04_retention_recruitment.R).
low_years  <- c(2010, 2015, 2023)
high_years <- c(2000, 2016, 2022)

grid_low  <- agg_weighted_grid(filter(df_w, Year %in% low_years))
grid_high <- agg_weighted_grid(filter(df_w, Year %in% high_years))

max_common <- max(grid_low$Eggs_Weighted, grid_high$Eggs_Weighted, na.rm = TRUE)

plot_weighted_group <- function(df_grid) {
  ggplot(df_grid, aes(x = Longitude, y = Latitude, fill = Eggs_Weighted)) +
    geom_tile(colour = NA) +
    front_layers +
    geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "H", name = "Number of particles\n(weighted)",
                         limits = c(0, max_common), oob = scales::squish) +
    coord_sf(xlim = c(60, 85), ylim = c(-54, -44), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() + theme_paper() +
    theme(legend.position = "right")
}

save_fig("trajectories_low_retention.png",  plot_weighted_group(grid_low),
         width = 15, height = 10)
save_fig("trajectories_high_retention.png", plot_weighted_group(grid_high),
         width = 15, height = 10)
