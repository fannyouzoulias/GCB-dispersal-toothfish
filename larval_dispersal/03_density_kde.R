################################################################################
# 03_density_kde.R
#
# Kernel densities of the particle positions at the end of the 18-week
# simulation, weighted by the egg production of the release cell.
#
#   - one map per year, on a common colour scale (appendix figure)
#   - one map pooling all years, restricted to particles that ended on the
#     0-150 m shelf, i.e. the recruits (main-text figure)
#
# Requires: 00_setup.R, 01_load_trajectories.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

release_weeks <- 24:28
years         <- 2000:2023

## End-of-simulation positions -------------------------------------------------
final_positions <- df_traj %>%
  mutate(Week_release = Week - 1) %>%
  filter(Week_release %in% release_weeks) %>%
  group_by(Unique_ID) %>%
  filter(Date == max(Date, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(dplyr::select(df_release, Unique_ID, Eggs_Released), by = "Unique_ID") %>%
  filter(!is.na(Longitude), !is.na(Latitude), !is.na(Eggs_Released))

## Yearly densities on a common scale ------------------------------------------
# The whole dispersal domain, no depth filter: this shows where the larvae end
# up, whether or not they land on the shelf.
x_range <- c(60, 85)
y_range <- c(-52, -44)

kde_by_year <- lapply(years, function(y) {
  d <- filter(final_positions, Year == y)
  if (nrow(d) == 0) return(NULL)
  k <- kde2d_weighted(d$Longitude, d$Latitude, w = d$Eggs_Released,
                      n = 300, lims = c(x_range, y_range))
  expand.grid(Longitude = k$x, Latitude = k$y) %>%
    mutate(Density = as.vector(k$z), Year = y)
})
names(kde_by_year) <- years
kde_by_year <- compact(kde_by_year)

# Common colour limits, so the panels can be compared year to year.
fill_limits <- c(0, max(sapply(kde_by_year,
                               function(d) max(log(d$Density + 1), na.rm = TRUE))))

map_kde_year <- function(df_kde) {
  ggplot() +
    geom_raster(data = df_kde,
                aes(x = Longitude, y = Latitude, fill = log(Density + 1))) +
    scale_fill_gradientn(colours = density_palette, limits = fill_limits,
                         name = "log(Density + 1)",
                         guide = guide_colorbar(barwidth  = unit(15, "cm"),
                                                barheight = unit(0.5, "cm"))) +
    geom_sf(data = recruitment_zones, colour = "darkgreen", fill = NA,
            linewidth = 1, inherit.aes = FALSE) +
    geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
    coord_sf(xlim = c(62, 75), ylim = c(-51, -46)) +
    labs(x = "Longitude", y = "Latitude", title = unique(df_kde$Year)) +
    theme_bw() + theme_paper() +
    theme(axis.text = element_text(size = 14),
          plot.margin = margin(0, 0, 0, 0))
}

panel <- wrap_plots(lapply(kde_by_year, map_kde_year), ncol = 5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_fig("density_kde_by_year.png", panel, width = 30, height = 20)

## Pooled density of the recruits ----------------------------------------------
# Recruits = particles whose final position lies on the 0-150 m shelf. The map is
# zoomed on the northern plateau, where essentially all of them are found.
recruits <- filter(final_positions, Depth >= -150, Depth <= 0)

kde_all <- kde2d_weighted(recruits$Longitude, recruits$Latitude,
                          w = recruits$Eggs_Released,
                          n = 400, lims = c(68, 71, -50, -48))

df_kde <- expand.grid(Longitude = kde_all$x, Latitude = kde_all$y) %>%
  mutate(Density = as.vector(kde_all$z))

p_recruits <- ggplot() +
  geom_raster(data = df_kde,
              aes(x = Longitude, y = Latitude, fill = log(Density + 1))) +
  scale_fill_gradientn(
    colours = density_palette,
    name    = "log\n(Density + 1)",
    limits  = c(0, log(max(df_kde$Density, na.rm = TRUE) + 1)),
    guide   = guide_colorbar(direction = "vertical", title.position = "top",
                             label.position = "right",
                             barwidth = unit(0.5, "cm"), barheight = unit(5, "cm"))
  ) +
  geom_sf(data = recruitment_zones, colour = "darkgreen", fill = NA,
          linewidth = 1, inherit.aes = FALSE) +
  geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
  coord_sf(xlim = c(67.9, 71.1), ylim = c(-50.1, -47.9), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + theme_paper() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10),
        axis.title   = element_text(size = 18),
        axis.text    = element_text(size = 12),
        panel.grid   = element_blank())

save_fig("density_kde_recruits_all_years.png", p_recruits, width = 10, height = 7)
