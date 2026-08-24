################################################################################
# 00_setup.R
#
# Packages, plotting theme and background layers shared by the scripts that
# model female abundance, female length and fecundity, and combine them into
# gridded egg production over the spawning area.
# Source this first, then the numbered scripts in order.
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

## Packages --------------------------------------------------------------------
pkgs <- c(
  "tidyverse",      # dplyr, ggplot2, tidyr, readr, purrr, stringr
  "lubridate",
  "sf",             # vector spatial data
  "sdmTMB",         # spatio-temporal models of abundance and length
  "ggnewscale",     # several colour scales on the same map
  "patchwork",      # side-by-side diagnostics
  "rnaturalearth",  # coastlines (avoids redistributing a shapefile)
  "here",           # project-root-relative paths
  "conflicted"
)
invisible(lapply(pkgs, library, character.only = TRUE))

conflicted::conflicts_prefer(
  dplyr::select, dplyr::filter, dplyr::mutate, dplyr::arrange,
  dplyr::summarise, dplyr::first, dplyr::rename, dplyr::count,
  ggplot2::annotate, ggplot2::ggsave, lubridate::year, base::intersect,
  .quiet = TRUE
)

## Paths -----------------------------------------------------------------------
source(here::here("config.R"))   # defines data_dir, fig_dir, dpath()

## Plotting theme --------------------------------------------------------------
theme_paper <- function() {
  theme_bw() +
    theme(
      axis.title       = element_text(colour = "black", size = 20),
      axis.text        = element_text(size = 20),
      legend.title     = element_text(colour = "black", size = 20,
                                      margin = margin(b = 20)),
      legend.text      = element_text(colour = "black", size = 20),
      plot.title       = element_text(colour = "black", size = 22, face = "bold"),
      strip.text       = element_text(size = 18),
      strip.background = element_rect(fill = "lightgrey", colour = "black"),
      panel.spacing    = unit(1, "lines"),
      plot.margin      = margin(t = 10, r = 30, b = 10, l = 10)
    )
}

## Background layers -----------------------------------------------------------
# Detailed coastline, from the archive. The models use a generalised version of
# the same data as their barrier (study_domain.geojson); this one is for maps.
land_sf <- st_read(dpath("coastline.geojson"), quiet = TRUE)

# Isobaths (150, 500, 700, 1500, 1800 m), derived from GEBCO.
isobaths_sf <- st_read(dpath("isobaths_kerguelen.geojson"), quiet = TRUE) %>%
  mutate(level = as.numeric(as.character(level)))

# Drawn light to dark with increasing depth, as reusable layers.
iso_layers <- lapply(
  list(c(-500, "#C6DBEF"), c(-700, "#9ECAE1"),
       c(-1500, "#6BAED6"), c(-1800, "#2171B5")),
  function(x) geom_sf(data = filter(isobaths_sf, level == as.numeric(x[1])),
                      colour = x[2], linewidth = 0.2, inherit.aes = FALSE)
)

# Spawning areas over which egg production is summed.
spawning_zones    <- st_read(dpath("zones_spawning.geojson"),    quiet = TRUE)

# Recruitment sectors, drawn on the final egg-production map to place the
# spawning hotspot relative to where the larvae have to end up.
recruitment_zones <- st_read(dpath("zones_recruitment.geojson"), quiet = TRUE)

# Release area: the single polygon the particles are seeded in. The models are
# predicted over this area, not over the three spawning sectors, which cover a
# slightly different set of cells.
release_area <- st_read(dpath("release_area.geojson"), quiet = TRUE)

## Base map --------------------------------------------------------------------
base_map <- ggplot() +
  geom_sf(data = filter(isobaths_sf, level %in% c(-500, -1500, -1800)),
          colour = "grey80", linewidth = 0.2) +
  geom_sf(data = land_sf, fill = "darkgrey", colour = "black") +
  coord_sf(xlim = c(64, 75), ylim = c(-52, -46), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.grid = element_blank())

## Helper ----------------------------------------------------------------------
save_fig <- function(name, plot, width = 10, height = 7, dpi = 300) {
  ggsave(file.path(fig_dir, name), plot,
         width = width, height = height, dpi = dpi)
}
