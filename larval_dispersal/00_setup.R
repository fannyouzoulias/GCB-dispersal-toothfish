################################################################################
# 00_setup.R  --  larval dispersal
#
# Packages, plotting theme and the background layers (coastlines, isobaths,
# recruitment sectors) shared by all the scripts of this folder.
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
  "raster",         # seabed depth sampled along the trajectories
  "MASS",           # bandwidth.nrd, for the weighted kernel densities
  "FNN",            # fast nearest-neighbour lookup of the egg weights
  "patchwork", "cowplot",
  "rnaturalearth",  # coastlines (avoids redistributing a shapefile)
  "here",           # project-root-relative paths
  "conflicted"
)
invisible(lapply(pkgs, library, character.only = TRUE))

## Namespace conflicts ---------------------------------------------------------
# MASS and raster both mask dplyr verbs; settle them once for the whole chain.
conflicted::conflicts_prefer(
  dplyr::select, dplyr::filter, dplyr::mutate, dplyr::arrange,
  dplyr::summarise, dplyr::first, dplyr::rename, dplyr::count,
  ggplot2::annotate, ggplot2::theme_void, ggplot2::ggsave,
  lubridate::year, base::intersect,
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
      plot.subtitle    = element_text(colour = "black", size = 22),
      strip.text       = element_text(size = 18),
      strip.background = element_rect(fill = "lightgrey", colour = "black"),
      panel.spacing    = unit(1, "lines"),
      plot.margin      = margin(t = 10, r = 30, b = 10, l = 10)
    )
}

# Palette used for the particle-density maps.
density_palette <- c("white", "#FBE0DA", "#F6B8AF", "#F28F84",
                     "#ED675A", "#E83E2F", "#D93B2D")

## Background layers -----------------------------------------------------------
# Detailed coastline of the French islands, from the archive. Heard & McDonald
# are Australian and absent from that layer, so they come from Natural Earth,
# whose "medium" scale resolves them. coord_sf() clips both to each map extent.
land_sf <- rbind(
  st_read(dpath("coastline.geojson"), quiet = TRUE),
  st_sf(geometry = st_geometry(dplyr::filter(
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
    name == "Heard I. and McDonald Is.")))
)

# Recruitment sectors (north, east, south, west) and spawning areas.
recruitment_zones <- st_read(dpath("zones_recruitment.geojson"), quiet = TRUE)
spawning_zones    <- st_read(dpath("zones_spawning.geojson"),    quiet = TRUE)

# Mean positions of the Polar Front (PF) and Subantarctic Front (SAF), read from
# the same files that carry the yearly intensities.
front_df <- bind_rows(
  read_csv(dpath("front_intensity_PF.csv"), show_col_types = FALSE) %>%
    dplyr::select(lon = `Lon PF`, lat = `Lat PF`) %>%
    mutate(Front = "Polar Front"),
  read_csv(dpath("front_intensity_SAF.csv"), show_col_types = FALSE) %>%
    dplyr::select(lon = `Lon SAF`, lat = `Lat SAF`) %>%
    mutate(Front = "Subantarctic Front")
) %>%
  filter(!is.na(lon), !is.na(lat))


front_layers <- list(
  geom_path(data = front_df, aes(x = lon, y = lat, group = Front),
            colour = "white", linewidth = 1.4, inherit.aes = FALSE),
  geom_path(data = front_df, aes(x = lon, y = lat, linetype = Front, group = Front),
            colour = "black", linewidth = 0.7, inherit.aes = FALSE),
  scale_linetype_manual(
    name   = "Front",
    values = c("Polar Front" = "solid", "Subantarctic Front" = "dashed")
  )
)

## Helpers ---------------------------------------------------------------
# Weighted 2-D kernel density. Same as MASS::kde2d but every point carries a
# weight (here the egg production of its release cell).
kde2d_weighted <- function(x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) stop("x and y must have the same length")
  if (missing(w)) w <- rep(1, nx)
  if (missing(h)) h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
  gx <- seq(lims[1], lims[2], length.out = n)
  gy <- seq(lims[3], lims[4], length.out = n)
  h  <- h / 4
  ax <- outer(gx, x, "-") / h[1]
  ay <- outer(gy, y, "-") / h[2]
  z  <- (matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE) *
           matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx)) /
    (sum(w) * h[1] * h[2])
  list(x = gx, y = gy, z = z)
}

# Save a figure to outputs/figures with the dimensions used in the paper.
save_fig <- function(name, plot, width = 10, height = 7, dpi = 300) {
  ggsave(file.path(fig_dir, name), plot,
         width = width, height = height, dpi = dpi)
}
