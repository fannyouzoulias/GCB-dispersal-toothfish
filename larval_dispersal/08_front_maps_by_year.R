################################################################################
# 08_front_maps_by_year.R
#
# The fronts of 07_trajectory_maps_by_year.R, on their own.
#
# 09 draws them over the egg-weighted density of particle passages; the colours
# of that density are what the eye goes to first, and the three lines have to
# fight the viridis "H" scale to be read at all. Here the same lines are drawn
# on a plain background, so that the geometry can be compared year to year and
# line to line:
#
#   FIGURE A  fronts_by_year.png
#             the three lines of 09, no density: the annual northern limit of
#             Winter Water (black), the climatological PF of Park & Durand
#             (2019) (black dashed) and the annual PF streamline (blue).
#
#   FIGURE B  fronts_by_year_jet_window_vs_fullyear.png
#             the ADT streamline alone, in its two versions: contoured on the
#             ADT averaged over the ADVECTION WINDOW (blue, weeks 23 to 31 + 18,
#             the published choice) and on the ADT averaged over the WHOLE
#             CALENDAR YEAR (orange, the sensitivity test). Both come from
#             advection/PF_position_park.py, the second with --fullyear.
#             Figure B is skipped, with a message, when it has not been run.
#
# Needs no trajectories: everything is read from the fronts written by
# advection/ (`park_dir` and `ww_dir` in config.R). Both plates keep
# the panel layout, the extent and the theme of 09 so they can be laid beside it.
#
# Author: Fanny Ouzoulias
# Date:   2026-09-21
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

## ---------------------------------------------------------------------------
## PART 0. Paths and inputs (as 09)
## ---------------------------------------------------------------------------
# The fronts come from advection/; their folders are set in config.R.
dir_pos <- ww_dir
dir_jet <- park_dir
years <- 2000:2023

# Files of the data archive, through `data_dir` of config.R.
pick <- function(f) dpath(f)

if (!exists("out_dir")) out_dir <- "outputs"
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)

if (!exists("theme_paper")) theme_paper <- function() {
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

if (!exists("land_sf")) {
  land_sf <- rbind(
    st_read(pick("coastline.geojson"), quiet = TRUE),
    st_sf(geometry = st_geometry(dplyr::filter(
      rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
      name == "Heard I. and McDonald Is.")))
  )
}

## ---------------------------------------------------------------------------
## PART 1. The lines
## ---------------------------------------------------------------------------
# Annual northern limit of Winter Water (A. Nalivaev), as read by 09.
pf_year <- map_dfr(years, function(y) {
  suppressMessages(read_csv(
    file.path(dir_pos, sprintf("position_PF_%d_week_25-29.csv", y)),
    show_col_types = FALSE)) %>%
    dplyr::select(lon = lon_PF, lat = lat_PF) %>%
    mutate(Year = y)
})

# Climatological PF: the same file and point order as the published index.
pf_clim <- read_csv(pick("front_intensity_PF.csv"), show_col_types = FALSE) %>%
  dplyr::select(lon = `Lon PF`, lat = `Lat PF`) %>%
  filter(!is.na(lon), !is.na(lat))

# ADT streamline: the advection window, and the whole calendar year if
# PF_position_park.py has also been run with --fullyear.
read_jet <- function(dir, pat) map_dfr(years, function(y) {
  read_csv(file.path(dir, sprintf(pat, y)), show_col_types = FALSE) %>%
    dplyr::select(lon, lat) %>% mutate(Year = y)
})

jet_year <- read_jet(dir_jet, "pf_park_%d.csv")

dir_full <- file.path(dir_jet, "fullyear")
has_full <- all(file.exists(file.path(dir_full, sprintf("pf_park_%d.csv", years))))
if (has_full) jet_full <- read_jet(dir_full, "pf_park_%d.csv")

## ---------------------------------------------------------------------------
## PART 2. A plate of 24 panels, shared by both figures
## ---------------------------------------------------------------------------
to_lab <- function(df) mutate(df, panel = factor(Year, levels = years))
pf_year  <- to_lab(pf_year)
jet_year <- to_lab(jet_year)
if (has_full) jet_full <- to_lab(jet_full)

# Colour and line type share one legend: same name, same breaks, same guide.
guide_front <- function(...) guide_legend(order = 1, direction = "vertical",
                                          keywidth = unit(2.5, "cm"), ...)

plate <- function(layers, key, col, lty) {
  ggplot() +
    layers +
    geom_sf(data = land_sf, fill = "darkgrey", colour = NA, inherit.aes = FALSE) +
    scale_colour_manual(name = NULL, values = col, breaks = unname(key),
                        guide = guide_front(override.aes = list(linewidth = 1.2))) +
    scale_linetype_manual(name = NULL, values = lty, breaks = unname(key),
                          guide = guide_front()) +
    facet_wrap(~ panel, ncol = 4) +
    scale_x_continuous(breaks = seq(60, 80, by = 5)) +
    coord_sf(xlim = c(60, 85), ylim = c(-54, -44), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() + theme_paper() +
    theme(legend.position = "bottom", legend.box = "vertical",
          panel.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.spacing.x  = unit(1.8, "lines"))
}

## ---------------------------------------------------------------------------
## FIGURE A. The three fronts, no density
## ---------------------------------------------------------------------------
keyA <- c(annual = "Annual northern limit of Winter Water",
          clim   = "Climatological Polar Front (Park & Durand, 2019)",
          jet    = "Annual PF streamline")
colA <- set_names(c("black", "black", "#0072B2"), keyA)
ltyA <- set_names(c("solid", "42", "solid"), keyA)

# No white halo here: on a white panel it would only thicken the lines.
pA <- plate(list(
  geom_path(data = pf_clim, aes(lon, lat, linetype = keyA[["clim"]],
                                colour = keyA[["clim"]]), linewidth = 0.6),
  geom_path(data = jet_year, aes(lon, lat, group = Year, linetype = keyA[["jet"]],
                                 colour = keyA[["jet"]]), linewidth = 0.9),
  geom_path(data = pf_year, aes(lon, lat, group = Year, linetype = keyA[["annual"]],
                                colour = keyA[["annual"]]), linewidth = 0.6)),
  keyA, colA, ltyA)

ggsave(file.path(rev_fig, "fronts_by_year.png"), pA,
       width = 17, height = 18, dpi = 250, limitsize = FALSE)
message("Written: ", normalizePath(file.path(rev_fig, "fronts_by_year.png"),
                                   mustWork = FALSE))

## ---------------------------------------------------------------------------
## FIGURE B. The streamline: advection window against the whole year
## ---------------------------------------------------------------------------
if (!has_full) {
  message("Figure B skipped: run advection/PF_position_park.py --fullyear first ",
          "(expected in ", dir_full, ")")
} else {
  keyB <- c(win  = "Annual PF streamline (advection window, weeks 23-31 + 18)",
            full = "Annual PF streamline (whole calendar year)")
  colB <- set_names(c("#0072B2", "#D55E00"), keyB)
  ltyB <- set_names(c("solid", "solid"), keyB)

  pB <- plate(list(
    geom_path(data = jet_full, aes(lon, lat, group = Year, linetype = keyB[["full"]],
                                   colour = keyB[["full"]]), linewidth = 0.9),
    geom_path(data = jet_year, aes(lon, lat, group = Year, linetype = keyB[["win"]],
                                   colour = keyB[["win"]]), linewidth = 0.9)),
    keyB, colB, ltyB)

  ggsave(file.path(rev_fig, "fronts_by_year_jet_window_vs_fullyear.png"), pB,
         width = 17, height = 18, dpi = 250, limitsize = FALSE)
  message("Written: ", normalizePath(
    file.path(rev_fig, "fronts_by_year_jet_window_vs_fullyear.png"), mustWork = FALSE))
}
