################################################################################
# 07_trajectory_maps_by_year.R
#
# Figure 5 of the paper, split year by year, with the Polar Front of that same
# year drawn on top.
#
# 02_trajectory_maps.R pools all years into one density of particle passages and
# overlays the CLIMATOLOGICAL front. Here the same density is computed per year
# and each panel carries the front position actually observed that year
# (A. Nalivaev, northern edge of the Winter Water probability of presence from
# release week 25 to week 29 + 18, i.e. the dispersal window despite the
# "week_25-29" file names), so the dispersal pattern can be read against a front
# that moves. The climatological PF (Park & Durand 2019), along which the
# published intensity index is sampled, is drawn in every panel for reference
# (SHOW_STATIC). The Polar Front jet of each year can be added in blue
# (SHOW_PARK_JET): the ADT streamline of Park et al. (2019), with their
# constraint that it round Kerguelen from the south on the 500-1000 m escarpment
# (advection/PF_position_park.py).
#
# Note that this is a HYDROGRAPHIC boundary, while the particles are advected by
# surface geostrophic velocities: the two need not coincide, and in 2012, 2013,
# 2017 and 2023 the contour is dragged well south of the jet (see below).
#
# Requires: 00_setup.R and 01_load_trajectories.R -- or nothing at all, since the
#           script reloads what it needs from outputs/ and from the data archive.
#
# Author: Fanny Ouzoulias
# Date:   2026-09-09
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

## ---------------------------------------------------------------------------
## PART 0. Paths and inputs
## ---------------------------------------------------------------------------
# Which front contour to draw:
#   "sent"    the files received (northernmost point topping 12 contiguous cells
#             above 0.8, then a 2 deg running mean)
#   "azarian" the literal definition of Azarian et al. 2024 (northernmost point
#             above 0.8), recomputed by 10_pf_position_azarian.R -- smoothed
#             copy, because the raw edge steps by several degrees between
#             neighbouring longitudes
# Each switch can be set before sourcing this script; these are the defaults.
if (!exists("PF_SOURCE")) PF_SOURCE <- "sent"

# Draw the annual Winter Water edge (black) of PF_SOURCE.
if (!exists("SHOW_WW")) SHOW_WW <- TRUE

# Also draw the climatological PF (dashed) in every panel.
if (!exists("SHOW_STATIC")) SHOW_STATIC <- TRUE

# Also draw the PF jet of each year (blue), from advection/PF_position_park.py.
# Its folder is `park_dir`, set in config.R.
if (!exists("SHOW_PARK_JET")) SHOW_PARK_JET <- TRUE

# Legend label of that jet.
if (!exists("JET_LABEL")) JET_LABEL <- "Annual PF streamline"

# Write the intensity and mean latitude of the jet in the sector of the PF
# index (67-72 E, 51-48 S, as 05) in each panel title. Needs SHOW_PARK_JET.
if (!exists("SHOW_PF_STATS")) SHOW_PF_STATS <- FALSE

dir_pos <- switch(PF_SOURCE,
  sent    = ww_dir,
  azarian = file.path("outputs", "revision", "position_PF_azarian"),
  stop("PF_SOURCE must be \"sent\" or \"azarian\"")
)
pf_file <- switch(PF_SOURCE,
  sent    = "position_PF_%d_week_25-29.csv",
  azarian = "position_PF_azarian_smoothed_%d_week_25-29.csv"
)
# Files of the data archive, through `data_dir` of config.R.
pick <- function(f) dpath(f)

if (!exists("out_dir")) out_dir <- "outputs"
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)

# Same theme as the published figures; redefined here only so the script can run
# without 00_setup.R having been sourced.
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

## Background layers, if 00_setup.R was not sourced -----------------------------
if (!exists("land_sf")) {
  land_sf <- rbind(
    st_read(pick("coastline.geojson"), quiet = TRUE),
    st_sf(geometry = st_geometry(dplyr::filter(
      rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
      name == "Heard I. and McDonald Is.")))
  )
}
## Trajectories ----------------------------------------------------------------
# 01_load_trajectories.R saves both tables; reload them rather than rebuilding.
if (!exists("df_traj"))
  df_traj <- readRDS(file.path(out_dir, "trajectories.rds"))
if (!exists("df_release"))
  df_release <- readRDS(file.path(out_dir, "release_positions.rds"))

## ---------------------------------------------------------------------------
## PART 1. Weighted density of particle passages, by year
## ---------------------------------------------------------------------------
release_weeks <- 24:28          # spawning peak, as in 02 and 04

df_w <- df_traj %>%
  mutate(Week_release = Week - 1) %>%
  filter(Week_release %in% release_weeks) %>%
  left_join(dplyr::select(df_release, Unique_ID, Eggs_Released), by = "Unique_ID") %>%
  filter(!is.na(Eggs_Released))

# Twice the cell size of 02_trajectory_maps.R (~10 km rather than ~5 km): the
# panels of a 24-facet plate are far too small for 5 km cells to be visible, and
# the coarser grid keeps the figure a readable size.
res_lon <- 0.13
res_lat <- 0.09
x0 <- floor(min(df_w$Longitude, na.rm = TRUE))
y0 <- floor(min(df_w$Latitude,  na.rm = TRUE))

# One year at a time: the de-duplication below is the memory-hungry step, and
# 24 years of trajectories at once do not need to be held in one intermediate.
grid_by_year <- map_dfr(sort(unique(df_w$Year)), function(y) {
  out <- df_w %>%
    filter(Year == y) %>%
    transmute(Week_release, Unique_ID, Eggs_Released,
              ix = as.integer((Longitude - x0) / res_lon),
              iy = as.integer((Latitude  - y0) / res_lat)) %>%
    # a particle crossing the same cell on several days counts once
    distinct(Week_release, Unique_ID, ix, iy, .keep_all = TRUE) %>%
    group_by(ix, iy) %>%
    summarise(Eggs_Weighted = sum(Eggs_Released, na.rm = TRUE), .groups = "drop") %>%
    mutate(Year      = y,
           Longitude = x0 + (ix + 0.5) * res_lon,
           Latitude  = y0 + (iy + 0.5) * res_lat)
  gc(verbose = FALSE)
  out
})

## ---------------------------------------------------------------------------
## PART 2. The front of each year
## ---------------------------------------------------------------------------
years <- sort(unique(grid_by_year$Year))

pf_year <- map_dfr(years, function(y) {
  f <- file.path(dir_pos, sprintf(pf_file, y))
  suppressMessages(read_csv(f, show_col_types = FALSE)) %>%
    dplyr::select(lon = lon_PF, lat = lat_PF) %>%
    mutate(Year = y)
})

# Same flag as explo/08_pf_interannual_position.R: mean latitude over the 67-72 E
# window, north or south of 51 S. In those four years the Winter Water field is
# patchier (8-12k cells above 0.8 against 14-17k, mean probability 0.47-0.53
# against 0.54-0.64 in this window), so the "12 contiguous cells" criterion has
# to walk much further south before it finds a solid block: the raw WW edge is
# only ~1.1 deg south of normal, the criterion turns that into ~2.4 deg.
pf_flag <- pf_year %>%
  filter(lon >= 67, lon <= 72) %>%
  group_by(Year) %>%
  summarise(lat_var = mean(lat), .groups = "drop") %>%
  mutate(position = if_else(lat_var >= -51, "typical", "displaced south"))

pf_year <- left_join(pf_year, dplyr::select(pf_flag, Year, position), by = "Year")

# Climatological PF: the same file, and the same point order, as the published
# index (05_front_indices.R) and Figure 5 (00_setup.R).
pf_clim <- read_csv(pick("front_intensity_PF.csv"), show_col_types = FALSE) %>%
  dplyr::select(lon = `Lon PF`, lat = `Lat PF`) %>%
  filter(!is.na(lon), !is.na(lat))

# PF jet of the year: the ADT streamline of Park et al. (2019), levelled, with
# the Kerguelen constraint applied every year (see advection/PF_position_park.py).
# Written in path order, west to east.
if (SHOW_PARK_JET) {
  jet_year <- map_dfr(years, function(y) {
    read_csv(file.path(park_dir, sprintf("pf_park_%d.csv", y)),
             show_col_types = FALSE) %>%
      dplyr::select(lon, lat, speed = mean_current_speed_cm_s) %>%
      mutate(Year = y)
  })
  jet_stats <- jet_year %>%
    filter(between(lon, 67, 72), between(lat, -51, -48)) %>%
    group_by(Year) %>%
    summarise(I = mean(speed, na.rm = TRUE), lat_mean = mean(lat), .groups = "drop")
}

front_key <- c(annual = "Annual northern limit of Winter Water",
               clim   = "Climatological Polar Front (Park & Durand, 2019)",
               jet    = JET_LABEL)
front_col <- set_names(c("black", "black", "#0072B2"), front_key)
front_lty <- set_names(c("solid", "42", "solid"), front_key)
# unnamed: ggplot takes the names of the breaks as labels
shown     <- unname(front_key[c(SHOW_WW, SHOW_STATIC, SHOW_PARK_JET)])
front_guide <- function(...) guide_legend(order = 1, direction = "vertical",
                                          keywidth = unit(2.5, "cm"), ...)

## ---------------------------------------------------------------------------
## PART 3. Panel labels
## ---------------------------------------------------------------------------
panel_lab <- if (SHOW_PARK_JET && SHOW_PF_STATS) {
  with(jet_stats[match(years, jet_stats$Year), ],
       sprintf("%d   %.1f cm/s, %.2f°S", years, I, -lat_mean))
} else as.character(years)
to_lab <- function(df) mutate(df, panel = factor(panel_lab[match(Year, years)],
                                                 levels = panel_lab))

grid_by_year <- to_lab(grid_by_year)
pf_year      <- to_lab(pf_year)
if (SHOW_PARK_JET) jet_year <- to_lab(jet_year)

## ---------------------------------------------------------------------------
## PART 4. The plate
## ---------------------------------------------------------------------------
# One colour scale for all panels so that years are comparable. The top of the
# scale is squished at a high quantile: a handful of cells just off the spawning
# grounds are an order of magnitude above the rest and would flatten everything.
cap <- as.numeric(quantile(grid_by_year$Eggs_Weighted, 0.995, na.rm = TRUE))

p_year <- ggplot() +
  geom_tile(data = grid_by_year,
            aes(x = Longitude, y = Latitude, fill = Eggs_Weighted),
            width = res_lon, height = res_lat) +
  # climatological front, dashed, no Year column so it repeats in every panel
  { if (SHOW_STATIC) list(
      geom_path(data = pf_clim, aes(lon, lat), colour = "white", linewidth = 1.1),
      geom_path(data = pf_clim, aes(lon, lat, linetype = front_key[["clim"]],
                                    colour = front_key[["clim"]]), linewidth = 0.5)) } +
  # PF jet of the year, blue, under the Winter Water edge
  { if (SHOW_PARK_JET) list(
      geom_path(data = jet_year, aes(lon, lat, group = Year),
                colour = "white", linewidth = 1.6),
      geom_path(data = jet_year, aes(lon, lat, group = Year, linetype = front_key[["jet"]],
                                     colour = front_key[["jet"]]), linewidth = 0.85)) } +
  # front of the year, with a white halo so it reads over the colours
  { if (SHOW_WW) list(
      geom_path(data = pf_year, aes(lon, lat, group = Year),
                colour = "white", linewidth = 1.1),
      geom_path(data = pf_year, aes(lon, lat, group = Year, linetype = front_key[["annual"]],
                                    colour = front_key[["annual"]]), linewidth = 0.45)) } +
  # colour and line type share one legend: same name, breaks and guide
  scale_colour_manual(name = NULL, values = front_col, breaks = shown,
                      guide = front_guide(override.aes = list(linewidth = 1.2))) +
  scale_linetype_manual(name = NULL, values = front_lty, breaks = shown,
                        guide = front_guide()) +
  geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 4) +
  scale_fill_viridis_c(
    option = "H", limits = c(0, cap), oob = scales::squish,
    name   = "Number of particles (weighted)",
    guide  = guide_colorbar(direction = "horizontal", title.position = "top",
                            barwidth = unit(14, "cm"), barheight = unit(0.6, "cm"))
  ) +
  # no 85 E label: at the panel edge it ran into the 60 E of the next panel
  scale_x_continuous(breaks = seq(60, 80, by = 5)) +
  coord_sf(xlim = c(60, 85), ylim = c(-54, -44), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + theme_paper() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.spacing.x = unit(1.8, "lines")) +
  # the longer titles of SHOW_PF_STATS need a smaller strip font
  { if (SHOW_PF_STATS) theme(strip.text = element_text(size = 15)) }

fig_name <- paste0("trajectories_by_year_with",
                   if (SHOW_WW) paste0("_front_", PF_SOURCE) else "",
                   if (SHOW_STATIC) (if (SHOW_WW) "_and_static" else "_static") else "",
                   if (SHOW_PARK_JET) "_park_jet" else "", ".png")
ggsave(file.path(rev_fig, fig_name), p_year,
       width = 17, height = 18, dpi = 250, limitsize = FALSE)

message("Written: ", normalizePath(file.path(rev_fig, fig_name), mustWork = FALSE))
