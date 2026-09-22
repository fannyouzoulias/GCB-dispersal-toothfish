################################################################################
# 05_front_indices.R
#
# Annual intensity indices of the Polar Front (PF) and the Subantarctic Front
# (SAF), from the geostrophic current speed sampled along the position of each
# front (DUACS altimetry).
#
# THE PF INDEX IS THE CORRECTION MADE AT REVISION. The manuscript sampled the
# speed along the CLIMATOLOGICAL front, one contour fixed for all years: that
# index measures how fast the water runs where the front is ON AVERAGE, not
# where it was that year. The index used here is the same speed sampled along
# the front OF THAT YEAR, the ADT contour of Park et al. (2019) with their
# Kerguelen constraint, written by advection/PF_position_park.py.
#
#   PF_park_cm_s   annual ADT contour       the index, used everywhere
#   PF_mean_cm_s   climatological contour   the manuscript index, kept so that
#                                           10_pf_park_intensity_glm.R can show
#                                           what the correction changes
#
# Both are averaged over the same sector, 67-72 E and lat >= -51, so the PATH
# is the only difference between them.
#
#   - map of the front positions and of the sector they are averaged over (appendix)
#   - annual PF index, from the front of each year
#   - annual SAF index, averaged over 63-73 E
#
# Writes front_indices_annual.csv, used by 06_retention_front_glm.R.
#
# Requires: 00_setup.R, and the annual PF contours in `park_dir`
#           (advection/PF_position_park.py)
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19, PF index corrected 2026-09-22
################################################################################

col_pf  <- "#4575b4"
col_saf <- "#1b9e77"

# The sector the PF index is averaged over. The latitude cut keeps the branch
# that runs along the northern plateau: the annual contour loops south of 51 S
# west of the islands in some years.
LON_BAND <- c(67, 72)
LAT_MIN  <- -51

## Front points, long format ---------------------------------------------------
# Each row is a point along the climatological front; the columns named
# 2000...2023 hold the mean geostrophic speed at that point for that year.
pf_pts <- read_csv(dpath("front_intensity_PF.csv"), show_col_types = FALSE) %>%
  rename(lon = `Lon PF`, lat = `Lat PF`)

saf_pts <- read_csv(dpath("front_intensity_SAF.csv"), show_col_types = FALSE) %>%
  rename(lon = `Lon SAF`, lat = `Lat SAF`)

to_long <- function(df, value_name) {
  df %>%
    pivot_longer(cols = matches("^[0-9]{4}$"),
                 names_to = "Year", values_to = value_name) %>%
    mutate(Year = as.integer(Year)) %>%
    filter(!is.na(.data[[value_name]]))
}

pf_long  <- to_long(pf_pts,  "intensity_pf")
saf_long <- to_long(saf_pts, "intensity_saf")

## The front of each year ------------------------------------------------------
# pf_park_<year>.csv holds the ADT contour of that year and, at each of its
# points, the mean current speed over the advection window. Reading it is all
# the index needs: the speed is sampled along the line by the Python script.
years_pf <- sort(unique(pf_long$Year))
f_park   <- file.path(park_dir, sprintf("pf_park_%d.csv", years_pf))

if (!all(file.exists(f_park)))
  stop("missing the annual PF contours in\n  ", park_dir,
       "\n  -> run advection/PF_position_park.py, then set `park_dir` in config.R")

pf_park_pts <- purrr::map2_dfr(f_park, years_pf, function(f, y) {
  read_csv(f, show_col_types = FALSE) %>%
    transmute(Year = y, lon, lat, speed = mean_current_speed_cm_s)
})

## Appendix figure: front positions and averaging sectors ----------------------
# The two boxes are the sectors the yearly speeds are averaged over: the PF
# where it runs along the northern plateau, the SAF further north. Grey lines
# are the annual PF contours the index is sampled along, brown points the
# climatological front of the manuscript.
box_pf  <- tibble(xmin = LON_BAND[1], xmax = LON_BAND[2], ymin = LAT_MIN, ymax = -48)
box_saf <- tibble(xmin = 63, xmax = 73, ymin = -47, ymax = -44)

p_fronts <- ggplot() +
  geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
  geom_path(data = pf_park_pts, aes(lon, lat, group = Year),
            colour = "grey55", linewidth = 0.3, alpha = 0.7) +
  geom_point(data = pf_long,  aes(lon, lat), size = 0.8, alpha = 0.9,
             colour = "brown4") +
  geom_point(data = saf_long, aes(lon, lat), size = 0.8, alpha = 0.9,
             colour = "blue3") +
  geom_rect(data = box_saf, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "blue3", alpha = 0.2, colour = "blue3", linewidth = 1) +
  geom_rect(data = box_pf, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "brown4", alpha = 0.2, colour = "brown4", linewidth = 1) +
  coord_sf(xlim = c(60, 85), ylim = c(-54, -44), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + theme_paper()

save_fig("fronts_PF_SAF.png", p_fronts, width = 8, height = 6)

## Annual PF index, along the front of that year (the index) -------------------
pf_park_index <- pf_park_pts %>%
  filter(lon >= LON_BAND[1], lon <= LON_BAND[2], lat >= LAT_MIN) %>%
  group_by(Year) %>%
  summarise(PF_park_cm_s = mean(speed, na.rm = TRUE), .groups = "drop")

## Annual PF index, along the climatological front (the manuscript) ------------
pf_index <- pf_long %>%
  filter(lon >= LON_BAND[1], lon <= LON_BAND[2], lat >= LAT_MIN) %>%
  group_by(Year) %>%
  summarise(PF_mean_cm_s = mean(intensity_pf, na.rm = TRUE), .groups = "drop")

## Annual SAF index ------------------------------------------------------------
saf_index <- read_csv(dpath("mean_intensity_SAF_63_73E.csv"), show_col_types = FALSE) %>%
  dplyr::select(matches("^[0-9]{4}$")) %>%
  pivot_longer(cols = everything(), names_to = "Year", values_to = "SAF_mean_cm_s") %>%
  mutate(Year = as.integer(Year))

front_indices <- purrr::reduce(list(pf_park_index, pf_index, saf_index),
                               full_join, by = "Year") %>%
  arrange(Year)

message(sprintf(
  "PF index: front of the year %.2f cm/s, climatological %.2f, r = %.3f",
  mean(front_indices$PF_park_cm_s, na.rm = TRUE),
  mean(front_indices$PF_mean_cm_s, na.rm = TRUE),
  cor(front_indices$PF_park_cm_s, front_indices$PF_mean_cm_s, use = "complete.obs")))

readr::write_csv(front_indices, file.path(out_dir, "front_indices_annual.csv"))

## Time series of the indices --------------------------------------------------
plot_index <- function(df, y, ylab, colour) {
  ggplot(df, aes(x = Year, y = .data[[y]])) +
    geom_point(size = 2, alpha = 0.7, colour = colour) +
    geom_line(colour = colour, linewidth = 1) +
    geom_hline(yintercept = mean(df[[y]], na.rm = TRUE),
               colour = "grey40", linetype = "dashed") +
    scale_x_continuous(breaks = sort(unique(df$Year))) +
    labs(x = "Year", y = ylab) +
    theme_bw() + theme_paper() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(margin = margin(t = 15)))
}

save_fig("pf_mean_intensity.png",
         plot_index(front_indices, "PF_park_cm_s",
                    "Polar Front mean intensity (cm/s)\n(67-72 E)", col_pf),
         width = 10, height = 6)

# The manuscript index, for the comparison drawn by 10_pf_park_intensity_glm.R.
save_fig("pf_mean_intensity_climatological.png",
         plot_index(front_indices, "PF_mean_cm_s",
                    "Polar Front mean intensity (cm/s)\n(67-72 E, climatological contour)",
                    col_pf),
         width = 10, height = 6)

save_fig("saf_mean_intensity.png",
         plot_index(front_indices, "SAF_mean_cm_s",
                    "Subantarctic Front \nmean intensity (cm/s)\n(63-73 E)", col_saf),
         width = 10, height = 6)
