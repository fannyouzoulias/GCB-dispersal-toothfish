################################################################################
# 05_front_indices.R
#
# Annual intensity indices of the Polar Front (PF) and the Subantarctic Front
# (SAF), from the geostrophic current speed sampled along the mean position of
# each front (DUACS altimetry).
#
#   - map of the front positions and of the sectors used to average them (appendix)
#   - annual PF index, averaged over 67-72 E and south of 51 S
#   - annual SAF index, averaged over 63-73 E
#
# Writes front_indices_annual.csv, used by 06_retention_front_glm.R.
#
# Requires: 00_setup.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

col_pf  <- "#4575b4"
col_saf <- "#1b9e77"

## Front points, long format ---------------------------------------------------
# Each row is a point along the mean front; the columns named 2000...2023 hold
# the mean geostrophic speed at that point for that year.
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

## Appendix figure: front positions and averaging sectors ----------------------
# The two boxes are the sectors over which the yearly speeds are averaged: the
# PF where it runs along the northern plateau, the SAF further north.
box_pf  <- tibble(xmin = 67, xmax = 72, ymin = -51, ymax = -48)
box_saf <- tibble(xmin = 63, xmax = 73, ymin = -47, ymax = -44)

p_fronts <- ggplot() +
  geom_sf(data = land_sf, fill = "darkgrey", inherit.aes = FALSE) +
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

## Annual PF index -------------------------------------------------------------
pf_index <- pf_long %>%
  filter(lon >= 67, lon <= 72, lat >= -51) %>%
  group_by(Year) %>%
  summarise(PF_mean_cm_s = mean(intensity_pf, na.rm = TRUE), .groups = "drop")

## Annual SAF index ------------------------------------------------------------
saf_index <- read_csv(dpath("mean_intensity_SAF_63_73E.csv"), show_col_types = FALSE) %>%
  dplyr::select(matches("^[0-9]{4}$")) %>%
  pivot_longer(cols = everything(), names_to = "Year", values_to = "SAF_mean_cm_s") %>%
  mutate(Year = as.integer(Year))

front_indices <- full_join(pf_index, saf_index, by = "Year") %>% arrange(Year)

readr::write_csv(front_indices, file.path(out_dir, "front_indices_annual.csv"))

## Time series of the two indices ----------------------------------------------
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
         plot_index(front_indices, "PF_mean_cm_s",
                    "Polar Front mean intensity (cm/s)\n(67-72 E)", col_pf),
         width = 10, height = 6)

save_fig("saf_mean_intensity.png",
         plot_index(front_indices, "SAF_mean_cm_s",
                    "Subantarctic Front \nmean intensity (cm/s)\n(63-73 E)", col_saf),
         width = 10, height = 6)
