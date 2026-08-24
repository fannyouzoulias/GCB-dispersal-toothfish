################################################################################
# 04_retention_recruitment.R
#
# Larval retention: the proportion of the eggs released that ends the 18-week pelagic
# phase inside a recruitment area. Everything is weighted by the number of eggs released
#
#   - retention by year and release week, with the annual series (main figure)
#   - inter-weekly variability over release weeks 22-30 (appendix)
#   - contribution of each recruitment area, and of each spawning area (main figure)
#   - sensitivity to simulation duration and to the depth range defining a
#     suitable recruitment area (appendix)
#
# Writes retention_annual.csv and recruited_annual.csv, both used by
# 06_retention_front_glm.R.
#
# Requires: 00_setup.R, 01_load_trajectories.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

release_weeks <- 24:28

## Recruitment area ------------------------------------------------------------
# A particle recruits if its final position falls inside any of the recruitment
# sectors (north, east, south, west) or the Skiff bank
recruitment_all <- recruitment_zones %>%
  st_make_valid() %>%
  summarise(geometry = st_union(geometry))

## Release and end-of-simulation positions -------------------------------------
df_all <- df_traj %>% mutate(Week_release = Week - 1)

# release_positions.csv already carries Week_release (= Week - 1, the spawning
# week) and Eggs_Released, one row per particle.
releases <- df_release %>%
  dplyr::select(Unique_ID, Year, Week_release, Eggs_Released, Release_Date = Date)

final_positions <- df_all %>%
  group_by(Unique_ID) %>%
  filter(Date == max(Date, na.rm = TRUE)) %>%
  ungroup()

# Denominator: all eggs released, by year and release week.
released_by_week <- releases %>%
  group_by(Year, Week_release) %>%
  summarise(Total_Eggs_Released = sum(Eggs_Released, na.rm = TRUE), .groups = "drop")

## Retention by year and release week ------------------------------------------
recruited_by_week <- final_positions %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_join(recruitment_all, join = st_within, left = FALSE) %>%
  st_drop_geometry() %>%
  left_join(dplyr::select(releases, Unique_ID, Eggs_Released), by = "Unique_ID") %>%
  group_by(Year, Week_release) %>%
  summarise(Recruited_Weighted = sum(Eggs_Released, na.rm = TRUE), .groups = "drop")

retention_week <- recruited_by_week %>%
  right_join(released_by_week, by = c("Year", "Week_release")) %>%
  mutate(Recruited_Weighted = replace_na(Recruited_Weighted, 0),
         Success_Pct = 100 * Recruited_Weighted / Total_Eggs_Released)

## Main figure: annual series over the year x week heatmap ---------------------
ret_24_28 <- filter(retention_week, Week_release %in% release_weeks)

retention_annual <- ret_24_28 %>%
  group_by(Year) %>%
  summarise(mean_ret = mean(Success_Pct, na.rm = TRUE),
            min_ret  = min(Success_Pct,  na.rm = TRUE),
            max_ret  = max(Success_Pct,  na.rm = TRUE),
            .groups  = "drop") %>%
  arrange(Year)

readr::write_csv(retention_annual, file.path(out_dir, "retention_annual.csv"))

# Upper panel: one point per year, with a locally smoothed trend.
p_annual <- ggplot(retention_annual, aes(x = Year, y = mean_ret)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = mean(retention_annual$mean_ret, na.rm = TRUE),
             linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.2, se = FALSE,
              colour = "darkred", linewidth = 1) +
  scale_x_continuous(breaks = seq(2000, 2023, 1)) +
  labs(y = "Larval retention\nrate (%)") +
  theme_bw() + theme_paper() +
  # the years are labelled once, under the heatmap below
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Lower panel: retention by year and release week, with marginal means.
# "Mean" appears both as an extra column (mean over years) and an extra row
# (mean over weeks).
heat_data <- bind_rows(
  mutate(ret_24_28, Year = as.character(Year),
         Week_release = as.character(Week_release)),
  ret_24_28 %>% group_by(Week_release) %>%
    summarise(Success_Pct = mean(Success_Pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(Year = "Mean", Week_release = as.character(Week_release)),
  ret_24_28 %>% group_by(Year) %>%
    summarise(Success_Pct = mean(Success_Pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(Week_release = "Mean", Year = as.character(Year))
)

year_levels <- c(sort(setdiff(unique(heat_data$Year), "Mean")), "Mean")
heat_data <- heat_data %>%
  mutate(Year = factor(Year, levels = year_levels),
         Week_release = factor(Week_release,
                               levels = c("Mean", as.character(release_weeks))))

p_heat <- ggplot(heat_data, aes(x = Year, y = Week_release, fill = Success_Pct)) +
  geom_tile(colour = "white") +
  scale_y_discrete(limits = rev(levels(heat_data$Week_release))) +
  scale_fill_gradient2(low = "deepskyblue4", mid = "white", high = "red",
                       midpoint = mean(heat_data$Success_Pct, na.rm = TRUE),
                       name = "Larval retention \nrate (%)") +
  geom_text(data = filter(heat_data, Year == "Mean" | Week_release == "Mean"),
            aes(label = paste0(round(Success_Pct, 1), "%")), size = 3) +
  labs(x = "Year", y = "Week of release") +
  theme_bw() + theme_paper() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 15)),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = unit(2, "cm"))

save_fig("retention_annual_and_weekly.png",
         p_annual / p_heat + plot_layout(heights = c(1, 2)),
         width = 15, height = 11)

## Appendix: inter-weekly variability, release weeks 22-30 ---------------------
data_22_30 <- retention_week %>%
  filter(Week_release %in% 22:30) %>%
  mutate(Week_release = factor(Week_release, levels = 22:30))

p_box_week <- ggplot(data_22_30, aes(x = Week_release, y = Success_Pct)) +
  geom_boxplot(width = 0.6) +
  geom_hline(yintercept = mean(data_22_30$Success_Pct, na.rm = TRUE),
             linetype = "dashed", colour = "darkred") +
  labs(x = "Week of release", y = "Larval retention \nrate (%)") +
  theme_bw() + theme_paper() +
  theme(panel.grid = element_blank())

save_fig("retention_by_release_week.png", p_box_week, width = 8, height = 4)

## Where do the recruits come from, and where do they end up? -----------------
# Each particle is tagged with the spawning sector of its release position and
# the recruitment sector of its final position.
tag_zone <- function(df, zones) {
  pts  <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)
  hits <- st_intersects(pts, zones)
  out  <- rep(NA_character_, nrow(df))
  ok   <- lengths(hits) > 0
  out[ok] <- vapply(hits[ok], function(h) zones$zone[h[1]], character(1))
  out
}

df_init <- df_all %>%
  filter(Week_release %in% release_weeks) %>%
  group_by(Unique_ID) %>%
  filter(Date == min(Date, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(dplyr::select(releases, Unique_ID, Eggs_Released), by = "Unique_ID")
df_init$Zone_spawning <- tag_zone(df_init, spawning_zones)

df_end <- filter(final_positions, Week_release %in% release_weeks)
df_end$Zone_recruitment <- tag_zone(df_end, recruitment_zones)

df_join <- df_end %>%
  dplyr::select(Unique_ID, Year, Zone_recruitment) %>%
  inner_join(dplyr::select(df_init, Unique_ID, Year, Eggs_Released, Zone_spawning),
             by = c("Unique_ID", "Year")) %>%
  filter(!is.na(Zone_spawning), !is.na(Eggs_Released))

released_year <- df_init %>%
  group_by(Year) %>%
  summarise(Released_total = sum(Eggs_Released, na.rm = TRUE), .groups = "drop")

zone_share <- function(df, zone_col) {
  df %>%
    filter(!is.na(Zone_recruitment)) %>%
    group_by(Year, Zone = .data[[zone_col]]) %>%
    summarise(Eggs = sum(Eggs_Released, na.rm = TRUE), .groups = "drop") %>%
    left_join(released_year, by = "Year") %>%
    mutate(Percent = 100 * Eggs / Released_total,
           Percent = if_else(is.finite(Percent), Percent, 0),
           Year = factor(as.integer(Year)))
}

# Annual number of retained eggs, summed over the four main recruitment sectors.
# The Skiff bank is excluded.
recruited_annual <- df_join %>%
  filter(!is.na(Zone_recruitment), Zone_recruitment != "skiff") %>%
  group_by(Year) %>%
  summarise(Recruited_total = sum(Eggs_Released, na.rm = TRUE), .groups = "drop") %>%
  left_join(released_year, by = "Year")

readr::write_csv(recruited_annual, file.path(out_dir, "recruited_annual.csv"))

# Order fixes both the colours and the order of the legend entries.
cols_recruitment <- c(east = "cadetblue3", north = "darkseagreen3",
                      south = "lightsalmon", west = "mistyrose3",
                      skiff = "lightgoldenrod1")
labels_recruitment <- c(north = "north-west", east = "north-east",
                        south = "south", west = "west", skiff = "skiff")
cols_spawning <- c(north = "#1f78b4", south = "palegreen3", west = "lightcoral")

zone_lines <- function(df, cols, labels, legend_title) {
  ggplot(df, aes(x = Year, y = Percent, colour = Zone, group = Zone)) +
    geom_line(linewidth = 2) +
    scale_colour_manual(values = cols, labels = labels, breaks = names(cols),
                        drop = FALSE, name = legend_title) +
    # two-line label: on one line it is nearly as tall as the panel and gets
    # clipped once room is left above for the A / B labels
    labs(x = "Year", y = "Larval retention\nrate (%)") +
    theme_bw() + theme_paper() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 30, r = 30, b = 10, l = 10))
}

# Panel A: share of the eggs released that reaches each recruitment sector.
# Panel B: origin of the recruited particles.
panel_recruitment <- zone_lines(zone_share(df_join, "Zone_recruitment"),
                                cols_recruitment, labels_recruitment,
                                "Recruitment zone")
panel_spawning    <- zone_lines(zone_share(df_join, "Zone_spawning"),
                                cols_spawning, names(cols_spawning), "Spawning zone")

save_fig("retention_by_zone.png",
         cowplot::plot_grid(panel_recruitment, panel_spawning, nrow = 2,
                            align = "h", labels = c("A", "B"),
                            label_x = 0.01, label_y = 1,
                            hjust = 0, vjust = 1, label_size = 22),
         width = 15, height = 8)

## Sensitivity 1: simulation duration ------------------------------------------
# Same computation, truncating the trajectories at 10, 15 and 18 weeks.
sim_weeks <- c(10, 15, 18)

retention_for_duration <- function(n_weeks) {
  df_all %>%
    inner_join(dplyr::select(releases, Unique_ID, Release_Date), by = "Unique_ID") %>%
    filter(Date >= Release_Date, Date <= Release_Date + weeks(n_weeks)) %>%
    group_by(Unique_ID, Year, Week_release) %>%
    slice_max(Date, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(Week_release %in% release_weeks) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_join(recruitment_all, join = st_within, left = FALSE) %>%
    st_drop_geometry() %>%
    left_join(dplyr::select(releases, Unique_ID, Eggs_Released), by = "Unique_ID") %>%
    group_by(Year, Week_release) %>%
    summarise(Recruited_Weighted = sum(Eggs_Released, na.rm = TRUE), .groups = "drop") %>%
    right_join(filter(released_by_week, Week_release %in% release_weeks),
               by = c("Year", "Week_release")) %>%
    mutate(Recruited_Weighted = replace_na(Recruited_Weighted, 0)) %>%
    group_by(Year) %>%
    summarise(Success_Pct = 100 * mean(Recruited_Weighted / Total_Eggs_Released),
              .groups = "drop") %>%
    mutate(Scenario = paste0(n_weeks, " weeks"))
}

res_duration <- map_dfr(sim_weeks, retention_for_duration) %>%
  mutate(Scenario = factor(Scenario, levels = paste0(sim_weeks, " weeks")))

p_duration <- ggplot(res_duration, aes(factor(Year), Success_Pct,
                                       colour = Scenario, group = Scenario)) +
  geom_line() +
  scale_colour_manual(name = "Simulation time",
                      values = c("10 weeks" = "#1f78b4", "15 weeks" = "#33a02c",
                                 "18 weeks" = "#e31a1c")) +
  labs(x = "Year", y = "Larval recruitment success (%)") +
  theme_bw() + theme_paper() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig("sensitivity_simulation_duration.png", p_duration, width = 15, height = 5)

## Sensitivity 2: depth range of a suitable recruitment area -------------------
# The 150 m scenario is the one used in the paper
depths <- c(150, 200, 300, 500, 700, 1000)

retention_for_depth <- function(dmax) {
  end_24_28 <- filter(final_positions, Week_release %in% release_weeks)

  recruited <- if (dmax == 150) {
    end_24_28 %>%
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
      st_join(recruitment_all, join = st_within, left = FALSE) %>%
      st_drop_geometry()
  } else {
    filter(end_24_28, abs(Depth) <= dmax)
  }

  recruited %>%
    left_join(dplyr::select(releases, Unique_ID, Eggs_Released), by = "Unique_ID") %>%
    filter(!is.na(Eggs_Released)) %>%
    group_by(Year, Week_release) %>%
    summarise(Recruited_Weighted = sum(Eggs_Released, na.rm = TRUE), .groups = "drop") %>%
    right_join(filter(released_by_week, Week_release %in% release_weeks),
               by = c("Year", "Week_release")) %>%
    mutate(Recruited_Weighted = replace_na(Recruited_Weighted, 0)) %>%
    group_by(Year) %>%
    summarise(Success_Pct = 100 * mean(Recruited_Weighted / Total_Eggs_Released),
              .groups = "drop") %>%
    mutate(Scenario = paste0("0-", dmax, " m"))
}

res_depth <- map_dfr(depths, retention_for_depth) %>%
  mutate(Scenario = factor(Scenario, levels = paste0("0-", depths, " m")))

p_depth <- ggplot(res_depth, aes(factor(Year), Success_Pct,
                                 colour = Scenario, group = Scenario)) +
  geom_line() + geom_point(size = 1) +
  labs(x = "Year", y = "Larval recruitment success (%)") +
  theme_bw() + theme_paper() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

save_fig("sensitivity_recruitment_depth.png", p_depth, width = 15, height = 5)
