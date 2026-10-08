################################################################################
# 14_metrics_correlation.R
#
# The annual metrics quoted across the Results, in one table, with the two
# supplementary figures built from it:
#
#   PART 2  pathway indices: share of the particles passing through the NW, NE
#           and southern corridors, and trapped in the western recirculation
#           cell (unweighted particle counts, as 13)
#   PART 3  retention by recruitment area and contribution of each spawning
#           area (weighted by egg production, same rules as 04)
#   PART 4  connectivity: where the eggs of each spawning area end up, over
#           2000-2023 and (4b) before and after 2010
#   PART 5 Pearson correlations among all annual metrics, 2000-2023
#
# Every correlation given in the text is one cell of the matrix of PART 5: PF-
# associated jet intensity against the corridors and against retention in the
# NW area, trapping against the southern corridor, SAM against the PF-
# associated jet intensity and against retention.
#
# Requires: 00_setup.R, the trajectories of 01_load_trajectories.R (read from
#           outputs/ if they are not in memory), the CSVs of 04 and 05, and the
#           SAM series in `sam_dir` (config.R). Needs data.table.
#
# Author: Fanny Ouzoulias
# Date:   2026-10-06
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(data.table)
})

## ---------------------------------------------------------------------------
## PART 0. Settings
## ---------------------------------------------------------------------------
if (!exists("out_dir")) out_dir <- "outputs"
rev_dir <- file.path(out_dir, "revision")
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)

release_weeks <- 24:28

# Boxes, in degrees: lon_min, lon_max, lat_min, lat_max (Supplementary Methods)
box_nw     <- c(68.5, 69.5, -49.25, -47)   # NW corridor, northern shelf pathway
box_ne     <- c(69.5, 71,   -49.25, -47)   # NE corridor, northern shelf pathway
box_south  <- c(69,   71,   -51,    -49.5) # southern PF-associated corridor
box_recirc <- c(63,   67,   -50,    -48)   # western recirculation cell

in_box <- function(lon, lat, b)
  lon >= b[1] & lon <= b[2] & lat >= b[3] & lat <= b[4]

## ---------------------------------------------------------------------------
## PART 1. Trajectories
## ---------------------------------------------------------------------------
if (!exists("df_traj")) df_traj <- readRDS(file.path(out_dir, "trajectories.rds"))
if (!exists("df_release")) df_release <- readRDS(file.path(out_dir, "release_positions.rds"))
if (!exists("recruitment_zones"))
  recruitment_zones <- st_read(dpath("zones_recruitment.geojson"), quiet = TRUE)
if (!exists("spawning_zones"))
  spawning_zones <- st_read(dpath("zones_spawning.geojson"), quiet = TRUE)

DT <- as.data.table(df_traj)[, .(Unique_ID, Year, Date, Longitude, Latitude, Week)]
DT <- DT[(Week - 1) %in% release_weeks]            # Week_release = Week - 1, as in 04
setkey(DT, Unique_ID, Date)

## ---------------------------------------------------------------------------
## PART 2. Pathway indices
## ---------------------------------------------------------------------------
# A corridor is crossed when any position of the particle falls in its box.
# Trapped: once inside the cell, every later position is inside it too (13).
flags <- DT[, {
  in_recirc <- in_box(Longitude, Latitude, box_recirc)
  first_in  <- match(TRUE, in_recirc)
  .(nw      = any(in_box(Longitude, Latitude, box_nw)),
    ne      = any(in_box(Longitude, Latitude, box_ne)),
    south   = any(in_box(Longitude, Latitude, box_south)),
    trapped = !is.na(first_in) && all(in_recirc[first_in:.N]))
}, by = .(Year, Unique_ID)]

pathways <- flags[, .(Pct_cross_NW    = 100 * mean(nw),
                      Pct_cross_NE    = 100 * mean(ne),
                      Pct_cross_south = 100 * mean(south),
                      Pct_trapped     = 100 * mean(trapped)), by = Year]

## ---------------------------------------------------------------------------
## PART 3. Retention by recruitment area and by spawning area
## ---------------------------------------------------------------------------
# Same rules as 04: a particle is tagged with the spawning area of its release
# position and the recruitment area of its final position; the particles
# released outside the three spawning areas are left out of the shares, and
# every share is a percentage of ALL the eggs released that year.
tag_zone <- function(lon, lat, zones) {
  pts  <- st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326)
  hits <- st_intersects(pts, zones)
  out  <- rep(NA_character_, length(lon))
  ok   <- lengths(hits) > 0
  out[ok] <- vapply(hits[ok], function(h) zones$zone[h[1]], character(1))
  out
}

final <- DT[, .SD[.N], by = Unique_ID]
rm(DT); gc(verbose = FALSE)
final[, Zone_recruitment := tag_zone(Longitude, Latitude, recruitment_zones)]

rel <- as.data.table(df_release)[Week_release %in% release_weeks & !is.na(Eggs_Released)]
rel[, Zone_spawning := tag_zone(Longitude, Latitude, spawning_zones)]

parts <- merge(final[, .(Unique_ID, Zone_recruitment)],
               rel[, .(Unique_ID, Year, Eggs_Released, Zone_spawning)],
               by = "Unique_ID")
released_year <- parts[, .(Released_total = sum(Eggs_Released)), by = Year]
parts <- parts[!is.na(Zone_spawning)]

share_by <- function(col, prefix) {
  s <- parts[!is.na(Zone_recruitment), .(Eggs = sum(Eggs_Released)),
             by = c("Year", col)]
  s <- merge(s, released_year, by = "Year")[, Percent := 100 * Eggs / Released_total]
  s <- dcast(s, as.formula(paste("Year ~", col)), value.var = "Percent", fill = 0)
  setnames(s, -1, paste0(prefix, names(s)[-1]))
  s
}
by_recruitment <- share_by("Zone_recruitment", "Ret_")     # Ret_north = NW, Ret_east = NE
by_spawning    <- share_by("Zone_spawning",    "From_")

## ---------------------------------------------------------------------------
## PART 4. Connectivity between spawning and recruitment areas
## ---------------------------------------------------------------------------
# Of the eggs released in each spawning area, the share ending in each
# recruitment area: pooled over 2000-2023 (the release is the same every year,
# so this is also the mean of the annual shares), with the range of the annual
# totals.
lab_rec <- c(north = "NW", east = "NE", south = "South", west = "West", skiff = "Skiff")
lab_spw <- c(north = "North", west = "West", south = "South")
cols_rec <- c(NW = "darkseagreen3", NE = "cadetblue3", South = "lightsalmon",
              West = "mistyrose3", Skiff = "lightgoldenrod1")   # as in 04

spawned <- parts[, .(Released = sum(Eggs_Released)), by = Zone_spawning]
spawned[, Share_of_eggs := 100 * Released / sum(released_year$Released_total)]

connectivity <- parts[!is.na(Zone_recruitment), .(Eggs = sum(Eggs_Released)),
                      by = .(Zone_spawning, Zone_recruitment)]
connectivity <- merge(connectivity, spawned, by = "Zone_spawning")
connectivity[, Percent := 100 * Eggs / Released]

annual_total <- merge(
  parts[, .(Released = sum(Eggs_Released)), by = .(Year, Zone_spawning)],
  parts[!is.na(Zone_recruitment), .(Retained = sum(Eggs_Released)),
        by = .(Year, Zone_spawning)],
  by = c("Year", "Zone_spawning"), all.x = TRUE)
annual_total[is.na(Retained), Retained := 0]
total_range <- annual_total[, .(Total_min = min(100 * Retained / Released),
                                Total_max = max(100 * Retained / Released)),
                            by = Zone_spawning]

connectivity <- merge(connectivity, total_range, by = "Zone_spawning")
connectivity[, Total := sum(Percent), by = Zone_spawning]
fwrite(connectivity[, .(Zone_spawning, Zone_recruitment, Percent, Total, Total_min,
                        Total_max, Share_of_eggs)][order(Zone_spawning, -Percent)],
       file.path(rev_dir, "connectivity_spawning_recruitment.csv"))

cat("\n======== CONNECTIVITY (% of the eggs of each spawning area) ========\n")
print(as.data.frame(dcast(connectivity, Zone_spawning + Share_of_eggs + Total ~ Zone_recruitment,
                          value.var = "Percent")), digits = 3)

conn_plot <- connectivity %>%
  as_tibble() %>%
  mutate(Recruitment = factor(lab_rec[Zone_recruitment], levels = rev(lab_rec)),
         Spawning = factor(
           sprintf("%s\n(%.0f%% of eggs)", lab_spw[Zone_spawning], Share_of_eggs),
           levels = sprintf("%s\n(%.0f%% of eggs)", lab_spw[names(lab_spw)],
                            spawned$Share_of_eggs[match(names(lab_spw), spawned$Zone_spawning)])))

p_conn <- ggplot(conn_plot, aes(x = Percent, y = fct_rev(Spawning))) +
  geom_col(aes(fill = Recruitment), colour = "white", linewidth = 0.8, width = 0.6) +
  # segments under 3 % are left unlabelled (kept in the data so that the
  # others stay stacked at the right place)
  geom_text(aes(label = if_else(Percent >= 3, sprintf("%.0f", Percent), ""),
                group = Recruitment),
            position = position_stack(vjust = 0.5), size = 6, colour = "grey15") +
  # total retained, with the range of the annual totals over 2000-2023
  geom_text(data = distinct(conn_plot, Spawning, Total, Total_min, Total_max),
            aes(x = Total, y = fct_rev(Spawning),
                label = sprintf("%.0f%%  (%.0f–%.0f)", Total, Total_min, Total_max)),
            inherit.aes = FALSE, hjust = -0.08, size = 6, colour = "grey15") +
  scale_fill_manual(values = cols_rec, breaks = unname(lab_rec),
                    name = "Recruitment area") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = "Eggs retained (% of the eggs released in the spawning area)",
       y = "Spawning area") +
  theme_bw() + theme_paper() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_text(margin = margin(r = 15)))

ggsave(file.path(rev_fig, "connectivity_spawning_recruitment.png"), p_conn,
       width = 12, height = 6, dpi = 300)

## ---------------------------------------------------------------------------
## PART 5. Annual table and correlations
## ---------------------------------------------------------------------------
sam <- read_csv(file.path(sam_dir, "SAM_annual_2000_2023.csv"),
                show_col_types = FALSE) %>%
  rename(Year = year) %>%
  left_join(read_csv(file.path(sam_dir, "SAM_seasonal_2000_2023.csv"),
                     show_col_types = FALSE) %>%
              rename(Year = year) %>%
              filter(season %in% c("JJA", "SON")) %>%
              pivot_wider(names_from = season, values_from = SAM_season,
                          names_prefix = "SAM_"),
            by = "Year")

dat <- read_csv(file.path(out_dir, "front_indices_annual.csv"),
                show_col_types = FALSE) %>%
  dplyr::select(Year, PF_park_cm_s, SAF_park_cm_s) %>%
  left_join(sam, by = "Year") %>%
  left_join(as_tibble(pathways), by = "Year") %>%
  left_join(read_csv(file.path(out_dir, "recruited_annual.csv"),
                     show_col_types = FALSE) %>%
              transmute(Year, logR = log(Recruited_total)), by = "Year") %>%
  left_join(as_tibble(by_recruitment), by = "Year") %>%
  left_join(as_tibble(by_spawning), by = "Year") %>%
  filter(Year %in% 2000:2023)

# Order and labels of the matrix: circulation, pathways, retention, areas.
metric_labels <- c(
  PF_park_cm_s    = "PF-associated jet intensity",
  SAF_park_cm_s   = "SAF intensity",
  SAM_annual      = "SAM, annual",
  SAM_JJA         = "SAM, winter (JJA)",
  SAM_SON         = "SAM, spring (SON)",
  Pct_cross_south = "Southern corridor (%)",
  Pct_cross_NW    = "NW corridor (%)",
  Pct_cross_NE    = "NE corridor (%)",
  Pct_trapped     = "Trapped in western cell (%)",
  logR            = "Larval retention (log R)",
  Ret_north       = "Retained in NW area (%)",
  Ret_east        = "Retained in NE area (%)",
  Ret_south       = "Retained in South area (%)",
  Ret_west        = "Retained in West area (%)",
  Ret_skiff       = "Retained on Skiff Bank (%)",
  From_west       = "From West spawning area (%)",
  From_north      = "From North spawning area (%)",
  From_south      = "From South spawning area (%)"
)
metrics <- names(metric_labels)
stopifnot(nrow(dat) == 24, all(metrics %in% names(dat)),
          !anyNA(dplyr::select(dat, all_of(metrics))))

write_csv(dat, file.path(rev_dir, "annual_metrics.csv"))

cor_tbl <- expand_grid(i = seq_along(metrics), j = seq_along(metrics)) %>%
  filter(i > j) %>%
  mutate(x = metrics[j], y = metrics[i]) %>%
  mutate(test = map2(x, y, ~ cor.test(dat[[.x]], dat[[.y]])),
         r = map_dbl(test, ~ unname(.x$estimate)),
         p = map_dbl(test, "p.value"),
         n = nrow(dat)) %>%
  dplyr::select(-test)

write_csv(dplyr::select(cor_tbl, x, y, r, p, n),
          file.path(rev_dir, "metrics_correlation.csv"))

cat("\n======== CORRELATIONS QUOTED IN THE TEXT ========\n")
quoted <- tribble(~x, ~y,
  "PF_park_cm_s", "Pct_cross_south", "PF_park_cm_s", "Pct_cross_NW",
  "Pct_cross_south", "Pct_trapped",  "PF_park_cm_s", "Ret_north",
  "PF_park_cm_s", "logR",            "PF_park_cm_s", "SAM_JJA",
  "PF_park_cm_s", "SAM_SON",         "SAM_annual",   "logR")
print(as.data.frame(inner_join(quoted, cor_tbl, by = c("x", "y")) %>%
                      dplyr::select(x, y, r, p)), digits = 2)
cat(sprintf("\n%d of the %d pairs have p < 0.05 (no correction for multiple testing).\n",
            sum(cor_tbl$p < 0.05), nrow(cor_tbl)))

# Lower triangle; r in every cell, in bold where p < 0.05.
p_cor <- cor_tbl %>%
  mutate(X = factor(metric_labels[x], levels = metric_labels),
         Y = factor(metric_labels[y], levels = rev(metric_labels))) %>%
  ggplot(aes(X, Y, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r),
                fontface = if_else(p < 0.05, "bold", "plain")),
            colour = "black", size = 4.2) +
  # diverging, neutral at zero, readable with a colour-vision deficiency; the
  # two ends are kept light enough for the black figures to stay legible
  scale_fill_gradient2(low = "#4393C3", mid = "#F7F7F7", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.85, 0.8),
        legend.key.height = unit(1.2, "cm"))

ggsave(file.path(rev_fig, "metrics_correlation.png"), p_cor,
       width = 13, height = 12, dpi = 300)
