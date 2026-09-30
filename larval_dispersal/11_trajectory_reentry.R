################################################################################
# 11_trajectory_reentry.R
#
# Why does retention increase with the pelagic larval duration?
#
# Reviewer question: "Is it possible that this increased retention rate at
# longer time-scales is because particles have more time to exit the domain of
# interest and be transported back into it?"
#
# The answer needs a decomposition, not just a re-entry count. At any duration
# T, a particle sitting inside recruitment habitat got there in one of two ways:
#   FIRST STAY  it arrived at some point before T and has not left since;
#   RETURNED    it arrived, left, and came back -- the reviewer's mechanism.
# Retention at T is the sum of the two, so the rise of retention between 10 and
# 18 weeks can be split into the part due to particles still arriving for the
# first time and the part due to particles coming back. That split is the answer
# to the question, and it is what reentry_retention_decomposition.csv holds.
#
# Definitions, all on particles released in weeks 24-28, 2000-2023:
#   - a SPELL is a maximal run of consecutive days spent inside any recruitment
#     area (the union of the four sectors and the Skiff bank, as in
#     04_retention_recruitment.R, without a depth filter);
#   - a particle RE-ENTERS when it has two or more spells: it entered, left
#     every recruitment area, and entered one again before the end of the run;
#   - the amplitude of an excursion is measured three ways, see below.
#
# The per-particle spells are saved to reentry_spells.rds. They are the
# expensive part (40 million point-in-polygon tests), and every table below is
# derived from them, so re-running an analysis on a different duration, weighting
# or subset does not need the trajectories again.
#
# Three readings of "distance travelled outside recruitment habitat", over the
# days spent outside between the first and the last spell:
#   dist_path_km     length of the track swum outside, cumulated along the
#                    trajectory, maximised over the excursions
#   dist_exit_km     greatest displacement from the point where it left
#   dist_habitat_km  greatest distance to the nearest recruitment area
# The appendix sentence quotes a median of ~19 km with 74% under 50 km and 16%
# beyond 100 km; dist_habitat_km is far smaller than that (the union of the
# sectors is one 220 x 180 km block, so anything that steps outside stays close
# to it), so the published figures describe the amplitude of the excursion, not
# the distance to the habitat.
#
# Writes, in outputs/revision/, reentry_spells.rds,
# reentry_retention_decomposition.csv, reentry_by_duration.csv,
# reentry_summary.csv, reentry_distance_definitions.csv, reentry_by_year.csv
# and reentry_particles.csv; figures go to outputs/figures/revision/.
#
# Requires: 00_setup.R, and the trajectories of 01_load_trajectories.R (read
#           from outputs/ if they are not in memory). Needs the data.table and
#           terra packages on top of those of 00_setup.R.
#
# Author: Fanny Ouzoulias
# Date:   2026-09-09
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(data.table)     # terra is only called with terra::, never attached
})

# data.table exports shift() and between(); this script means the data.table
# ones (lag within a group).
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate,
                             data.table::shift, data.table::between,
                             .quiet = TRUE)

## ---------------------------------------------------------------------------
## PART 0. Inputs
## ---------------------------------------------------------------------------
if (!exists("out_dir")) out_dir <- "outputs"
rev_dir <- file.path(out_dir, "revision")
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)
if (!exists("theme_paper")) theme_paper <- function() theme_bw(base_size = 13)

release_weeks <- 24:28          # spawning peak, as in 02 and 04
sim_weeks     <- c(10, 15, 18)  # the durations of the sensitivity analysis in 04
f_spells      <- file.path(rev_dir, "reentry_spells.rds")

if (!exists("recruitment_zones"))
  recruitment_zones <- st_read(dpath("zones_recruitment.geojson"), quiet = TRUE)

zones_valid     <- st_make_valid(recruitment_zones)
recruitment_all <- st_sf(geometry = st_union(st_geometry(zones_valid)))
bb <- st_bbox(recruitment_all)

if (!exists("df_release")) df_release <- readRDS(file.path(out_dir, "release_positions.rds"))

## ---------------------------------------------------------------------------
## PART 1. Helpers
## ---------------------------------------------------------------------------
# One point-in-polygon pass against the individual sectors gives both flags:
# inside any recruitment area, and inside one other than the Skiff bank. Only
# the points inside the bounding box are tested; the rest are outside by
# construction.
#
# terra rather than sf here: sf would build one POINT geometry per position and
# test it with the spherical s2 engine, which for tens of millions of positions
# costs far more than the test itself. terra keeps the coordinates in a matrix
# and uses prepared planar geometries -- the same answer for a topological test
# at this scale, about an order of magnitude faster.
v_zones  <- terra::vect(zones_valid)
not_skiff <- zones_valid$zone != "skiff"

inside_flags <- function(lon, lat) {
  in_zone <- logical(length(lon)); in_nsk <- logical(length(lon))
  cand <- which(lon >= bb[["xmin"]] & lon <= bb[["xmax"]] &
                lat >= bb[["ymin"]] & lat <= bb[["ymax"]])
  if (!length(cand)) return(list(in_zone = in_zone, in_nsk = in_nsk))
  pts <- terra::vect(cbind(lon[cand], lat[cand]), type = "points",
                     crs = "EPSG:4326")
  # "within", not "intersects": 04_retention_recruitment.R tests the final
  # positions with st_within, which excludes points lying exactly on a boundary.
  m <- terra::relate(pts, v_zones, "within")       # points x sectors
  in_zone[cand] <- rowSums(m) > 0
  in_nsk[cand]  <- rowSums(m[, not_skiff, drop = FALSE]) > 0
  list(in_zone = in_zone, in_nsk = in_nsk)
}

dist_to_habitat_km <- function(lon, lat) {
  pts <- terra::vect(cbind(lon, lat), type = "points", crs = "EPSG:4326")
  as.numeric(terra::distance(pts, terra::vect(recruitment_all))) / 1000
}

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6378.137; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
       cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

## ---------------------------------------------------------------------------
## PART 2. Spells inside recruitment habitat, one year at a time
## ---------------------------------------------------------------------------
# 40 million daily positions over 24 years. This is the only pass over the
# trajectories; its output is cached.

spells_year <- function(y) {
  dt <- DT[.(y)]                       # keyed subset, already ordered
  fl <- inside_flags(dt$Longitude, dt$Latitude)
  dt[, `:=`(in_zone = fl$in_zone, in_nsk = fl$in_nsk)]
  dt[, day := as.numeric(Date - min(Date)), by = .(Unique_ID, Year)]

  # spells = maximal runs of consecutive days inside habitat
  dt[, runid := rleid(in_zone), by = .(Unique_ID, Year)]
  sp <- dt[in_zone == TRUE, .(start_day = min(day), end_day = max(day)),
           by = .(Unique_ID, Year, runid)]
  setorder(sp, Unique_ID, Year, start_day)
  sp[, spell_id := seq_len(.N), by = .(Unique_ID, Year)]
  sp[, runid := NULL]

  # per-particle constants, including the particles that never entered.
  # A spell starts on a day inside habitat whose previous day was outside, so
  # counting those starts is the same as counting runs, without calling rle()
  # once per particle.
  info <- dt[, .(t_end = max(day),
                 Week_release = Week_release[1],
                 n_spells_nsk = sum(in_nsk & !shift(in_nsk, fill = FALSE))),
             by = .(Unique_ID, Year)]

  # excursions = days outside between the first and the last spell
  dt[, `:=`(in_before = cumsum(in_zone),
            in_after  = rev(cumsum(rev(in_zone)))), by = .(Unique_ID, Year)]
  dt[, excursion := !in_zone & in_before > 0 & in_after > 0]

  reent <- sp[, .N, by = .(Unique_ID, Year)][N >= 2, .(Unique_ID, Year)]
  sub   <- dt[reent, on = .(Unique_ID, Year)]

  if (nrow(sub)) {
    sub[, xrun := rleid(excursion), by = .(Unique_ID, Year)]
    sub[, `:=`(prev_lon = shift(Longitude), prev_lat = shift(Latitude)),
        by = .(Unique_ID, Year)]
    ex <- sub[excursion == TRUE]
    exits <- ex[, .(exit_lon = prev_lon[1], exit_lat = prev_lat[1]),
                by = .(Unique_ID, Year, xrun)]
    ex <- exits[ex, on = .(Unique_ID, Year, xrun)]

    ex[, dist_habitat := dist_to_habitat_km(Longitude, Latitude)]
    ex[, dist_exit    := haversine_km(exit_lon, exit_lat, Longitude, Latitude)]
    ex[, step := c(0, haversine_km(Longitude[-.N], Latitude[-.N],
                                   Longitude[-1],  Latitude[-1])),
       by = .(Unique_ID, Year, xrun)]
    ex[, dist_path := cumsum(step), by = .(Unique_ID, Year, xrun)]

    # One row per EXCURSION. The appendix sentence is phrased per travel ("the
    # travels ... have a median amplitude"), not per particle, and a particle
    # that leaves three times contributes three travels, so the two units give
    # different medians; both are reported.
    per_excursion <- ex[, .(n_days          = .N,
                            dist_path_km    = max(dist_path),
                            dist_exit_km    = max(dist_exit),
                            dist_habitat_km = max(dist_habitat)),
                        by = .(Unique_ID, Year, xrun)]

    dmax <- per_excursion[, .(n_out           = sum(n_days),
                              n_excursions    = .N,
                              dist_path_km    = max(dist_path_km),
                              dist_exit_km    = max(dist_exit_km),
                              dist_habitat_km = max(dist_habitat_km)),
                          by = .(Unique_ID, Year)]
    info <- dmax[info, on = .(Unique_ID, Year)]
    rm(sub, ex, exits, dmax)
  } else {
    per_excursion <- data.table()
    info[, `:=`(n_out = NA_integer_, n_excursions = NA_integer_,
                dist_path_km = NA_real_, dist_exit_km = NA_real_,
                dist_habitat_km = NA_real_)]
  }

  rm(dt); gc(verbose = FALSE)
  message("  ", y, ": ", nrow(info), " particles, ",
          nrow(reent), " re-entering")
  list(spells = sp, info = info, excursions = per_excursion)
}

if (file.exists(f_spells)) {
  message("Reusing ", f_spells)
  cached <- readRDS(f_spells)
  spells <- cached$spells; info <- cached$info; excursions <- cached$excursions
} else {
  if (!exists("df_traj")) df_traj <- readRDS(file.path(out_dir, "trajectories.rds"))

  # Convert and subset ONCE, keyed on Year, rather than re-converting the whole
  # 72-million-row table inside every yearly call. Only the columns the analysis
  # needs are carried over.
  DT <- as.data.table(df_traj)[, .(Unique_ID, Date, Longitude, Latitude,
                                   Week, Year)]
  DT[, Week_release := Week - 1]
  DT <- DT[Week_release %in% release_weeks]
  DT[, Week := NULL]
  setkey(DT, Year, Unique_ID, Date)
  gc(verbose = FALSE)

  years <- sort(unique(DT$Year))
  message("Spell extraction over ", length(years), " years:")
  res    <- lapply(years, spells_year)
  spells     <- rbindlist(lapply(res, `[[`, "spells"))
  info       <- rbindlist(lapply(res, `[[`, "info"))
  excursions <- rbindlist(lapply(res, `[[`, "excursions"), fill = TRUE)
  saveRDS(list(spells = spells, info = info, excursions = excursions), f_spells)
  rm(res, DT); gc(verbose = FALSE)
}

# Egg production weight of each particle, so the decomposition can be read on
# the same scale as the retention of 04_retention_recruitment.R.
eggs <- as.data.table(df_release)[, .(Unique_ID, Year, Eggs_Released)]
info <- eggs[info, on = .(Unique_ID, Year)]
info[is.na(Eggs_Released), Eggs_Released := 0]

stopifnot(!anyDuplicated(info[, .(Unique_ID, Year)]))
years <- sort(unique(info$Year))

## ---------------------------------------------------------------------------
## PART 3. Why retention rises with duration
## ---------------------------------------------------------------------------
# At each cutoff, a particle inside habitat is on its FIRST spell (it arrived
# and has not left) or on a later one (it left and came back).

cutoffs <- seq(7, 18 * 7, by = 7)
n_tot   <- nrow(info)
eggs_tot <- sum(info$Eggs_Released)

decompose <- function(T_days) {
  ins <- spells[start_day <= T_days & end_day >= T_days]     # one spell at most
  ins <- info[, .(Unique_ID, Year, Eggs_Released)][ins, on = .(Unique_ID, Year)]
  ins[, returned := spell_id >= 2]
  data.table(
    weeks          = T_days / 7,
    day            = T_days,
    pct_inside     = 100 * nrow(ins) / n_tot,
    pct_first_stay = 100 * sum(!ins$returned) / n_tot,
    pct_returned   = 100 * sum(ins$returned)  / n_tot,
    wpct_inside     = 100 * sum(ins$Eggs_Released) / eggs_tot,
    wpct_first_stay = 100 * sum(ins$Eggs_Released[!ins$returned]) / eggs_tot,
    wpct_returned   = 100 * sum(ins$Eggs_Released[ins$returned])  / eggs_tot
  )
}

decomp <- rbindlist(lapply(cutoffs, decompose))
fwrite(decomp, file.path(rev_dir, "reentry_retention_decomposition.csv"))

cat("\n============ WHY RETENTION RISES WITH DURATION ============\n")
cat("Share of released particles inside recruitment habitat at each duration,\n")
cat("split by how they got there (unweighted, then weighted by egg production).\n\n")
print(as.data.frame(decomp[weeks %in% sim_weeks]), digits = 3)

gain <- function(v1, v2, col_all, col_ret) {
  a <- decomp[weeks == v1]; b <- decomp[weeks == v2]
  d_all <- b[[col_all]] - a[[col_all]]
  d_ret <- b[[col_ret]] - a[[col_ret]]
  c(total = d_all, returned = d_ret, share = 100 * d_ret / d_all)
}

g_u <- gain(10, 18, "pct_inside",  "pct_returned")
g_w <- gain(10, 18, "wpct_inside", "wpct_returned")

cat(sprintf("\nFrom 10 to 18 weeks, retention rises by %.2f points (unweighted);\n",
            g_u["total"]))
cat(sprintf("  %.2f of those points are particles that had left and come back,\n",
            g_u["returned"]))
cat(sprintf("  i.e. %.0f %% of the increase. The rest is particles still arriving\n",
            g_u["share"]))
cat("  in recruitment habitat for the first time.\n")
cat(sprintf("\nWeighted by egg production: +%.2f points, of which %.2f from returns (%.0f %%).\n",
            g_w["total"], g_w["returned"], g_w["share"]))

## ---------------------------------------------------------------------------
## PART 4. Re-entry counts and excursion amplitude
## ---------------------------------------------------------------------------
n_spells <- spells[, .(n_spells = .N), by = .(Unique_ID, Year)]
info <- n_spells[info, on = .(Unique_ID, Year)]
info[is.na(n_spells), n_spells := 0L]
info[, `:=`(entered = n_spells >= 1, reentered = n_spells >= 2)]

re <- info[reentered == TRUE]
fwrite(re, file.path(rev_dir, "reentry_particles.csv"))

# Re-entry as a function of the duration allowed: a particle counts as
# re-entering by T if two spells have begun by then.
reentry_by_duration <- rbindlist(lapply(cutoffs, function(T_days) {
  k <- spells[start_day <= T_days, .N, by = .(Unique_ID, Year)][N >= 2, .N]
  data.table(weeks = T_days / 7, pct_reentered = 100 * k / n_tot)
}))
fwrite(reentry_by_duration, file.path(rev_dir, "reentry_by_duration.csv"))

dist_stats <- function(x, label) tibble(
  definition       = label,
  median_km        = median(x, na.rm = TRUE),
  mean_km          = mean(x, na.rm = TRUE),
  pct_within_50km  = 100 * mean(x <= 50,  na.rm = TRUE),
  pct_beyond_100km = 100 * mean(x > 100, na.rm = TRUE)
)

distance_tbl <- bind_rows(
  dist_stats(re$dist_path_km,    "per particle: path travelled outside"),
  dist_stats(re$dist_exit_km,    "per particle: displacement from the exit point"),
  dist_stats(re$dist_habitat_km, "per particle: distance to nearest area"),
  dist_stats(excursions$dist_path_km,    "per travel: path travelled outside"),
  dist_stats(excursions$dist_exit_km,    "per travel: displacement from the exit point"),
  dist_stats(excursions$dist_habitat_km, "per travel: distance to nearest area")
)
write_csv(distance_tbl, file.path(rev_dir, "reentry_distance_definitions.csv"))

summary_tbl <- tibble(
  n_particles            = n_tot,
  n_entered              = sum(info$entered),
  n_reentered            = nrow(re),
  pct_reentered          = 100 * nrow(re) / n_tot,
  pct_of_entrants        = 100 * nrow(re) / sum(info$entered),
  median_days_outside    = median(re$n_out, na.rm = TRUE),
  pct_reentered_no_skiff = 100 * mean(info$n_spells_nsk >= 2)
)
write_csv(summary_tbl, file.path(rev_dir, "reentry_summary.csv"))

cat("\n================ RE-ENTRY ================\n")
cat(sprintf("Particles released in weeks %d-%d:          %s\n",
            min(release_weeks), max(release_weeks), format(n_tot, big.mark = " ")))
cat(sprintf("  entered a recruitment area at least once: %s (%.1f %%)\n",
            format(summary_tbl$n_entered, big.mark = " "),
            100 * summary_tbl$n_entered / n_tot))
cat(sprintf("  entered, left, and re-entered:            %s (%.1f %% of released, %.1f %% of entrants)\n",
            format(summary_tbl$n_reentered, big.mark = " "),
            summary_tbl$pct_reentered, summary_tbl$pct_of_entrants))
cat(sprintf("  without the Skiff bank:                   %.1f %%\n",
            summary_tbl$pct_reentered_no_skiff))
cat("\nAmplitude of the excursions, under three readings of the sentence:\n")
print(as.data.frame(distance_tbl), digits = 3)

by_year <- info[, .(n_particles = .N,
                    pct_reentered = 100 * mean(reentered),
                    median_dist_path_km = median(dist_path_km, na.rm = TRUE)),
                by = Year][order(Year)]
fwrite(by_year, file.path(rev_dir, "reentry_by_year.csv"))

## ---------------------------------------------------------------------------
## PART 5. Supplementary figures
## ---------------------------------------------------------------------------
# The one that answers the reviewer: retention against duration, split by how
# the particles got there.
p_decomp <- decomp %>%
  dplyr::select(weeks, `Arrived and stayed` = pct_first_stay,
                `Left and came back` = pct_returned) %>%
  pivot_longer(-weeks, names_to = "route", values_to = "pct") %>%
  mutate(route = factor(route, levels = c("Left and came back",
                                          "Arrived and stayed"))) %>%
  ggplot(aes(weeks, pct, fill = route)) +
  geom_area(colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = sim_weeks, linetype = 3, colour = "grey30") +
  scale_fill_manual(values = c("Arrived and stayed" = "grey65",
                               "Left and came back" = "#b2182b"), name = NULL) +
  scale_x_continuous(breaks = seq(2, 18, 2)) +
  labs(x = "Simulation duration (weeks)",
       y = "Particles inside \nrecruitment habitat (%)") +
  theme_bw() + theme_paper() +
  theme(panel.grid = element_blank(), legend.position = "top")

p_reentry <- ggplot(reentry_by_duration, aes(weeks, pct_reentered)) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_x_continuous(breaks = seq(2, 18, 2)) +
  labs(x = "Simulation duration (weeks)", y = "Particles having re-entered (%)") +
  theme_bw() + theme_paper() + theme(panel.grid = element_blank())

p_dist <- ggplot(re[is.finite(dist_path_km) & dist_path_km > 0], aes(dist_path_km)) +
  geom_histogram(bins = 60, fill = "grey30", colour = NA) +
  geom_vline(xintercept = median(re$dist_path_km, na.rm = TRUE),
             linetype = 2, colour = "#b2182b", linewidth = 1) +
  scale_x_log10(breaks = c(1, 5, 10, 25, 50, 100, 250, 500)) +
  labs(x = "Distance travelled outside recruitment \nhabitat before re-entry (km)",
       y = "Number of particles") +
  theme_bw() + theme_paper() + theme(panel.grid = element_blank())

ggsave(file.path(rev_fig, "reentry_retention_decomposition.png"), p_decomp,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(rev_fig, "reentry_by_duration.png"), p_reentry,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(rev_fig, "reentry_distance_distribution.png"), p_dist,
       width = 9, height = 6, dpi = 300)
