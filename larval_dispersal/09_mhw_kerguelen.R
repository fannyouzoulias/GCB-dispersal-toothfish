################################################################################
# 09_mhw_kerguelen.R
#
# Marine heatwaves at 200 m around Kerguelen, 2000-2023.
#
# The four years in which the Winter Water front is displaced south (2012, 2013,
# 2017, 2023) could simply be warm years. This script tests that directly: it
# applies the marine-heatwave definition of Hobday et al. (2016) to the
# temperature at 200 m -- the depth that carries the Winter Water signature the
# front is read from, so a subsurface heatwave is the quantity that matters, not
# a surface one.
#
# Definition: above the seasonally varying 90th percentile for at least 5
# consecutive days, percentile built on an 11-day window and smoothed over 31
# days, gaps of at most 2 days merged, no detrending. Climatology 1993-2022, the
# thirty years GLORYS12 covers; it therefore contains 2012, 2013 and 2017, which
# makes the test conservative for those years.
#
# The series is the area-weighted daily mean over 65-72 E, 55-45 S.
#
# Reads the GLORYS12 daily potential temperature interpolated to 200 m and
# averaged to 0.25 degrees by advection/download_glorys_T200_kerguelen.py
# (33 MB; that script is the only way to rebuild the NetCDF). Its folder is
# `t200_dir`, set in config.R.
#
# Writes
#   outputs/revision/mhw200_annual_2000_2023.csv    days, events, cumulative and
#                                                   maximum anomaly, per year
#   outputs/figures/revision/mhw200_chronology.png  the 24-year chronology
#
# Author: Fanny Ouzoulias
# Date:   2026-09-21
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(terra)
  library(heatwaveR)
})

## ---------------------------------------------------------------------------
## PART 0. Paths and settings
## ---------------------------------------------------------------------------
data_file <- file.path(t200_dir, "glorys12_T200m_1993_2023_Kerguelen_025deg.nc")
data_var  <- "temp"
var_label <- "T 200 m"

if (!exists("out_dir")) out_dir <- "outputs"
rev_dir <- file.path(out_dir, "revision")
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)

analysis_start <- 2000
analysis_end   <- 2023
clim_start     <- "1993-01-01"
clim_end       <- "2022-12-31"
PCTILE <- 90; MIN_DURATION <- 5; MAX_GAP <- 2

if (!file.exists(data_file))
  stop("missing ", data_file,
       "\n  -> rebuild it with advection/download_glorys_T200_kerguelen.py")

## ---------------------------------------------------------------------------
## PART 1. The regional daily series
## ---------------------------------------------------------------------------
# ncdf4 then terra: reading thousands of layers through GDAL directly is far
# slower than building the raster from the array in memory.
load_grid <- function(file, var) {
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc))
  lon <- as.numeric(ncdf4::ncvar_get(nc, "lon"))
  lat <- as.numeric(ncdf4::ncvar_get(nc, "lat"))
  time_raw <- ncdf4::ncvar_get(nc, "time")
  time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  arr <- ncdf4::ncvar_get(nc, var, collapse_degen = FALSE)   # [lon, lat, time]

  stopifnot(grepl("^days since", time_units))
  dates <- as.Date(sub("^days since ([0-9-]+).*", "\\1", time_units)) + floor(time_raw)

  arr <- aperm(arr, c(2, 1, 3))[order(lat, decreasing = TRUE), order(lon), , drop = FALSE]
  res_deg <- abs(diff(sort(lon))[1])
  r <- rast(arr, crs = "EPSG:4326",
            extent = ext(min(lon) - res_deg / 2, max(lon) + res_deg / 2,
                         min(lat) - res_deg / 2, max(lat) + res_deg / 2))
  time(r) <- dates
  list(r = r, dates = dates)
}

g <- load_grid(data_file, data_var)
keep <- g$dates >= as.Date(clim_start) &
        g$dates <= as.Date(sprintf("%d-12-31", analysis_end))
grd <- g$r[[which(keep)]]
dates <- g$dates[keep]
rm(g)

# Area-weighted mean: cells shrink poleward, so a plain mean would over-weight
# the south of the box.
area <- values(cellSize(grd[[1]], unit = "km"))[, 1]
v <- values(grd)
ts_var <- tibble(t = dates,
                 temp = as.numeric(colSums(v * area, na.rm = TRUE) /
                                     colSums((!is.na(v)) * area))) %>%
  arrange(t)
stopifnot(all(diff(ts_var$t) == 1))
message(sprintf("series: %s to %s, %d days, %.2f to %.2f degC",
                min(ts_var$t), max(ts_var$t), nrow(ts_var),
                min(ts_var$temp), max(ts_var$temp)))

## ---------------------------------------------------------------------------
## PART 2. Detection (Hobday et al. 2016)
## ---------------------------------------------------------------------------
clim <- ts2clm(data = ts_var, x = t, y = temp,
               climatologyPeriod = c(clim_start, clim_end),
               pctile = PCTILE, windowHalfWidth = 5,
               smoothPercentile = TRUE, smoothPercentileWidth = 31)
mhw <- detect_event(clim, x = t, y = temp, minDuration = MIN_DURATION,
                    joinAcrossGaps = TRUE, maxGap = MAX_GAP)

annual <- mhw$climatology %>%
  mutate(year = year(t), is_mhw = as.logical(event), anom = temp - seas) %>%
  filter(year >= analysis_start, year <= analysis_end) %>%
  group_by(year) %>%
  summarise(
    MHW_days      = sum(is_mhw, na.rm = TRUE),
    MHW_events    = if (any(is_mhw, na.rm = TRUE)) n_distinct(event_no[is_mhw], na.rm = TRUE) else 0L,
    MHW_cum_anom  = if (any(is_mhw, na.rm = TRUE)) sum(anom[is_mhw], na.rm = TRUE) else 0,
    MHW_max_anom  = if (any(is_mhw, na.rm = TRUE)) max(anom[is_mhw], na.rm = TRUE) else NA_real_,
    mean_anomaly  = mean(anom, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(pct_MHW_days = 100 * percent_rank(MHW_days),
         pct_MHW_cum  = 100 * percent_rank(MHW_cum_anom))

write_csv(annual, file.path(rev_dir, "mhw200_annual_2000_2023.csv"))
print(annual %>% filter(MHW_days > 0), n = 30)

## ---------------------------------------------------------------------------
## PART 3. The chronology, 2000-2023
## ---------------------------------------------------------------------------
# Four years per row, six rows. The year is written inside the panel because a
# facet strip per year would eat a third of the height.
years_all <- analysis_start:analysis_end
bloc_of <- function(y) {
  b <- analysis_start + 4 * ((y - analysis_start) %/% 4)
  sprintf("%d-%d", b, pmin(b + 3, analysis_end))
}
bloc_levels <- unique(bloc_of(years_all))

chrono <- mhw$climatology %>%
  filter(t >= as.Date(sprintf("%d-01-01", analysis_start)),
         t <= as.Date(sprintf("%d-12-31", analysis_end))) %>%
  mutate(mhw_top = if_else(as.logical(event), pmax(temp, thresh), thresh),
         panel   = factor(bloc_of(year(t)), levels = bloc_levels))

y_lab <- max(chrono$temp, na.rm = TRUE) + 0.12 * diff(range(chrono$temp, na.rm = TRUE))
year_labels <- tibble(year = years_all) %>%
  mutate(t = as.Date(sprintf("%d-07-02", year)), y = y_lab,
         panel = factor(bloc_of(year), levels = bloc_levels))

col_lines <- set_names(c("black", "#1fcc1f", "#1f1fff"),
                       c(var_label, "90th percentile", "Climatology"))

p <- ggplot(chrono, aes(t)) +
  geom_ribbon(aes(ymin = thresh, ymax = mhw_top, fill = "Marine heatwave")) +
  geom_line(aes(y = seas,   colour = "Climatology"), linewidth = 0.7) +
  geom_line(aes(y = thresh, colour = "90th percentile"), linewidth = 0.7) +
  geom_line(aes(y = temp,   colour = var_label), linewidth = 0.4) +
  # every year in the same style: the four displaced years are not singled out
  geom_text(data = year_labels, aes(t, y, label = year), inherit.aes = FALSE,
            colour = "grey25", size = 5.4) +
  scale_colour_manual(values = col_lines, breaks = names(col_lines), name = NULL) +
  scale_fill_manual(values = c("Marine heatwave" = "#f4a3a3"), name = NULL) +
  scale_x_date(breaks = as.Date(sprintf("%d-01-01", years_all)),
               labels = NULL, expand = expansion(mult = 0.003)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(linewidth = 2.2)),
         fill = guide_legend(order = 2, nrow = 1)) +
  facet_wrap(~ panel, scales = "free_x", ncol = 1) +
  labs(x = NULL, y = "Temperature at 200 m (\u00b0C)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(size = 17),
        legend.key.width = unit(2.2, "cm"),
        legend.background = element_rect(fill = "white", colour = NA),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, margin = margin(r = 10)),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.4),
        strip.text = element_blank(),
        strip.background = element_blank())

ggsave(file.path(rev_fig, "mhw200_chronology.png"), p,
       width = 15, height = 14, dpi = 250)
message("Written: ", normalizePath(file.path(rev_fig, "mhw200_chronology.png"),
                                   mustWork = FALSE))
