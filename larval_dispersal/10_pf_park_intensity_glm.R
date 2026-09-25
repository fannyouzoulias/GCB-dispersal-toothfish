################################################################################
# 10_pf_park_intensity_glm.R
#
# WHY THE PF INDEX WAS CHANGED AT REVISION. 05_front_indices.R and
# 06_retention_front_glm.R now sample the front intensity along the front OF
# THAT YEAR, the ADT contour of Park et al. (2019); the manuscript sampled it
# along the climatological contour, fixed for all years. This script is the
# justification: it refits the model of 06 on three sampling paths and shows
# what the change does to the result.
#
# The three indices are sampled in the SAME yearly velocity field
# (advection/average_velocity_field.ipynb, release week 23 to week 31 + 18),
# so the sampling PATH is the only thing that differs between them:
#
#   I_static  mean speed along the climatological PF   (the manuscript index)
#   I_ww      mean speed along that year's Winter Water edge
#   I_park    mean speed along that year's ADT contour (the index now used)
#
# All three are restricted to 67-72 E. The 51-48 S cut of the manuscript
# index (the sector of the Methods, and the box of 05) picks the Kerguelen
# branch of a contour that wanders, and so it applies to the two CONTOURS, not
# to the Winter Water edge:
#
#   I_static  cut: the climatological contour loops
#   I_park    cut: the ADT contour loops down to 52-54 S west of the islands in
#                  2000, 2009 and 2020 -- 40 % of its points in the band in
#                  2020, which would otherwise be averaged in
#   I_ww      no cut: a northern limit, single-valued per meridian, and in 2023
#                  it lies entirely south of -51 in this band, so a cut would
#                  delete the year outright
#
# The other version of each is written out beside it.
#
# Then, for each index, the model of 06_retention_front_glm.R:
#     log R ~ z(I),  Gaussian GLM, R the annual number of retained eggs.
#
# Two consistency checks are run on the way. I_static must reproduce
# PF_mean_cm_s, and I_park must agree with PF_park_cm_s: 05 takes the speed
# along the annual contour from the Python script itself rather than
# re-sampling the velocity field here, and the two must not disagree.
#
# Requires: 00_setup.R, the CSVs of 04_retention_recruitment.R and
#           05_front_indices.R, and the front positions written by advection/
#           (park_dir, ww_int_dir, uv_dir in config.R)
#
# Author: Fanny Ouzoulias
# Date:   2026-09-21
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(ncdf4)
})

## ---------------------------------------------------------------------------
## PART 0. Paths and settings (as 15)
## ---------------------------------------------------------------------------
# All three fronts, and the field they are sampled in, come from advection/.
# Their folders are set in config.R.
dir_fields <- uv_dir       # average_velocity_field.ipynb
dir_sent   <- ww_int_dir   # surface_currents_at_PF_location.ipynb
dir_park   <- park_dir     # PF_position_park.py

if (!exists("out_dir")) out_dir <- "outputs"
rev_dir <- file.path(out_dir, "revision")
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)

if (!exists("theme_paper")) theme_paper <- function() theme_bw(base_size = 13)

years       <- 2000:2023
south_years <- c(2012, 2013, 2017, 2023)
LON_MIN <- 67; LON_MAX <- 72; LAT_MIN <- -51; LAT_MAX <- -48
cols <- c(other = "#2166ac", displaced = "#b2182b")

## ---------------------------------------------------------------------------
## PART 1. The three indices, in one and the same velocity field
## ---------------------------------------------------------------------------
read_field <- function(y) {
  nc <- nc_open(file.path(dir_fields, sprintf("%d_uv_norm_field.nc", y)))
  on.exit(nc_close(nc))
  list(lon = as.numeric(ncvar_get(nc, "lon")), lat = as.numeric(ncvar_get(nc, "lat")),
       speed = ncvar_get(nc, "Geostrophic velocity field norm"))   # [lon, lat]
}
fields <- set_names(map(years, read_field), years)

# Nearest grid cell, as for the published index.
sample_field <- function(f, lon, lat) {
  i <- round((lon - f$lon[1]) / diff(f$lon[1:2])) + 1
  j <- round((lat - f$lat[1]) / diff(f$lat[1:2])) + 1
  ok <- !is.na(lat) & i >= 1 & i <= length(f$lon) & j >= 1 & j <= length(f$lat)
  out <- rep(NA_real_, length(lon)); out[ok] <- f$speed[cbind(i[ok], j[ok])]
  out
}

# One index from any table of front points, with and without the latitude cut.
index_along <- function(tab, y) {
  b <- filter(tab, lon >= LON_MIN, lon <= LON_MAX)
  s <- sample_field(fields[[as.character(y)]], b$lon, b$lat)
  bc <- b$lat >= LAT_MIN & b$lat <= LAT_MAX
  tibble(Year = y,
         I              = mean(s[bc], na.rm = TRUE),
         I_nocut        = mean(s,     na.rm = TRUE),
         lat_mean       = mean(b$lat[bc], na.rm = TRUE),
         lat_mean_nocut = mean(b$lat,     na.rm = TRUE),
         n_points       = sum(bc),
         n_points_all   = nrow(b))
}

# (a) the published static index: the climatological contour
pf_clim <- read_csv(dpath("front_intensity_PF.csv"), show_col_types = FALSE) %>%
  transmute(lon = `Lon PF`, lat = `Lat PF`) %>%
  filter(!is.na(lon), !is.na(lat))

static <- map_dfr(years, ~ index_along(pf_clim, .x)) %>%
  transmute(Year, I_static = I)

# The published index is reproduced by this sampling: check before going on.
written <- read_csv(file.path(out_dir, "front_indices_annual.csv"),
                    show_col_types = FALSE) %>%
  transmute(Year, PF_mean_cm_s, PF_park_cm_s)
chk <- left_join(static, written, by = "Year")
stopifnot(max(abs(chk$I_static - chk$PF_mean_cm_s)) < 1e-4)
message("climatological index reproduced from the fields (max |diff| < 1e-4 cm/s)")

# (b) the annual Winter Water edge, as in 15
ww <- map_dfr(years, function(y)
  read_csv(file.path(dir_sent, sprintf("intensity_PF_%d_week_25-29.csv", y)),
           show_col_types = FALSE, name_repair = "unique_quiet",
           col_select = c(lon_PF, lat_PF)) %>%
    transmute(lon = lon_PF, lat = lat_PF) %>%
    index_along(y)) %>%
  transmute(Year, I_ww = I_nocut, I_ww_cut = I, lat_ww = lat_mean_nocut)

# (c) the annual ADT contour of Park et al. (2019), the new path
park <- map_dfr(years, function(y)
  read_csv(file.path(dir_park, sprintf("pf_park_%d.csv", y)),
           show_col_types = FALSE, col_select = c(lon, lat)) %>%
    index_along(y)) %>%
  transmute(Year, I_park = I, I_park_nocut = I_nocut, lat_park = lat_mean,
            n_park = n_points)

# 05 takes the speed straight from PF_position_park.py; here it is re-sampled
# in the LAMTA-averaged field. Same path, two velocity products: they must
# agree, or the index of 06 and the index compared here are not the same thing.
cmp <- left_join(park, written, by = "Year")
message(sprintf(
  "annual-contour index: 05 vs re-sampled here, r = %.3f, max |diff| = %.2f cm/s",
  cor(cmp$I_park, cmp$PF_park_cm_s), max(abs(cmp$I_park - cmp$PF_park_cm_s))))
stopifnot(cor(cmp$I_park, cmp$PF_park_cm_s) > 0.98)

rec <- read_csv(file.path(out_dir, "recruited_annual.csv"), show_col_types = FALSE) %>%
  transmute(Year, Recruited_total, logR = log(Recruited_total))

dat <- reduce(list(static, ww, park, rec), left_join, by = "Year") %>%
  mutate(group = if_else(Year %in% south_years, "displaced", "other"))
stopifnot(nrow(dat) == 24, !anyNA(dplyr::select(dat, I_static, I_ww, I_park, logR)))

write_csv(dat, file.path(rev_dir, "pf_park_intensity_annual.csv"))

## ---------------------------------------------------------------------------
## PART 2. Do the three indices agree?
## ---------------------------------------------------------------------------
idx <- c(I_static = "Climatological PF (manuscript)",
         I_ww     = "Annual Winter Water edge",
         I_park   = "Annual PF streamline (Park)")

agreement <- map_dfr(c("I_ww", "I_park"), function(v) {
  pe <- cor.test(dat$I_static, dat[[v]])
  sp <- suppressWarnings(cor.test(dat$I_static, dat[[v]], method = "spearman"))
  o  <- dat$group == "other"
  tibble(index = idx[[v]],
         mean_cm_s = mean(dat[[v]]), sd_cm_s = sd(dat[[v]]),
         pearson_r = unname(pe$estimate), p_pearson = pe$p.value,
         spearman_rho = unname(sp$estimate),
         R2_with_static = unname(pe$estimate)^2,
         pearson_r_without_4_years = cor(dat$I_static[o], dat[[v]][o]))
})
write_csv(agreement, file.path(rev_dir, "pf_park_intensity_agreement.csv"))
print(agreement)

## ---------------------------------------------------------------------------
## PART 3. The model of 06, refitted on each index
## ---------------------------------------------------------------------------
fit_one <- function(v) {
  d <- mutate(dat, z = as.numeric(scale(.data[[v]])))
  m <- glm(logR ~ z, data = d, family = gaussian())
  s <- summary(m)$coefficients
  r2 <- 1 - deviance(m) / m$null.deviance
  tibble(index = idx[[v]], variable = v,
         slope = s["z", "Estimate"], se = s["z", "Std. Error"],
         t = s["z", "t value"], p = s["z", "Pr(>|t|)"],
         R2 = r2,
         R2_adj = 1 - (1 - r2) * ((nrow(d) - 1) / (nrow(d) - length(coef(m)))),
         pct_per_sd = 100 * (exp(s["z", "Estimate"]) - 1))
}

models <- map_dfr(names(idx), fit_one)
write_csv(models, file.path(rev_dir, "pf_park_intensity_glm.csv"))
print(models)

## ---------------------------------------------------------------------------
## PART 4. Figure
## ---------------------------------------------------------------------------
long <- dat %>%
  dplyr::select(Year, group, logR, Recruited_total, all_of(names(idx))) %>%
  pivot_longer(all_of(names(idx)), names_to = "variable", values_to = "I") %>%
  group_by(variable) %>% mutate(z = as.numeric(scale(I))) %>% ungroup() %>%
  mutate(panel = factor(idx[variable], levels = unname(idx)))

lab <- models %>%
  mutate(panel = factor(index, levels = unname(idx)),
         text  = sprintf("slope %+.3f (p = %.3f)\nR2 = %.2f", slope, p, R2))

p_fit <- ggplot(long, aes(z, Recruited_total)) +
  geom_smooth(method = "glm", formula = y ~ x,
              method.args = list(family = gaussian(link = "log")),
              colour = "grey20", fill = "grey80", linewidth = 0.7) +
  geom_point(aes(colour = group), size = 2.6) +
  geom_text_repel(data = filter(long, group == "displaced"), aes(label = Year),
                  colour = cols[["displaced"]], size = 3.2, seed = 1,
                  min.segment.length = 0) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = text), hjust = -0.08,
            vjust = 1.15, size = 3.4, inherit.aes = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "Mean surface geostrophic current speed along the path (standardised)",
       y = "Number of larvae retained") +
  theme_paper()

p_scatter <- ggplot(dat, aes(I_static, I_park)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey60") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey20", fill = "grey80",
              linewidth = 0.7) +
  geom_point(aes(colour = group), size = 2.8) +
  geom_text_repel(data = filter(dat, group == "displaced"), aes(label = Year),
                  colour = cols[["displaced"]], size = 3.2, seed = 1,
                  min.segment.length = 0) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "Intensity along the climatological PF (cm/s)",
       y = "PF-associated jet intensity\n(annual PF streamline, cm/s)") +
  theme_paper()

p <- p_fit / p_scatter + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave(file.path(rev_fig, "pf_park_intensity_glm.png"), p,
       width = 12, height = 8.5, dpi = 250)
message("Written: ", normalizePath(file.path(rev_fig, "pf_park_intensity_glm.png"),
                                   mustWork = FALSE))
