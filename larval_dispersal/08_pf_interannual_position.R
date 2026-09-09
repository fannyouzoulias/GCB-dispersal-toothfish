################################################################################
# 08_pf_interannual_position.R
#
# Exploratory re-analysis for the GCB revision.
#
# 05_front_indices.R samples the geostrophic speed along the CLIMATOLOGICAL
# (time-mean) position of the Polar Front (Park & Durand), a contour that does
# not move between years. A. Nalivaev has since produced, for every year
# 2000-2023:
#   - the position of the PF that year, defined as the NORTHERN edge of the
#     Winter Water probability of presence (> 0.8) over weeks 25-29 -- more
#     precisely the northernmost point topping a run of 12 contiguous latitude
#     cells above the threshold -- smoothed with a 24-point (2 deg) running
#     mean in longitude;
#   - the DUACS 1/8 deg geostrophic speed averaged over the advection window of
#     that year (Thursday of week 25 to Thursday of week 29 + 18 weeks),
#     sampled AT that year's front position.
#
# That temporal window is the right one: 04_retention_recruitment.R keeps
# spawning weeks 24-28 (`release_weeks <- 24:28`), i.e. release weeks 25-29, so
# weeks 25-29 plus the 18-week PLD covers exactly the dispersal of the particles
# that both annual responses are computed from. Position, velocity field and
# retention therefore share one window.
#
# Two questions:
#   Q1. Does log(retention) ~ PF intensity give the same answer when the front
#       is allowed to move, using the same spatial window as in the paper
#       (67-72 E, and the paper's additional lat >= -51 filter)?
#   Q2. Is the mean latitude of the PF correlated with the SAM?
#   Q3. Given the answer to Q1, which model should the paper report?
#
# Nothing here overwrites the published outputs: everything is written to
# outputs/revision/ and outputs/figures/revision/.
#
# Requires: 00_setup.R  (for theme_paper(); the script is otherwise standalone)
#
# Author: Fanny Ouzoulias
# Date:   2026-09-04
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

## ---------------------------------------------------------------------------
## PART 0. Paths
## ---------------------------------------------------------------------------
# The interannual front files and the SAM series are NOT in the SEANOE archive,
# so they are pointed at explicitly here rather than through dpath().

pf_root  <- "P:/MNHN/PhD/dispersion_larvaire/position_PF"
dir_int  <- file.path(pf_root, "intensity_PF_interannual")  # position + speed
dir_pos  <- file.path(pf_root, "position_PF_interannual")   # position only
sam_dir  <- "P:/MNHN/PhD/donnees_env/env_data/SAM"

# Static (published) inputs. Uses the archive if config.R found it, otherwise
# the working copy of the SEANOE deposit.
seanoe   <- "P:/MNHN/PhD/dispersion_larvaire/simus_dispersion/seanoe"
pick <- function(f) if (exists("data_dir") && file.exists(file.path(data_dir, f)))
  file.path(data_dir, f) else file.path(seanoe, f)

f_pf_static <- pick("front_intensity_PF.csv")

# Annual retention comes from 04_retention_recruitment.R.
if (!exists("out_dir")) out_dir <- "outputs"
f_recruited <- file.path(out_dir, "recruited_annual.csv")

rev_dir <- file.path(out_dir, "revision")
rev_fig <- file.path(out_dir, "figures", "revision")
dir.create(rev_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rev_fig, showWarnings = FALSE, recursive = TRUE)

if (!exists("theme_paper")) theme_paper <- function() theme_bw(base_size = 13)
save_rev <- function(name, plot, width = 9, height = 6)
  ggsave(file.path(rev_fig, name), plot, width = width, height = height, dpi = 300)

## Spatial window of the paper --------------------------------------------------
LON_MIN <- 67; LON_MAX <- 72   # 05_front_indices.R: filter(lon >= 67, lon <= 72)
LAT_MIN <- -51                 # 05_front_indices.R: filter(lat >= -51)

## ---------------------------------------------------------------------------
## PART 1. Read the year-by-year front
## ---------------------------------------------------------------------------
years <- 2000:2023

pf_var <- map_dfr(years, function(y) {
  read_csv(file.path(dir_int, sprintf("intensity_PF_%d_week_25-29.csv", y)),
           show_col_types = FALSE) %>%
    dplyr::select(lon = lon_PF, lat = lat_PF,
                  intensity = `surface currents at PF(cm/s)`) %>%
    mutate(Year = y)
})

message("Interannual front: ", nrow(pf_var), " points, ",
        n_distinct(pf_var$Year), " years, lon ",
        round(min(pf_var$lon), 2), "-", round(max(pf_var$lon), 2))

## The paper's static front, same file that 05_front_indices.R reads -----------
pf_static <- read_csv(f_pf_static, show_col_types = FALSE) %>%
  rename(lon = `Lon PF`, lat = `Lat PF`) %>%
  pivot_longer(matches("^[0-9]{4}$"), names_to = "Year", values_to = "intensity") %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(intensity))

## ---------------------------------------------------------------------------
## PART 2. How far does the front move relative to the paper's window?
## ---------------------------------------------------------------------------
# The static contour meanders, so at a given longitude it can carry several
# latitudes: there, `lat >= -51` selects the northern branch. The interannual
# contour is single valued in longitude (one point per 1/12 deg), so the same
# filter no longer selects a branch -- it drops whole years in which the front
# sat south of 51 S. This table quantifies that.

window_check <- pf_var %>%
  filter(lon >= LON_MIN, lon <= LON_MAX) %>%
  group_by(Year) %>%
  summarise(n_points   = n(),
            lat_mean   = mean(lat),
            lat_min    = min(lat),
            lat_max    = max(lat),
            n_north_51 = sum(lat >= LAT_MIN),
            .groups = "drop")

write_csv(window_check, file.path(rev_dir, "pf_window_check.csv"))
print(window_check, n = Inf)

dropped <- window_check$Year[window_check$n_north_51 == 0]
if (length(dropped))
  message("!! Years with NO front point north of ", LAT_MIN, " deg in ",
          LON_MIN, "-", LON_MAX, " E: ", paste(dropped, collapse = ", "),
          "\n   -> the paper's lat filter cannot be applied verbatim to the ",
          "moving front.")

## ---------------------------------------------------------------------------
## PART 3. Annual indices
## ---------------------------------------------------------------------------
# Three indices, all averaged over the SAME longitude window as the paper:
#   I_static  : paper's index (static contour, lat >= -51)          [reference]
#   I_var     : moving front, all points in 67-72 E                 [main test]
#   I_var_n51 : moving front, plus the paper's lat >= -51 filter    [sensitivity]

idx_static <- pf_static %>%
  filter(lon >= LON_MIN, lon <= LON_MAX, lat >= LAT_MIN) %>%
  group_by(Year) %>%
  summarise(I_static = mean(intensity, na.rm = TRUE), .groups = "drop")

idx_var <- pf_var %>%
  filter(lon >= LON_MIN, lon <= LON_MAX) %>%
  group_by(Year) %>%
  summarise(I_var    = mean(intensity, na.rm = TRUE),
            lat_var  = mean(lat),
            .groups  = "drop")

idx_var_n51 <- pf_var %>%
  filter(lon >= LON_MIN, lon <= LON_MAX, lat >= LAT_MIN) %>%
  group_by(Year) %>%
  summarise(I_var_n51 = mean(intensity, na.rm = TRUE), .groups = "drop")

# Whole-domain latitude, for the SAM analysis (PART 6): the SAM is a
# hemispheric mode, so the front position is also averaged over all longitudes
# available (61-79 E), not only over the retention window.
lat_domain <- pf_var %>%
  group_by(Year) %>%
  summarise(lat_domain = mean(lat),
            I_domain   = mean(intensity, na.rm = TRUE), .groups = "drop")

indices <- reduce(list(idx_static, idx_var, idx_var_n51, lat_domain),
                  full_join, by = "Year") %>% arrange(Year)

write_csv(indices, file.path(rev_dir, "pf_indices_static_vs_varying.csv"))

## Do the two intensity indices agree? ------------------------------------------
cmp <- indices %>% filter(!is.na(I_static), !is.na(I_var))
r_pearson  <- cor(cmp$I_static, cmp$I_var)
r_spearman <- cor(cmp$I_static, cmp$I_var, method = "spearman")
message(sprintf("I_static vs I_var: r = %.3f (Pearson), rho = %.3f (Spearman), n = %d",
                r_pearson, r_spearman, nrow(cmp)))

p_cmp_ts <- indices %>%
  dplyr::select(Year, I_static, I_var) %>%
  pivot_longer(-Year, names_to = "index", values_to = "cm_s") %>%
  mutate(index = recode(index,
                        I_static = "Static (climatological) front, paper",
                        I_var    = "Interannual front position")) %>%
  ggplot(aes(Year, cm_s, colour = index)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = c("#b2182b", "#2166ac"), name = NULL) +
  labs(y = "PF intensity (cm/s), 67-72 E", x = NULL) +
  theme_paper() + theme(legend.position = "top", panel.grid = element_blank())

p_cmp_sc <- ggplot(cmp, aes(I_static, I_var)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", alpha = 0.15) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(aes(label = Year), size = 3, seed = 1) +
  labs(x = "Static index (cm/s)", y = "Interannual-position index (cm/s)",
       subtitle = sprintf("r = %.2f", r_pearson)) +
  theme_paper() + theme(panel.grid = element_blank())

save_rev("pf_index_static_vs_varying.png", p_cmp_ts, width = 10, height = 5)
save_rev("pf_index_static_vs_varying_scatter.png", p_cmp_sc, width = 7, height = 6)

## ---------------------------------------------------------------------------
## PART 3b. The front position is bimodal -- and the excursions are real
## ---------------------------------------------------------------------------
# In 67-72 E the detected front sits near 49.5 S in 20 years out of 24 and jumps
# to 51.4-52.2 S in 2012, 2013, 2017 and 2023: a 2.5 deg step, not a gradual
# shift. Those same four years carry the four lowest sampled intensities, so at
# first sight they look like a failure of the WW-edge detector.
#
# They are not. Their latitude is corroborated by the independent Mercator front
# position at 70 E used in 07_front_validation.R (Azarian et al. 2024): over
# 2000-2019 the two latitude series correlate, and 2017, 2013 and 2012 are among
# the most southern years of that independent record. The southern excursions
# are a real feature of the front, so these years must NOT be dropped -- the
# 20-year fits below are a sensitivity, not the result.
#
# What the low intensities then mean is that a displaced front is also a weaker
# one AT ITS OWN POSITION, which is why the single-index model of PART 4 fails:
# I_var confounds how fast the jet is with where it is. PART 4c separates them.

indices <- indices %>%
  mutate(branch = if_else(lat_var >= LAT_MIN, "northern", "southern"))

message("Southern-displacement years: ",
        paste(indices$Year[indices$branch == "southern"], collapse = ", "))

north <- indices %>% filter(branch == "northern")
message(sprintf("I_static vs I_var, northern years only: r = %.3f (n = %d)",
                cor(north$I_static, north$I_var), nrow(north)))

## Independent check of the front latitude --------------------------------------
# pf_index_validation.csv is written by 07_front_validation.R and holds, per
# year 2000-2019, the mean latitude of the PF at 70 E from the Mercator
# reanalysis. It never saw the WW detector, so it is an independent verdict on
# the southern excursions. 2020-2023 are outside that record: the most extreme
# year, 2023, therefore remains unverified.
f_val <- file.path(out_dir, "pf_index_validation.csv")
if (file.exists(f_val)) {
  val <- read_csv(f_val, show_col_types = FALSE) %>%
    dplyr::select(Year, lat_mercator = lat_mean)
  vv <- indices %>% inner_join(val, by = "Year")
  ct <- cor.test(vv$lat_var, vv$lat_mercator)
  cs <- suppressWarnings(cor.test(vv$lat_var, vv$lat_mercator, method = "spearman"))
  message(sprintf(
    "WW-detector latitude vs Mercator latitude at 70 E (n = %d): r = %.3f (p = %.4f), rho = %.3f (p = %.4f)",
    nrow(vv), ct$estimate, ct$p.value, cs$estimate, cs$p.value))
  write_csv(vv %>% dplyr::select(Year, lat_var, lat_mercator, I_static, I_var, branch),
            file.path(rev_dir, "pf_latitude_vs_mercator.csv"))

  p_val <- ggplot(vv, aes(lat_mercator, lat_var, colour = branch)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", alpha = 0.15) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = Year), size = 3, seed = 1,
                             show.legend = FALSE) +
    scale_colour_manual(values = c(northern = "#2166ac", southern = "#b2182b"),
                        name = NULL) +
    labs(x = "PF latitude at 70 E, Mercator (deg N)",
         y = "PF latitude 67-72 E, WW detector (deg N)",
         subtitle = sprintf("r = %.2f, p = %.3f (2000-2019)",
                            ct$estimate, ct$p.value)) +
    theme_paper() + theme(panel.grid = element_blank())
  save_rev("pf_latitude_vs_mercator.png", p_val, width = 8, height = 6)
} else {
  message("pf_index_validation.csv not found: run 07_front_validation.R first ",
          "to check the southern excursions against Mercator.")
}

p_branch <- ggplot(indices, aes(lat_var, I_var, colour = branch)) +
  geom_hline(yintercept = mean(indices$I_var), linetype = 3, colour = "grey50") +
  geom_vline(xintercept = LAT_MIN, linetype = 2, colour = "grey40") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = Year), size = 3, seed = 1,
                           show.legend = FALSE) +
  scale_colour_manual(values = c(northern = "#2166ac", southern = "#b2182b"),
                      name = "Detected branch") +
  labs(x = "Mean PF latitude, 67-72 E (deg N)",
       y = "Sampled PF intensity (cm/s)") +
  theme_paper() + theme(panel.grid = element_blank())

save_rev("pf_latitude_vs_intensity_branch.png", p_branch, width = 8, height = 6)

## ---------------------------------------------------------------------------
## PART 4. Q1 -- log(retention) ~ PF intensity
## ---------------------------------------------------------------------------
rec <- read_csv(f_recruited, show_col_types = FALSE) %>%
  dplyr::select(Year, Recruited_total)

# Retention (% of released eggs still on the shelf) is the alternative response
# the paper reports alongside recruitment; it is used as a robustness check.
ret <- read_csv(file.path(out_dir, "retention_annual.csv"), show_col_types = FALSE) %>%
  dplyr::select(Year, mean_ret)

dat <- rec %>%
  left_join(indices, by = "Year") %>%
  left_join(ret, by = "Year") %>%
  mutate(logR  = log(Recruited_total),
         south = as.integer(branch == "southern"))

# Same model as 06_retention_front_glm.R: Gaussian GLM of log R on the
# standardized index.
fit_index <- function(df, index, label, response = "logR") {
  d <- df %>% filter(!is.na(.data[[index]]), !is.na(.data[[response]])) %>%
    mutate(z = as.numeric(scale(.data[[index]])),
           y = .data[[response]])
  m  <- glm(y ~ z, data = d, family = gaussian())
  s  <- summary(m)$coefficients
  r2 <- 1 - deviance(m) / m$null.deviance
  tibble(
    model     = label,
    n         = nrow(d),
    intercept = s["(Intercept)", "Estimate"],
    slope     = s["z", "Estimate"],
    se        = s["z", "Std. Error"],
    t         = s["z", "t value"],
    p         = s["z", "Pr(>|t|)"],
    R2        = r2,
    R2_adj    = 1 - (1 - r2) * ((nrow(d) - 1) / (nrow(d) - 2)),
    # effect size on the natural scale: multiplicative change in retention
    # for a one-SD increase in the index
    ratio_per_SD = exp(s["z", "Estimate"])
  )
}

dat_n <- dat %>% filter(branch == "northern")

glm_compare <- bind_rows(
  fit_index(dat,   "I_static",  "Paper: static front, 67-72 E, lat >= -51"),
  fit_index(dat,   "I_var",     "Interannual front, 67-72 E (all lat)"),
  fit_index(dat,   "I_var_n51", "Interannual front, 67-72 E, lat >= -51"),
  fit_index(dat,   "I_domain",  "Interannual front, whole domain 61-79 E"),
  # like-for-like on the 20 years where the detector stayed on the northern branch
  fit_index(dat_n, "I_static",  "Static front, northern-branch years only"),
  fit_index(dat_n, "I_var",     "Interannual front, northern-branch years only"),
  # same two models with retention (%) instead of recruitment as the response
  fit_index(dat %>% mutate(log_ret = log(mean_ret)), "I_static",
            "Response log(retention %): static front", response = "log_ret"),
  fit_index(dat %>% mutate(log_ret = log(mean_ret)), "I_var",
            "Response log(retention %): interannual front", response = "log_ret")
)

write_csv(glm_compare, file.path(rev_dir, "glm_retention_pf_comparison.csv"))
print(as.data.frame(glm_compare), digits = 3)

## Main figure: the paper's panel, redrawn with the moving front ----------------
# Two panels: all 24 years, and the 20 years on the northern branch.
effect_panel <- function(d, title) {
  d <- d %>% filter(!is.na(I_var)) %>% mutate(z = as.numeric(scale(I_var)))
  m <- glm(logR ~ z, data = d, family = gaussian())
  g <- tibble(z = seq(min(d$z), max(d$z), length.out = 200))
  pr <- predict(m, newdata = g, se.fit = TRUE)
  g <- g %>% mutate(fit = exp(pr$fit),
                    lo  = exp(pr$fit - 1.96 * pr$se.fit),
                    hi  = exp(pr$fit + 1.96 * pr$se.fit))
  s <- summary(m)$coefficients
  ggplot(d, aes(z, Recruited_total)) +
    geom_ribbon(data = g, aes(z, ymin = lo, ymax = hi), alpha = 0.2,
                inherit.aes = FALSE) +
    geom_line(data = g, aes(z, fit), linewidth = 1.1, inherit.aes = FALSE) +
    geom_point(aes(colour = branch), size = 2.8) +
    scale_colour_manual(values = c(northern = "#2166ac", southern = "#b2182b"),
                        name = "Detected branch") +
    labs(x = "Polar Front intensity (interannual position, standardized)",
         y = "Number of \nlarvae retained",
         subtitle = sprintf("%s: slope = %.3f, p = %.3g", title,
                            s["z", "Estimate"], s["z", "Pr(>|t|)"])) +
    theme_paper() + theme(panel.grid = element_blank(), legend.position = "top")
}

p_effect <- effect_panel(dat, "All 24 years") |
  effect_panel(dat_n, "Northern-branch years")

save_rev("pf_effect_on_retention_interannual.png", p_effect, width = 14, height = 6)

## The plain answer, on the retention rate ---------------------------------------
# log(retention %) against the intensity recomputed on the yearly front
# positions, 67-72 E. One predictor, 24 years, nothing else in the model: this
# is the direct replacement of the published panel.
d_plain <- dat %>% mutate(z = as.numeric(scale(I_var)), y = log(mean_ret))
m_plain <- glm(y ~ z, data = d_plain, family = gaussian())
s_plain <- summary(m_plain)$coefficients
g <- tibble(z = seq(min(d_plain$z), max(d_plain$z), length.out = 200))
pr <- predict(m_plain, newdata = g, se.fit = TRUE)
g  <- g %>% mutate(fit = exp(pr$fit),
                   lo = exp(pr$fit - 1.96 * pr$se.fit),
                   hi = exp(pr$fit + 1.96 * pr$se.fit))

p_plain <- ggplot(d_plain, aes(z, mean_ret)) +
  geom_ribbon(data = g, aes(z, ymin = lo, ymax = hi), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = g, aes(z, fit), linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(aes(colour = branch), size = 2.8) +
  ggrepel::geom_text_repel(aes(label = Year), size = 3, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c(northern = "#2166ac", southern = "#b2182b"),
                      name = "Front position") +
  labs(x = "PF intensity, 67-72 E, on the yearly front position (standardized)",
       y = "Larval retention 
rate (%)",
       subtitle = sprintf("slope = %.3f, p = %.2f, R2 = %.3f",
                          s_plain[2, 1], s_plain[2, 4],
                          1 - deviance(m_plain) / m_plain$null.deviance)) +
  theme_paper() + theme(panel.grid = element_blank(), legend.position = "top")

save_rev("retention_vs_pf_intensity_interannual.png", p_plain, width = 9, height = 6)

# The annual index next to the two responses, so the series can be read directly.
write_csv(dat %>% dplyr::select(Year, I_var, lat_var, I_static, mean_ret,
                                Recruited_total, branch),
          file.path(rev_dir, "annual_index_and_retention.csv"))

## ---------------------------------------------------------------------------
## PART 4c. Which model should the paper report?
## ---------------------------------------------------------------------------
# The moving-front index confounds two things: how fast the jet is (intensity)
# and where it is (latitude). PART 3b showed the two are strongly linked --
# displaced years are also weak years -- so a single-index model cannot separate
# them. The candidate models below do, on all 24 years, with no exclusion.
#
#   M1  logR ~ I_static                 the published model
#   M2  logR ~ I_var                    the naive replacement
#   M3  logR ~ lat_var                  position alone
#   M4  logR ~ I_var + lat_var          intensity and position separated
#   M5  logR ~ I_static + lat_var       does position add to the published index?
#   M6  logR ~ I_static + I_var         which index carries the signal?
#   M4b logR ~ I_var + south(0/1)       is the position effect a gradient or a
#                                       contrast between the displaced years and
#                                       the rest?

model_set <- list(
  M1 = logR ~ scale(I_static),
  M2 = logR ~ scale(I_var),
  M3 = logR ~ scale(lat_var),
  M4 = logR ~ scale(I_var) + scale(lat_var),
  M5 = logR ~ scale(I_static) + scale(lat_var),
  M6 = logR ~ scale(I_static) + scale(I_var),
  M4b = logR ~ scale(I_var) + south
)

model_table <- imap_dfr(model_set, function(f, nm) {
  m <- lm(f, data = dat)
  s <- summary(m)
  vif <- if (length(coef(m)) > 2) {
    max(diag(solve(cor(model.matrix(m)[, -1, drop = FALSE])))) } else NA_real_
  as.data.frame(s$coefficients) %>%
    tibble::rownames_to_column("term") %>%
    filter(term != "(Intercept)") %>%
    transmute(model = nm, formula = deparse1(f), term,
              estimate = Estimate, se = `Std. Error`, p = `Pr(>|t|)`,
              R2 = s$r.squared, R2_adj = s$adj.r.squared,
              AIC = AIC(m), max_VIF = vif)
})

write_csv(model_table, file.path(rev_dir, "model_selection.csv"))

# I_var and lat_var correlate strongly, so M4 has a suppression flavour: the
# intensity coefficient moves from near zero (M2) to clearly negative once
# position is held constant. Leave-one-out shows whether that rests on one year.
message(sprintf("cor(I_var, lat_var) = %.3f | cor(I_static, lat_var) = %.3f",
                cor(dat$I_var, dat$lat_var), cor(dat$I_static, dat$lat_var)))

loo <- map_dfr(dat$Year, function(y) {
  s <- summary(lm(logR ~ scale(I_var) + scale(lat_var),
                  data = filter(dat, Year != y)))$coefficients
  tibble(dropped = y, b_intensity = s[2, 1], p_intensity = s[2, 4],
         b_latitude = s[3, 1], p_latitude = s[3, 4])
})
write_csv(loo, file.path(rev_dir, "model_M4_leave_one_out.csv"))
message(sprintf("M4 leave-one-out: worst p(intensity) = %.4f (dropping %d), worst p(latitude) = %.4f",
                max(loo$p_intensity), loo$dropped[which.max(loo$p_intensity)],
                max(loo$p_latitude)))

# Is the position term a gradient, or just the four displaced years sitting
# lower? Refit within the 20 non-displaced years: if latitude still matters
# there, it is a gradient; if it does not, the effect is a regime contrast and
# M4b is the honest description of it.
s_north <- summary(lm(logR ~ scale(I_var) + scale(lat_var),
                      data = filter(dat, branch == "northern")))$coefficients
message(sprintf(
  "Within the 20 non-displaced years: intensity %.3f (p = %.4f), latitude %.3f (p = %.3f)",
  s_north[2, 1], s_north[2, 4], s_north[3, 1], s_north[3, 4]))

# Influence of each year on the PUBLISHED model, to check that the disputed
# years are not the ones holding it up.
m1 <- lm(logR ~ scale(I_static), data = dat)
influence_m1 <- tibble(Year = dat$Year, cooks_d = cooks.distance(m1),
                       leverage = hatvalues(m1)) %>% arrange(desc(cooks_d))
write_csv(influence_m1, file.path(rev_dir, "model_M1_influence.csv"))

## Partial-effect figure for M4 --------------------------------------------------
m4 <- lm(logR ~ scale(I_var) + scale(lat_var), data = dat)
partial <- bind_rows(
  tibble(effect = "PF intensity (position held constant)",
         x = as.numeric(scale(dat$I_var)),
         y = residuals(lm(logR ~ scale(lat_var), data = dat)),
         Year = dat$Year, branch = dat$branch),
  tibble(effect = "PF latitude (intensity held constant)",
         x = as.numeric(scale(dat$lat_var)),
         y = residuals(lm(logR ~ scale(I_var), data = dat)),
         Year = dat$Year, branch = dat$branch)
)

p_partial <- ggplot(partial, aes(x, y)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", alpha = 0.15) +
  geom_point(aes(colour = branch), size = 2.8) +
  facet_wrap(~ effect, scales = "free_x") +
  scale_colour_manual(values = c(northern = "#2166ac", southern = "#b2182b"),
                      name = "Front position") +
  labs(x = "Standardized predictor", y = "Partial residual of log R") +
  theme_paper() + theme(panel.grid = element_blank(), legend.position = "top")

save_rev("model_M4_partial_effects.png", p_partial, width = 12, height = 6)

## ---------------------------------------------------------------------------
## PART 5. Sensitivity to the longitude window
## ---------------------------------------------------------------------------
# The 67-72 E window was chosen on the static contour. Sliding a window of the
# same width along the domain shows whether the retention signal is specific to
# that sector or holds more broadly.
slide <- expand_grid(lo = seq(61, 74, by = 1),
                     subset = c("all years", "northern-branch years")) %>%
  pmap_dfr(function(lo, subset) {
    hi <- lo + (LON_MAX - LON_MIN)
    idx <- pf_var %>% filter(lon >= lo, lon <= hi) %>%
      group_by(Year) %>% summarise(I = mean(intensity, na.rm = TRUE),
                                   lat = mean(lat), .groups = "drop")
    d <- rec %>% left_join(idx, by = "Year") %>%
      left_join(dplyr::select(indices, Year, branch), by = "Year") %>%
      filter(subset == "all years" | branch == "northern") %>%
      mutate(logR = log(Recruited_total), z = as.numeric(scale(I)))
    m <- glm(logR ~ z, data = d, family = gaussian())
    s <- summary(m)$coefficients
    tibble(lon_min = lo, lon_max = hi, subset = subset, n = nrow(d),
           slope = s["z", "Estimate"], se = s["z", "Std. Error"],
           p = s["z", "Pr(>|t|)"],
           R2 = 1 - deviance(m) / m$null.deviance)
  })

write_csv(slide, file.path(rev_dir, "glm_longitude_window_sensitivity.csv"))

p_slide <- ggplot(slide, aes((lon_min + lon_max) / 2, slope)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_ribbon(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se),
              alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(aes(colour = p < 0.05), size = 2.6) +
  facet_wrap(~ subset) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#b2182b"),
                      name = "p < 0.05") +
  annotate("rect", xmin = (LON_MIN + LON_MAX) / 2 - 0.5,
           xmax = (LON_MIN + LON_MAX) / 2 + 0.5, ymin = -Inf, ymax = Inf,
           alpha = 0.10, fill = "#2166ac") +
  labs(x = paste0("Centre of the ", LON_MAX - LON_MIN, " deg longitude window"),
       y = "Slope of logR ~ z(I_PF)") +
  theme_paper() + theme(panel.grid = element_blank())

save_rev("glm_longitude_window_sensitivity.png", p_slide, width = 12, height = 6)

## ---------------------------------------------------------------------------
## PART 6. Q2 -- mean latitude of the PF vs SAM
## ---------------------------------------------------------------------------
# Three SAM series are available. The Jun-Dec mean is the one that matches the
# advection window (release weeks 25-29 plus 18 weeks), the annual mean is the
# conventional one, and DJF is the summer mode.
sam <- read_csv(file.path(sam_dir, "SAM_annual_2000_2023.csv"),
                show_col_types = FALSE) %>%
  rename(Year = year) %>%
  full_join(read_csv(file.path(sam_dir, "SAM_JunDec_2000_2023.csv"),
                     show_col_types = FALSE) %>% rename(Year = year),
            by = "Year") %>%
  full_join(read_csv(file.path(sam_dir, "SAM_seasonal_2000_2023.csv"),
                     show_col_types = FALSE) %>%
              rename(Year = year) %>%
              filter(season %in% c("DJF", "JJA", "SON")) %>%
              pivot_wider(names_from = season, values_from = SAM_season,
                          names_prefix = "SAM_"),
            by = "Year") %>%
  arrange(Year)

sam_dat <- indices %>%
  left_join(sam, by = "Year") %>%
  left_join(rec, by = "Year") %>%
  mutate(logR = log(Recruited_total))

write_csv(sam_dat, file.path(rev_dir, "pf_latitude_sam.csv"))

## Correlation table ------------------------------------------------------------
sam_vars <- c("SAM_annual", "SAM_JunDec", "SAM_DJF", "SAM_JJA", "SAM_SON")
pf_vars  <- c("lat_var", "lat_domain", "I_var", "I_domain", "I_static", "logR")

# The latitude of the detected front is bimodal (PART 3b), which by itself can
# create or destroy a Pearson correlation, so every pair is also computed on the
# 20 northern-branch years alone.
cor_tab <- expand_grid(pf = pf_vars, sam = sam_vars,
                       subset = c("all years", "northern-branch years")) %>%
  pmap_dfr(function(pf, sam, subset) {
    d <- if (subset == "all years") sam_dat else filter(sam_dat, branch == "northern")
    x <- d[[pf]]; y <- d[[sam]]
    ok <- is.finite(x) & is.finite(y)
    ct <- cor.test(x[ok], y[ok])
    cs <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
    tibble(pf_variable = pf, sam_index = sam, subset = subset, n = sum(ok),
           r = unname(ct$estimate), p = ct$p.value,
           rho = unname(cs$estimate), p_spearman = cs$p.value)
  }) %>%
  arrange(p)

write_csv(cor_tab, file.path(rev_dir, "correlations_pf_sam.csv"))
print(as.data.frame(cor_tab), digits = 3)

## Figures -----------------------------------------------------------------------
# Time series of the PF latitude and of the SAM, on twin scaled axes.
p_lat_ts <- sam_dat %>%
  dplyr::select(Year, lat_var, lat_domain, SAM_JunDec) %>%
  pivot_longer(-Year) %>%
  group_by(name) %>% mutate(z = as.numeric(scale(value))) %>% ungroup() %>%
  mutate(name = recode(name,
                       lat_var    = "PF latitude, 67-72 E",
                       lat_domain = "PF latitude, 61-79 E",
                       SAM_JunDec = "SAM (Jun-Dec)")) %>%
  ggplot(aes(Year, z, colour = name)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = c("#4575b4", "#91bfdb", "#d73027"), name = NULL) +
  labs(y = "Standardized anomaly", x = NULL) +
  theme_paper() + theme(legend.position = "top", panel.grid = element_blank())

save_rev("pf_latitude_sam_timeseries.png", p_lat_ts, width = 10, height = 5)

# Scatter, PF latitude against each SAM index.
p_lat_sc <- sam_dat %>%
  dplyr::select(Year, lat_var, lat_domain, all_of(sam_vars)) %>%
  pivot_longer(all_of(sam_vars), names_to = "sam_index", values_to = "SAM") %>%
  pivot_longer(c(lat_var, lat_domain), names_to = "window", values_to = "lat") %>%
  mutate(window = recode(window, lat_var = "67-72 E", lat_domain = "61-79 E")) %>%
  filter(is.finite(SAM), is.finite(lat)) %>%
  ggplot(aes(SAM, lat, colour = window)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15) +
  geom_point(size = 2) +
  facet_wrap(~ sam_index, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c("67-72 E" = "#4575b4", "61-79 E" = "#f46d43"),
                      name = "Longitude window") +
  labs(x = "SAM index", y = "Mean PF latitude (deg N)") +
  theme_paper() + theme(legend.position = "top", panel.grid = element_blank(),
                        axis.text = element_text(size = 11),
                        axis.title = element_text(size = 14),
                        strip.text = element_text(size = 11))

save_rev("pf_latitude_vs_sam.png", p_lat_sc, width = 13, height = 5)

## Map of the yearly front positions --------------------------------------------
p_map <- ggplot(pf_var, aes(lon, lat, group = Year, colour = Year)) +
  annotate("rect", xmin = LON_MIN, xmax = LON_MAX, ymin = -53, ymax = -47,
           fill = NA, colour = "black", linetype = 2) +
  geom_path(linewidth = 0.7, alpha = 0.85) +
  scale_colour_viridis_c(option = "plasma") +
  labs(x = "Longitude", y = "Latitude",
       subtitle = "Yearly PF position (WW edge, weeks 25-29); dashed box = paper window") +
  theme_paper() + theme(panel.grid = element_blank())

save_rev("pf_position_by_year.png", p_map, width = 10, height = 6)

## ---------------------------------------------------------------------------
## PART 7. Summary printed to the console
## ---------------------------------------------------------------------------
cat("\n================ SUMMARY ================\n")
cat("\nQ1. log(retention) ~ PF intensity\n")
print(as.data.frame(glm_compare %>%
  dplyr::select(model, n, slope, se, p, R2, ratio_per_SD)), digits = 3)

cat("\nQ2. PF latitude vs SAM (strongest first)\n")
print(as.data.frame(cor_tab %>% filter(pf_variable %in% c("lat_var", "lat_domain"))),
      digits = 3)

cat("\nOutputs in ", normalizePath(rev_dir, mustWork = FALSE), "\n", sep = "")
cat("Figures in  ", normalizePath(rev_fig, mustWork = FALSE), "\n", sep = "")
