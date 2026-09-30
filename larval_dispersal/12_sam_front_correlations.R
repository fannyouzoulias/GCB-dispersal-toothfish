################################################################################
# 12_sam_front_correlations.R
#
# Southern Annular Mode, front intensity and larval retention: the numbers of
# the SAM paragraph of the Results.
#
#   PART 2  SAM (annual and each season) against PF-associated jet intensity
#           and SAF intensity, Pearson correlations over 2000-2023
#   PART 3  pooled model: standardized PF-associated jet intensity against the
#           seasonal SAM index and season, one row per year x season
#   PART 4  SAM and PF-associated jet intensity against larval retention, on
#           two readings of retention:
#             logR      log of the annual number of retained eggs, the response
#                       of the GLM of 06 (so the PF test here has the p-value
#                       of Table S3)
#             ret_pct   annual retention rate (% of the eggs released), the
#                       series of the main retention figure of 04
#
# SAM: the monthly station-based index of Marshall (2003)
# (https://legacy.bas.ac.uk/met/gjma/sam.html), averaged by calendar year
# (SAM_annual) and by season. DJF of year y is December y-1 to February y, the
# austral summer before that year's spawning; MAM, JJA and SON are within y.
#
# The front indices are those of 05 (PF_park_cm_s along the annual PF
# streamline, SAF_mean_cm_s), and retention is summed over all suitable
# recruitment areas, Skiff Bank included, as in 04.
#
# The pooled model treats the four seasons of a year as separate rows, so its
# 96 rows are not independent (the response is the same annual PF value four
# times): read its p-value as indicative only.
#
# Writes, in outputs/revision/, sam_front_correlations.csv and
# sam_pf_pooled_model.csv.
#
# Requires: 00_setup.R, the CSVs of 04_retention_recruitment.R and
#           05_front_indices.R, and the SAM series in `sam_dir` (config.R)
#
# Author: Fanny Ouzoulias
# Date:   2026-09-28
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---------------------------------------------------------------------------
## PART 0. Paths
## ---------------------------------------------------------------------------
if (!exists("out_dir")) out_dir <- "outputs"
rev_dir <- file.path(out_dir, "revision")
dir.create(rev_dir, showWarnings = FALSE, recursive = TRUE)

seasons <- c("DJF", "MAM", "JJA", "SON")

## ---------------------------------------------------------------------------
## PART 1. Annual table
## ---------------------------------------------------------------------------
sam <- read_csv(file.path(sam_dir, "SAM_annual_2000_2023.csv"),
                show_col_types = FALSE) %>%
  rename(Year = year) %>%
  left_join(read_csv(file.path(sam_dir, "SAM_seasonal_2000_2023.csv"),
                     show_col_types = FALSE) %>%
              rename(Year = year) %>%
              filter(season %in% seasons) %>%
              pivot_wider(names_from = season, values_from = SAM_season,
                          names_prefix = "SAM_"),
            by = "Year")

dat <- read_csv(file.path(out_dir, "front_indices_annual.csv"),
                show_col_types = FALSE) %>%
  left_join(sam, by = "Year") %>%
  left_join(read_csv(file.path(out_dir, "recruited_annual.csv"),
                     show_col_types = FALSE), by = "Year") %>%
  left_join(read_csv(file.path(out_dir, "retention_annual.csv"),
                     show_col_types = FALSE) %>%
              dplyr::select(Year, ret_pct = mean_ret), by = "Year") %>%
  mutate(logR = log(Recruited_total)) %>%
  filter(Year %in% 2000:2023)

stopifnot(nrow(dat) == 24, !anyNA(dplyr::select(dat, PF_park_cm_s, SAF_mean_cm_s,
                                                SAM_annual, logR, ret_pct)))

cor_row <- function(x, y) {
  ct <- cor.test(dat[[x]], dat[[y]])
  tibble(x = x, y = y, n = sum(complete.cases(dat[[x]], dat[[y]])),
         r = unname(ct$estimate), p = ct$p.value)
}

## ---------------------------------------------------------------------------
## PART 2. SAM against front intensity
## ---------------------------------------------------------------------------
sam_vars <- c("SAM_annual", paste0("SAM_", seasons))

cor_fronts <- expand_grid(x = sam_vars, y = c("PF_park_cm_s", "SAF_mean_cm_s")) %>%
  pmap_dfr(cor_row) %>%
  mutate(block = "SAM vs front intensity")

## ---------------------------------------------------------------------------
## PART 3. Pooled model: PF intensity ~ seasonal SAM + season
## ---------------------------------------------------------------------------
pooled <- dat %>%
  mutate(z_PF = as.numeric(scale(PF_park_cm_s))) %>%
  dplyr::select(Year, z_PF, all_of(paste0("SAM_", seasons))) %>%
  pivot_longer(starts_with("SAM_"), names_to = "season", values_to = "SAM",
               names_prefix = "SAM_") %>%
  mutate(season = factor(season, levels = seasons))

m_pooled <- lm(z_PF ~ SAM + season, data = pooled)
pooled_tbl <- summary(m_pooled)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  as_tibble()
write_csv(pooled_tbl, file.path(rev_dir, "sam_pf_pooled_model.csv"))

## ---------------------------------------------------------------------------
## PART 4. SAM and PF intensity against retention
## ---------------------------------------------------------------------------
cor_ret <- expand_grid(x = c("SAM_annual", "PF_park_cm_s"),
                       y = c("logR", "ret_pct")) %>%
  pmap_dfr(cor_row) %>%
  mutate(block = "vs larval retention")

cor_all <- bind_rows(cor_fronts, cor_ret)
write_csv(cor_all, file.path(rev_dir, "sam_front_correlations.csv"))

## ---------------------------------------------------------------------------
## PART 5. Summary
## ---------------------------------------------------------------------------
cat("\n================ SAM, FRONTS AND RETENTION (2000-2023) ================\n")
print(as.data.frame(cor_all %>% mutate(r = round(r, 2), p = signif(p, 2))))
cat("\nPooled model: z(PF) ~ SAM + season\n")
print(as.data.frame(pooled_tbl), digits = 3)
