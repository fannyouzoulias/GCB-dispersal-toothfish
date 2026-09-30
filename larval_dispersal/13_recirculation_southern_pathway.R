################################################################################
# 13_recirculation_southern_pathway.R
#
# Trapping in the western recirculation cell against transport along the
# southern PF-associated pathway, year by year (reviewer point: the statistics
# of that relationship, now reported in the Results).
#
# Two annual indices, on the particles released in weeks 24-28:
#   Pct_trapped      % of released particles that enter the western
#                    recirculation cell (63-67 E, 50-48 S) and never leave it
#                    again before the end of the 18-week run
#   Pct_cross_south  % of released particles that pass at least once through
#                    the southern corridor (69-71 E, 51-49.5 S), the pathway
#                    along the PF south of the islands
# N_cross_south is the same count as a number; the number of particles
# released is the same every year, so the two give the same correlation.
#
# Both are particle counts, NOT weighted by egg production (unlike the
# retention of 04): they describe where the flow takes the particles, not how
# many eggs reach the nurseries.
#
# Then: Pearson (and Spearman) correlation between the two over 2000-2023.
#
# Writes outputs/revision/recirculation_southern_pathway_annual.csv,
#        outputs/revision/recirculation_southern_pathway_cor.csv and a figure.
#
# Requires: 00_setup.R, and the trajectories of 01_load_trajectories.R (read
#           from outputs/ if they are not in memory). Needs data.table.
#
# Author: Fanny Ouzoulias
# Date:   2026-09-28
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
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
if (!exists("theme_paper")) theme_paper <- function() theme_bw(base_size = 13)

release_weeks <- 24:28

# Boxes, in degrees: lon_min, lon_max, lat_min, lat_max
box_recirc <- c(63, 67, -50, -48)       # western recirculation cell
box_south  <- c(69, 71, -51, -49.5)     # southern PF-associated corridor

in_box <- function(lon, lat, b)
  lon >= b[1] & lon <= b[2] & lat >= b[3] & lat <= b[4]

## ---------------------------------------------------------------------------
## PART 1. Per-particle flags
## ---------------------------------------------------------------------------
if (!exists("df_traj")) df_traj <- readRDS(file.path(out_dir, "trajectories.rds"))

DT <- as.data.table(df_traj)[, .(Unique_ID, Year, Date, Longitude, Latitude, Week)]
DT <- DT[(Week - 1) %in% release_weeks]            # Week_release = Week - 1, as in 04
setkey(DT, Year, Unique_ID, Date)

DT[, `:=`(in_recirc = in_box(Longitude, Latitude, box_recirc),
          in_south  = in_box(Longitude, Latitude, box_south))]

# Trapped: once inside the cell, every later position is inside it too, i.e.
# the last position outside comes before the first position inside.
flags <- DT[, {
  first_in <- match(TRUE, in_recirc)
  .(trapped = !is.na(first_in) && all(in_recirc[first_in:.N]),
    south   = any(in_south))
}, by = .(Year, Unique_ID)]
rm(DT); gc(verbose = FALSE)

## ---------------------------------------------------------------------------
## PART 2. Annual indices and correlation
## ---------------------------------------------------------------------------
annual <- flags[, .(N_total         = .N,
                    N_trapped       = sum(trapped),
                    N_cross_south   = sum(south)), by = Year][order(Year)]
annual[, `:=`(Pct_trapped     = 100 * N_trapped / N_total,
              Pct_cross_south = 100 * N_cross_south / N_total)]
fwrite(annual, file.path(rev_dir, "recirculation_southern_pathway_annual.csv"))

ct <- cor.test(annual$Pct_trapped, annual$Pct_cross_south)
cs <- suppressWarnings(cor.test(annual$Pct_trapped, annual$Pct_cross_south,
                                method = "spearman", exact = FALSE))
cor_tbl <- tibble(n = nrow(annual),
                  pearson_r = unname(ct$estimate), pearson_p = ct$p.value,
                  ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
                  spearman_rho = unname(cs$estimate), spearman_p = cs$p.value)
write_csv(cor_tbl, file.path(rev_dir, "recirculation_southern_pathway_cor.csv"))

cat("\n======== WESTERN RECIRCULATION vs SOUTHERN PF PATHWAY ========\n")
print(as.data.frame(annual), digits = 3)
cat(sprintf("\nPearson r = %.2f (95%% CI %.2f to %.2f), p = %.2g, n = %d\n",
            cor_tbl$pearson_r, cor_tbl$ci_low, cor_tbl$ci_high,
            cor_tbl$pearson_p, cor_tbl$n))
cat(sprintf("Spearman rho = %.2f, p = %.2g\n", cor_tbl$spearman_rho, cor_tbl$spearman_p))

## ---------------------------------------------------------------------------
## PART 3. Figure
## ---------------------------------------------------------------------------
p <- ggplot(annual, aes(Pct_trapped, Pct_cross_south)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey40",
              fill = "grey85", linewidth = 0.8) +
  geom_point(size = 3, colour = "#6a3d9a") +
  geom_text(aes(label = Year), size = 3.2, vjust = -0.9, colour = "grey20") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 5,
           label = sprintf("r = %.2f, p = %.1g", cor_tbl$pearson_r, cor_tbl$pearson_p)) +
  labs(x = "Particles trapped in the western\nrecirculation cell (%)",
       y = "Particles advected along the southern\nPF-associated pathway (%)") +
  theme_bw() + theme_paper() + theme(panel.grid = element_blank())

ggsave(file.path(rev_fig, "recirculation_vs_southern_pathway.png"), p,
       width = 9, height = 7, dpi = 300)
