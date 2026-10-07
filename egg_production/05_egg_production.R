################################################################################
# 05_egg_production.R
#
# Gridded egg production over the spawning area, combining the three models:
#
#     eggs per cell = N females (04) x eggs per female (02 length -> 03 fecundity)
#
# The three sources of uncertainty are propagated by multiplying the posterior
# draws cell by cell and draw by draw, rather than multiplying the means. The
# resulting grid is what weights the particles released in the Lagrangian
# simulations.
#
# Requires: 00_setup.R, 02_length_model.R, 03_fecundity_length.R,
#           04_abundance_model.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

length_pred    <- readRDS(file.path(out_dir, "length_predictions.rds"))
abundance_pred <- readRDS(file.path(out_dir, "abundance_predictions.rds"))
fit_eggs       <- readRDS(file.path(out_dir, "fecundity_model.rds"))

## Eggs per female, draw by draw -----------------------------------------------
# The fecundity relation is applied to each length draw, so a cell where length
# is poorly resolved ends up with a wide egg distribution.
a <- coef(fit_eggs)["a"]
b <- coef(fit_eggs)["b"]

LT_sims        <- exp(length_pred$logLT_sims)   # cells x draws
eggs_per_female <- a * exp(b * LT_sims)

## Total egg production per cell -----------------------------------------------
# Both matrices are [cells x draws] over the same grid, so the product is
# element-wise. Summarising after the product, not before, is what keeps the
# uncertainty honest.
N_sims <- as.matrix(abundance_pred$N_sims)
stopifnot(identical(dim(N_sims), dim(eggs_per_female)))

total_eggs_sims <- N_sims * eggs_per_female

grid_final <- abundance_pred$grid %>%
  mutate(LT_mean          = rowMeans(LT_sims),
         LT_sd            = apply(LT_sims, 1, sd),
         eggs_female_mean = rowMeans(eggs_per_female),
         total_oeufs_mean = rowMeans(total_eggs_sims),
         total_oeufs_sd   = apply(total_eggs_sims, 1, sd))

saveRDS(grid_final, file.path(out_dir, "egg_production_grid.rds"))
readr::write_csv(
  dplyr::select(grid_final, LON, LAT, DEPTH, pente,
                N_mean, total_oeufs_mean, total_oeufs_sd),
  file.path(out_dir, "egg_production_grid.csv")
)

## Maps ------------------------------------------------------------------------
# Panel D of the egg-production figure: the gridded field over the release area,
# without the recruitment sectors.
map_grid <- function(fill, name) {
  base_map +
    geom_tile(data = filter(grid_final, LON < 69),
              aes(x = LON, y = LAT, fill = {{ fill }})) +
    scale_fill_viridis_c(option = "inferno", name = name) +
    iso_layers +
    geom_sf(data = land_sf, fill = "darkgrey", colour = "darkgrey") +
    coord_sf(xlim = c(62, 70), ylim = c(-51.5, -47.9), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_paper() +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))
}

save_fig("egg_production_mean.png",
         map_grid(log(total_oeufs_mean), "Number of eggs (log)   "),
         width = 15, height = 10)

save_fig("egg_production_sd.png",
         map_grid(total_oeufs_sd, "sd   "),
         width = 15, height = 10)

## Egg production with the recruitment sectors --------------------------------
# The panel of the paper: the spawning hotspot together with the areas the
# larvae have to reach, which is what makes the retention results readable.
cols_recruitment <- c(east = "cadetblue3", north = "darkseagreen3",
                      west = "lightsalmon", south = "mistyrose3",
                      skiff = "lightgoldenrod1")
labels_recruitment <- c(east = "north-east", north = "north-west",
                        west = "west", south = "south", skiff = "skiff")

p_eggs_zones <- base_map +
  geom_tile(data = grid_final,
            aes(x = LON, y = LAT, fill = log(total_oeufs_mean))) +
  scale_fill_viridis_c(
    option = "inferno", name = "Number of eggs (log)",
    guide = guide_colorbar(barwidth = unit(5, "cm"), barheight = unit(0.6, "cm"),
                           title.position = "left")
  ) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = st_make_valid(recruitment_zones), aes(fill = zone),
          linewidth = 0.5, alpha = 0.55, inherit.aes = FALSE) +
  scale_fill_manual(
    name = "Recruitment\nzones", values = cols_recruitment,
    labels = labels_recruitment, breaks = names(labels_recruitment),
    guide = guide_legend(ncol = 2, byrow = FALSE, title.position = "left",
                         override.aes = list(alpha = 0.55))
  ) +
  iso_layers +
  geom_sf(data = land_sf, fill = "darkgrey", colour = "darkgrey") +
  coord_sf(xlim = c(62, 70), ylim = c(-51.5, -47.9), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_paper() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.7, "cm"),
        legend.spacing.x = unit(0.8, "cm"))

save_fig("egg_production_with_recruitment_zones.png", p_eggs_zones,
         width = 15, height = 10)

## Four-panel figure -----------------------------------------------------------
# A) female abundance (04), B) female length (02), C) fecundity-length relation
# (03), D) egg production. The first three plots are those built by 02-04, so
# this block only runs when they are in memory (i.e. the whole chain was run).
if (all(sapply(c("p_abundance", "p_length", "p_fecundity"), exists))) {
  # at half the width a longitude label every degree no longer fits, and the
  # 51.5 S label runs into the 62 E one at the corner
  panel_axes <- function(p) p +
    scale_x_continuous(breaks = seq(62, 70, by = 2)) +
    scale_y_continuous(breaks = seq(-51, -48, by = 0.5))

  # panel C takes the shape of the maps and a legend at the bottom like theirs,
  # so that the four panels have the same size and line up in rows and columns
  map_ratio <- diff(c(-51.5, -47.9)) /
    (diff(c(62, 70)) * cos(mean(c(-51.5, -47.9)) * pi / 180))
  p_fec_panel <- p_fecundity +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme(aspect.ratio = map_ratio, legend.position = "bottom")

  p_panels <- (panel_axes(p_abundance) | panel_axes(p_length)) /
    (p_fec_panel |
       panel_axes(map_grid(log(total_oeufs_mean), "Number of eggs (log)   "))) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag          = element_text(size = 30, face = "bold"),
          plot.tag.location = "panel",
          plot.tag.position = c(0.07, 0.92))

  save_fig("egg_production_panels.png", p_panels, width = 24, height = 17)
}

## Total over the spawning area ------------------------------------------------
total_by_draw <- colSums(total_eggs_sims)
message("Total egg production: mean = ", signif(mean(total_by_draw), 3),
        ", sd = ", signif(sd(total_by_draw), 3))
