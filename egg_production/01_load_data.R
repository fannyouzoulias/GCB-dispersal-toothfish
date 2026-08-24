################################################################################
# 01_load_data.R
#
# Loads the two fishery-observer tables the abundance and length models are
# fitted to. Both are in the data archive.
#
# They are extracts of the Pecheker database (MNHN / TAAF), already restricted
# to stage-4 (ripe) females sampled on the Kerguelen Plateau in weeks 22 to 30
# of 2020 to 2024, with the counts raised to the total catch of each station and
# the seabed slope attached. The raw extracts are not redistributed; the
# selection and the raising are documented in the data README.
#
# Produces: spawner_counts  (one row per longline station)
#           spawner_lengths (one row per measured fish)
#
# Requires: 00_setup.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-24
################################################################################

# nb_corr is the number of mature females at the station, raised to its total
# catch: only a fraction of the fish landed is measured, so the measured fish
# give the composition rather than the numbers. HAM_VIR is the effort, in hooks
# hauled. Position and depth are the mid-point of the longline.
spawner_counts <- readr::read_csv(dpath("spawner_counts.csv"),
                                  show_col_types = FALSE)

spawner_lengths <- readr::read_csv(dpath("spawner_lengths.csv"),
                                   show_col_types = FALSE)

## Coordinate offset -----------------------------------------------------------
# The two tables above are archived with their coordinates centred on the mean
# station position, so that the deposit alone does not circulate georeferenced
# longline positions. The offset lives here, with the code, and is added back
# before anything else runs. Everything downstream sees ordinary lon/lat.
#
# This is a separation of the two halves, not anonymisation: the relative
# geometry, the depths and the slopes are unchanged, and this file restores the
# absolute positions.
lon0 <- 69.5825763403
lat0 <- -48.3116530992

restore_coords <- function(df) {
  df$LON <- df$LON + lon0
  df$LAT <- df$LAT + lat0
  df
}

spawner_counts  <- restore_coords(spawner_counts)
spawner_lengths <- restore_coords(spawner_lengths)

stopifnot(
  all(c("LON", "LAT", "DEPTH", "pente", "nb_corr", "HAM_VIR") %in%
        names(spawner_counts)),
  all(c("LON", "LAT", "DEPTH", "pente", "LT") %in% names(spawner_lengths))
)
