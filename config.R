################################################################################
# config.R
#
# Single place where the paths used by every script are defined.
# Edit `data_dir` to point to the data archive, then run scripts from the
# repository root (e.g. via the .Rproj file).
#
# The data are NOT in this repository. They are archived separately on SEANOE:
#   [SEANOE DOI to be added]
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

## Root of the downloaded data archive -----------------------------------------
data_dir <- "data"

## Where figures and intermediate objects are written --------------------------
out_dir <- "outputs"
fig_dir <- file.path(out_dir, "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Helper: build a path inside the data archive --------------------------------
dpath <- function(...) file.path(data_dir, ...)

if (!dir.exists(data_dir)) {
  warning("data_dir does not exist: ", normalizePath(data_dir, mustWork = FALSE),
          "\nDownload the data archive and set `data_dir` in config.R.")
}

## Front positions written by advection/ ---------------------------------------
# Not part of the data archive: these are rebuilt, year by year, by the scripts
# of advection/. Edit if you ran them somewhere else.
#
#   park_dir    the `output/` folder of advection/PF_position_park.py, i.e. its
#               ROOT + "/output". Holds pf_park_<year>.csv, the annual PF
#               contour and the current speed along it: THE FRONT INTENSITY
#               INDEX used by 05 and 06.
#   ww_dir      the Winter Water front of each year, written by
#               advection/PF_position_from_ww_prob_presence.ipynb.
#   ww_int_dir  that same front with the current speed sampled along it,
#               written by advection/surface_currents_at_PF_location.ipynb.
#   uv_dir      the mean velocity field of each advection window, written by
#               advection/average_velocity_field.ipynb. Only 10 reads it: it is
#               the one field all three sampling paths are compared in.
#   t200_dir    the GLORYS12 temperature at 200 m, written by
#               advection/download_glorys_T200_kerguelen.py (its OUT_DIR).
#               Only 09 reads it.
front_dir  <- "fronts"
park_dir   <- file.path(front_dir, "park")
ww_dir     <- file.path(front_dir, "winter_water")
ww_int_dir <- file.path(front_dir, "winter_water_intensity")
uv_dir     <- file.path(front_dir, "velocity_fields")
t200_dir   <- file.path(front_dir, "glorys_T200")

## Local overrides -------------------------------------------------------------
# config.local.R is gitignored. Put YOUR paths there, with the same variable
# names as above, and they take precedence over the placeholders. That way no
# tracked file has to carry a machine-specific path.
local_cfg <- if (requireNamespace("here", quietly = TRUE))
  here::here("config.local.R") else "config.local.R"
if (file.exists(local_cfg)) source(local_cfg)
rm(local_cfg)
