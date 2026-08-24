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
