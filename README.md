# Larval dispersal and retention of Patagonian toothfish on the Kerguelen Plateau

R code supporting:

> Ouzoulias F. et al. *Climate-driven oceanographic changes shape larval dispersal success in the Southern Ocean]*. **Global Change Biology**.

The study combines (i) a spatially explicit model of egg production on the
northern Kerguelen Plateau and (ii) Lagrangian simulations of the 18-week
pelagic larval phase (2000-2023), to quantify larval retention on the shelf and
its relationship with Polar Front intensity.

## Contents

The repository follows the paper from the fishery data to the figures: the
egg-production models, the Lagrangian simulation they feed, and the analysis of
its output. The two R chains can also be run independently of each other.

```
egg_production/     spatial models of female abundance, length and fecundity
                    -> gridded egg production used to weight the released particles
advection/          the Lagrangian simulation itself (Python notebook)
larval_dispersal/   trajectories, retention metrics, front-intensity indices
                    and the retention ~ front model
config.R            the ONE file to edit: path to the data archive
outputs/figures/    where figures are written
```

### `egg_production/`

| Script | What it does |
|---|---|
| `00_setup.R` | packages, plotting theme, base map |
| `01_load_data.R` | loads the observer counts and lengths from the archive |
| `02_length_model.R` | spatial model of female length (sdmTMB, barrier mesh) |
| `03_fecundity_length.R` | fecundity-length relationship (exponential, nls) |
| `04_abundance_model.R` | spatial model of female abundance (sdmTMB, Tweedie) |
| `05_egg_production.R` | combines the three models into gridded egg production |

### `advection/`

`Forward_advection_eggs.ipynb` releases particles on a 0.05 degree grid inside
the release area, once a week on the Thursday of weeks 23 to 31, and advects each
release forward for 18 weeks on DUACS surface geostrophic velocities, with a
fourth-order Runge-Kutta scheme at a 6-hour timestep. It writes one CSV of daily
positions per release week and year: the input of `larval_dispersal/`.

Written by A. Nalivaev. It needs the **LAMTA** software (Rousselet et al. 2025),
which is not a pip package, plus `numpy`, `pandas` and `geopandas`.

### `larval_dispersal/`

| Script | What it does |
|---|---|
| `00_setup.R` | packages, theme, base map, sector polygons, front positions |
| `01_load_trajectories.R` | builds the trajectories from the raw advection output, adds seabed depth, applies beaching, weights each particle by egg production |
| `02_trajectory_maps.R` | dispersal pathways, weighted by egg production |
| `03_density_kde.R` | kernel densities of end-of-simulation positions |
| `04_retention_recruitment.R` | retention by sector, spawning-to-recruitment contributions, sensitivity analyses |
| `05_front_indices.R` | annual Polar Front / Subantarctic Front intensity indices |
| `06_retention_front_glm.R` | Gaussian GLM of log recruitment against front intensity |

Scripts are numbered in run order and are meant to be sourced in sequence within
a chain: `00_setup.R` first, then the others. `06_retention_front_glm.R` needs
the annual tables written by `04` and `05`.

`01_load_trajectories.R` is the slow one: it reads about 900 MB of raw advection
output, so it takes several minutes and a few GB of memory. It saves its two
tables to `outputs/`, so a later session can skip it:

```r
source("larval_dispersal/00_setup.R")
df_traj    <- readRDS(file.path(out_dir, "trajectories.rds"))
df_release <- readRDS(file.path(out_dir, "release_positions.rds"))
```

## Data

**No data are stored in this repository.** The datasets needed to run the code
(the raw output of the advection model, gridded egg production, front-intensity
indices, sector polygons, isobaths, coastline) are archived on SEANOE. The
deposit is cited in the associated publication.

Download the archive, then set `data_dir` in `config.R` to point at it. Every
script builds its paths from there.

The archive also holds the fishery-observer extracts that the egg-production
models are fitted to (station counts raised to the total catch, individual
lengths, fecundity samples), so both chains can be run from it.

One category of input is *not* redistributed: **third-party products** (DUACS
geostrophic velocities, GEBCO bathymetry, Park & Durand front climatologies)
are cited in the paper and in the data archive, and must be downloaded from
their own repositories. Only two scripts read them directly: the advection
notebook needs the DUACS velocity fields, and `01_load_trajectories.R` needs the
GEBCO tile covering 55-89 E / 40-59 S to sample the seabed depth along the
trajectories. The `egg_production/` chain needs neither: the slope it uses is
already attached to the archived tables.

## Figures

Figures are written to `outputs/figures/`. Colours, panel sizes and label
positions may differ slightly from the published versions, which received minor
cosmetic editing; multi-panel plates were assembled in a vector graphics editor
from the individual panels produced here.


## Requirements

R >= 4.3 for both R chains. Packages are listed at the top of each `00_setup.R`;
the two spatial models load `sdmTMBextra` themselves. The advection notebook is
Python 3.

- **sdmTMBextra** is not on CRAN: `pak::pak("pbs-assess/sdmTMBextra")`. It
  provides the barrier mesh used by both `egg_production/` models.
- Coastlines come from `coastline.geojson` in the data archive. Heard &
  McDonald are not in that layer and are taken from `rnaturalearth` at
  `scale = "medium"`, the finest scale available without the non-CRAN
  `rnaturalearthhires` package.
- **LAMTA** (Rousselet et al. 2025) is required by `advection/`, together with
  `numpy`, `pandas` and `geopandas`. It is not distributed on PyPI.

`sessionInfo.txt` records the exact environment the results were produced in.


## Contact

Fanny Ouzoulias - fanny.ouzoulias@mnhn.fr
Félix Massiot-Granier - felix.massiot-granier@mnhn.fr
