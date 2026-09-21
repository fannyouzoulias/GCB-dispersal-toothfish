# Larval dispersal and retention of Patagonian toothfish on the Kerguelen Plateau

R code supporting:

> Ouzoulias F. et al. *Climate-driven oceanographic changes shape larval dispersal success in the Southern Ocean*. **Global Change Biology**.

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
advection/          the Lagrangian simulation, and the front positions it is read against
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

**The front of each year.** The dispersal is read against two independent
definitions of the Polar Front: a hydrographic one, the northern limit of Winter
Water, and a dynamical one, a contour of sea surface topography. The first is
built by the notebooks below, in this order; the second is `PF_position_park.py`.

| Script | What it does |
|---|---|
| `PF_position_clean.ipynb` | flags, day by day, every grid point of GLORYS12V1 whose temperature minimum between 100 and 400 m is below 3 degrees and shallower than 350 m: the Winter Water signature of Azarian et al. (2024). Averaging those daily flags over the advection window gives a probability of Winter Water presence. One year per run (`year` is set by hand). Writes `WW_prob_presence_<year>_week_25-29.nc` |
| `PF_position_from_ww_prob_presence.ipynb` | turns that field into a line: per meridian, the northernmost cell above a probability of 0.8, kept only where the Winter Water block below it is 12 cells deep, so that detached patches further north are ignored, then smoothed with a 24-point running mean. Loops over 2000-2023 and writes `position_PF_<year>_week_25-29.csv`, the front used by `larval_dispersal/`. The `to_csv` call is commented out in the saved notebook |
| `average_velocity_field.ipynb` | averages the norm of the DUACS surface geostrophic velocity over the same window, one field per year, read through the LAMTA `loadCMEMSuv` function. Writes `<year>_uv_norm_field.nc` |
| `surface_currents_at_PF_location.ipynb` | samples that mean velocity field along the front line of each year: the surface proxy for front intensity the published index is built on |
| `PF_position_park.py` | the dynamical definition: the front as a single contour of dynamic topography, as Park et al. (2019) define it, with their constraint that it round Kerguelen from the south on the 500-1000 m escarpment. Self-contained, `--download` included. Writes one line per year, an annual table and a plate of 24 panels |

Note the file names. They all say `week_25-29`, but the window is the Thursday
of week 25 to the Thursday of week 29 + 18, i.e. the release weeks plus the
drift: roughly the end of June to the end of November.

`PF_position_park.py` runs on its own: it needs no other script and no module
of `front_position/`. Edit the `ROOT` path at the top, run
`python PF_position_park.py --download` once, download the Park & Durand front
file by hand from doi:10.17882/59800, then run the script. `--download` fetches
the CNES-CLS22 mean dynamic topography, the DUACS fields of each advection
window and the GLORYS12 bathymetry: a few hundred MB, a few minutes.


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

One category of input is *not* redistributed: **third-party products** (DUACS
geostrophic velocities, GEBCO bathymetry, Park & Durand front climatologies)
are cited in the paper and in the data archive, and must be downloaded from
their own repositories. Only two scripts read them directly: the advection
notebook needs the DUACS velocity fields, and `01_load_trajectories.R` needs the
GEBCO tile covering 55-89 E / 40-59 S to sample the seabed depth along the
trajectories. The `egg_production/` chain needs neither: the slope it uses is
already attached to the archived tables.


## Requirements

R >= 4.3 for both R chains. Packages are listed at the top of each `00_setup.R`;
the two spatial models load `sdmTMBextra` themselves. The advection notebook is
Python 3.

- **sdmTMBextra** is not on CRAN: `pak::pak("pbs-assess/sdmTMBextra")`. It
  provides the barrier mesh used by both `egg_production/` models.
- Coastlines come from `coastline.geojson` in the data archive. 
- **LAMTA** (Rousselet et al. 2025) is required by `advection/`, together with
  `numpy`, `pandas` and `geopandas`. It is not distributed on PyPI.

`sessionInfo.txt` records the exact environment the results were produced in.


## Contact

Fanny Ouzoulias - fanny.ouzoulias@mnhn.fr
Félix Massiot-Granier - felix.massiot-granier@mnhn.fr
Clara Péron - clara.peron@mnhn.fr
