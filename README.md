# Larval dispersal and retention of Patagonian toothfish on the Kerguelen Plateau

Code supporting:

> Ouzoulias F. et al. *Climate-driven oceanographic changes shape larval dispersal success in the Southern Ocean*. **Global Change Biology**.

## Overview

This repository contains the R and Python code used to model egg production and larval dispersal of Patagonian toothfish on the Kerguelen Plateau from 2000 to 2023.

The workflow combines:

1. spatial models of female abundance, length and fecundity;
2. Lagrangian simulations of the 18-week pelagic larval phase;
3. estimates of larval retention and recruitment;
4. analyses of their relationship with Polar Front intensity.

## Important note

One result differs slightly from the manuscript.

In this repository, Polar Front intensity is sampled along the front position of each year rather than along a fixed climatological front.

The relationship remains negative and significant:

- repository: slope = -0.167, p = 0.014
- manuscript: slope = -0.180, p = 0.007

See:

- `larval_dispersal/05_front_indices.R`
- `larval_dispersal/06_retention_front_glm.R`
- `larval_dispersal/10_pf_park_intensity_glm.R`

## Repository structure

```text
egg_production/     Spatial models of female abundance, length and fecundity
                    -> gridded egg production

advection/          Lagrangian simulations and reconstruction of annual fronts

larval_dispersal/   Trajectories, retention metrics, front indices
                    and retention-front analyses

config.R            Main configuration file for data paths

outputs/figures/    Output figures
```

The `egg_production/` and `larval_dispersal/` R workflows can be run independently.

## Workflow

### 1. Egg production

Run the scripts in `egg_production/` in numerical order:

```text
00_setup.R
01_load_data.R
02_length_model.R
03_fecundity_length.R
04_abundance_model.R
05_egg_production.R
```

These scripts estimate spatial patterns in female abundance, length and fecundity and combine them into gridded egg production.

### 2. Advection

`advection/Forward_advection_eggs.ipynb` releases particles weekly during spawning and tracks them for 18 weeks using DUACS surface geostrophic velocities.

The simulation uses:

- a 0.05° release grid;
- weekly releases during weeks 23-31;
- a 6-hour timestep;
- fourth-order Runge-Kutta integration.

The `advection/` folder also contains scripts used to reconstruct annual Polar Front positions and environmental fields.

`PF_position_park.py` reconstructs the annual Polar Front following Park et al. (2019) and provides the front used for the main intensity index.

Set `ROOT` in the script or define `PF_PARK_ROOT`, then run:

```bash
python PF_position_park.py --download
python PF_position_park.py
```

The Park & Durand front file must be downloaded separately from DOI `10.17882/59800`.

### 3. Larval dispersal

Run the scripts in `larval_dispersal/` in numerical order.

The main workflow:

1. builds and weights larval trajectories;
2. calculates dispersal and retention metrics;
3. derives annual front-intensity indices;
4. models retention as a function of Polar Front intensity.

Scripts `07` to `10` contain additional analyses added during manuscript revision.

`01_load_trajectories.R` is the most memory-intensive step. Once completed, its saved outputs can be loaded directly:

```r
source("larval_dispersal/00_setup.R")

df_traj <- readRDS(file.path(out_dir, "trajectories.rds"))
df_release <- readRDS(file.path(out_dir, "release_positions.rds"))
```

## Data

No data are stored directly in this repository.

The main datasets are archived on SEANOE and cited in the associated publication.

After downloading the archive, set `data_dir` in `config.R`.

`config.R` also contains paths to larger datasets reconstructed by the scripts in `advection/`, including annual front positions, velocity fields and temperature fields.

Some third-party products are not redistributed and must be obtained from their original sources, including:

- DUACS geostrophic velocities;
- GEBCO bathymetry;
- Park & Durand Polar Front climatologies.

## Requirements

### R

- R >= 4.3
- package requirements are listed in each `00_setup.R`

`sdmTMBextra` is not available on CRAN:

```r
pak::pak("pbs-assess/sdmTMBextra")
```

### Python

Python 3 is required for the advection workflow.

Main dependencies:

```text
numpy
pandas
geopandas
LAMTA
```

LAMTA is required for the advection simulations and is not distributed through PyPI.

`sessionInfo.txt` records the environment used to produce the results.

## Contact

Fanny Ouzoulias — fanny.ouzoulias@mnhn.fr  
Félix Massiot-Granier — felix.massiot-granier@mnhn.fr  
Clara Péron — clara.peron@mnhn.fr
