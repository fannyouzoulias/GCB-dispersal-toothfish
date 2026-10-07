"""
Annual SAF streamline north of Kerguelen.

The Subantarctic Front (SAF), by the method PF_position_park.py applies to the
Polar Front: for each year, the ADT streamline associated with the SAF (the
annual SAF streamline), from the climatological SAF of Park et al. (2019), and
the SAF intensity: the mean surface geostrophic current speed along that
streamline.

NO PATHWAY CONSTRAINT. The PF streamline is pinned to the escarpment south of
Kerguelen. The SAF runs in deep water about 350 km north of the islands, with
no comparable choke point, so the reference contour is kept as it is every
year.

Reads what PF_position_park.py has downloaded, from the same ROOT: run
`python PF_position_park.py --download` first.

  1. CNES-CLS22 MDT: `mdt`
  2. DUACS L4: `adt`, `ugos`, `vgos`, over each larval-advection window
  3. Park & Durand (2019) ACC fronts, doi:10.17882/59800
     -> ROOT/park_durand_2019_ACC_fronts.nc

Writes, in ROOT/output_saf/ (`saf_dir` of config.R):

  saf_park_<year>.csv                    the annual SAF streamline, with the
                                         current speed along it
  saf_park_annual.csv                    one row per year
  saf_park_climatology_<first>-<last>.csv  the same streamline, read from the
                                         mean of the annual ADT fields

"""

import os
from datetime import datetime
from pathlib import Path

import contourpy
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Same ROOT, and same DUACS cache, as PF_position_park.py.
ROOT = Path(os.environ.get("PF_PARK_ROOT", "pf_park"))

YEARS = range(2000, 2024)
WINDOW = (23, 31, 18)  # first release week, last release week, drift duration (weeks)

ACC_BAND = (-56.0, -46.0)       # regional level removed from daily ADT
SMOOTH = 1.5                    # grid cells
SPEED_BAND = (63.0, 73.0)       # sector used for the SAF intensity
SPEED_LAT = (-47.0, -44.0)

YEAR_DIR = Path(os.environ.get("PF_DUACS_DIR", ROOT / "duacs_yearly"))
OUT = ROOT / "output_saf"
MDT_FILE = ROOT / "cnes_cls22_mdt.nc"
PARK_FILE = ROOT / "park_durand_2019_ACC_fronts.nc"


def thursday(year, week):
    """Thursday of week `week`, matching the advection notebooks."""
    return datetime.strptime(f"{year}-W{week}-4", "%Y-W%W-%w")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def need(path, message):
    if not path.exists():
        raise SystemExit(f"Missing {path}\n  -> {message}")
    return path


def fill_land(field):
    """Fill NaNs from neighbouring ocean values so contours can cross islands."""
    valid = np.isfinite(field)
    out = field.copy()

    for sigma in (1, 2, 4, 8):
        weight = gaussian_filter(valid.astype(float), sigma)
        value = gaussian_filter(np.where(valid, field, 0.0), sigma)
        value = value / np.where(weight > 1e-6, weight, np.nan)
        out = np.where(np.isfinite(out), out, value)

    return out


def contour_line(field, value, lon, lat):
    """Return the contour branch spanning the largest longitude range."""
    segments = [
        seg
        for seg in contourpy.contour_generator(x=lon, y=lat, z=field).lines(value)
        if len(seg) > 5
    ]
    if not segments:
        raise RuntimeError(f"No contour found for ADT={value:.4f} m")

    line = max(segments, key=lambda seg: np.ptp(seg[:, 0]))
    return line if line[0, 0] <= line[-1, 0] else line[::-1]


def in_sector(line):
    """Points of a line inside the sector of the SAF intensity."""
    return (
        (line[:, 0] >= SPEED_BAND[0])
        & (line[:, 0] <= SPEED_BAND[1])
        & (line[:, 1] >= SPEED_LAT[0])
        & (line[:, 1] <= SPEED_LAT[1])
    )


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
OUT.mkdir(parents=True, exist_ok=True)

need(MDT_FILE, "run: python PF_position_park.py --download")
need(PARK_FILE, "download Park & Durand (2019), doi:10.17882/59800")

mdt = xr.open_dataset(MDT_FILE)["mdt"].squeeze(drop=True).load()
mdt = mdt.sortby("latitude").sortby("longitude")
lat = mdt.latitude.values.astype(float)
lon = mdt.longitude.values.astype(float)

# Published climatological SAF in the Kerguelen sector
park = xr.open_dataset(PARK_FILE)
lo, la = park.LonSAF.values, park.LatSAF.values
idx = np.where((lo >= lon[0]) & (lo <= lon[-1]))[0]
lon_park = lo[idx.min() : idx.max() + 1]
lat_park = la[idx.min() : idx.max() + 1]

# Transfer the published SAF pathway onto the current CNES-CLS22 MDT field.
mdt_interp = RegularGridInterpolator(
    (lat, lon), mdt.values.astype(float), bounds_error=False, fill_value=np.nan
)
SAF_VALUE = float(np.nanmedian(mdt_interp(np.c_[lat_park, lon_park])))

# Reference regional sea level, as for the PF.
acc = (lat >= ACC_BAND[0]) & (lat <= ACC_BAND[1])
REFERENCE_LEVEL = float(np.nanmean(mdt.values[acc]))

print(f"Reference SAF-associated MDT contour: {100 * SAF_VALUE:.2f} cm")

rows = []
# Window-mean fields of every year, for the climatology written after the loop.
stack_adt, stack_u, stack_v = [], [], []

for year in YEARS:
    path = need(YEAR_DIR / f"duacs_{year}.nc", "run: python PF_position_park.py --download")
    ds = xr.open_dataset(path).sortby("latitude").sortby("longitude")

    t0, t1 = thursday(year, WINDOW[0]), thursday(year, WINDOW[1] + WINDOW[2])
    w = ds.sel(time=slice(t0, t1)).load()
    ds.close()

    if not w.sizes["time"]:
        raise SystemExit(f"{year}: no DUACS data between {t0:%Y-%m-%d} and {t1:%Y-%m-%d}")

    # Level each daily ADT field, then average over the larval-advection window.
    adt = w.adt.values.astype(float)
    daily_level = np.nanmean(adt[:, acc, :], axis=(1, 2))
    adt -= (daily_level - REFERENCE_LEVEL)[:, None, None]

    valid = np.isfinite(adt)
    annual_adt = np.where(
        valid.any(axis=0),
        np.nansum(adt, axis=0) / np.maximum(valid.sum(axis=0), 1),
        np.nan,
    )

    field = gaussian_filter(fill_land(annual_adt), SMOOTH)
    line = contour_line(field, SAF_VALUE, lon, lat)

    # Strength of the time-mean geostrophic current over the same advection window.
    mean_u = w.ugos.mean("time").values.astype(float)
    mean_v = w.vgos.mean("time").values.astype(float)
    mean_speed = np.hypot(mean_u, mean_v) * 100.0  # cm/s

    stack_adt.append(annual_adt)
    stack_u.append(mean_u)
    stack_v.append(mean_v)

    speed_at = RegularGridInterpolator(
        (lat, lon), mean_speed, bounds_error=False, fill_value=np.nan
    )
    line_speed = speed_at(line[:, ::-1])
    in_band = in_sector(line)

    pd.DataFrame(
        {"lon": line[:, 0], "lat": line[:, 1], "mean_current_speed_cm_s": line_speed}
    ).to_csv(OUT / f"saf_park_{year}.csv", index=False)

    record = {
        "Year": year,
        "n_days": int(w.sizes["time"]),
        "adt_contour_m": SAF_VALUE,
        "saf_current_cm_s": float(np.nanmean(line_speed[in_band])),
        "lat_mean": float(np.mean(line[in_band, 1])),
        "lat_south": float(np.min(line[in_band, 1])),
        "lat_north": float(np.max(line[in_band, 1])),
    }
    rows.append(record)

    print(
        f"{year}: {record['n_days']:3d} days; "
        f"{-record['lat_mean']:.2f} S; "
        f"SAF intensity {record['saf_current_cm_s']:.1f} cm/s"
    )

    w.close()

summary = pd.DataFrame(rows).set_index("Year")
summary.to_csv(OUT / "saf_park_annual.csv")

print(
    f"\nSAF intensity, {SPEED_BAND[0]:.0f}-{SPEED_BAND[1]:.0f} E: "
    f"median {summary.saf_current_cm_s.median():.1f} cm/s, "
    f"range {summary.saf_current_cm_s.min():.1f}-{summary.saf_current_cm_s.max():.1f}"
)
print(
    f"Latitude of the streamline in that sector: "
    f"{-summary.lat_north.max():.2f} to {-summary.lat_south.min():.2f} S"
)


# ---------------------------------------------------------------------------
# Climatology
# ---------------------------------------------------------------------------
# The same streamline, read from the mean of the annual window-mean ADT fields,
# every year weighing the same. Drawn on the front map of larval_dispersal/05.
clim_adt = np.nanmean(np.stack(stack_adt), axis=0)
clim_line = contour_line(gaussian_filter(fill_land(clim_adt), SMOOTH), SAF_VALUE, lon, lat)

clim_speed = np.hypot(
    np.nanmean(np.stack(stack_u), axis=0), np.nanmean(np.stack(stack_v), axis=0)
) * 100.0
clim_speed_at = RegularGridInterpolator(
    (lat, lon), clim_speed, bounds_error=False, fill_value=np.nan
)

pd.DataFrame(
    {
        "lon": clim_line[:, 0],
        "lat": clim_line[:, 1],
        "mean_current_speed_cm_s": clim_speed_at(clim_line[:, ::-1]),
    }
).to_csv(OUT / f"saf_park_climatology_{YEARS[0]}-{YEARS[-1]}.csv", index=False)

print(f"\nWritten to {OUT}")
