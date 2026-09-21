"""
Annual PF-associated ADT streamline around Kerguelen.

The Polar Front (PF) is a hydrographic boundary. This script derives an annual
ADT streamline associated with the PF, using the climatological PF and the
Kerguelen pathway constraint described by Park et al. (2019).

Downloads with --download:
  1. CNES-CLS22 MDT: `mdt`
  2. DUACS L4: `adt`, `ugos`, `vgos`, only over each larval-advection window
  3. GLORYS12 static bathymetry: `deptho`

Manual input:
  Park & Durand (2019) ACC fronts, doi:10.17882/59800
  -> ROOT/park_durand_2019_ACC_fronts.nc

"""

import sys
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
ROOT = Path(r"P:\MNHN\PhD\dispersion_larvaire\position_PF\park_standalone")

# Optional figure inputs
COASTLINE = Path(
    r"P:\MNHN\PhD\dispersion_larvaire\simus_dispersion\seanoe\coastline.geojson"
)
LARVAL_PATHWAY = None  # folder containing larval_pathway_<year>.csv, or None

YEARS = range(2000, 2024)
WINDOW = (23, 31, 18)  # first release week, last release week, drift duration (weeks)

BOX = dict(
    minimum_longitude=55,
    maximum_longitude=85,
    minimum_latitude=-58,
    maximum_latitude=-42,
)

# Copernicus Marine datasets
MDT_DATASET = "cnes_obs-sl_glo_phy-mdt_my_0.125deg_P20Y"
ADT_DATASET = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D"
BATHY_DATASET = "cmems_mod_glo_phy_my_0.083deg_static"

ACC_BAND = (-56.0, -46.0)       # regional level removed from daily ADT
SMOOTH = 1.5                    # grid cells
SPEED_BAND = (67.0, 72.0)       # sector used for PF-associated current intensity
ISLAND_LON = 69.0
ISLAND_LAT = (-49.8, -48.6)
ISOBATHS = (500.0, 1000.0)
TOL = 0.125                     # one DUACS grid cell
STEP = 0.001                    # m
MAX_SHIFT = 0.05                # m
XLIM, YLIM = (60.0, 85.0), (-54.0, -44.0)
DISPLACED_YEARS = (2012, 2013, 2017, 2023)  # flagged in the figure

YEAR_DIR = ROOT / "duacs_yearly"
OUT = ROOT / "output"
MDT_FILE = ROOT / "cnes_cls22_mdt.nc"
BATHY_FILE = ROOT / "glorys12_deptho.nc"
PARK_FILE = ROOT / "park_durand_2019_ACC_fronts.nc"


def thursday(year, week):
    """Thursday of week `week`, matching the advection notebooks."""
    return datetime.strptime(f"{year}-W{week}-4", "%Y-W%W-%w")


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
def download():
    import time
    import copernicusmarine as cm

    YEAR_DIR.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    tic = time.time()

    def grab(path, **kwargs):
        if path.exists():
            print(f"{path.name}: already there")
            return

        for attempt in range(5):
            try:
                ds = cm.open_dataset(**kwargs).load()
                encoding = {
                    name: dict(zlib=True, complevel=4, dtype="float32")
                    for name in ds.data_vars
                }
                ds.to_netcdf(path, encoding=encoding)
                print(
                    f"{path.name}: {path.stat().st_size / 1e6:.1f} MB, "
                    f"{time.time() - tic:.0f} s",
                    flush=True,
                )
                return
            except Exception as exc:
                print(
                    f"{path.name}: attempt {attempt + 1} failed: {exc!r}"[:160],
                    flush=True,
                )
                time.sleep(10)

        raise RuntimeError(f"{path.name}: five failed download attempts")

    print("(1) CNES-CLS22 mean dynamic topography")
    grab(
        MDT_FILE,
        dataset_id=MDT_DATASET,
        variables=["mdt"],
        **BOX,
    )

    print("(2) DUACS annual advection windows")
    w0, w1, nadv = WINDOW
    for year in YEARS:
        t0 = thursday(year, w0)
        t1 = thursday(year, w1 + nadv)
        grab(
            YEAR_DIR / f"duacs_{year}.nc",
            dataset_id=ADT_DATASET,
            variables=["adt", "ugos", "vgos"],
            **BOX,
            start_datetime=t0.strftime("%Y-%m-%d"),
            end_datetime=t1.strftime("%Y-%m-%d"),
        )

    print("(3) GLORYS12 static bathymetry")
    grab(
        BATHY_FILE,
        dataset_id=BATHY_DATASET,
        variables=["deptho"],
        **BOX,
    )

    print(f"\nDownloads finished in {(time.time() - tic) / 60:.1f} min.")

    if not PARK_FILE.exists():
        print(
            "\nManual input still missing:\n"
            f"  {PARK_FILE}\n"
            "Download Park & Durand (2019) from doi:10.17882/59800 "
            "and place the NetCDF file there."
        )


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


def crossings(line, longitude):
    """Latitudes at which a line crosses a given meridian."""
    a, b = line[:-1], line[1:]
    idx = np.where((a[:, 0] - longitude) * (b[:, 0] - longitude) <= 0)[0]
    idx = idx[a[idx, 0] != b[idx, 0]]

    if not len(idx):
        return np.array([])

    frac = (longitude - a[idx, 0]) / (b[idx, 0] - a[idx, 0])
    return a[idx, 1] + frac * (b[idx, 1] - a[idx, 1])


def side_of_islands(lats):
    """Return south, north, both, or none for crossings at 69 E."""
    if not len(lats):
        return "none"
    if np.all(lats < ISLAND_LAT[0]):
        return "south"
    if np.all(lats > ISLAND_LAT[1]):
        return "north"
    return "both"


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
OUT.mkdir(parents=True, exist_ok=True)

if "--download" in sys.argv:
    download()
    sys.exit()

need(MDT_FILE, "run: python PF_position_park.py --download")
need(BATHY_FILE, "run: python PF_position_park.py --download")
need(PARK_FILE, "download Park & Durand (2019), doi:10.17882/59800")

mdt = xr.open_dataset(MDT_FILE)["mdt"].squeeze(drop=True).load()
mdt = mdt.sortby("latitude").sortby("longitude")
lat = mdt.latitude.values.astype(float)
lon = mdt.longitude.values.astype(float)

# Published climatological PF in the Kerguelen sector
park = xr.open_dataset(PARK_FILE)
lo, la = park.LonPF.values, park.LatPF.values
idx = np.where((lo >= lon[0] - 2) & (lo <= lon[-1] + 2))[0]
lon_park = lo[idx.min() : idx.max() + 1]
lat_park = la[idx.min() : idx.max() + 1]

# Transfer the published PF pathway onto the current CNES-CLS22 MDT field.
mdt_interp = RegularGridInterpolator(
    (lat, lon), mdt.values.astype(float), bounds_error=False, fill_value=np.nan
)
PF_VALUE = float(np.nanmedian(mdt_interp(np.c_[lat_park, lon_park])))

# Bathymetric constraint south of Kerguelen
bathy = xr.open_dataset(BATHY_FILE)["deptho"].squeeze(drop=True).load()
bathy = bathy.sortby("latitude").sortby("longitude")
depth_at = RegularGridInterpolator(
    (bathy.latitude.values, bathy.longitude.values),
    bathy.values,
    bounds_error=False,
    fill_value=np.nan,
)


def isobath_lat(depth):
    """Latitude of an isobath on 69 E, immediately south of Kerguelen."""
    section = (
        bathy.sel(longitude=ISLAND_LON, method="nearest")
        .sel(latitude=slice(-51.5, ISLAND_LAT[0]))
    )
    slat = section.latitude.values
    sdep = section.values
    idx = np.where((sdep[:-1] >= depth) & (sdep[1:] < depth))[0][-1]
    return float(
        slat[idx]
        + (depth - sdep[idx])
        / (sdep[idx + 1] - sdep[idx])
        * (slat[idx + 1] - slat[idx])
    )


LAT_500 = isobath_lat(ISOBATHS[0])
LAT_1000 = isobath_lat(ISOBATHS[1])


def on_escarpment(line):
    """
    Kerguelen constraint used when the reference contour rounds the islands
    from the north: all 69 E crossings must be south of the islands, and the
    nearest crossing must lie on the 500-1000 m escarpment (± one grid cell).
    """
    c = crossings(line, ISLAND_LON)
    return (
        side_of_islands(c) == "south"
        and LAT_1000 - TOL <= c.max() <= LAT_500 + TOL
    )


# Reference regional sea level. Subtracting a spatially uniform daily offset
# changes contour labels but not horizontal ADT gradients or geostrophic flow.
acc = (lat >= ACC_BAND[0]) & (lat <= ACC_BAND[1])
REFERENCE_LEVEL = float(np.nanmean(mdt.values[acc]))

print(f"Reference PF-associated MDT contour: {100 * PF_VALUE:.2f} cm")
print(
    f"69 E escarpment: 500 m at {LAT_500:.2f}, 1000 m at {LAT_1000:.2f}; "
    f"accepted {LAT_1000 - TOL:.2f} to {LAT_500 + TOL:.2f}"
)

w0, w1, nadv = WINDOW
rows = []
lines = {}

for year in YEARS:
    path = need(
        YEAR_DIR / f"duacs_{year}.nc",
        "run: python PF_position_park.py --download",
    )
    ds = xr.open_dataset(path).sortby("latitude").sortby("longitude")

    t0 = thursday(year, w0)
    t1 = thursday(year, w1 + nadv)
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

    fixed = contour_line(field, PF_VALUE, lon, lat)
    shift = 0.0

    # Keep the reference contour when it already rounds Kerguelen from the south.
    # Otherwise move to the nearest lower contour satisfying the Park et al.
    # topographic pathway constraint.
    if side_of_islands(crossings(fixed, ISLAND_LON)) != "south":
        for delta in np.arange(STEP, MAX_SHIFT + STEP / 2, STEP):
            candidate = contour_line(field, PF_VALUE - delta, lon, lat)
            if on_escarpment(candidate):
                shift = -delta
                break
        else:
            print(
                f"{year}: no contour within {100 * MAX_SHIFT:.0f} cm satisfies "
                "the Kerguelen constraint; keeping the reference contour"
            )

    line = contour_line(field, PF_VALUE + shift, lon, lat) if shift else fixed

    # Strength of the time-mean geostrophic current over the same advection window.
    mean_u = w.ugos.mean("time").values.astype(float)
    mean_v = w.vgos.mean("time").values.astype(float)
    mean_speed = np.hypot(mean_u, mean_v) * 100.0  # cm/s

    speed_at = RegularGridInterpolator(
        (lat, lon), mean_speed, bounds_error=False, fill_value=np.nan
    )
    line_speed = speed_at(line[:, ::-1])
    line_depth = depth_at(line[:, ::-1])

    in_band = (line[:, 0] >= SPEED_BAND[0]) & (line[:, 0] <= SPEED_BAND[1])

    pd.DataFrame(
        {
            "lon": line[:, 0],
            "lat": line[:, 1],
            "mean_current_speed_cm_s": line_speed,
            "depth_m": line_depth,
        }
    ).to_csv(OUT / f"pf_park_{year}.csv", index=False)

    record = {
        "Year": year,
        "n_days": int(w.sizes["time"]),
        "adt_contour_m": PF_VALUE + shift,
        "contour_shift_cm": 100 * shift,
        "side_reference_contour": side_of_islands(crossings(fixed, ISLAND_LON)),
        "side_selected_streamline": side_of_islands(crossings(line, ISLAND_LON)),
        "pf_associated_current_cm_s": float(np.nanmean(line_speed[in_band])),
        "share_shallower_500m": float(np.nanmean(line_depth[in_band] < 500)),
    }

    for x in (65.0, ISLAND_LON, 72.0):
        record[f"lat_{x:.0f}E"] = "/".join(f"{v:.2f}" for v in crossings(line, x))

    rows.append(record)
    lines[year] = (line, fixed if shift else None)

    print(
        f"{year}: {record['n_days']:3d} days; contour "
        f"{100 * record['adt_contour_m']:6.2f} cm "
        f"({record['contour_shift_cm']:+.1f}); "
        f"{record['side_selected_streamline']:5s} of Kerguelen; "
        f"PF-associated current {record['pf_associated_current_cm_s']:.1f} cm/s"
    )

    w.close()

summary = pd.DataFrame(rows).set_index("Year")
summary.to_csv(OUT / "pf_park_annual.csv")

moved = summary.contour_shift_cm != 0
print(
    f"\nContour adjusted in {moved.sum()} years: "
    + ", ".join(
        f"{year} ({summary.contour_shift_cm[year]:+.1f} cm)"
        for year in summary.index[moved]
    )
)
print(
    "South of Kerguelen after selection: "
    f"{(summary.side_selected_streamline == 'south').sum()} / {len(summary)}"
)
print(
    f"PF-associated current, {SPEED_BAND[0]:.0f}-{SPEED_BAND[1]:.0f} E: "
    f"median {summary.pf_associated_current_cm_s.median():.1f} cm/s, "
    f"range {summary.pf_associated_current_cm_s.min():.1f}-"
    f"{summary.pf_associated_current_cm_s.max():.1f}"
)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
if "--no-figure" in sys.argv:
    print(f"\nWritten to {OUT}")
    sys.exit()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INK = "#0b0b0b"
STREAM = "#0072B2"
PATH = "#d95f02"
FLAG = "#c0392b"

coast = None
if COASTLINE is not None and Path(COASTLINE).exists():
    try:
        import geopandas as gpd

        coast = gpd.read_file(COASTLINE)
    except ImportError:
        print("geopandas not installed: coastline omitted")

ncol = 4
nrow = int(np.ceil(len(list(YEARS)) / ncol))
fig, axes = plt.subplots(
    nrow,
    ncol,
    figsize=(3.9 * ncol, 2.25 * nrow),
    sharex=True,
    sharey=True,
    gridspec_kw=dict(hspace=0.28, wspace=0.06),
)

for ax, year in zip(axes.ravel(), YEARS):
    line, fixed = lines[year]

    ax.plot(lon_park, lat_park, color=INK, lw=0.8, ls="--", zorder=3)

    if fixed is not None:
        ax.plot(fixed[:, 0], fixed[:, 1], color=STREAM, lw=1.0, ls=":", zorder=3.5)

    ax.plot(line[:, 0], line[:, 1], color="white", lw=4.2, zorder=4)
    ax.plot(line[:, 0], line[:, 1], color=STREAM, lw=2.6, zorder=5)

    if LARVAL_PATHWAY is not None:
        f = Path(LARVAL_PATHWAY) / f"larval_pathway_{year}.csv"
        if f.exists():
            pathway = pd.read_csv(f)
            ax.plot(pathway.lon, pathway.lat, color=PATH, lw=1.1, zorder=6)

    if coast is not None:
        coast.plot(ax=ax, color="#6b6b68", linewidth=0, zorder=8)

    value = summary.loc[year, "pf_associated_current_cm_s"]
    displaced = year in DISPLACED_YEARS
    ax.set_title(
        f"{year}   {value:.1f} cm/s",
        fontsize=8.5,
        pad=3,
        color=FLAG if displaced else INK,
        fontweight="bold" if displaced else "normal",
    )

    for spine in ax.spines.values():
        spine.set_edgecolor(FLAG if displaced else "#c3c2b7")
        spine.set_linewidth(1.6 if displaced else 0.8)

    if summary.loc[year, "contour_shift_cm"]:
        ax.text(
            0.98,
            0.04,
            f"{100 * summary.loc[year, 'adt_contour_m']:.1f} cm",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.5,
            color=STREAM,
        )

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.tick_params(labelsize=7)

for ax in axes.ravel()[len(list(YEARS)) :]:
    ax.set_visible(False)

fig.subplots_adjust(top=0.94, bottom=0.075, left=0.04, right=0.98)
fig.suptitle(
    "Annual PF-associated ADT streamline around Kerguelen\n"
    f"(larval-advection window: weeks {w0}-{w1} + {nadv} weeks)",
    fontsize=12,
    y=0.99,
)

handles = [
    Line2D([], [], color=STREAM, lw=2.4, label="annual PF-associated ADT streamline"),
    Line2D(
        [],
        [],
        color=STREAM,
        lw=1.0,
        ls=":",
        label="reference contour when replaced by the Kerguelen constraint",
    ),
    Line2D(
        [],
        [],
        color=INK,
        lw=1.0,
        ls="--",
        label="climatological PF (Park & Durand, 2019)",
    ),
    Line2D([], [], color=FLAG, lw=1.6, label="the four displaced years"),
]
if LARVAL_PATHWAY is not None:
    handles.insert(
        2,
        Line2D([], [], color=PATH, lw=2.0, label="annual larval pathway"),
    )

fig.legend(
    handles=handles,
    loc="lower center",
    ncol=2,
    fontsize=8.5,
    frameon=False,
    bbox_to_anchor=(0.5, 0.0),
)

fig.savefig(OUT / "pf_park_by_year.png", dpi=160)
print(f"\nWritten to {OUT}")
