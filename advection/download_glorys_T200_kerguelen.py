# python advection/download_glorys_T200_kerguelen.py [year ...]
"""
Daily potential temperature at 200 m from GLORYS12 (GLOBAL_MULTIYEAR_PHY_001_030,
cmems_mod_glo_phy_my_0.083deg_P1D-m) over the Kerguelen PF/jet box used by
larval_dispersal/09_mhw_kerguelen.R (65-72 E, 55-45 S), 1993-2023.

Each year: thetao on the 186.13 and 222.48 m levels, linear interpolation to
200 m (NaN where either level is missing, i.e. bottom shallower than 222 m),
3x3 block mean from 1/12 deg to 0.25 deg (cells with < 50 % valid pixels set to
missing), cached, then concatenated into
OUT_DIR/glorys12_T200m_1993_2023_Kerguelen_025deg.nc (variable 'temp'), which
09_mhw_kerguelen.R reads as `t200_dir` (config.R).
"""

import os
import sys
import asyncio
import numpy as np
import xarray as xr
import copernicusmarine

if sys.platform.startswith("win"):
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

DATASET_ID = "cmems_mod_glo_phy_my_0.083deg_P1D-m"
LON_MIN, LON_MAX = 65.0, 72.0
LAT_MIN, LAT_MAX = -55.0, -45.0
DEPTH_TARGET = 200.0
YEARS = range(1993, 2024)
COARSEN = 3  # 1/12 deg -> 0.25 deg

# Where the NetCDF and its caches are written: 33 MB for the result, more for
# the yearly caches, so keep it OUTSIDE the repository. Set it here, or in the
# GLORYS_T200_DIR environment variable. This is `t200_dir` in config.R.
OUT_DIR = os.environ.get("GLORYS_T200_DIR", "glorys_T200")

data_dir = OUT_DIR
os.makedirs(data_dir, exist_ok=True)
raw_dir = os.path.join(data_dir, "glorys_raw_tmp")
cache_dir = os.path.join(data_dir, "glorys_T200_025deg_yearly")
out_nc = os.path.join(data_dir, "glorys12_T200m_1993_2023_Kerguelen_025deg.nc")
for d in (raw_dir, cache_dir):
    os.makedirs(d, exist_ok=True)


def year_file(y):
    return os.path.join(cache_dir, f"glorys_T200_025deg_{y}.nc")


def process_year(y):
    out = year_file(y)
    if os.path.exists(out):
        print(f"{y}: cached", flush=True)
        return
    raw = os.path.join(raw_dir, f"glorys_raw_{y}.nc")
    if not os.path.exists(raw):
        copernicusmarine.subset(
            dataset_id=DATASET_ID,
            variables=["thetao"],
            minimum_longitude=LON_MIN, maximum_longitude=LON_MAX,
            minimum_latitude=LAT_MIN, maximum_latitude=LAT_MAX,
            minimum_depth=180.0, maximum_depth=225.0,
            start_datetime=f"{y}-01-01T00:00:00",
            end_datetime=f"{y}-12-31T23:59:59",
            coordinates_selection_method="inside",
            output_directory=raw_dir,
            output_filename=os.path.basename(raw),
            disable_progress_bar=True,
        )
    with xr.open_dataset(raw) as ds:
        th = ds["thetao"].load()
    depths = th["depth"].values
    if len(depths) != 2 or not (depths[0] < DEPTH_TARGET < depths[1]):
        raise RuntimeError(f"{y}: unexpected depth levels {depths}")
    w = (DEPTH_TARGET - depths[0]) / (depths[1] - depths[0])
    t200 = th.isel(depth=0) * (1 - w) + th.isel(depth=1) * w    # NaN if either level is NaN
    t200 = t200.drop_vars("depth", errors="ignore")
    t200 = t200.sortby("latitude").sortby("longitude")

    coarse = t200.coarsen(latitude=COARSEN, longitude=COARSEN, boundary="trim")
    valid = t200.notnull().coarsen(latitude=COARSEN, longitude=COARSEN, boundary="trim").mean()
    temp = coarse.mean(skipna=True).where(valid >= 0.5)
    temp = temp.rename("temp").rename({"latitude": "lat", "longitude": "lon"})
    temp["lat"] = np.round(temp["lat"].values.astype("float64"), 4)
    temp["lon"] = np.round(temp["lon"].values.astype("float64"), 4)
    temp["time"] = temp["time"].dt.floor("D")
    temp.attrs = {"units": "degC",
                  "long_name": "Potential temperature at 200 m (GLORYS12, linear interp. 186-222 m, 3x3 block mean)"}
    temp.to_dataset().to_netcdf(out, encoding={"temp": {"dtype": "float32", "zlib": True,
                                                        "complevel": 4, "_FillValue": -9.99}})
    os.remove(raw)
    print(f"{y}: {temp.sizes['time']} days, grid {temp.sizes['lat']}x{temp.sizes['lon']}, "
          f"T200 {float(temp.min()):.2f}..{float(temp.max()):.2f} degC, "
          f"masked cells {int(temp.isel(time=0).isnull().sum())}", flush=True)


def assemble():
    files = [year_file(y) for y in YEARS]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        print(f"not assembled, {len(missing)} years missing", flush=True)
        return
    ds = xr.open_mfdataset(files, combine="by_coords").load()
    t = ds["time"].values
    gaps = np.diff(t).astype("timedelta64[D]").astype(int)
    assert (gaps == 1).all(), "gaps in the daily series"
    ds.attrs = {
        "title": "GLORYS12 daily potential temperature at 200 m, Kerguelen PF/jet box 65-72E 55-45S, 0.25 deg",
        "source": f"Copernicus Marine {DATASET_ID} (GLOBAL_MULTIYEAR_PHY_001_030)",
        "processing": ("linear interpolation between the 186.13 and 222.48 m levels; 3x3 block mean "
                       "of the 1/12 deg grid; cells with < 50 % valid pixels set to missing"),
    }
    ds.to_netcdf(out_nc, encoding={"temp": {"dtype": "float32", "zlib": True, "complevel": 4,
                                            "_FillValue": -9.99},
                                   "time": {"units": "days since 1800-01-01", "dtype": "float64"}})
    print(f"written {out_nc}: {ds.sizes['time']} days {str(t[0])[:10]} -> {str(t[-1])[:10]}, "
          f"{os.path.getsize(out_nc) / 1e6:.1f} MB", flush=True)


if __name__ == "__main__":
    years = [int(a) for a in sys.argv[1:]] or list(YEARS)
    for y in years:
        process_year(y)
    if not sys.argv[1:]:
        assemble()
