#!/usr/bin/env python3
"""
emit_ratio_check.py v0.2: do the standard Landsat / Sentinel-2 mineral band ratios agree with
NASA EMIT imaging-spectrometer mineral maps?

For each site the script automatically
  1. finds EMIT L2B mineral scenes (EMITL2BMIN) in the site's dry season (NASA Earthdata, free login)
  2. pairs each with the closest clear Sentinel-2 L2A and Landsat 8/9 C2 L2 acquisition (Planetary Computer)
  3. keeps the best N scene pairs with different dates (default 2) for repeatability
  4. downloads EMIT minerals + their fit uncertainty + cloud/water mask (+ optionally EMIT reflectance),
     orthorectifies them (GLT) and computes the spectral abundance of the 10 EMIT minerals
  5. reads Sentinel-2/Landsat onto the EMIT 60 m grid (area-averaged), masks cloud/shadow/water/vegetation
  6. tests every classic ratio:
       MAIN test     clear presence (abundance >= T, confident fit) vs clear absence (abundance = 0);
                     weak/uncertain detections are left out (grey zone). ROC AUC + spatial block-bootstrap CI.
       PURE test     pixels where the target makes up >= 70 % of the EMIT absorption vs absence
       SPECIFICITY   pure target vs pure other mineral classes (e.g. carbonate vs clay)
       SENSITIVITY   thresholds T, with/without the fit filter
       PHYSICS       the same ratios computed from EMIT's own reflectance resampled to Sentinel-2/Landsat
                     bands (same pixels, same time): separates "the ratio is physically wrong" from
                     "the multispectral data are different from EMIT"
  7. writes one report per scene and an index.html comparing all sites (opens automatically)

Usage (Windows CMD):
  python emit_ratio_check.py                                   (all preset sites, 2 scenes each)
  python emit_ratio_check.py --site ain_chair cuprite         (selected sites)
  python emit_ratio_check.py --list-sites
  python emit_ratio_check.py --site figuig --list-only        (show candidate scene pairs, no download)
  python emit_ratio_check.py --emit-reflectance none          (skip the physics control, much smaller download)

Credentials: NASA Earthdata account (https://urs.earthdata.nasa.gov/users/new) - put the user name (line 1)
and password (line 2) in earthdata_login.txt next to this script, or type them when asked.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import socket
import sys
import time
import webbrowser
from pathlib import Path

import numpy as np
import pandas as pd
def _use_bundled_proj_data():
    """Other software (PostGIS, old GDAL installs) can set PROJ_LIB/PROJ_DATA to an old proj.db,
    which breaks rasterio ("DATABASE.LAYOUT.VERSION.MINOR"). Use the data shipped inside rasterio."""
    import importlib.util
    spec = importlib.util.find_spec("rasterio")
    if spec and spec.submodule_search_locations:
        pkg = list(spec.submodule_search_locations)[0]
        proj, gdal = os.path.join(pkg, "proj_data"), os.path.join(pkg, "gdal_data")
        if os.path.exists(os.path.join(proj, "proj.db")):
            os.environ["PROJ_DATA"] = os.environ["PROJ_LIB"] = proj
        if os.path.isdir(gdal):
            os.environ["GDAL_DATA"] = gdal


_use_bundled_proj_data()
import rasterio
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject, transform_bounds
from rasterio.windows import from_bounds
from scipy.stats import spearmanr, rankdata
from shapely.geometry import box, shape, Polygon
from shapely.ops import unary_union
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

__version__ = "0.3.8"
HERE = os.path.dirname(os.path.abspath(__file__))
GROUPING_URL = ("https://raw.githubusercontent.com/emit-sds/emit-sds-l2b/develop/data/"
                "mineral_grouping_matrix_20230503.csv")
EMIT10 = ["Calcite", "Chlorite", "Dolomite", "Goethite", "Gypsum", "Hematite",
          "Illite+Muscovite", "Kaolinite", "Montmorillonite", "Vermiculite"]
GROUP_OF = {m: (1 if m in ("Goethite", "Hematite") else 2) for m in EMIT10}   # EMIT spectral group of each mineral
CLASSES = {
    "Carbonate": ["Calcite", "Dolomite"],
    "Clay (Al-OH)": ["Kaolinite", "Montmorillonite", "Illite+Muscovite"],
    "Fe-oxide": ["Hematite", "Goethite"],
    "Gypsum": ["Gypsum"],
    "Chlorite/Vermiculite": ["Chlorite", "Vermiculite"],
}
for _c, _m in CLASSES.items():
    GROUP_OF[_c] = GROUP_OF[_m[0]]
GROUP2_CLASSES = ["Carbonate", "Clay (Al-OH)", "Gypsum", "Chlorite/Vermiculite"]
ALL_TARGETS = list(dict.fromkeys(list(CLASSES) + EMIT10))

# Sites: box (W, S, E, N) and the months used (dry season, high sun, little vegetation).
SITES = {
    "ain_chair":     dict(bbox=[-2.75, 32.05, -2.30, 32.50], months=[4, 5, 6, 7, 8, 9, 10], country="Morocco",
                          note="Tamlelt inlier: limestones, dolomites, quartzites, schists (thesis site)"),
    "figuig":        dict(bbox=[-1.57, 32.17, -1.13, 32.52], months=[4, 5, 6, 7, 8, 9, 10], country="Morocco",
                          note="Figuig: Jurassic limestones/dolomites, marls, sands (centre 32.3436 N, 1.3474 W)"),
    "tiwit":         dict(bbox=[-6.05, 30.98, -5.62, 31.34], months=[4, 5, 6, 7, 8, 9, 10], country="Morocco",
                          note="Tiwit / Iknioun (Saghro): volcanics, alteration, Fe-oxides (centre 31.1606 N, 5.8334 W)"),
    "cuprite":       dict(bbox=[-117.32, 37.45, -117.12, 37.62], months=[5, 6, 7, 8, 9, 10], country="USA",
                          note="Cuprite, Nevada: kaolinite/alunite alteration, limestone, muscovite (benchmark site)"),
    "mountain_pass": dict(bbox=[-115.66, 35.38, -115.42, 35.58], months=[4, 5, 6, 7, 8, 9, 10], country="USA",
                          note="Mountain Pass, California: limestone vs dolomite, carbonatite (Rowan & Mars 2003)"),
    "rosh_pinah":    dict(bbox=[16.55, -28.15, 17.00, -27.75], months=[10, 11, 12, 1, 2, 3], country="Namibia",
                          note="Rosh Pinah, Gariep Belt: dolostones, schists, quartzites (hyper-arid; summer = dry season)"),
    # v0.3.8 candidates (chlorite / Mg-OH and Mg-smectite end-members) - screen them first with --screen
    "bou_azzer":     dict(bbox=[-6.85, 30.30, -6.40, 30.65], months=[3, 4, 5, 6, 7, 8, 9, 10, 11], country="Morocco",
                          note="Bou Azzer-El Graara inlier (Anti-Atlas): serpentinized Neoproterozoic ophiolite, "
                               "carbonate-chlorite-epidote alteration (chlorite / Mg-OH end-member)"),
    "jbel_rhassoul": dict(bbox=[-4.55, 32.80, -4.10, 33.10], months=[4, 5, 6, 7, 8, 9, 10], country="Morocco",
                          note="Missour basin, Jbel Rhassoul (Ghassoul): stevensite (Mg-smectite), secondary gypsum, "
                               "dolomite (lacustrine Miocene)"),
    # v0.3.1 candidates (calcite / gypsum end-members) - check them first with --screen
    "hamar_laghdad": dict(bbox=[-4.25, 31.25, -3.90, 31.52], months=[4, 5, 6, 7, 8, 9, 10], country="Morocco",
                          note="Hamar Laghdad / Kess-Kess, 18 km SE of Erfoud (31.38 N, 4.07 W): Devonian microsparitic "
                               "limestone mud mounds, crinoidal limestone, shales and marls (calcite end-member)"),
    "safi_gypsum":   dict(bbox=[-9.35, 32.05, -8.85, 32.50], months=[5, 6, 7, 8, 9, 10], country="Morocco",
                          note="Safi basin: Upper Jurassic evaporite formation (gypsum-carbonate, up to 200 m), "
                               "Lalla Fatna cliffs, Sidi Tiji gypsum quarries (gypsum end-member)"),
    "white_sands":   dict(bbox=[-106.50, 32.65, -106.15, 33.00], months=[4, 5, 6, 7, 8, 9, 10], country="USA",
                          note="White Sands, New Mexico: gypsum dune field, Alkali Flat, Lake Lucero (gypsum benchmark)",
                          emit_cloud_mask=False),   # EMIT's cloud test flags the bright gypsum dunes as cloud
}

RATIOS = {
    "sentinel-2": {
        "Ferric iron B4/B2": dict(num="B04", den="B02", claim="Fe-oxide"),
        "Clay/OH B11/B12": dict(num="B11", den="B12", claim="Clay (Al-OH)"),
        "Carbonate B12/B11": dict(num="B12", den="B11", claim="Carbonate"),
    },
    "landsat": {
        "Ferric iron B4/B2": dict(num="red", den="blue", claim="Fe-oxide"),
        "Clay/OH B6/B7": dict(num="swir16", den="swir22", claim="Clay (Al-OH)"),
        "Carbonate B7/B6": dict(num="swir22", den="swir16", claim="Carbonate"),
    },
}
# Band centre and full width at half maximum (nm), used to resample EMIT reflectance (flat-top response)
BANDPASS = {
    "sentinel-2": {"B02": (492, 66), "B03": (560, 36), "B04": (665, 31), "B08": (833, 106),
                   "B11": (1614, 91), "B12": (2202, 175)},
    "landsat": {"blue": (482, 60), "green": (561, 57), "red": (655, 37), "nir08": (865, 28),
                "swir16": (1609, 85), "swir22": (2201, 187)},
}
SIM_PREFIX = "EMIT-as-"


def base_sensor(s):
    return s.replace(SIM_PREFIX, "")


_ASCII = {"≥": ">=", "≤": "<=", "×": "x", "–": "-", "—": "-", "ρ": "rho", "→": "->"}


def log(msg):
    msg = str(msg)
    for k, v in _ASCII.items():
        msg = msg.replace(k, v)
    print(time.strftime("%H:%M:%S"), msg, flush=True)


for _stream in (sys.stdout, sys.stderr):      # Windows console/PowerShell redirection: never crash on a character
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass


# =============================================================================
# Grid
# =============================================================================
class Grid:
    def __init__(self, transform: Affine, width: int, height: int, crs="EPSG:4326"):
        self.transform, self.width, self.height, self.crs = transform, width, height, crs

    @property
    def bounds(self):
        t = self.transform
        return (t.c, t.f + self.height * t.e, t.c + self.width * t.a, t.f)


def write_tif(path, arr, grid, nodata=np.nan):
    with rasterio.open(path, "w", driver="GTiff", height=grid.height, width=grid.width, count=1,
                       dtype="float32", crs=grid.crs, transform=grid.transform, nodata=nodata,
                       compress="deflate") as dst:
        dst.write(np.asarray(arr, dtype=np.float32), 1)


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").lower()


# =============================================================================
# 1. NASA Earthdata login, EMIT search and download
# =============================================================================
CRED_FILE = os.path.join(HERE, "earthdata_login.txt")


def _try_login(strategy, persist=False):
    import earthaccess
    try:
        auth = earthaccess.login(strategy=strategy, persist=persist)
        return auth is not None and getattr(auth, "authenticated", False)
    except Exception as e:
        log(f"  Earthdata login ({strategy}) failed: {e}")
        return False


def earthdata_login():
    import earthaccess  # noqa: F401
    if os.environ.get("EARTHDATA_USERNAME") and os.environ.get("EARTHDATA_PASSWORD"):
        if _try_login("environment", persist=True):
            return
    if os.path.exists(CRED_FILE):
        lines = [l.strip() for l in open(CRED_FILE, encoding="utf-8-sig").read().splitlines() if l.strip()]
        if len(lines) >= 2:
            os.environ["EARTHDATA_USERNAME"], os.environ["EARTHDATA_PASSWORD"] = lines[0], lines[1]
            log(f"  using credentials from {os.path.basename(CRED_FILE)} (user '{lines[0]}')")
            if _try_login("environment", persist=True):
                return
            log("  ! the user name/password in earthdata_login.txt were refused")
    elif _try_login("netrc"):
        return
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        raise SystemExit("Earthdata login needed, but the output is redirected to a file, so the password prompt "
                         "would be invisible.\nPut earthdata_login.txt (line 1: user name, line 2: password) next to "
                         f"the script ({os.path.dirname(CRED_FILE) or '.'}) or run once without '> file'.")
    print("\nNASA Earthdata login (https://urs.earthdata.nasa.gov) - not the USGS EarthExplorer account.")
    print("The password WILL BE VISIBLE on screen while you type it.")
    for _ in range(3):
        user = input("Earthdata user name: ").strip()
        pwd = input("Earthdata password (visible): ").strip()
        os.environ["EARTHDATA_USERNAME"], os.environ["EARTHDATA_PASSWORD"] = user, pwd
        if _try_login("environment", persist=True):
            log("  login OK (stored for next runs)")
            return
        print("  Refused. Use the NASA Earthdata user name (not e-mail), and confirm the account via NASA's e-mail.")
    raise SystemExit("Earthdata login failed 3 times.")


def _granule_time(g):
    te = g["umm"].get("TemporalExtent", {})
    s = te.get("RangeDateTime", {}).get("BeginningDateTime") or te.get("SingleDateTime")
    t = pd.Timestamp(s)
    return t.tz_convert(None) if t.tzinfo else t


def _granule_footprint(g):
    try:
        pts = g["umm"]["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]["GPolygons"][0]["Boundary"]["Points"]
        return Polygon([(p["Longitude"], p["Latitude"]) for p in pts])
    except (KeyError, IndexError):
        return None


EMIT_VERSION = "002"   # set by --emit-version. v002 (Oct 2026): same IDs 1-294 as v001 + 18 new library entries
L2A_VERSION = "001"    # EMIT L2A reflectance / mask (shared by both mineral versions)

# v002 library entries 295-312 (emit-sds-l2b data/mineral_id_matrix_v6.00a6.csv) mapped like the v1 table does:
# value = fraction of the mineral in the library sample (from the sample title); jarosites/bloedite are not among the
# 10 EMIT minerals -> 0 (jarosite still counts as group-1 iron absorption for the Fe-oxide 'absent' rule).
V2_EXTRA = {295: {"Hematite": 1.0}, 296: {"Hematite": 0.03}, 297: {"Goethite": 1.0}, 298: {"Goethite": 0.5},
            299: {"Hematite": 0.1, "Goethite": 0.9}, 300: {"Hematite": 0.2, "Goethite": 0.8},
            301: {"Hematite": 0.3, "Goethite": 0.7}, 302: {"Hematite": 0.1, "Goethite": 0.4},
            303: {"Hematite": 0.4, "Goethite": 0.1}}
V2_GROUP = {i: 1 for i in range(295, 312)} | {312: 2}


def search_emit(bbox, start, end, max_cloud, months=None):
    import earthaccess
    res = earthaccess.search_data(short_name="EMITL2BMIN", version=EMIT_VERSION, bounding_box=tuple(bbox),
                                  temporal=(start, end), cloud_cover=(0, max_cloud), count=1000)
    if not res:
        res = [g for g in earthaccess.search_data(short_name="EMITL2BMIN", version=EMIT_VERSION, bounding_box=tuple(bbox),
                                                  temporal=(start, end), count=1000)
               if not (float(g["umm"].get("CloudCover", 0) or 0) > max_cloud)]
    aoi = box(*bbox)
    rows = []
    for g in res:
        if f"_L2B_MIN_{EMIT_VERSION}_" not in g["umm"]["GranuleUR"]:
            continue                                   # never mix product versions
        fp = _granule_footprint(g)
        if fp is None or not fp.is_valid:
            continue
        t = _granule_time(g)
        inter = fp.intersection(aoi)
        rows.append(dict(granule=g, id=g["umm"]["GranuleUR"], time=t, cloud=float(g["umm"].get("CloudCover", np.nan)),
                         overlap=inter.area / aoi.area, aoi_part=inter, in_season=(months is None or t.month in months)))
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["in_season", "cloud", "overlap"], ascending=[False, True, False]).reset_index(drop=True)
    return df


def _links(granule, tag):
    return [u for u in granule.data_links() if f"_{tag}_" in u.split("/")[-1] and u.endswith(".nc")]


def find_l2a(row):
    """EMIT L2A (reflectance + mask) granule of the same acquisition. v001 mineral names carry orbit/scene numbers,
    v002 names only the timestamp, so match on the timestamp."""
    import earthaccess
    ts = re.search(r"_(\d{8}T\d{6})", row["id"]).group(1)
    try:
        rfl = earthaccess.search_data(short_name="EMITL2ARFL", version=L2A_VERSION,
                                      granule_name=f"EMIT_L2A_RFL_{L2A_VERSION}_{ts}*")
    except Exception:
        rfl = []
    if not rfl:
        t = pd.Timestamp(row["time"])
        rfl = earthaccess.search_data(short_name="EMITL2ARFL", version=L2A_VERSION,
                                      temporal=((t - pd.Timedelta(seconds=30)).isoformat(),
                                                (t + pd.Timedelta(seconds=30)).isoformat()), count=20)
    return [g for g in rfl if ts in g["umm"]["GranuleUR"]] or list(rfl)[:1]


def download_emit(row, outdir, with_reflectance=False):
    """EMIT L2B mineral + mineral uncertainty, L2A mask (+ L2A reflectance if requested)."""
    import earthaccess
    os.makedirs(outdir, exist_ok=True)
    g = row["granule"]
    want = {"min": _links(g, "MIN")[:1], "unc": _links(g, "MINUNCERT")[:1]}
    if not want["min"]:
        raise RuntimeError(f"no EMIT_L2B_MIN file in granule {row['id']}")
    rfl = find_l2a(row)
    want["mask"] = [u for gr in rfl for u in _links(gr, "MASK")][:1]
    if with_reflectance:
        want["rfl"] = [u for gr in rfl for u in _links(gr, "RFL")][:1]
    urls = [u for v in want.values() for u in v]
    todo = [u for u in urls if not os.path.exists(os.path.join(outdir, u.split("/")[-1]))]
    if todo:
        log(f"  downloading {len(todo)} EMIT file(s){' (reflectance is ~1-3 GB)' if any('_RFL_' in u for u in todo) else ''} ...")
        earthaccess.download(todo, local_path=outdir)
    return {k: (os.path.join(outdir, v[0].split("/")[-1]) if v else None) for k, v in want.items()}


# =============================================================================
# 2. EMIT reading (GLT orthorectification, cropping, abundances, fit filter, reflectance)
# =============================================================================
def load_grouping(version=None):
    local = os.path.join(HERE, "mineral_grouping_matrix_20230503.csv")
    if not os.path.exists(local):
        import requests
        r = requests.get(GROUPING_URL, timeout=60); r.raise_for_status()
        open(local, "wb").write(r.content)
    g = pd.read_csv(local)
    lut = {}
    for grp in (1, 2):
        m = np.zeros((int(g.Index.max()) + 1, len(EMIT10)), np.float32)
        sub = g[g.Group == grp]
        vals = sub[EMIT10].to_numpy(dtype=float)
        vals[np.isnan(vals)] = 0
        vals[vals == -1] = 1
        m[sub.Index.to_numpy(int)] = vals
        lut[grp] = m
    if (version or EMIT_VERSION) == "002":
        for grp in (1, 2):
            ext = np.zeros((max(V2_GROUP) + 1, len(EMIT10)), np.float32)
            ext[:lut[grp].shape[0]] = lut[grp]
            for idx, minerals in V2_EXTRA.items():
                if V2_GROUP[idx] == grp:
                    for mname, v in minerals.items():
                        ext[idx, EMIT10.index(mname)] = v
            lut[grp] = ext
    return lut


class EmitGeometry:
    """Crop of the EMIT GLT to a site box; reads only the raw (swath) block that the crop needs."""

    def __init__(self, nc_path, bbox):
        import h5py
        with h5py.File(nc_path, "r") as f:
            gt = [float(v) for v in f.attrs["geotransform"]]
            gx = f["location/glt_x"][()].astype(np.int64)
            gy = f["location/glt_y"][()].astype(np.int64)
        c0 = max(int(math.floor((bbox[0] - gt[0]) / gt[1])), 0)
        c1 = min(int(math.ceil((bbox[2] - gt[0]) / gt[1])), gx.shape[1])
        r0 = max(int(math.floor((bbox[3] - gt[3]) / gt[5])), 0)
        r1 = min(int(math.ceil((bbox[1] - gt[3]) / gt[5])), gx.shape[0])
        if c1 <= c0 or r1 <= r0:
            raise RuntimeError("EMIT scene does not overlap the site box")
        gx, gy = gx[r0:r1, c0:c1], gy[r0:r1, c0:c1]
        self.ok = (gx > 0) & (gy > 0)
        if not self.ok.any():
            raise RuntimeError("EMIT scene has no data inside the site box")
        self.y0, self.y1 = int(gy[self.ok].min()) - 1, int(gy[self.ok].max())
        self.x0, self.x1 = int(gx[self.ok].min()) - 1, int(gx[self.ok].max())
        self.iy = gy - 1 - self.y0
        self.ix = gx - 1 - self.x0
        self.grid = Grid(Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3]) * Affine.translation(c0, r0),
                         c1 - c0, r1 - r0)

    def read(self, f, name, bands=None, fill=np.nan):
        d = f[name]
        a = d[self.y0:self.y1, self.x0:self.x1] if bands is None else d[self.y0:self.y1, self.x0:self.x1, bands]
        a = a.astype(np.float32)
        out = np.full(self.ok.shape + a.shape[2:], fill, np.float32)
        out[self.ok] = a[self.iy[self.ok], self.ix[self.ok]]
        return out


def load_emit(files, bbox, lut, fit_min=0.5, cloud_mask=True, info=None):
    """Returns grid, valid mask, abundances (incl. classes) from confident detections, and per-group
    'uncertain' masks (detections whose fit is below fit_min: kept out of the tests)."""
    import h5py
    if f"_L2B_MIN_{EMIT_VERSION}_" not in os.path.basename(files["min"]):
        raise RuntimeError(f"{os.path.basename(files['min'])} is not EMIT version {EMIT_VERSION}: "
                           f"its mineral IDs do not match the grouping table")
    geo = EmitGeometry(files["min"], bbox)
    with h5py.File(files["min"], "r") as f:
        o = {k: geo.read(f, k) for k in ("group_1_band_depth", "group_1_mineral_id",
                                         "group_2_band_depth", "group_2_mineral_id")}
    fit = {}
    if files.get("unc") and os.path.exists(files["unc"]):
        with h5py.File(files["unc"], "r") as f:
            for g in (1, 2):
                if f"group_{g}_fit" in f:
                    fit[g] = geo.read(f, f"group_{g}_fit")
        top = np.nanmax([np.nanmax(np.where(v > -9990, v, np.nan)) if np.isfinite(v).any() else 0 for v in fit.values()]
                        or [1.0])
        if top <= 0.51:                       # v001 stores fit scaled by 1/2 (NASA tutorial); detect, don't assume
            fit = {g: v * 2.0 for g, v in fit.items()}
    valid = geo.ok & (o["group_1_band_depth"] > -9990) & (o["group_2_band_depth"] > -9990)
    if files.get("mask") and os.path.exists(files["mask"]):
        with h5py.File(files["mask"], "r") as f:
            m = geo.read(f, "mask", bands=[0, 1, 2, 3, 4], fill=1)
        cloud = np.nansum(m[..., [0, 1, 4]], axis=-1) > 0          # cloud, cirrus, dilated cloud
        hard = (m[..., 2] > 0) | (m[..., 3] > 0)                  # water, spacecraft
        if info is not None:
            info["cloud_flagged"] = float((cloud & valid).sum() / max(valid.sum(), 1))
        valid &= ~hard
        if cloud_mask:
            valid &= ~cloud
    abund = np.zeros(o["group_1_band_depth"].shape + (len(EMIT10),), np.float32)
    uncertain, bd_conf = {}, {}
    for grp in (1, 2):
        ids = np.nan_to_num(o[f"group_{grp}_mineral_id"], nan=0).astype(np.int64)
        ids[(ids < 0) | (ids >= lut[grp].shape[0])] = 0
        bd = np.clip(np.nan_to_num(o[f"group_{grp}_band_depth"], nan=0), 0, None)
        if grp in fit:
            low = (ids > 0) & ~(np.nan_to_num(fit[grp], nan=0) >= fit_min)
        else:
            low = np.zeros(ids.shape, bool)
        uncertain[grp] = low & valid
        bd = np.where(low, 0, bd)
        bd_conf[grp] = bd
        abund += lut[grp][ids] * bd[..., None]
    abund[~valid] = np.nan
    emit = {name: abund[..., i] for i, name in enumerate(EMIT10)}
    for cname, members in CLASSES.items():
        emit[cname] = np.sum([emit[m] for m in members], axis=0)
    total2 = np.nansum([emit[c] for c in GROUP2_CLASSES], axis=0)
    emit["_total_group2"] = np.where(valid, total2, np.nan)
    for grp in (1, 2):   # confident band depth of ANY mineral of the group (group 1 = all Fe3+/Fe2+ absorptions)
        emit[f"_bd_group{grp}"] = np.where(valid, bd_conf[grp], np.nan)
    if info is not None:   # share of detections that use v002-only library entries (IDs > 294)
        for grp in (1, 2):
            ids = np.nan_to_num(o[f"group_{grp}_mineral_id"], nan=0)[valid]
            det = ids > 0
            info[f"group{grp}_new_id_share"] = float((ids > 294).sum() / max(det.sum(), 1))
    return geo, valid, emit, uncertain, bool(fit)


def emit_simulated_bands(files, geo, sensor):
    """Resample EMIT L2A reflectance to the sensor's bands (flat-top band responses)."""
    import h5py
    with h5py.File(files["rfl"], "r") as f:
        wl = f["sensor_band_parameters/wavelengths"][()].astype(float)
        good = f["sensor_band_parameters/good_wavelengths"][()].astype(bool) \
            if "sensor_band_parameters/good_wavelengths" in f else np.ones(wl.size, bool)
        weights, need = {}, set()
        for b, (c, fw) in BANDPASS[sensor].items():
            w = np.exp(-math.log(2) * np.abs(2 * (wl - c) / fw) ** 8) * good
            w[w < 1e-3] = 0
            weights[b] = w
            need |= set(np.flatnonzero(w))
        idx = sorted(need)
        cube = geo.read(f, "reflectance", bands=idx)
    cube[cube < -9990] = np.nan
    out = {}
    for b, w in weights.items():
        ww = w[idx]
        out[b] = (cube * ww).sum(-1) / ww.sum()
    return out


# =============================================================================
# 3. Sentinel-2 / Landsat
# =============================================================================
def _stac():
    import pystac_client, planetary_computer
    return pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",
                                     modifier=planetary_computer.sign_inplace)


def find_multispectral(sensor, aoi_poly, t0, max_days, max_cloud):
    cat = _stac()
    coll = "sentinel-2-l2a" if sensor == "sentinel-2" else "landsat-c2-l2"
    q = {"eo:cloud_cover": {"lt": max_cloud}}
    if sensor == "landsat":
        q["platform"] = {"in": ["landsat-8", "landsat-9"]}
    rng = f"{(t0 - pd.Timedelta(days=max_days)):%Y-%m-%d}/{(t0 + pd.Timedelta(days=max_days)):%Y-%m-%d}"
    items = list(cat.search(collections=[coll], intersects=aoi_poly.__geo_interface__, datetime=rng, query=q).items())
    by_day = {}
    for it in items:
        by_day.setdefault(it.datetime.date(), []).append(it)
    best = None
    for day, its in by_day.items():
        cover = aoi_poly.intersection(unary_union([shape(i.geometry) for i in its])).area / aoi_poly.area
        if cover < 0.9:
            continue
        gap = abs((pd.Timestamp(day) - t0.normalize()).days)
        cloud = float(np.mean([i.properties.get("eo:cloud_cover", 0) for i in its]))
        if best is None or (gap, cloud) < best[0]:
            best = ((gap, cloud), day, its, cover)
    if best is None:
        return None
    (gap, cloud), day, its, cover = best
    return dict(sensor=sensor, date=pd.Timestamp(day), gap_days=gap, cloud=cloud, cover=cover, items=its,
                ids=[i.id for i in its])


def _resign(href):
    """Fresh Planetary Computer token (tokens expire; a long run can outlive them)."""
    if not href.startswith("http"):
        return href
    try:
        import planetary_computer
        return planetary_computer.sign_url(href.split("?")[0])
    except Exception:
        return href


def _read_to_grid(href, grid, resampling, transform_fn=None, fill=np.nan, tries=4):
    for k in range(tries):
        try:
            return _read_to_grid_once(href, grid, resampling, transform_fn, fill)
        except (rasterio.errors.RasterioIOError, OSError, socket.timeout) as e:
            if k == tries - 1:
                raise
            wait = 15 * (k + 1)
            log(f"    read failed ({type(e).__name__}); new token, retry {k + 2}/{tries} in {wait} s")
            time.sleep(wait)
            href = _resign(href)


def _read_to_grid_once(href, grid, resampling, transform_fn=None, fill=np.nan):
    dst = np.full((grid.height, grid.width), fill, dtype=np.float32)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="3",
                      GDAL_HTTP_TIMEOUT="90", GDAL_HTTP_CONNECTTIMEOUT="30", GDAL_HTTP_LOW_SPEED_TIME="60",
                      GDAL_HTTP_LOW_SPEED_LIMIT="1"), \
            rasterio.open(href) as src:
        b = transform_bounds(grid.crs, src.crs, *grid.bounds, densify_pts=21)
        win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
        a = src.read(1, window=win, boundless=True, fill_value=src.nodata if src.nodata is not None else 0)
        wtr = src.window_transform(win)
        a = transform_fn(a) if transform_fn else a.astype(np.float32)
        reproject(a, dst, src_transform=wtr, src_crs=src.crs, dst_transform=grid.transform, dst_crs=grid.crs,
                  src_nodata=np.nan, dst_nodata=fill, resampling=resampling)
    return dst


def read_multispectral(match, grid, bright_ok=False):
    """bright_ok: for bright evaporite sites (white gypsum is labelled cloud/snow by SCL and CFMask), keep only
    no-data, saturation, shadow and water as 'bad'. Only used on sites flagged emit_cloud_mask=False."""
    sensor = match["sensor"]
    bands = list(BANDPASS[sensor])
    out = {b: np.full((grid.height, grid.width), np.nan, np.float32) for b in bands}
    badfrac = np.full((grid.height, grid.width), np.nan, np.float32)
    for it in match["items"]:
        if sensor == "sentinel-2":
            pb = float(str(it.properties.get("s2:processing_baseline", "0")) or 0)
            off = 1000.0 if pb >= 4.0 else 0.0

            def conv(a, off=off):
                a = a.astype(np.float32); bad = a == 0
                r = (a - off) / 10000.0; r[bad] = np.nan
                return r

            def badm(a):
                a = a.astype(np.float32)
                keep = [4, 5, 7, 8, 9, 10, 11] if bright_ok else [4, 5, 7]
                b = (~np.isin(a, keep)).astype(np.float32); b[a == 0] = np.nan
                return b
            qa_key = "SCL"
        else:
            def conv(a):
                a = a.astype(np.float32); bad = a == 0
                r = a * 0.0000275 - 0.2; r[bad] = np.nan
                return r

            def badm(a):
                a = a.astype(np.int64)
                flags = ((1 << 4) | (1 << 7)) if bright_ok else ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 7))
                b = ((a & flags) > 0).astype(np.float32); b[(a & 1) > 0] = np.nan
                return b
            qa_key = "qa_pixel"
        for b in bands:
            arr = _read_to_grid(it.assets[b].href, grid, Resampling.average, transform_fn=conv)
            fill = np.isnan(out[b]) & np.isfinite(arr)
            out[b][fill] = arr[fill]
        bf = _read_to_grid(it.assets[qa_key].href, grid, Resampling.average, transform_fn=badm)
        fill = np.isnan(badfrac) & np.isfinite(bf)
        badfrac[fill] = bf[fill]
    valid = np.isfinite(badfrac) & (badfrac < 0.01)
    for b in bands:
        valid &= np.isfinite(out[b]) & (out[b] > 0.005)
    return out, valid


# =============================================================================
# 4. Statistics
# =============================================================================
def auc_score(score, label):
    n1 = int(label.sum()); n0 = int(label.size - n1)
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = rankdata(score)
    return float((ranks[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def block_bootstrap_auc(score, label, block_id, reps=200, seed=0):
    rng = np.random.default_rng(seed)
    ub, inv, cnt = np.unique(block_id, return_inverse=True, return_counts=True)
    members = np.split(np.argsort(inv, kind="stable"), np.cumsum(cnt)[:-1])
    vals = []
    for _ in range(reps):
        ix = np.concatenate([members[k] for k in rng.integers(0, len(ub), len(ub))])
        a = auc_score(score[ix], label[ix])
        if np.isfinite(a):
            vals.append(a)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if len(vals) > 20 else (np.nan, np.nan)


def roc_curve(score, label, n=200):
    thr = np.quantile(score, np.linspace(1, 0, n))
    P, N = label.sum(), (~label).sum()
    return (np.r_[0, [(score[~label] >= t).sum() / N for t in thr], 1],
            np.r_[0, [(score[label] >= t).sum() / P for t in thr], 1])


def verdict(auc, lo, hi, n_pos, n_neg):
    if n_pos < 200 or n_neg < 200:
        return "TOO FEW PIXELS"
    if not np.isfinite(auc):
        return "N/A"
    if auc >= 0.70 and (not np.isfinite(lo) or lo >= 0.60):
        return "WORKS"
    if auc <= 0.30 and (not np.isfinite(hi) or hi <= 0.40):
        return "INVERTED"
    if auc >= 0.60:
        return "WEAK"
    if auc <= 0.40:
        return "WEAK INVERTED"
    return "FAILS"


def threshold_for(target, args, T=None):
    """Presence threshold. Group-1 (Fe) band depths are ~5-10x smaller than group-2 ones, so they get their own."""
    if T is not None:
        return T
    return args.min_abund_fe if GROUP_OF[target] == 1 else args.min_abund


def labels_for(emit, uncertain, target, T, purity=None):
    """+1 clear presence (abund >= T, optionally pure), 0 clear absence, -1 excluded.
    Absence: group-2 targets -> target abundance == 0 (other SWIR minerals allowed: realistic contrast);
             group-1 (Fe) targets -> NO group-1 absorption at all (v0.2 used 'hematite+goethite == 0', which
             counted pixels dominated by other ferric minerals as 'absent' and inverted the Fe test)."""
    a = emit[target]
    g = GROUP_OF[target]
    lab = np.full(a.shape, -1, np.int8)
    unc = uncertain.get(g, np.zeros(a.shape, bool))
    pos = (a >= T) & ~unc
    if purity is not None:
        tot = emit["_total_group2"] if g == 2 else emit["_bd_group1"]
        pos &= (a / np.where(tot > 0, tot, np.nan)) >= purity
    lab[pos] = 1
    neg = (emit["_bd_group1"] == 0) if g == 1 else (a == 0)
    lab[neg & ~unc] = 0
    return lab


def labels_quantile(emit, uncertain, target, idx, q=0.2):
    """Threshold-free contrast: top q of the target abundance (and > 0) vs bottom q, among the sampled pixels."""
    rr, cc = idx
    a = np.nan_to_num(emit[target][rr, cc], nan=-1)
    unc = uncertain.get(GROUP_OF[target], np.zeros(emit[target].shape, bool))[rr, cc]
    lab = np.full(a.shape, -1, np.int8)
    ok = (a >= 0) & ~unc
    if GROUP_OF[target] == 1:   # other iron absorbers (ferrihydrite, jarosite, Fe2+ silicates...) are not 'low Fe'
        ok &= ~((a == 0) & (np.nan_to_num(emit["_bd_group1"][rr, cc]) > 0))
    if ok.sum() < 10:
        return lab
    hi, lo = np.quantile(a[ok], 1 - q), np.quantile(a[ok], q)
    pos = ok & (a >= hi) & (a > 0)
    lab[pos] = 1
    lab[ok & (a <= lo) & ~pos] = 0
    return lab


def dominant_group2(emit, T, purity):
    """Index into GROUP2_CLASSES of the pixel's pure dominant class, -1 if none."""
    st = np.nan_to_num(np.stack([emit[c] for c in GROUP2_CLASSES], -1), nan=0)
    tot = st.sum(-1)
    k = st.argmax(-1)
    frac = np.take_along_axis(st, k[..., None], -1)[..., 0] / np.where(tot > 0, tot, np.nan)
    return np.where((st.max(-1) >= T) & (frac >= purity), k, -1)


def analyse(site, scene, grid, emit_valid, emit, uncertain, ms, ms_valid, sensor, outdir, args, seed=0):
    rat_def = RATIOS[base_sensor(sensor)]
    b = BANDPASS[base_sensor(sensor)]
    nir, red = [k for k in b][3], [k for k in b][2]
    ndvi = (ms[nir] - ms[red]) / (ms[nir] + ms[red])
    bare = emit_valid & ms_valid & (ndvi < args.ndvi_max)
    n_bare = int(bare.sum())
    log(f"  {sensor}: {n_bare} bare, clear pixels shared with EMIT")
    if n_bare < 1000:
        raise RuntimeError(f"only {n_bare} usable pixels")
    rr, cc = np.nonzero(bare)
    rng = np.random.default_rng(seed)
    if rr.size > args.max_pixels:
        k = rng.choice(rr.size, args.max_pixels, replace=False); rr, cc = rr[k], cc[k]
    block = (rr // 16) * 100000 + (cc // 16)
    rat = {name: (ms[r["num"]] / ms[r["den"]]) for name, r in rat_def.items()}
    if not sensor.startswith(SIM_PREFIX):
        for name, a in rat.items():
            write_tif(os.path.join(outdir, f"ratio_{slug(sensor)}_{slug(name)}.tif"), np.where(bare, a, np.nan), grid)
    P = args.purity
    dom = dominant_group2(emit, args.min_abund, P)
    main_rows, sens_rows, mat_rows = [], [], []
    for name, r in rat_def.items():
        x = rat[name][rr, cc].astype(np.float64)
        ok = np.isfinite(x)
        # full ratio x mineral matrix (main test definition)
        for cname in ALL_TARGETS:
            lab = labels_for(emit, uncertain, cname, threshold_for(cname, args))[rr, cc]
            sel = ok & (lab >= 0)
            mat_rows.append(dict(site=site, scene=scene, sensor=sensor, ratio=name, mineral=cname,
                                 auc=auc_score(x[sel], lab[sel] == 1), n_pos=int((lab[sel] == 1).sum()),
                                 n_neg=int((lab[sel] == 0).sum())))
        c = r["claim"]
        T = threshold_for(c, args)
        # MAIN
        lab = labels_for(emit, uncertain, c, T)[rr, cc]
        sel = ok & (lab >= 0)
        xs, ys, bs = x[sel], lab[sel] == 1, block[sel]
        a = auc_score(xs, ys)
        lo, hi = block_bootstrap_auc(xs, ys, bs, seed=seed)
        rho = spearmanr(x[ok], np.nan_to_num(emit[c][rr, cc][ok])).correlation
        top = xs >= np.quantile(xs, 0.9) if xs.size else np.zeros(0, bool)
        lift = ys[top].mean() / ys.mean() if ys.size and ys.mean() > 0 else np.nan
        # QUANTILE CONTRAST (threshold-free): top 20 % vs bottom 20 % of EMIT abundance
        labq = labels_quantile(emit, uncertain, c, (rr, cc), args.quantile)
        selq = ok & (labq >= 0)
        xq, yq = x[selq], labq[selq] == 1
        aq = auc_score(xq, yq)
        qlo, qhi = block_bootstrap_auc(xq, yq, block[selq], seed=seed + 1)
        # PURE
        labp = labels_for(emit, uncertain, c, T, purity=P)[rr, cc]
        selp = ok & (labp >= 0)
        a_pure = auc_score(x[selp], labp[selp] == 1)
        # SPECIFICITY (group-2 classes): pure target vs pure other classes
        spec = np.nan; spec_n = 0
        if c in GROUP2_CLASSES:
            d = dom[rr, cc]
            tgt = ok & (d == GROUP2_CLASSES.index(c)); oth = ok & (d >= 0) & (d != GROUP2_CLASSES.index(c))
            if tgt.sum() >= 50 and oth.sum() >= 50:
                spec = auc_score(x[tgt | oth], tgt[tgt | oth]); spec_n = int(oth.sum())
        main_rows.append(dict(site=site, scene=scene, sensor=sensor, ratio=name, claimed_target=c,
                              presence_threshold=T,
                              n_pos=int(ys.sum()), n_neg=int((~ys).sum()), auc=a, auc_ci95_low=lo, auc_ci95_high=hi,
                              spearman_rho=float(rho), top10_lift=float(lift), auc_pure=a_pure,
                              n_pure=int((labp[selp] == 1).sum()), specificity_auc=spec, n_other_pure=spec_n,
                              verdict=verdict(a, lo, hi, int(ys.sum()), int((~ys).sum())),
                              auc_quantile=aq, auc_q_ci95_low=qlo, auc_q_ci95_high=qhi,
                              n_pos_q=int(yq.sum()), n_neg_q=int((~yq).sum()),
                              verdict_quantile=verdict(aq, qlo, qhi, int(yq.sum()), int((~yq).sum()))))
        # SENSITIVITY
        for Ts in (args.sens_thresholds_fe if GROUP_OF[c] == 1 else args.sens_thresholds):
            for pur in (None, P):
                labs = labels_for(emit, uncertain, c, Ts, purity=pur)[rr, cc]
                s2 = ok & (labs >= 0)
                sens_rows.append(dict(site=site, scene=scene, sensor=sensor, ratio=name, claimed_target=c,
                                      threshold=Ts, pure_only=pur is not None, auc=auc_score(x[s2], labs[s2] == 1),
                                      n_pos=int((labs[s2] == 1).sum()), n_neg=int((labs[s2] == 0).sum())))
    return (pd.DataFrame(main_rows), pd.DataFrame(sens_rows), pd.DataFrame(mat_rows), rat, bare, dom, (rr, cc))


# =============================================================================
# 5. Figures and reports
# =============================================================================
def _stretch(a, p=(2, 98)):
    v = a[np.isfinite(a)]
    if v.size == 0:
        return a
    lo, hi = np.percentile(v, p)
    return np.clip((a - lo) / (hi - lo + 1e-12), 0, 1)


def make_figures(site, sensor, grid, emit, ms, rat, bare, dom, sample, res, mat, outdir, T):
    ext = [grid.bounds[0], grid.bounds[2], grid.bounds[1], grid.bounds[3]]
    tag = slug(sensor)
    figs = {}
    bnames = list(BANDPASS[base_sensor(sensor)])
    rgb = np.nan_to_num(np.dstack([_stretch(ms[bnames[2]]), _stretch(ms[bnames[1]]), _stretch(ms[bnames[0]])]))
    fig, ax = plt.subplots(2, 4, figsize=(20, 9.5))
    ax[0, 0].imshow(rgb, extent=ext); ax[0, 0].set_title(f"{sensor} true colour")
    for k, c in enumerate(["Carbonate", "Clay (Al-OH)", "Fe-oxide"]):
        a = np.where(bare, emit[c], np.nan)
        vmax = np.nanpercentile(a[a > 0], 98) if np.any(a > 0) else 1
        im = ax[0, k + 1].imshow(a, extent=ext, cmap="viridis", vmin=0, vmax=vmax)
        ax[0, k + 1].set_title(f"EMIT {c} (spectral abundance)"); fig.colorbar(im, ax=ax[0, k + 1], shrink=0.7)
    ax[1, 0].imshow(np.where(bare, 1, 0), extent=ext, cmap="gray"); ax[1, 0].set_title("Pixels used (white)")
    for k, (name, a) in enumerate(rat.items()):
        im = ax[1, k + 1].imshow(np.where(bare, a, np.nan), extent=ext, cmap="magma",
                                 vmin=np.nanpercentile(a[bare], 2), vmax=np.nanpercentile(a[bare], 98))
        ax[1, k + 1].set_title(f"{sensor}: {name}"); fig.colorbar(im, ax=ax[1, k + 1], shrink=0.7)
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{site}: EMIT reference vs {sensor} ratios", fontsize=14)
    fig.tight_layout(); p = os.path.join(outdir, f"maps_{tag}.png"); fig.savefig(p, dpi=100); plt.close(fig)
    figs["maps"] = p
    rr, cc = sample
    fig, ax = plt.subplots(1, len(rat), figsize=(6 * len(rat), 4.8))
    labels_g2 = ["no pure class"] + GROUP2_CLASSES
    for k, (name, a) in enumerate(rat.items()):
        x = a[rr, cc]
        if "Ferric" in name:
            fe, g1 = emit["Fe-oxide"][rr, cc], emit["_bd_group1"][rr, cc]
            q = np.nanquantile(fe, [0.25, 0.5, 0.75]) if np.isfinite(fe).any() else [0, 0, 0]
            groups = [x[g1 == 0], x[(fe > 0) & (fe < q[0])], x[(fe >= q[0]) & (fe < q[1])],
                      x[(fe >= q[1]) & (fe < q[2])], x[fe >= q[2]]]
            lab = ["no Fe absorption", "Fe-oxide Q1", "Q2", "Q3", "Q4 (most)"]
        else:
            d = dom[rr, cc]
            groups = [x[d == i - 1] for i in range(len(labels_g2))]; lab = labels_g2
        keep = [i for i, gq in enumerate(groups) if np.isfinite(gq).sum() >= 30]
        ax[k].boxplot([groups[i][np.isfinite(groups[i])] for i in keep], showfliers=False)
        ax[k].set_xticks(range(1, len(keep) + 1))
        ax[k].set_xticklabels([f"{lab[i]}\n(n={np.isfinite(groups[i]).sum()})" for i in keep], fontsize=8)
        ax[k].set_title(f"{name} (claimed: high for {RATIOS[base_sensor(sensor)][name]['claim']})"); ax[k].grid(alpha=0.3)
    fig.suptitle(f"{site} – {sensor}: ratio by EMIT class (pure pixels: class ≥ 70 % of the absorption)")
    fig.tight_layout(); p = os.path.join(outdir, f"boxplots_{tag}.png"); fig.savefig(p, dpi=100); plt.close(fig)
    figs["box"] = p
    m = mat.pivot(index="ratio", columns="mineral", values="auc")[ALL_TARGETS]
    npos = mat.pivot(index="ratio", columns="mineral", values="n_pos")[ALL_TARGETS]
    fig, ax = plt.subplots(figsize=(16, 2.6 + 0.5 * len(m)))
    im = ax.imshow(m.values, cmap="RdBu_r", vmin=0.2, vmax=0.8, aspect="auto")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v, n = m.values[i, j], npos.values[i, j]
            ax.text(j, i, "–" if (not np.isfinite(v) or n < 50) else f"{v:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if (np.isfinite(v) and abs(v - 0.5) > 0.2) else "black")
    fig.colorbar(im, ax=ax, label="AUC (0.5 = no skill; <0.5 = inverse)")
    ax.set_title(f"{site} – {sensor}: AUC per EMIT mineral (presence ≥ {T} vs absence; '–' = fewer than 50 pixels)")
    fig.tight_layout(); p = os.path.join(outdir, f"auc_matrix_{tag}.png"); fig.savefig(p, dpi=100); plt.close(fig)
    figs["matrix"] = p
    return figs


VCOL = {"WORKS": "#1a9850", "WEAK": "#fdae61", "FAILS": "#d73027", "WEAK INVERTED": "#9e6fb8",
        "INVERTED": "#7b3294", "TOO FEW PIXELS": "#999999", "N/A": "#999999"}
CSS = """<style>body{font-family:Segoe UI,Arial,sans-serif;max-width:1500px;margin:20px auto;color:#222}
table{border-collapse:collapse;margin:10px 0}td,th{border:1px solid #ccc;padding:4px 8px;font-size:13px}
th{background:#eee}img{max-width:100%;margin:8px 0;border:1px solid #ddd}.note{color:#555;font-size:13px}</style>"""
LEGEND = """<p class="note">AUC = probability that a pixel where EMIT clearly finds the claimed mineral class has a higher ratio
than a pixel where EMIT finds none of it (0.5 = no skill, 1 = perfect, &lt;0.5 = the ratio works the opposite way).
Weak or uncertain EMIT detections (below the abundance threshold, or fit below the minimum) are left out.
Presence threshold: 0.05 for SWIR classes, 0.01 for Fe-oxide (group-1 band depths are much smaller); Fe-oxide
"absent" = no group-1 (iron) absorption at all. Top/bottom 20% = threshold-free contrast: the 20% of pixels with the
most EMIT abundance vs the 20% with the least. 95% CI: spatial block bootstrap (~1 km blocks). Pure = target ≥ 70% of EMIT absorption. Specificity = pure target vs
pure other classes. EMIT-as-… rows = the same ratios computed from EMIT's own reflectance resampled to the sensor bands
(same pixels and time): if they agree with the real sensor, the result reflects the band physics, not data quality.</p>"""


def _fmt(v, f="{:.3f}"):
    return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f.format(v)


def results_table(df):
    h = ["<table><tr>" + "".join(f"<th>{c}</th>" for c in
         ["Sensor", "Ratio", "Claimed", "Present / absent px", "AUC [95% CI]", "AUC pure", "Specificity",
          "Spearman ρ", "Top-10% lift", "Verdict", "AUC top/bottom 20% [95% CI]", "Verdict (top/bottom)"]) + "</tr>"]
    for r in df.itertuples():
        h.append(f"<tr><td>{r.sensor}</td><td>{r.ratio}</td><td>{r.claimed_target}</td><td>{r.n_pos} / {r.n_neg}</td>"
                 f"<td>{_fmt(r.auc)} [{_fmt(r.auc_ci95_low)}–{_fmt(r.auc_ci95_high)}]</td><td>{_fmt(r.auc_pure)}</td>"
                 f"<td>{_fmt(r.specificity_auc)}</td><td>{_fmt(r.spearman_rho, '{:+.3f}')}</td>"
                 f"<td>{_fmt(r.top10_lift, '{:.2f}')}</td>"
                 f"<td style='background:{VCOL.get(r.verdict, '#999')};color:white;font-weight:bold'>{r.verdict}</td>"
                 f"<td>{_fmt(r.auc_quantile)} [{_fmt(r.auc_q_ci95_low)}–{_fmt(r.auc_q_ci95_high)}]</td>"
                 f"<td style='background:{VCOL.get(r.verdict_quantile, '#999')};color:white;font-weight:bold'>"
                 f"{r.verdict_quantile}</td></tr>")
    return "\n".join(h + ["</table>"])


def write_scene_report(site, scene, meta, res, sens, figs, outdir):
    meta_rows = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in meta.items())
    piv = sens.pivot_table(index=["sensor", "ratio"], columns=["pure_only", "threshold"], values="auc")
    img = "".join(f"<h3>{s}</h3>" + "".join(f"<img src='{os.path.basename(p)}'>" for p in fd.values())
                  for s, fd in figs.items())
    html = (f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{site} {scene}</title>{CSS}</head><body>"
            f"<p><a href='../../index.html'>← all sites</a></p><h1>{site} – {scene}</h1>{LEGEND}"
            f"<h2>Result</h2>{results_table(res)}<h2>Sensitivity (AUC by abundance threshold; pure_only=True: pure pixels)</h2>"
            f"{piv.to_html(float_format=lambda v: f'{v:.3f}')}<h2>Scenes</h2><table>{meta_rows}</table>{img}"
            f"<p class='note'>emit_ratio_check.py v{__version__}</p></body></html>")
    p = os.path.join(outdir, "report.html")
    open(p, "w", encoding="utf-8").write(html)
    return p


def write_index(out, allres, allsens, links, args, mat=None):
    df = allres.copy()
    df["case"] = df.site + " " + df.scene
    # dot plot: AUC with CI per case, one panel per (claimed ratio), colours by sensor
    ratios = list(dict.fromkeys(df.claimed_target))
    sensors = list(dict.fromkeys(df.sensor))
    cols = {s: c for s, c in zip(sensors, ["#2c7fb8", "#d95f0e", "#7fcdbb", "#fdae6b"])}
    cases = list(dict.fromkeys(df.case))
    fig, ax = plt.subplots(1, len(ratios), figsize=(6 * len(ratios), 0.45 * len(cases) + 2), sharey=True)
    ax = np.atleast_1d(ax)
    for k, c in enumerate(ratios):
        for j, s in enumerate(sensors):
            d = df[(df.claimed_target == c) & (df.sensor == s)].set_index("case").reindex(cases)
            y = np.arange(len(cases)) + (j - (len(sensors) - 1) / 2) * 0.18
            ax[k].errorbar(d.auc, y, xerr=[d.auc - d.auc_ci95_low, d.auc_ci95_high - d.auc], fmt="o", ms=4,
                           color=cols[s], label=s, capsize=2)
            ax[k].plot(d.auc_quantile, y, "D", ms=5, mfc="none", color=cols[s],
                       label="top/bottom 20%" if j == 0 else None)
        ax[k].axvline(0.5, color="k", lw=1); ax[k].axvspan(0.7, 1, color="#1a9850", alpha=0.08)
        ax[k].axvspan(0, 0.3, color="#7b3294", alpha=0.08); ax[k].set_xlim(0, 1)
        name = df[df.claimed_target == c].ratio.iloc[0]
        ax[k].set_title(f"'{c}' ratio ({name.split()[-1]} on S2)"); ax[k].set_xlabel("AUC vs EMIT")
        ax[k].grid(alpha=0.3)
    ax[0].set_yticks(range(len(cases))); ax[0].set_yticklabels(cases, fontsize=8)
    ax[0].set_ylim(len(cases) - 0.5, -0.5)
    ax[-1].legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(out, "summary_auc.png"), dpi=110); plt.close(fig)
    # summary statistics per sensor x ratio
    g = df.groupby(["sensor", "ratio", "claimed_target"])
    summ = g.agg(cases=("auc", "size"), median_auc=("auc", "median"), min_auc=("auc", "min"), max_auc=("auc", "max"),
                 median_auc_top_bottom=("auc_quantile", "median"), median_rho=("spearman_rho", "median"),
                 median_auc_pure=("auc_pure", "median"), median_specificity=("specificity_auc", "median")).reset_index()
    vc = df.groupby(["sensor", "ratio"]).verdict_quantile.value_counts().unstack(fill_value=0)
    vc.columns = [f"top/bottom: {c}" for c in vc.columns]
    summ = summ.merge(vc.reset_index(), on=["sensor", "ratio"], how="left")
    summ.to_csv(os.path.join(out, "summary_by_ratio.csv"), index=False)
    piv = df.pivot_table(index=["sensor", "ratio"], columns="case", values="auc")
    mat_html = ""
    if mat is not None and len(mat):
        m = mat[~mat.sensor.str.startswith(SIM_PREFIX) & (mat.n_pos >= 50) & (mat.n_neg >= 50)]
        if len(m):
            mm = m.pivot_table(index=["sensor", "ratio"], columns="mineral", values="auc", aggfunc="median")
            mm = mm[[c for c in ALL_TARGETS if c in mm.columns]]
            nn = m.pivot_table(index=["sensor", "ratio"], columns="mineral", values="auc", aggfunc="size")
            nn = nn[mm.columns]
            mm.to_csv(os.path.join(out, "summary_mineral_matrix.csv"))
            mat_html = (f"<h2>Which EMIT mineral does each ratio actually follow? (median AUC across scenes)</h2>"
                        f"<p class='note'>Cells need ≥ 50 present and ≥ 50 absent pixels in a scene; "
                        f"the second table gives the number of scenes behind each median.</p>"
                        f"{mm.to_html(float_format=lambda v: f'{v:.2f}')}"
                        f"<details><summary>scenes per cell</summary>{nn.to_html()}</details>")
    link_html = "".join(f"<li><a href='{os.path.relpath(p, out)}'>{n}</a></li>" for n, p in links)
    html = (f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>EMIT ratio check – all sites</title>{CSS}</head><body>"
            f"<h1>Do Landsat / Sentinel-2 mineral ratios agree with EMIT? – {len(cases)} scene(s), "
            f"{df.site.nunique()} site(s)</h1>{LEGEND}<img src='summary_auc.png'>"
            f"<h2>Summary per ratio</h2>{summ.to_html(index=False, float_format=lambda v: f'{v:.3f}')}"
            f"{mat_html}<h2>AUC per scene</h2>{piv.to_html(float_format=lambda v: f'{v:.3f}')}"
            f"<h2>All results</h2>{results_table(df)}<h2>Scene reports</h2><ul>{link_html}</ul>"
            f"<p class='note'>Settings: presence threshold {args.min_abund} (Fe-oxide {args.min_abund_fe}), "
            f"EMIT mineral product v{EMIT_VERSION}, top/bottom share {args.quantile}, purity {args.purity}, fit ≥ {args.fit_min}, "
            f"NDVI &lt; {args.ndvi_max}, max {args.max_days} days between EMIT and multispectral. "
            f"emit_ratio_check.py v{__version__}</p></body></html>")
    p = os.path.join(out, "index.html")
    open(p, "w", encoding="utf-8").write(html)
    return p


# =============================================================================
# 6. Site workflow
# =============================================================================
def choose_scenes(name, cfg, args):
    cand = search_emit(cfg["bbox"], args.start, args.end, args.emit_max_cloud, cfg.get("months"))
    if not len(cand):
        log(f"  ! no EMIT mineral scene over {name} with ≤{args.emit_max_cloud}% cloud")
        return [], pd.DataFrame()
    forced = args.forced_dates.get(name)
    n_want = args.scenes_per_site
    if forced:                                       # --scene-dates: exactly these EMIT dates, no other filter
        cand = cand[cand.time.dt.strftime("%Y-%m-%d").isin(forced) & (cand.overlap >= args.min_overlap)]
        cand = (cand.assign(_day=cand.time.dt.date).sort_values("overlap", ascending=False)
                .drop_duplicates("_day").drop(columns="_day"))           # best-covering granule of each date
        n_want = len(forced)
        missing = sorted(set(forced) - set(cand.time.dt.strftime("%Y-%m-%d")))
        if missing:
            log(f"  ! no EMIT v{EMIT_VERSION} scene on {', '.join(missing)} over {name}")
        log(f"  using the requested date(s) {', '.join(forced)}; pairing with Sentinel-2/Landsat ...")
    else:
        cand = cand[cand.overlap >= args.min_overlap]
        if args.season_only:
            cand = cand[cand.in_season]
        log(f"  {len(cand)} EMIT scene(s) with ≥{args.min_overlap:.0%} overlap"
            f"{' in the dry season' if args.season_only else ''}; pairing with Sentinel-2/Landsat ...")
    picked, tried = [], []
    for _, row in cand.head(args.max_candidates).iterrows():
        if len(picked) >= n_want and not args.list_only:
            break
        if any(abs((row["time"] - p[0]["time"]).days) < 30 for p in picked):
            continue
        matches = {}
        for s in args.sensors:
            try:
                m = find_multispectral(s, row["aoi_part"], row["time"], args.max_days, args.ms_max_cloud)
            except Exception as e:
                log(f"  ! {s} search failed: {e}"); m = None
            if m:
                matches[s] = m
        tried.append(dict(emit=row["id"], time=str(row["time"]), cloud=row["cloud"], overlap=round(row["overlap"], 2),
                          in_season=row["in_season"],
                          **{f"{s}": (f"{matches[s]['date']:%Y-%m-%d} ({matches[s]['gap_days']} d)" if s in matches else "-")
                             for s in args.sensors}))
        if "sentinel-2" in matches or (matches and "sentinel-2" not in args.sensors):
            picked.append((row, matches))
    return picked, pd.DataFrame(tried)


def run_scene(name, cfg, row, matches, args, lut):
    scene = f"{row['time']:%Y%m%d}"
    outdir = os.path.join(args.out, name, scene)
    os.makedirs(outdir, exist_ok=True)
    done = [os.path.join(outdir, f) for f in ("results.csv", "sensitivity.csv", "auc_matrix.csv", "report.html", "results.json")]
    if args.resume and all(os.path.exists(f) for f in done):
        try:
            meta = json.load(open(done[4]))["meta"]
            rd = lambda f: pd.read_csv(f, dtype={"scene": str})
            res = rd(done[0])
            want_cm = "on" if (cfg.get("emit_cloud_mask", True) and name not in args.no_emit_cloud_mask) else "off"
            if (meta.get("EMIT granule") == row["id"] and "auc_quantile" in res.columns
                    and meta.get("EMIT cloud mask", "on").split()[0] == want_cm):
                log(f"  --- {name} {scene}: already done with {row['id']} - reusing (use --no-resume to redo)")
                return res, rd(done[1]), rd(done[2]), done[3]
        except Exception:
            pass
    want_rfl = args.emit_reflectance == "all" or (args.emit_reflectance == "first" and getattr(args, "_first", True))
    log(f"  --- {name} {scene}: EMIT {row['id']} (cloud {row['cloud']:.0f}%, overlap {row['overlap']:.0%})")
    files = download_emit(row, args.cache, with_reflectance=want_rfl)
    use_cm = cfg.get("emit_cloud_mask", True) and name not in args.no_emit_cloud_mask
    info = {}
    geo, emit_valid, emit, uncertain, has_fit = load_emit(files, cfg["bbox"], lut, args.fit_min, use_cm, info)
    log(f"  EMIT cloud mask {'ON' if use_cm else 'OFF (site setting)'}: {info.get('cloud_flagged', 0):.0%} of pixels "
        f"flagged as cloud" + (f"; v002-only library IDs: group 1 {info.get('group1_new_id_share', 0):.0%}, "
                               f"group 2 {info.get('group2_new_id_share', 0):.0%} of detections"
                               if EMIT_VERSION == "002" else ""))
    grid = geo.grid
    log(f"  EMIT grid {grid.width}×{grid.height} px; valid {emit_valid.mean():.0%}; fit filter "
        f"{'on' if has_fit else 'OFF (no uncertainty file)'}; uncertain group-2 detections "
        f"{uncertain[2][emit_valid].mean():.1%}")
    for c in ["Carbonate", "Clay (Al-OH)", "Fe-oxide", "Gypsum"]:
        write_tif(os.path.join(outdir, f"emit_{slug(c)}.tif"), emit[c], grid)
    meta = {"Site": f"{name} {cfg['bbox']} – {cfg['note']}", "EMIT granule": row["id"],
            "EMIT mineral product": f"EMITL2BMIN v{EMIT_VERSION}",
            "EMIT cloud mask": "on" if use_cm else "off (bright surface flagged as cloud)",
            "EMIT time (UTC)": str(row["time"]), "EMIT cloud %": f"{row['cloud']:.0f}",
            "Site covered by EMIT": f"{row['overlap']:.0%}", "EMIT fit filter": f"fit ≥ {args.fit_min}" if has_fit else "not available",
            "Presence threshold / purity": f"{args.min_abund} (Fe-oxide: {args.min_abund_fe}) / {args.purity}",
            "Quantile contrast": f"top {args.quantile:.0%} vs bottom {args.quantile:.0%} of EMIT abundance",
            "NDVI threshold": args.ndvi_max}
    inputs = []
    for s, m in matches.items():
        ms, ms_valid = read_multispectral(m, grid, bright_ok=not use_cm)
        meta[f"{s} scenes"] = ", ".join(m["ids"])
        meta[f"{s} date / gap"] = f"{m['date']:%Y-%m-%d} / {m['gap_days']} days"
        inputs.append((s, ms, ms_valid))
    if files.get("rfl") and os.path.exists(files["rfl"]):
        for s in list(matches):
            log(f"  simulating {s} bands from EMIT reflectance")
            sim = emit_simulated_bands(files, geo, s)
            sv = np.all([np.isfinite(v) & (v > 0.005) for v in sim.values()], axis=0)
            inputs.append((SIM_PREFIX + s, sim, sv))
        meta["Physics control"] = "EMIT L2A reflectance resampled to sensor bands (flat-top band responses)"
    res_all, sens_all, mat_all, figs = [], [], [], {}
    for s, ms, ms_valid in inputs:
        try:
            res, sens, mat, rat, bare, dom, sample = analyse(name, scene, grid, emit_valid, emit, uncertain,
                                                             ms, ms_valid, s, outdir, args)
        except RuntimeError as e:
            log(f"  ! {s}: {e}"); continue
        figs[s] = make_figures(name, s, grid, emit, ms, rat, bare, dom, sample, res, mat, outdir,
                               f"{args.min_abund} (Fe: {args.min_abund_fe})")
        res_all.append(res); sens_all.append(sens); mat_all.append(mat)
    if not res_all:
        return None
    res, sens, mat = (pd.concat(x, ignore_index=True) for x in (res_all, sens_all, mat_all))
    res.to_csv(os.path.join(outdir, "results.csv"), index=False)
    sens.to_csv(os.path.join(outdir, "sensitivity.csv"), index=False)
    mat.to_csv(os.path.join(outdir, "auc_matrix.csv"), index=False)
    json.dump(dict(meta=meta, results=res.to_dict("records")), open(os.path.join(outdir, "results.json"), "w"),
              indent=2, default=str)
    rep = write_scene_report(name, scene, meta, res, sens, figs, outdir)
    print(res[["sensor", "ratio", "n_pos", "n_neg", "auc", "auc_quantile", "spearman_rho",
               "specificity_auc", "verdict", "verdict_quantile"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    return res, sens, mat, rep


SCREEN_MINERALS = ["Calcite", "Dolomite", "Gypsum", "Kaolinite", "Illite+Muscovite", "Montmorillonite",
                   "Chlorite", "Fe-oxide"]


def screen_sites(jobs, args, lut):
    """Cheap check before a full run: EMIT mineral files only (no reflectance, no Sentinel-2/Landsat).
    Counts confident and pure pixels per mineral, so a candidate site is kept only if EMIT really sees the target."""
    rows = []
    for name, cfg in jobs:
        log(f"=== screening {name} {cfg['bbox']}")
        try:
            cand = search_emit(cfg["bbox"], args.start, args.end, args.emit_max_cloud, cfg.get("months"))
            if len(cand):
                cand = cand[cand.overlap >= args.min_overlap]
                if args.season_only:
                    cand = cand[cand.in_season]
            if not len(cand):
                log(f"  ! no usable EMIT scene (try --any-season / --emit-max-cloud 30)"); continue
            for _, row in cand.head(args.screen_scenes).iterrows():
                files = download_emit(row, args.cache, with_reflectance=False)
                use_cm = cfg.get("emit_cloud_mask", True) and name not in args.no_emit_cloud_mask
                geo, valid, emit, unc, _ = load_emit(files, cfg["bbox"], lut, args.fit_min, use_cm)
                r = dict(site=name, scene=f"{row['time']:%Y-%m-%d}", cloud=round(row["cloud"]), overlap=round(row["overlap"], 2),
                         valid_px=int(valid.sum()))
                for m in SCREEN_MINERALS:
                    T = threshold_for(m, args)
                    a = np.nan_to_num(emit[m])
                    tot = emit["_bd_group1"] if GROUP_OF[m] == 1 else emit["_total_group2"]
                    pure = (a >= T) & (a / np.where(np.nan_to_num(tot) > 0, tot, np.nan) >= args.purity)
                    r[f"{m} present"] = int(((a >= T) & valid).sum())
                    r[f"{m} pure"] = int((pure & valid).sum())
                rows.append(r)
                print("  " + ", ".join(f"{m}: {r[m + ' present']} ({r[m + ' pure']} pure)" for m in SCREEN_MINERALS))
        except Exception as e:
            log(f"  ! {name} failed: {type(e).__name__}: {e}")
    if rows:
        df = pd.DataFrame(rows)
        p = os.path.join(args.out, "screen_sites.csv")
        df.to_csv(p, index=False)
        show = ["site", "scene", "valid_px"] + [f"{m} pure" for m in SCREEN_MINERALS]
        print("\nPure EMIT pixels per mineral (a site is useful for a mineral with >= 200 pure pixels):")
        print(df[show].to_string(index=False))
        log(f"screening table: {p}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", nargs="*", default=[], help="preset site(s); default: all; see --list-sites")
    ap.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"))
    ap.add_argument("--name", default="custom_site")
    ap.add_argument("--list-sites", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "emit_check_results"))
    ap.add_argument("--start", default="2022-08-01"); ap.add_argument("--end", default=f"{dt.date.today()}")
    ap.add_argument("--scenes-per-site", type=int, default=2)
    ap.add_argument("--season-only", action="store_true", default=True)
    ap.add_argument("--any-season", dest="season_only", action="store_false", help="also use scenes outside the dry season")
    ap.add_argument("--emit-max-cloud", type=float, default=15)
    ap.add_argument("--min-overlap", type=float, default=0.3)
    ap.add_argument("--sensors", nargs="+", default=["sentinel-2", "landsat"], choices=["sentinel-2", "landsat"])
    ap.add_argument("--max-days", type=int, default=30)
    ap.add_argument("--ms-max-cloud", type=float, default=10)
    ap.add_argument("--max-candidates", type=int, default=20)
    ap.add_argument("--ndvi-max", type=float, default=0.25)
    ap.add_argument("--min-abund", type=float, default=0.05, help="EMIT spectral abundance for 'clearly present'")
    ap.add_argument("--min-abund-fe", type=float, default=0.01,
                    help="presence threshold for Fe-oxide (group-1 band depths are much smaller)")
    ap.add_argument("--quantile", type=float, default=0.2, help="top/bottom share for the threshold-free contrast")
    ap.add_argument("--purity", type=float, default=0.7, help="share of EMIT absorption for 'pure' pixels")
    ap.add_argument("--fit-min", type=float, default=0.5, help="minimum EMIT fit (0-1) for a detection to be trusted")
    ap.add_argument("--sens-thresholds", nargs="+", type=float, default=[0.02, 0.05, 0.10])
    ap.add_argument("--sens-thresholds-fe", nargs="+", type=float, default=[0.005, 0.01, 0.02])
    ap.add_argument("--emit-reflectance", default="first", choices=["none", "first", "all"],
                    help="download EMIT reflectance (1-3 GB/scene) for the physics control: first scene per site, all, none")
    ap.add_argument("--max-pixels", type=int, default=150_000)
    ap.add_argument("--list-only", action="store_true")
    ap.add_argument("--no-emit-cloud-mask", nargs="*", default=[], metavar="SITE",
                    help="sites where EMIT's cloud flags are ignored (bright evaporites are flagged as cloud); "
                         "white_sands is set by default")
    ap.add_argument("--emit-version", choices=["001", "002"], default="001",
                    help="EMIT L2B mineral product version (002 = current, reprocessed Oct 2026; 001 = original)")
    ap.add_argument("--cache", default=None, help="folder for downloaded EMIT files (default: <out>/_emit_cache); "
                    "point several result folders to the same cache to avoid downloading twice")
    ap.add_argument("--scene-dates", nargs="+", default=[], metavar="SITE=YYYY-MM-DD[,YYYY-MM-DD]",
                    help="use exactly these EMIT dates for a site, e.g. white_sands=2023-04-20 (other sites unchanged)")
    ap.add_argument("--no-resume", dest="resume", action="store_false", default=True,
                    help="recompute scenes that already have results (default: reuse them if the EMIT granule matches)")
    ap.add_argument("--net-timeout", type=float, default=180, help="seconds before a stalled download/search is aborted")
    ap.add_argument("--screen", action="store_true",
                    help="only count EMIT minerals per site (small downloads, no Sentinel-2/Landsat): use before adding a site")
    ap.add_argument("--screen-scenes", type=int, default=2, help="EMIT scenes per site checked by --screen")
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args(argv)

    if args.list_sites:
        for k, v in SITES.items():
            print(f"{k:14s} {v['country']:8s} {str(v['bbox']):38s} months {v['months']}\n{'':14s} {v['note']}")
        return
    unknown = [s for s in args.site if s not in SITES]
    if unknown:
        raise SystemExit(f"unknown site(s) {unknown}; use --list-sites or --bbox")
    jobs = [(s, SITES[s]) for s in args.site]
    if args.bbox:
        jobs.append((args.name, dict(bbox=list(args.bbox), months=None, country="", note="custom box")))
    if not jobs:
        jobs = list(SITES.items())
        log(f"no site given: running all {len(jobs)} preset sites")
    global EMIT_VERSION
    EMIT_VERSION = args.emit_version
    args.cache = args.cache or os.path.join(args.out, "_emit_cache")
    args.forced_dates = {}
    for item in args.scene_dates:
        if "=" not in item:
            ap.error(f"--scene-dates expects SITE=YYYY-MM-DD, got {item!r}")
        site, dates = item.split("=", 1)
        args.forced_dates[site.strip()] = [d.strip() for d in dates.split(",") if d.strip()]
    os.makedirs(args.out, exist_ok=True)
    socket.setdefaulttimeout(args.net_timeout)        # no network call may hang forever
    earthdata_login()
    lut = load_grouping(EMIT_VERSION)
    if args.screen:
        screen_sites(jobs, args, lut)
        log("done"); return
    allres, allsens, allmat, links = [], [], [], []
    for name, cfg in jobs:
        log(f"=== {name} ({cfg.get('country', '')}) {cfg['bbox']}")
        try:
            picked, tried = choose_scenes(name, cfg, args)
            os.makedirs(os.path.join(args.out, name), exist_ok=True)
            if len(tried):
                tried.to_csv(os.path.join(args.out, name, "candidates.csv"), index=False)
                print(tried.to_string(index=False))
            if args.list_only:
                continue
            if not picked:
                log(f"  ! {name}: no usable EMIT + Sentinel-2 pair (try --any-season, --max-days 60, --emit-max-cloud 30)")
                continue
            for k, (row, matches) in enumerate(picked):
                args._first = (k == 0)
                try:
                    out = run_scene(name, cfg, row, matches, args, lut)
                except Exception as e:
                    log(f"  ! {name} {row['time']:%Y%m%d} failed: {type(e).__name__}: {e}")
                    out = None
                if out:
                    allres.append(out[0]); allsens.append(out[1]); allmat.append(out[2]); links.append((f"{name} {row['time']:%Y-%m-%d}", out[3]))
        except Exception as e:
            log(f"  ! {name} failed: {type(e).__name__}: {e}")
    if allres:
        R, S = pd.concat(allres, ignore_index=True), pd.concat(allsens, ignore_index=True)
        R.to_csv(os.path.join(args.out, "summary_all_results.csv"), index=False)
        S.to_csv(os.path.join(args.out, "summary_all_sensitivity.csv"), index=False)
        M = pd.concat(allmat, ignore_index=True)
        M.to_csv(os.path.join(args.out, "summary_all_matrix.csv"), index=False)
        idx = write_index(args.out, R, S, links, args, M)
        print("\nMedian AUC across scenes:")
        print(R.groupby(["sensor", "ratio"]).agg(n=("auc", "size"), median_auc=("auc", "median"),
                                                  median_auc_q=("auc_quantile", "median"),
                                                  median_rho=("spearman_rho", "median"))
              .to_string(float_format=lambda v: f"{v:.3f}"))
        log(f"index: {idx}")
        if not args.no_open:
            webbrowser.open(Path(idx).resolve().as_uri())
    log("done")


if __name__ == "__main__":
    main()
