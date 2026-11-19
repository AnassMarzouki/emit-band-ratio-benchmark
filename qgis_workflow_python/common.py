"""Shared helpers for the step scripts. Nothing to edit here."""
import os
import sys
import json
import time

import numpy as np

import config as C


def _use_bundled_proj_data():
    """Other software (here PostGIS) can set PROJ_LIB / PROJ_DATA to an old proj.db, which breaks
    rasterio. Point PROJ and GDAL to the data files shipped inside the rasterio package instead.
    Only affects this Python process."""
    import importlib.util
    spec = importlib.util.find_spec("rasterio")
    if not spec or not spec.submodule_search_locations:
        return
    pkg = list(spec.submodule_search_locations)[0]
    proj = os.path.join(pkg, "proj_data")
    gdal = os.path.join(pkg, "gdal_data")
    if os.path.exists(os.path.join(proj, "proj.db")):
        os.environ["PROJ_DATA"] = proj
        os.environ["PROJ_LIB"] = proj
    if os.path.isdir(gdal):
        os.environ["GDAL_DATA"] = gdal
    try:                                   # tell an already-loaded PROJ as well
        from rasterio.env import set_proj_data_search_path
        set_proj_data_search_path(proj)
    except Exception:
        try:
            from rasterio._env import set_proj_data_search_path
            set_proj_data_search_path(proj)
        except Exception:
            pass


_use_bundled_proj_data()
COMMON_VERSION = "v6 (PROJ fix)"

for _s in (sys.stdout, sys.stderr):          # Windows consoles: never crash on a character
    try:
        _s.reconfigure(errors="replace")
    except Exception:
        pass


def log(msg):
    print(time.strftime("%H:%M:%S"), str(msg).encode("ascii", "replace").decode(), flush=True)


def ensure_dirs():
    for d in (C.DIR_EMIT, C.DIR_S2, C.DIR_LS, C.DIR_WORK, C.DIR_OUT):
        os.makedirs(d, exist_ok=True)


def fail(msg):
    log("ERROR: " + msg)
    sys.exit(1)


# ------------------------------------------------------------------ credentials
CREDENTIAL_FILES = ("earthdata_login.txt", "earthdata_login.txt.txt", "credentials.txt", "login.txt")


def _parse_credentials(text):
    """Accepts either two plain lines (user, password) or 'user = ...' / 'password = ...' lines."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    kv = {}
    for ln in lines:
        cut = [i for i in (ln.find("="), ln.find(":")) if i > 0]
        if cut:
            i = min(cut)
            kv[ln[:i].strip().lower()] = ln[i + 1:].strip().strip('"').strip("'")
    user = next((kv[k] for k in ("user", "username", "login", "user name", "earthdata_user") if k in kv), None)
    pw = next((kv[k] for k in ("password", "pass", "pwd", "earthdata_password") if k in kv), None)
    if user and pw:
        return user, pw
    if len(lines) >= 2:
        return lines[0], lines[1]
    return None


def earthdata_credentials():
    """Reads the NASA Earthdata login from a text file next to the scripts. Never asks, never prints it."""
    for name in CREDENTIAL_FILES:
        p = os.path.join(C.HERE, name)
        if os.path.exists(p):
            cred = _parse_credentials(open(p, encoding="utf-8-sig").read())
            if cred:
                log(f"Earthdata login read from {name}")
                return cred
            fail(f"{name} was found but could not be read. Put the user name on line 1 and the "
                 f"password on line 2 (or 'user = ...' and 'password = ...').")
    fail(f"No login file in {C.HERE}. Create earthdata_login.txt there with your NASA Earthdata "
         f"user name on line 1 and password on line 2.")


# ------------------------------------------------------------------ rasters on the EMIT grid
def grid_path():
    return os.path.join(C.DIR_WORK, "grid.json")


def save_grid(transform, width, height):
    json.dump({"transform": list(transform)[:6], "width": int(width), "height": int(height),
               "crs": "EPSG:4326"}, open(grid_path(), "w"), indent=1)


def load_grid():
    from affine import Affine
    if not os.path.exists(grid_path()):
        fail("EMIT grid not found. Run step2_emit_abundance.py first.")
    g = json.load(open(grid_path()))
    return Affine(*g["transform"]), g["width"], g["height"]


def write_tif(name, arr, folder=None):
    import rasterio
    tr, w, h = load_grid()
    path = os.path.join(folder or C.DIR_WORK, name if name.endswith(".tif") else name + ".tif")
    a = np.asarray(arr, np.float32)
    with rasterio.open(path, "w", driver="GTiff", width=w, height=h, count=1, dtype="float32",
                       crs="EPSG:4326", transform=tr, nodata=np.nan, compress="deflate") as ds:
        ds.write(a, 1)
    return path


def read_tif(name, folder=None):
    import rasterio
    path = os.path.join(folder or C.DIR_WORK, name if name.endswith(".tif") else name + ".tif")
    if not os.path.exists(path):
        fail(f"{os.path.basename(path)} is missing. Run the previous step first.")
    with rasterio.open(path) as ds:
        a = ds.read(1).astype(np.float32)
        if ds.nodata is not None and not np.isnan(ds.nodata):
            a[a == ds.nodata] = np.nan
    return a


# ------------------------------------------------------------------ statistics
def auc_quantile(abund, ratio, q=None):
    """Top q of abundance (and > 0) vs bottom q, Mann-Whitney AUC with average ranks for ties.
    Same rule as the spreadsheet. Returns (auc, n_high, n_low)."""
    import pandas as pd
    q = C.QUANTILE if q is None else q
    a = np.asarray(abund, float)
    r = np.asarray(ratio, float)
    ok = np.isfinite(a) & np.isfinite(r)
    a, r = a[ok], r[ok]
    if a.size < 10:
        return np.nan, 0, 0
    hi, lo = np.quantile(a, 1 - q), np.quantile(a, q)
    pos = (a >= hi) & (a > 0)
    neg = (a <= lo) & ~pos
    n1, n0 = int(pos.sum()), int(neg.sum())
    if n1 == 0 or n0 == 0:
        return np.nan, n1, n0
    ranks = pd.Series(np.concatenate([r[pos], r[neg]])).rank(method="average").to_numpy()
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), n1, n0


def verdict(auc):
    if not np.isfinite(auc):
        return "NOT ENOUGH DATA"
    if auc >= 0.7:
        return "WORKS"
    if auc <= 0.3:
        return "INVERTED"
    if auc >= 0.6:
        return "WEAK"
    if auc <= 0.4:
        return "WEAK INVERTED"
    return "NO SKILL"
