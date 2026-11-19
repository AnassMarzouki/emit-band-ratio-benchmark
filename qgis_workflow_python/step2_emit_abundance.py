"""
STEP 2 - EMIT abundance maps (manual Steps 2 "split the bands" and 4).

Run:  python step2_emit_abundance.py
Reads the EMIT "_min" ENVI file(s) from step 1, keeps the scene that covers the box best,
and writes to runs/<site>_<date>/work/:
  group_1_band_depth, group_1_mineral_id, group_2_band_depth, group_2_mineral_id
  emit_valid                    1 = usable EMIT pixel
  abund_carbonate, abund_clay, abund_fe, abund_calcite, abund_dolomite
  fe_other                      1 = iron absorption that is not hematite/goethite
This grid (about 60 m, EPSG:4326) is the reference grid for every later step.
"""
import os
import re
import glob

import numpy as np
import pandas as pd

import config as C
from common import log, fail, ensure_dirs, save_grid, write_tif

BANDS = ["group_1_band_depth", "group_1_mineral_id", "group_2_band_depth", "group_2_mineral_id"]


def open_min(path_img):
    """Returns dict band_name -> array, plus transform, for the part of the file inside BBOX."""
    import rasterio
    from rasterio.windows import from_bounds
    with rasterio.open(path_img) as ds:
        if ds.count < 4:
            fail(f"{os.path.basename(path_img)} has {ds.count} bands, expected 4")
        desc = [(d or "").lower() for d in ds.descriptions]
        order = {}
        for b in BANDS:
            hit = [i for i, d in enumerate(desc) if b in d and "uncert" not in d]
            if hit:
                order[b] = hit[0]
        if len(order) < 4:                        # no usable names: the order you saw in QGIS
            order = {b: i for i, b in enumerate(BANDS)}
            log(f"  band names not found in the header, assuming order {BANDS}")
        degrees_like = abs(ds.transform.a) < 0.01 and -180 <= ds.bounds.left <= 180 and -90 <= ds.bounds.top <= 90
        if (ds.crs is not None and not ds.crs.is_geographic) or (ds.crs is None and not degrees_like):
            fail(f"{os.path.basename(path_img)} is not on a latitude/longitude grid "
                 f"(was 'Orthorectify' ticked in AppEEARS?)")
        win = from_bounds(*C.BBOX, transform=ds.transform).round_offsets().round_lengths()
        win = win.intersection(rasterio.windows.Window(0, 0, ds.width, ds.height))
        tr = ds.window_transform(win)
        out = {}
        for b, i in order.items():
            a = ds.read(i + 1, window=win).astype(np.float32)
            nod = ds.nodatavals[i]
            if nod is not None and np.isfinite(nod):
                a[a == nod] = np.nan
            a[a < -9990] = np.nan
            out[b] = a
    check_bands(out, os.path.basename(path_img))
    return out, tr


def check_bands(d, name):
    """Band depths must lie in 0-1 and mineral IDs must be whole numbers 0-294; stop if the order is wrong."""
    msgs = []
    for b, a in d.items():
        v = a[np.isfinite(a)]
        if v.size == 0:
            continue
        lo, hi = float(v.min()), float(v.max())
        if "depth" in b and (hi > 1.0 or lo < -0.01):
            msgs.append(f"{b} ranges {lo:.2f}..{hi:.2f}, expected 0..1")
        if "mineral_id" in b and (hi > 400 or np.any(np.abs(v - np.round(v)) > 1e-3)):
            msgs.append(f"{b} ranges {lo:.2f}..{hi:.2f} or is not whole numbers")
    if msgs:
        fail(f"{name}: the bands are not in the expected order - " + "; ".join(msgs))


def lut(csv, column):
    t = pd.read_csv(os.path.join(C.HERE, csv))
    m = np.full(int(t["id"].max()) + 1, np.nan, np.float32)
    m[t["id"].to_numpy(int)] = t[column].to_numpy(float)
    return m


def apply_lut(ids, m):
    out = np.full(ids.shape, np.nan, np.float32)
    ok = np.isfinite(ids) & (ids >= 0) & (ids < m.size)
    out[ok] = m[ids[ok].astype(int)]
    return out


def main():
    import common
    log(f"common.py {getattr(common, 'COMMON_VERSION', 'OLD - replace common.py')}, "
        f"PROJ data: {os.environ.get('PROJ_DATA', 'not set')}")
    ensure_dirs()
    files = sorted(glob.glob(os.path.join(C.DIR_EMIT, "*_min.img")))
    files = [f for f in files if "MINUNCERT" not in os.path.basename(f).upper()]
    if C.EMIT_SCENE:
        files = [f for f in files if C.EMIT_SCENE in os.path.basename(f)]
    if not files:
        fail(f"no EMIT '_min.img' file in {C.DIR_EMIT}. Run step 1 (or copy your AppEEARS files there).")
    best, best_n = None, -1
    for f in files:
        d, tr = open_min(f)
        n = int(np.sum(np.isfinite(d["group_2_band_depth"]) & (d["group_2_band_depth"] >= 0)))
        log(f"{os.path.basename(f)}: {n} valid pixels in the box")
        if n > best_n:
            best, best_n, data, transform = f, n, d, tr
    ts = re.search(r"\d{8}T\d{6}", os.path.basename(best))
    log(f"using {os.path.basename(best)} (scene {ts.group(0) if ts else '?'})")
    if best_n < 1000:
        fail("fewer than 1000 valid EMIT pixels in the box: wrong scene or box?")

    h, w = data["group_2_band_depth"].shape
    save_grid(transform, w, h)
    for b in BANDS:
        write_tif(b, data[b])

    g1, g2 = data["group_1_band_depth"], data["group_2_band_depth"]
    valid = np.isfinite(g1) & np.isfinite(g2) & (g1 >= 0) & (g2 >= 0)
    write_tif("emit_valid", valid.astype(np.float32))
    log(f"EMIT grid {w} x {h} pixels, {valid.mean():.1%} usable")

    fr = {
        "carbonate": apply_lut(data["group_2_mineral_id"], lut("emit_group2_reclass.csv", "carbonate")),
        "clay": apply_lut(data["group_2_mineral_id"], lut("emit_group2_reclass.csv", "clay_aloh")),
        "calcite": apply_lut(data["group_2_mineral_id"], lut("emit_group2_reclass.csv", "calcite")),
        "dolomite": apply_lut(data["group_2_mineral_id"], lut("emit_group2_reclass.csv", "dolomite")),
        "fe": apply_lut(data["group_1_mineral_id"], lut("emit_group1_reclass.csv", "fe_oxide")),
    }
    for name, f in fr.items():
        depth = g1 if name == "fe" else g2
        ab = np.where(valid, f * np.clip(depth, 0, None), np.nan)
        write_tif(f"abund_{name}", ab)
        log(f"abund_{name}: max {np.nanmax(ab):.3f}, {np.nanmean(ab > 0):.1%} of usable pixels > 0")
    ab_fe = np.where(valid, fr["fe"] * np.clip(g1, 0, None), np.nan)
    fe_other = np.where(valid, ((ab_fe == 0) & (g1 > 0)).astype(np.float32), np.nan)
    write_tif("fe_other", fe_other)
    log("STEP 2 DONE. Next: python step3_download_s2_landsat.py")


if __name__ == "__main__":
    main()
