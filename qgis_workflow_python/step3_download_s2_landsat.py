"""
STEP 3 - Download Sentinel-2 L2A and Landsat 8/9 Collection 2 Level-2 for the box (manual Step 3).

Run:  python step3_download_s2_landsat.py
Source: Microsoft Planetary Computer (free, no account needed). Only the part of each
image that covers the box is downloaded, at full resolution and in its own UTM projection.
Picks the acquisition day closest to the EMIT day (within MAX_DAYS, cloud cover <= 20%).
Writes clipped GeoTIFFs and an items.json (scene names, dates, processing baseline)
into runs/<site>_<date>/s2/ and .../landsat/.
"""
import os
import json
import time
from datetime import datetime, timedelta

import numpy as np

import config as C
from common import log, fail, ensure_dirs

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
SENSORS = {
    "s2": dict(collection="sentinel-2-l2a", assets=["B02", "B04", "B8A", "B11", "B12", "SCL"], folder=C.DIR_S2),
    "landsat": dict(collection="landsat-c2-l2", assets=["blue", "red", "nir08", "swir16", "swir22", "qa_pixel"],
                    folder=C.DIR_LS),
}
MARGIN_M = 500


def emit_points(n=3000):
    """Lon/lat of usable EMIT pixels (from step 2), or a regular grid over the box if step 2 has not run."""
    import os
    from common import load_grid, read_tif
    if os.path.exists(os.path.join(C.DIR_WORK, "grid.json")) and os.path.exists(os.path.join(C.DIR_WORK, "emit_valid.tif")):
        tr, w, h = load_grid()
        rows, cols = np.nonzero(read_tif("emit_valid") == 1)
        if rows.size:
            pick = np.random.default_rng(0).choice(rows.size, size=min(n, rows.size), replace=False)
            return np.c_[tr.c + (cols[pick] + 0.5) * tr.a, tr.f + (rows[pick] + 0.5) * tr.e]
    w, s, e, n_ = C.BBOX
    gx, gy = np.meshgrid(np.linspace(w, e, 60), np.linspace(s, n_, 60))
    return np.c_[gx.ravel(), gy.ravel()]


def coverage(items, pts):
    """Share of the points that fall inside at least one scene footprint."""
    from matplotlib.path import Path
    inside = np.zeros(len(pts), bool)
    for it in items:
        g = it.geometry or {}
        polys = [g["coordinates"]] if g.get("type") == "Polygon" else g.get("coordinates", [])
        for poly in polys:
            inside |= Path(np.asarray(poly[0])[:, :2]).contains_points(pts)
    return float(inside.mean())


def search(key):
    import pystac_client
    import planetary_computer
    cat = pystac_client.Client.open(STAC, modifier=planetary_computer.sign_inplace)
    d0 = datetime.strptime(C.DATE, "%Y-%m-%d").date()
    forced = getattr(C, "S2_DATE" if key == "s2" else "LS_DATE", "") or ""
    if forced:
        f = datetime.strptime(forced, "%Y-%m-%d").date()
        rng = f"{f:%Y-%m-%d}/{f:%Y-%m-%d}"
    else:
        rng = f"{(d0 - timedelta(days=C.MAX_DAYS)):%Y-%m-%d}/{(d0 + timedelta(days=C.MAX_DAYS)):%Y-%m-%d}"
    query = {"eo:cloud_cover": {"lt": 20}}
    if key == "landsat":
        query["platform"] = {"in": ["landsat-8", "landsat-9"]}
    items = list(cat.search(collections=[SENSORS[key]["collection"]], bbox=C.BBOX, datetime=rng,
                            query=query).items())
    if not items:
        fail(f"no {key} image {'on ' + forced if forced else f'within {C.MAX_DAYS} days of {C.DATE}'} with < 20% cloud")
    groups = {}                               # one group = one day (and one platform for Landsat)
    for it in items:
        k = (it.datetime.date(), it.properties.get("platform") if key == "landsat" else "")
        groups.setdefault(k, []).append(it)
    pts = emit_points()
    cand = []
    for (day, plat), its in groups.items():
        cov = coverage(its, pts)
        cand.append((day, plat, its, cov))
        log(f"  {key} candidate {day} {plat}: {len(its)} scene(s), covers {cov:.0%} of the usable EMIT area")
    good = [c for c in cand if c[3] >= C.MIN_COVERAGE]
    if good:
        day, plat, chosen, cov = min(good, key=lambda c: (abs((c[0] - d0).days), -c[3], c[0]))
    else:
        day, plat, chosen, cov = max(cand, key=lambda c: (c[3], -abs((c[0] - d0).days)))
        log(f"WARNING: no {key} day covers {C.MIN_COVERAGE:.0%} of the area; best is {cov:.0%}. "
            f"Results will use only that part of the box.")
    log(f"{key}: using {day} ({(day - d0).days:+d} days from EMIT, covers {cov:.0%}): "
        + ", ".join(it.id for it in chosen))
    return chosen


def clip(href, out):
    import rasterio
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds, Window
    for attempt in range(4):
        try:
            with rasterio.Env(GDAL_HTTP_TIMEOUT=120, GDAL_HTTP_MAX_RETRY=4, GDAL_HTTP_RETRY_DELAY=5):
                with rasterio.open(href) as ds:
                    b = transform_bounds("EPSG:4326", ds.crs, *C.BBOX, densify_pts=21)
                    b = (b[0] - MARGIN_M, b[1] - MARGIN_M, b[2] + MARGIN_M, b[3] + MARGIN_M)
                    win = from_bounds(*b, transform=ds.transform).round_offsets().round_lengths()
                    win = win.intersection(Window(0, 0, ds.width, ds.height))
                    a = ds.read(1, window=win)
                    prof = dict(driver="GTiff", width=a.shape[1], height=a.shape[0], count=1, dtype=a.dtype,
                                crs=ds.crs, transform=ds.window_transform(win), nodata=ds.nodata,
                                compress="deflate")
            with rasterio.open(out + ".part.tif", "w", **prof) as o:
                o.write(a, 1)
            os.replace(out + ".part.tif", out)
            return True
        except Exception as e:                # outside this tile, or a network hiccup
            if "intersect" in str(e).lower():
                return False
            log(f"    retry {attempt + 1}/3 after: {type(e).__name__}: {str(e)[:120]}")
            time.sleep(10 * (attempt + 1))
    fail(f"could not read {href[:80]}...")


def main():
    ensure_dirs()
    for key, s in SENSORS.items():
        meta_path = os.path.join(s["folder"], "items.json")
        if os.path.exists(meta_path):
            log(f"{key}: already downloaded ({meta_path}) - delete the folder to download again")
            continue
        items = search(key)
        meta = []
        for it in items:
            files = {}
            for a in s["assets"]:
                if a not in it.assets:
                    fail(f"{it.id} has no asset '{a}'")
                out = os.path.join(s["folder"], f"{it.id}_{a}.tif")
                if os.path.exists(out) or clip(it.assets[a].href, out):
                    files[a] = os.path.basename(out)
                    log(f"  {os.path.basename(out)}")
            if len(files) == len(s["assets"]):
                meta.append({"id": it.id, "datetime": it.datetime.isoformat(),
                             "platform": it.properties.get("platform"),
                             "processing_baseline": it.properties.get("s2:processing_baseline"),
                             "files": files})
        if not meta:
            fail(f"{key}: none of the scenes overlaps the box")
        json.dump(meta, open(meta_path, "w"), indent=1)
    log("STEP 3 DONE. Next: python step4_align_and_ratios.py")


if __name__ == "__main__":
    main()
