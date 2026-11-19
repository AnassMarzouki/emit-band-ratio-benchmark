"""
STEP 4 - Masks, averaging onto the EMIT grid, reflectance, ratios and NDVI (manual Steps 5 and 6).

Run:  python step4_align_and_ratios.py
For each sensor:
  1. bad-pixel mask at full resolution
     Sentinel-2: SCL not in {4 vegetation, 5 bare soil, 7 unclassified}
     Landsat:    QA_PIXEL bits 1-5 or 7 set (dilated cloud, cirrus, cloud, shadow, snow, water)
  2. every band and the mask averaged onto the EMIT grid ("Average" resampling)
  3. reflectance: Sentinel-2 (DN - 1000) / 10000 for processing baseline >= 4.0, else DN / 10000
                  Landsat DN * 0.0000275 - 0.2
  4. ratios: fe = red/blue, clay = SWIR1/SWIR2, carb = SWIR2/SWIR1; NDVI
  5. usable pixel: EMIT valid AND < 1% bad AND NDVI < 0.25 AND all bands > 0
Writes to runs/<site>_<date>/work/: s2_*.tif and ls_*.tif.
"""
import os
import json

import numpy as np

import config as C
from common import log, fail, ensure_dirs, load_grid, write_tif, read_tif

BANDS = {
    "s2": {"blue": "B02", "red": "B04", "nir": "B8A", "swir1": "B11", "swir2": "B12", "qa": "SCL"},
    "ls": {"blue": "blue", "red": "red", "nir": "nir08", "swir1": "swir16", "swir2": "swir22", "qa": "qa_pixel"},
}
FOLDERS = {"s2": C.DIR_S2, "ls": C.DIR_LS}


def to_grid(arr, src_tr, src_crs):
    from rasterio.warp import reproject, Resampling
    tr, w, h = load_grid()
    dst = np.full((h, w), np.nan, np.float32)
    reproject(source=arr.astype(np.float32), destination=dst, src_transform=src_tr, src_crs=src_crs,
              dst_transform=tr, dst_crs="EPSG:4326", resampling=Resampling.average,
              src_nodata=np.nan, dst_nodata=np.nan)
    return dst


def read_clip(path):
    import rasterio
    with rasterio.open(path) as ds:
        return ds.read(1), ds.transform, ds.crs, ds.nodata


def bad_mask(key, qa):
    qa = qa.astype(np.int64)
    if key == "s2":
        bad = (~np.isin(qa, [4, 5, 7])).astype(np.float32)
        bad[qa == 0] = np.nan                                  # outside the image
    else:
        bad = ((qa & 0b10111110) > 0).astype(np.float32)       # bits 1,2,3,4,5,7 = 190
        bad[(qa & 1) > 0] = np.nan                             # bit 0 = fill
    return bad


def process(key):
    meta_path = os.path.join(FOLDERS[key], "items.json")
    if not os.path.exists(meta_path):
        fail(f"{meta_path} missing. Run step3_download_s2_landsat.py first.")
    items = json.load(open(meta_path))
    tr, w, h = load_grid()
    out = {name: np.full((h, w), np.nan, np.float32) for name in list(BANDS[key]) + ["bad"]}
    for it in items:
        offset = 0.0
        if key == "s2":
            try:
                offset = 1000.0 if float(it.get("processing_baseline") or 0) >= 4.0 else 0.0
            except ValueError:
                offset = 1000.0
        for name, asset in BANDS[key].items():
            a, s_tr, s_crs, nod = read_clip(os.path.join(FOLDERS[key], it["files"][asset]))
            if name == "qa":
                g = to_grid(bad_mask(key, a), s_tr, s_crs)
                tgt = "bad"
            else:
                a = a.astype(np.float32)
                a[a == (nod if nod is not None else 0)] = np.nan
                a[a == 0] = np.nan
                g = to_grid(a, s_tr, s_crs)
                g = (g - offset) / 10000.0 if key == "s2" else g * 0.0000275 - 0.2
                tgt = name
            fill = np.isnan(out[tgt]) & np.isfinite(g)       # several scenes: first one wins
            out[tgt][fill] = g[fill]
        log(f"{key}: {it['id']} averaged onto the EMIT grid")

    b, r, n, s1, s2 = out["blue"], out["red"], out["nir"], out["swir1"], out["swir2"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = {"fe": r / b, "clay": s1 / s2, "carb": s2 / s1, "ndvi": (n - r) / (n + r)}
    clear = np.isfinite(out["bad"]) & (out["bad"] < C.BAD_FRACTION_MAX)
    positive = np.all([np.isfinite(x) & (x > 0) for x in (b, r, n, s1, s2)], axis=0)
    emit_valid = read_tif("emit_valid") == 1
    ok = emit_valid & clear & positive & (ratios["ndvi"] < C.NDVI_MAX)

    for name in ("blue", "red", "nir", "swir1", "swir2"):
        write_tif(f"{key}_refl_{name}", out[name])
    write_tif(f"{key}_bad_al", out["bad"])
    write_tif(f"{key}_clear", clear.astype(np.float32))
    for name, a in ratios.items():
        write_tif(f"{key}_{name}", a)
    write_tif(f"{key}_ok", ok.astype(np.float32))
    log(f"{key}: bad fraction max {np.nanmax(out['bad']):.3f}; clear {clear.mean():.1%}; "
        f"usable (EMIT valid, clear, NDVI<{C.NDVI_MAX}) {ok.mean():.1%} of the box")
    log(f"{key}: median ratios on usable pixels - fe {np.nanmedian(ratios['fe'][ok]):.3f}, "
        f"clay {np.nanmedian(ratios['clay'][ok]):.3f}, carb {np.nanmedian(ratios['carb'][ok]):.3f}")


def main():
    ensure_dirs()
    for key in ("s2", "ls"):
        process(key)
    log("STEP 4 DONE. Next: python step5_sample_and_auc.py")


if __name__ == "__main__":
    main()
