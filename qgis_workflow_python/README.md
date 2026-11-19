# EMIT band-ratio test, automated (the manual QGIS tutorial as six scripts)

Each script does one part of the manual tutorial and tells you which script to run next.
Everything is written to `runs\<site>_<date>\`. You end with two AUC workbooks
(`auc_S2_...xlsx`, `auc_LS_...xlsx`) in the same format as the one you filled in by hand.

| Script | Manual step | What it does | Time |
| --- | --- | --- | --- |
| `step1_download_emit.py` | 2 | AppEEARS request for EMIT minerals (ENVI, orthorectified), waits, downloads | 10-60 min (waiting) |
| `step2_emit_abundance.py` | 2 + 4 | Picks the best scene, splits the bands, abundance maps, EMIT valid mask | < 1 min |
| `step3_download_s2_landsat.py` | 3 | Sentinel-2 L2A and Landsat C2 L2 for the box (Planetary Computer, no account) | 2-5 min |
| `step4_align_and_ratios.py` | 5 + 6 | Cloud masks, averaging onto the EMIT grid, reflectance, ratios, NDVI, usable pixels | < 1 min |
| `step5_sample_and_auc.py` | 7 | 8000 random pixels, exports, AUC workbooks + Python cross-check | < 1 min |
| `step6_maps.py` (optional) | 9 | Four-panel comparison maps (PNG) | < 1 min |

## One-time setup (Windows, no installer needed)

Open a Command Prompt in this folder (for example `C:\emit_auto`) and run:

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

The virtual environment (`.venv`) keeps these packages separate from QGIS and from your
other Python work, so the NumPy conflict you hit in QGIS cannot happen here.
Next time, only run `.venv\Scripts\activate` before the scripts.

Then create `earthdata_login.txt` in this folder (next to the scripts): line 1 your NASA
Earthdata user name, line 2 your password. Step 1 reads it automatically and never asks or prints it.
`user = ...` / `password = ...` lines also work, and the names `credentials.txt` or `login.txt` are
accepted too. Keep this file private (not in OneDrive or a shared folder).
No Copernicus or USGS account is needed.

## Optional shortcuts for Hamar Laghdad

- **Reuse your manual EMIT download.** Copy the `..._006_..._min.img` and `.hdr` files (and the
  `_005` ones if you like) into `runs\hamar_laghdad_20230425\emit\`. Step 1 then does nothing.
- **Reuse your AppEEARS settings.** Copy your `hamar-emit-20230425-request.json` next to the scripts
  (any file name containing `request` and ending in `.json` is picked up). Step 1 then submits exactly
  the same kind of request, with only the box and date changed.

## Run

```
python step1_download_emit.py
python step2_emit_abundance.py
python step3_download_s2_landsat.py
python step4_align_and_ratios.py
python step5_sample_and_auc.py
python step6_maps.py
```

Each step checks that the previous one finished and stops with a clear message if not.
Step 1 can be interrupted while it waits: run it again and it resumes the same AppEEARS request.

## What to send back

From `runs\hamar_laghdad_20230425\results\`: `auc_S2_hamar_laghdad_20230425.xlsx`,
`auc_LS_hamar_laghdad_20230425.xlsx` and `auc_results.csv` (and the two `maps_*.png` if you made them).

## Another site or date

Edit `config.py`: `SITE`, `BBOX` (west, south, east, north), `DATE`. Set `REFERENCE_AUC = {}`
unless you have script values for that scene. Nothing else changes.

## Differences from the manual run (on purpose, all small)

- Sentinel-2 blue and red are averaged from the 10 m bands (the manual run used the 20 m copies).
- Landsat gets its own QA_PIXEL cloud mask (the manual run reused the Sentinel-2 mask).
- `fe_other` is used: pixels with an iron absorption that is not hematite or goethite are left out
  of the iron test (`EXCLUDE_OTHER_IRON = True` in `config.py`).
- Points are random EMIT pixels drawn with a fixed seed, so every run uses the same points.
- No EMIT cloud mask and no fit filter, exactly as in the manual run.
