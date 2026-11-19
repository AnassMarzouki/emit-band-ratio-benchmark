"""
Settings for the automated version of the manual QGIS test.
Edit this file only. Every step script reads it.
"""
import os

# ---------------------------------------------------------------- site and date
# Pick ONE run by its name, then run steps 1-6. Each run writes to its own folder runs/<site>_<date>/.
RUN = "tiwit_1"

# Every site and EMIT date of the paper (v001 run). "ref" = the script's AUC (top/bottom 20%) for that scene.
# Boxes: west, south, east, north (degrees). Hamar Laghdad uses the tutorial box (east edge -3.96).
RUNS = {
    "hamar_laghdad_1": dict(site="hamar_laghdad", bbox=(-4.25, 31.25, -3.96, 31.52), date="2023-04-25", max_days=10,  # Morocco, calcite (Kess-Kess mud mounds)
                     ref={("Landsat", "Carbonate"): 0.05, ("Landsat", "Clay"): 0.28, ("Landsat", "Fe-oxide"): 0.96, ("Sentinel-2", "Carbonate"): 0.06, ("Sentinel-2", "Clay"): 0.29, ("Sentinel-2", "Fe-oxide"): 0.96}),
    "hamar_laghdad_2": dict(site="hamar_laghdad", bbox=(-4.25, 31.25, -3.96, 31.52), date="2023-10-09", max_days=10,  # Morocco, calcite (Kess-Kess mud mounds)
                     ref={("Landsat", "Carbonate"): 0.08, ("Landsat", "Clay"): 0.14, ("Landsat", "Fe-oxide"): 0.96, ("Sentinel-2", "Carbonate"): 0.09, ("Sentinel-2", "Clay"): 0.15, ("Sentinel-2", "Fe-oxide"): 0.96}),
    "figuig_1": dict(site="figuig", bbox=(-1.57, 32.17, -1.13, 32.52), date="2023-06-21", max_days=10,  # Morocco, calcite + iron oxides
                     ref={("Landsat", "Carbonate"): 0.26, ("Landsat", "Clay"): 0.34, ("Landsat", "Fe-oxide"): 0.97, ("Sentinel-2", "Carbonate"): 0.28, ("Sentinel-2", "Clay"): 0.36, ("Sentinel-2", "Fe-oxide"): 0.98}),
    "figuig_2": dict(site="figuig", bbox=(-1.57, 32.17, -1.13, 32.52), date="2024-06-09", max_days=10,  # Morocco, calcite + iron oxides
                     ref={("Landsat", "Carbonate"): 0.33, ("Landsat", "Clay"): 0.46, ("Landsat", "Fe-oxide"): 0.98, ("Sentinel-2", "Carbonate"): 0.46, ("Sentinel-2", "Clay"): 0.54, ("Sentinel-2", "Fe-oxide"): 0.98}),
    "tiwit_1": dict(site="tiwit", bbox=(-6.05, 30.98, -5.62, 31.34), date="2023-04-01", max_days=10,  # Morocco, Al-OH clays + iron (Tiouit, Jbel Saghro)
                     ref={("Landsat", "Carbonate"): 0.4, ("Landsat", "Clay"): 0.57, ("Landsat", "Fe-oxide"): 0.75, ("Sentinel-2", "Carbonate"): 0.4, ("Sentinel-2", "Clay"): 0.55, ("Sentinel-2", "Fe-oxide"): 0.85}),
    "tiwit_2": dict(site="tiwit", bbox=(-6.05, 30.98, -5.62, 31.34), date="2023-06-24", max_days=10,  # Morocco, Al-OH clays + iron (Tiouit, Jbel Saghro)
                     ref={("Landsat", "Carbonate"): 0.4, ("Landsat", "Clay"): 0.44, ("Landsat", "Fe-oxide"): 0.78, ("Sentinel-2", "Carbonate"): 0.41, ("Sentinel-2", "Clay"): 0.44, ("Sentinel-2", "Fe-oxide"): 0.92}),
    "cuprite_1": dict(site="cuprite", bbox=(-117.32, 37.45, -117.12, 37.62), date="2024-07-27", max_days=10,  # USA, Al-OH clays (reference site)
                     ref={("Landsat", "Carbonate"): 0.5, ("Landsat", "Clay"): 0.78, ("Landsat", "Fe-oxide"): 0.57, ("Sentinel-2", "Carbonate"): 0.51, ("Sentinel-2", "Clay"): 0.78, ("Sentinel-2", "Fe-oxide"): 0.66}),
    "cuprite_2": dict(site="cuprite", bbox=(-117.32, 37.45, -117.12, 37.62), date="2026-06-16", max_days=10,  # USA, Al-OH clays (reference site)
                     ref={("Landsat", "Carbonate"): 0.57, ("Landsat", "Clay"): 0.82, ("Landsat", "Fe-oxide"): 0.73, ("Sentinel-2", "Carbonate"): 0.57, ("Sentinel-2", "Clay"): 0.81, ("Sentinel-2", "Fe-oxide"): 0.82}),
    "mountain_pass_1": dict(site="mountain_pass", bbox=(-115.66, 35.38, -115.42, 35.58), date="2024-04-16", max_days=10,  # USA, dolomite vs limestone
                     ref={("Landsat", "Carbonate"): 0.56, ("Landsat", "Clay"): 0.73, ("Landsat", "Fe-oxide"): 0.96, ("Sentinel-2", "Carbonate"): 0.56, ("Sentinel-2", "Clay"): 0.7, ("Sentinel-2", "Fe-oxide"): 0.98}),
    "mountain_pass_2": dict(site="mountain_pass", bbox=(-115.66, 35.38, -115.42, 35.58), date="2025-06-23", max_days=10,  # USA, dolomite vs limestone
                     ref={("Landsat", "Carbonate"): 0.44, ("Landsat", "Clay"): 0.59, ("Landsat", "Fe-oxide"): 0.95, ("Sentinel-2", "Carbonate"): 0.37, ("Sentinel-2", "Clay"): 0.48, ("Sentinel-2", "Fe-oxide"): 0.95}),
    "rosh_pinah_1": dict(site="rosh_pinah", bbox=(16.55, -28.15, 17.0, -27.75), date="2023-01-07", max_days=10,  # Namibia, dolostone
                     ref={("Landsat", "Carbonate"): 0.46, ("Landsat", "Clay"): 0.42, ("Landsat", "Fe-oxide"): 0.83, ("Sentinel-2", "Carbonate"): 0.39, ("Sentinel-2", "Clay"): 0.38, ("Sentinel-2", "Fe-oxide"): 0.87}),
    "rosh_pinah_2": dict(site="rosh_pinah", bbox=(16.55, -28.15, 17.0, -27.75), date="2025-01-14", max_days=10,  # Namibia, dolostone
                     ref={("Landsat", "Carbonate"): 0.3, ("Landsat", "Clay"): 0.44, ("Landsat", "Fe-oxide"): 0.96, ("Sentinel-2", "Carbonate"): 0.29, ("Sentinel-2", "Clay"): 0.42, ("Sentinel-2", "Fe-oxide"): 0.97}),
    # New Moroccan sites, dates from the screen (15 Nov 2026). "ref" stays empty until the script has run them.
    "bou_azzer_1": dict(site="bou_azzer", bbox=(-6.85, 30.30, -6.40, 30.65), date="2023-04-01", max_days=10,
                        ref={}),  # Morocco, chlorite + illite + kaolinite (EMIT overlap 100%)
    "bou_azzer_2": dict(site="bou_azzer", bbox=(-6.85, 30.30, -6.40, 30.65), date="2023-10-01", max_days=10,
                        ref={}),  # Morocco, chlorite + illite (EMIT overlap 77%)
    "jbel_rhassoul_1": dict(site="jbel_rhassoul", bbox=(-4.55, 32.80, -4.10, 33.10), date="2023-04-01", max_days=10,
                            ref={}),  # Morocco, calcite + dolomite (EMIT overlap 90%)
    "jbel_rhassoul_2": dict(site="jbel_rhassoul", bbox=(-4.55, 32.80, -4.10, 33.10), date="2022-08-27", max_days=10,
                            ref={}),  # Morocco, calcite + dolomite (EMIT overlap 84%)
}

if RUN not in RUNS:
    raise SystemExit(f"config.py: RUN = '{RUN}' is not in RUNS. Choose one of: {', '.join(RUNS)}")
_r = RUNS[RUN]
if _r["date"].startswith("YYYY"):
    raise SystemExit(f"config.py: run '{RUN}' has no date yet. Screen the site first and fill in its EMIT date.")
SITE = _r["site"]
BBOX = tuple(_r["bbox"])
DATE = _r["date"]
MAX_DAYS = _r["max_days"]
# Force a Sentinel-2 or Landsat acquisition day ("" = automatic: closest day covering >= MIN_COVERAGE)
S2_DATE = _r.get("s2_date", "")
LS_DATE = _r.get("ls_date", "")
MIN_COVERAGE = 0.90     # share of the usable EMIT area a Sentinel-2/Landsat day must cover

# Optional: keep only the EMIT scene whose name contains this timestamp
# (e.g. "20230425T112640"). Leave "" to pick the scene that covers the box best.
EMIT_SCENE = ""

# ---------------------------------------------------------------- folders
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "runs", f"{SITE}_{DATE.replace('-', '')}")
DIR_EMIT = os.path.join(ROOT, "emit")        # AppEEARS downloads
DIR_S2 = os.path.join(ROOT, "s2")            # Sentinel-2 clips (native resolution)
DIR_LS = os.path.join(ROOT, "landsat")       # Landsat clips (native resolution)
DIR_WORK = os.path.join(ROOT, "work")        # rasters on the EMIT grid
DIR_OUT = os.path.join(ROOT, "results")      # samples, AUC workbooks, maps

# ---------------------------------------------------------------- credentials
# Text file with the NASA Earthdata user name on line 1 and password on line 2.
# Keep it private: do not put it in a shared or synced folder.
EARTHDATA_LOGIN_FILE = os.path.join(HERE, "earthdata_login.txt")

# Optional: the "...-request.json" that AppEEARS gave you for a request that worked.
# If this file exists, step 1 reuses its exact settings and only changes name, dates and box.
APPEEARS_TEMPLATE = os.path.join(HERE, "appeears_request_template.json")

# ---------------------------------------------------------------- method (same as the tutorial)
NDVI_MAX = 0.25          # bare ground only
BAD_FRACTION_MAX = 0.01  # an EMIT pixel is usable if < 1% of it is cloud/shadow/water
N_POINTS = 8000          # random pixels drawn in the box (one per EMIT pixel, max 9999)
SEED = 42                # same seed = same points every run
QUANTILE = 0.2           # top 20% vs bottom 20% of EMIT abundance
EXCLUDE_OTHER_IRON = True  # drop pixels with an iron absorption that is not hematite/goethite

# AppEEARS polling
POLL_SECONDS = 60
MAX_WAIT_HOURS = 4

# Script values for the selected run (from RUNS), shown next to your results in the summary.
REFERENCE_AUC = _r["ref"]
