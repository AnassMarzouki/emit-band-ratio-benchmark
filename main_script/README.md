# Main processing script: `emit_ratio_check.py` (v0.3.8)

This is the script that produced the paper's results (Tables 3–6, Figs. 5, 7 and 8). For each site it:

1. searches the EMIT L2B mineral scenes (NASA Earthdata) and pairs each with the closest Sentinel-2 L2A and
   Landsat 8/9 C2 L2 acquisition (Microsoft Planetary Computer) within 30 days, with < 10 % cloud and ≥ 90 % cover;
2. downloads the EMIT minerals, their fit uncertainty, the L2A mask and, for the first scene of each site, the L2A
   reflectance; orthorectifies them with the GLT and computes the abundance of the EMIT mineral classes;
3. averages Sentinel-2 / Landsat onto the EMIT 60 m grid, masks cloud, shadow, water and vegetation (NDVI < 0.25);
4. computes the AUC of each ratio against EMIT, with 95 % spatial block-bootstrap confidence intervals
   (16 × 16-pixel blocks, about 1 km; 200 draws), the per-mineral AUC matrix and the physics control
   (bands simulated from EMIT reflectance);
5. writes `results.csv`, `report.html` per scene and `index.html` + `summary_*.csv` for all sites.

## Setup

```
python -m venv .venv
.venv\Scripts\activate          # Windows  (Linux/macOS: source .venv/bin/activate)
pip install -r requirements.txt
```

Put `earthdata_login.txt` next to the script: NASA Earthdata user name on line 1, password on line 2
(https://urs.earthdata.nasa.gov). It is listed in `.gitignore`: never commit it.
NASA's mineral grouping table (`mineral_grouping_matrix_20230503.csv`) is downloaded automatically on the first run.

## Reproducing the paper

Main run (EMIT v001, the 16 scenes of Table 2):

```
python emit_ratio_check.py --site hamar_laghdad figuig jbel_rhassoul tiwit bou_azzer cuprite mountain_pass rosh_pinah ^
  --scene-dates hamar_laghdad=2023-04-25,2023-10-09 figuig=2023-06-21,2024-06-09 jbel_rhassoul=2022-08-27,2023-04-01 ^
  tiwit=2023-04-01,2023-06-24 bou_azzer=2023-04-01,2023-10-01 cuprite=2024-07-27,2026-06-16 ^
  mountain_pass=2024-04-16,2025-06-23 rosh_pinah=2023-01-07,2025-01-14 ^
  --emit-version 001 --out results_v001
```

(`^` continues a command on the next line in Windows CMD; use `\` on Linux/macOS.)

Robustness run with EMIT v002 (Table 6, left): the same command with `--emit-version 002 --emit-reflectance none
--out results_v002`, restricted to the sites and dates for which v002 scenes exist.

Useful options: `--list-sites`, `--list-only` (show candidate scene pairs without downloading), `--screen`
(count EMIT mineral pixels before adding a site), `--emit-reflectance none` (skip the physics control and its
1–3 GB downloads), `--cache` (share downloaded EMIT files between result folders).

Settings used in the paper (the defaults; they are printed in every `index.html`): presence threshold 0.05
(Fe-oxide 0.01), top/bottom share 0.2, purity 0.7, fit ≥ 0.5, NDVI < 0.25, at most 30 days between EMIT and the
satellite images, at most 150 000 usable pixels per scene (random subsample, fixed seed).

The preset list also contains candidate sites that were screened but not used in the paper; they are only
processed when named with `--site` or when the script is run without `--site`.
