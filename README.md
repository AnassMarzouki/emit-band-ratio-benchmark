# Do multispectral mineral band ratios map what they claim?

Code for the paper *"Do multispectral mineral band ratios map what they claim? A benchmark of Sentinel-2 and
Landsat 8/9 against EMIT"* (Marzouki, submitted to *Earth Science Informatics*).

The code tests three standard band ratios against NASA EMIT Level-2B mineral maps, pixel by pixel:

| Ratio | Sentinel-2 | Landsat 8/9 | Claimed target |
| --- | --- | --- | --- |
| Ferric iron | B4/B2 | B4/B2 | hematite, goethite |
| "Clay" | B11/B12 | B6/B7 | Al-OH clays |
| "Carbonate" | B12/B11 | B7/B6 | calcite, dolomite |

Agreement is measured by the AUC between the top and bottom 20 % of EMIT target abundance, over 16 EMIT scenes at
eight dryland sites (Morocco, USA, Namibia).

## Repository content

| Folder | Content |
| --- | --- |
| `main_script/` | `emit_ratio_check.py` (v0.3.8), the script that produced the paper's results: scene search and download, masking, AUC with block-bootstrap confidence intervals, per-mineral analysis, physics control, HTML reports. |
| `qgis_workflow_python/` | The manual QGIS workflow automated as six Python scripts (`step1` … `step6`), with its own README. Used as the independent re-implementation in Section 4.6 of the paper. |
| `figures/` | `plot_library_spectra.py`: Fig. 3, the USGS Spectral Library v7 spectra with the Sentinel-2 and Landsat bands. |
| `data/` | `site_boxes.geojson`: the study boxes, ready to load in QGIS. |

## Requirements

- Python 3.10 or newer
- A free NASA Earthdata account (https://urs.earthdata.nasa.gov) to download EMIT data.
  Sentinel-2 and Landsat are read from Microsoft Planetary Computer and need no account.

```
python -m venv .venv
.venv\Scripts\activate          # Windows  (Linux/macOS: source .venv/bin/activate)
pip install -r main_script/requirements.txt              # main script
pip install -r qgis_workflow_python/requirements.txt     # automated QGIS workflow
```

## Earthdata login

Create a text file named `earthdata_login.txt` next to the scripts: your Earthdata user name on line 1 and your
password on line 2. The scripts read it automatically and never print it.
**This file is listed in `.gitignore`. Never commit it, and do not keep it in a shared or synced folder.**

## Reproducing the paper

1. **Main results** (Tables 3–6, Figs. 5, 7, 8): see `main_script/README.md` for the exact command.
2. **Independent re-implementation** (Table 6, Fig. S1): in `qgis_workflow_python/`, set `RUN` in `config.py` to
   one of the listed site/date runs, then run `step1` to `step6` in order (see `qgis_workflow_python/README.md`).
3. **Fig. 3**: download `ASCIIdata_splib07a.zip` from the USGS Spectral Library Version 7
   (https://doi.org/10.5066/F7RR1WDJ), put it next to `figures/plot_library_spectra.py` and run the script.

Note: the Hamar Laghdad run in `qgis_workflow_python/config.py` uses the tutorial box (eastern edge 3.96° W);
the paper's main run uses 3.90° W (Table 2).

## Data sources

- EMIT L2B mineral identification v001 (https://doi.org/10.5067/EMIT/EMITL2BMIN.001) and v002
  (https://doi.org/10.5067/EMIT/EMITL2BMIN.002); EMIT L2A surface reflectance v001
  (https://doi.org/10.5067/EMIT/EMITL2ARFL.001), from NASA LP DAAC.
- Sentinel-2 Level-2A and Landsat 8/9 Collection 2 Level-2, through Microsoft Planetary Computer.
- `emit_group1_reclass.csv` / `emit_group2_reclass.csv` are derived from NASA's EMIT mineral grouping table
  (`mineral_grouping_matrix_20230503.csv`).

## Citation

If you use this code, please cite the paper (reference to be added on publication).

## License

MIT License (see `LICENSE`).
