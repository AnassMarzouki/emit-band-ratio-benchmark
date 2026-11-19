# Main processing script

Add `emit_ratio_check.py` (version 0.3.8, the version used for the paper) to this folder.

It reads the NASA Earthdata login from `earthdata_login.txt` next to the script, searches and downloads the EMIT
scenes, pairs them with Sentinel-2 and Landsat images from Microsoft Planetary Computer, and writes the result
tables and an HTML report per run.

Settings used in the paper (also printed in every report): presence threshold 0.05 (Fe-oxide 0.01), top/bottom
share 0.2, purity 0.7, fit ≥ 0.5, NDVI < 0.25, at most 30 days between EMIT and the satellite images, EMIT
mineral product v001 (v002 for the robustness test).

[Add the exact command lines used for the 16 scenes and for the v002 run.]
