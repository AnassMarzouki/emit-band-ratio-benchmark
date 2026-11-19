"""
STEP 5 - Random pixels, sampling and AUC workbooks (manual Step 7).

Run:  python step5_sample_and_auc.py
Draws N_POINTS random EMIT pixels (one point per pixel, fixed SEED), reads every layer at
those pixels and writes to runs/<site>_<date>/results/:
  samples_all.csv                 every drawn pixel with all layers
  s2_all.xlsx, ls_all.xlsx        usable pixels per sensor (same as the QGIS exports)
  auc_S2_<site>_<date>.xlsx       the spreadsheet workbooks, filled (Excel computes on opening)
  auc_LS_<site>_<date>.xlsx
  auc_results.csv                 the same AUCs computed in Python, as a cross-check
"""
import os
import shutil

import numpy as np
import pandas as pd

import config as C
from common import log, fail, ensure_dirs, load_grid, read_tif, auc_quantile, verdict

LAYERS = {
    "emit_carb": "abund_carbonate", "emit_clay": "abund_clay", "emit_fe": "abund_fe",
    "emit_calcite": "abund_calcite", "emit_dolomite": "abund_dolomite", "fe_other": "fe_other",
    "emit_valid": "emit_valid",
    "s2_carb": "s2_carb", "s2_clay": "s2_clay", "s2_fe": "s2_fe", "s2_ndvi": "s2_ndvi", "s2_ok": "s2_ok",
    "ls_carb": "ls_carb", "ls_clay": "ls_clay", "ls_fe": "ls_fe", "ls_ndvi": "ls_ndvi", "ls_ok": "ls_ok",
}
TESTS = [  # sheet, EMIT column, ratio suffix, label
    ("Fe-oxide", "emit_fe", "fe", {"s2": "Iron B4/B2", "ls": "Iron B4/B2"}),
    ("Clay", "emit_clay", "clay", {"s2": "Clay B11/B12", "ls": "Clay B6/B7"}),
    ("Carbonate", "emit_carb", "carb", {"s2": "Carbonate B12/B11", "ls": "Carbonate B7/B6"}),
]
SENSOR_NAME = {"s2": "Sentinel-2", "ls": "Landsat"}
MAX_ROWS = 10000   # the workbook template holds rows 2..10001


def sample():
    tr, w, h = load_grid()
    rng = np.random.default_rng(C.SEED)
    n = min(C.N_POINTS, w * h, MAX_ROWS - 1)
    flat = np.sort(rng.choice(w * h, size=n, replace=False))
    rows, cols = np.divmod(flat, w)
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    df = pd.DataFrame({"id": np.arange(1, n + 1), "row": rows, "col": cols, "lon": lon, "lat": lat})
    for col, layer in LAYERS.items():
        df[col] = read_tif(layer)[rows, cols]
    return df


def selection(df, key, test_suffix):
    m = (df[f"{key}_ok"] == 1) & df[f"{key}_{test_suffix}"].notna()
    if test_suffix == "fe" and C.EXCLUDE_OTHER_IRON:
        m &= df["fe_other"] == 0
    return df[m]


def fill_workbook(df, key, out_path, results):
    from openpyxl import load_workbook
    from openpyxl.styles import Font
    template = os.path.join(C.HERE, "auc_calculator.xlsx")
    if not os.path.exists(template):
        fail("auc_calculator.xlsx is missing next to the scripts")
    shutil.copy(template, out_path)
    wb = load_workbook(out_path)
    F, Bd = Font(name="Arial"), Font(name="Arial", bold=True)
    for sheet, emit_col, suf, _ in TESTS:
        ws = wb[sheet]
        ws["A2"].comment = None
        for r in range(2, MAX_ROWS + 2):
            ws.cell(r, 1).value = None
            ws.cell(r, 2).value = None
        sel = selection(df, key, suf)
        for i, (a, b) in enumerate(zip(sel[emit_col], sel[f"{key}_{suf}"]), start=2):
            ws.cell(i, 1, float(a))
            ws.cell(i, 2, float(b))
        ws["H14"] = "Source"
        ws["H14"].font = Bd
        ws["I14"] = f"A = {emit_col}, B = {key}_{suf} from step5_sample_and_auc.py ({len(sel)} pixels)"
        ws["I14"].font = F
    sm = wb.create_sheet("Summary", 0)
    head = ["Ratio", "AUC (Excel)", "Verdict (Excel)", "AUC (Python check)", "Script reference", "Pixels high / low"]
    for j, t in enumerate(head, 1):
        sm.cell(1, j, t).font = Bd
    for i, (sheet, _, suf, label) in enumerate(TESTS, start=2):
        res = results[(key, sheet)]
        ref = C.REFERENCE_AUC.get((SENSOR_NAME[key], sheet))
        vals = [label[key], f"='{sheet}'!I9", f"='{sheet}'!I10", round(res[0], 4), ref,
                f"='{sheet}'!I6&\" / \"&'{sheet}'!I7"]
        for j, v in enumerate(vals, 1):
            c = sm.cell(i, j, v)
            c.font = F
            if j in (2, 4):
                c.number_format = "0.000"
    sm["A6"] = (f"{SENSOR_NAME[key]}, {C.SITE}, {C.DATE}. Top {int(C.QUANTILE * 100)}% vs bottom "
                f"{int(C.QUANTILE * 100)}% of EMIT abundance. Excel and Python columns should agree.")
    sm["A6"].font = F
    for col, wd in zip("ABCDEF", [22, 14, 18, 20, 18, 18]):
        sm.column_dimensions[col].width = wd
    wb.calculation.fullCalcOnLoad = True
    wb.save(out_path)


def main():
    ensure_dirs()
    df = sample()
    df.to_csv(os.path.join(C.DIR_OUT, "samples_all.csv"), index=False)
    log(f"{len(df)} random pixels drawn; usable: Sentinel-2 {int((df.s2_ok == 1).sum())}, "
        f"Landsat {int((df.ls_ok == 1).sum())}")
    keep = ["id", "lon", "lat", "emit_carb", "emit_clay", "emit_fe", "fe_other"]
    for key in ("s2", "ls"):
        m = df[f"{key}_ok"] == 1
        df.loc[m, keep + [f"{key}_carb", f"{key}_clay", f"{key}_fe"]].to_excel(
            os.path.join(C.DIR_OUT, f"{key}_all.xlsx"), index=False)

    results, rows = {}, []
    for key in ("s2", "ls"):
        for sheet, emit_col, suf, label in TESTS:
            sel = selection(df, key, suf)
            auc, n1, n0 = auc_quantile(sel[emit_col], sel[f"{key}_{suf}"])
            results[(key, sheet)] = (auc, n1, n0)
            ref = C.REFERENCE_AUC.get((SENSOR_NAME[key], sheet))
            rows.append({"sensor": SENSOR_NAME[key], "ratio": label[key], "auc": round(auc, 4),
                         "verdict": verdict(auc), "n_high": n1, "n_low": n0, "n_pixels": len(sel),
                         "script_reference": ref})
            log(f"{SENSOR_NAME[key]:10s} {label[key]:18s} AUC {auc:.3f}  {verdict(auc):13s} "
                f"(high {n1}, low {n0})" + (f"   script: {ref:.2f}" if ref is not None else ""))
    pd.DataFrame(rows).to_csv(os.path.join(C.DIR_OUT, "auc_results.csv"), index=False)

    tag = f"{C.SITE}_{C.DATE.replace('-', '')}"
    for key, short in (("s2", "S2"), ("ls", "LS")):
        p = os.path.join(C.DIR_OUT, f"auc_{short}_{tag}.xlsx")
        fill_workbook(df, key, p, results)
        log(f"wrote {os.path.basename(p)}")
    log("STEP 5 DONE. Send me the two auc_*.xlsx files and auc_results.csv. "
        "Optional: python step6_maps.py")


if __name__ == "__main__":
    main()
