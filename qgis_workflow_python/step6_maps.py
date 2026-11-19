"""
STEP 6 (optional) - Maps of EMIT abundance against a band ratio (paper figures and check maps).

Run (after step 5):
  python step6_maps.py                          carbonate figure for Sentinel-2 and Landsat + check maps
  python step6_maps.py --target fe --sensor s2  iron figure, Sentinel-2 only
  python step6_maps.py --target all             carbonate, clay and iron figures
  python step6_maps.py --no-overview            skip the 2 x 2 check maps

Writes to runs/<site>_<date>/results/:
  fig_<target>_<S2|LS>_<site>_<date>.png (300 dpi) and .pdf   three-panel figure:
      (a) EMIT abundance of the target
      (b) the band ratio, with the outline of the pixels in the top 20 % of EMIT abundance
      (c) the ratio of the top-20 % and bottom-20 % EMIT pixels, with the AUC
  fig_<target>_<S2|LS>_<site>_<date>_caption.txt            draft caption: dates, day gap, colour stretches, AUC, n
  maps_S2.png, maps_LS.png (300 dpi)                          2 x 2 check maps: carbonate and iron, both as above

Both maps of a figure show exactly the pixels of the AUC test: mid grey = excluded (cloud or shadow,
NDVI >= NDVI_MAX, invalid EMIT, and for iron other Fe absorptions); light grey = no EMIT data.
Panel (c) and the AUC use every usable pixel in the box, as the processing script behind the paper
tables does (Table S1). With --sample they use the random pixels of step 5 instead (= auc_results.csv).
"""
import argparse
import json
import os
from datetime import date

import numpy as np

import config as C
from common import log, ensure_dirs, load_grid, read_tif, auc_quantile, verdict
from step5_sample_and_auc import sample, selection

# ---------------------------------------------------------------- what each target means
TARGETS = {
    "carbonate": dict(emit_layer="abund_carbonate", emit_col="emit_carb", suffix="carb",
                      emit_title="EMIT carbonate abundance", ratio_title="‘Carbonate’ ratio",
                      bands={"s2": "B12/B11", "ls": "B7/B6"}, cmap="YlGnBu", outline="#d7301f"),
    "clay": dict(emit_layer="abund_clay", emit_col="emit_clay", suffix="clay",
                 emit_title="EMIT Al-OH clay abundance", ratio_title="‘Clay’ ratio",
                 bands={"s2": "B11/B12", "ls": "B6/B7"}, cmap="YlGnBu", outline="#d7301f"),
    "fe": dict(emit_layer="abund_fe", emit_col="emit_fe", suffix="fe",
               emit_title="EMIT iron-oxide abundance", ratio_title="Iron ratio",
               bands={"s2": "B4/B2", "ls": "B4/B2"}, cmap="OrRd", outline="#2166ac"),
}
SENSORS = {"s2": dict(name="Sentinel-2", short="S2", color="#2a78d6", folder=C.DIR_S2),
           "ls": dict(name="Landsat 8/9", short="LS", color="#eb6834", folder=C.DIR_LS)}
GREY_MASKED, GREY_NODATA, INK, MUTED = "#a9a9a9", "#e3e3e3", "#222222", "#666666"
ABUND_UNIT = "band depth × group fraction"


def style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": "#999999",
                         "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.labelcolor": INK, "axes.titlesize": 8.5, "savefig.dpi": 300,
                         "savefig.bbox": "tight"})


def site_name():
    return C.SITE.replace("_", " ").title()


def scene_days(key):
    """Acquisition day(s) of the Sentinel-2 or Landsat scene(s) and the gap to the EMIT date in days."""
    path = os.path.join(SENSORS[key]["folder"], "items.json")
    if not os.path.exists(path):
        return []
    d0 = date.fromisoformat(C.DATE)
    days = sorted({date.fromisoformat(it["datetime"][:10]) for it in json.load(open(path))})
    return [(d.isoformat(), (d - d0).days) for d in days]


def days_text(key):
    days = scene_days(key)
    if not days:
        return "date unknown (items.json missing)"
    return ", ".join(f"{d} ({g:+d} d)" for d, g in days)


OUTLINE_WINDOW = 5   # pixels; the outline is drawn around areas where most pixels of a 5 x 5 window are top-20 %


def majority(mask, n):
    """Share of True pixels in an n x n window (edges padded with False), for a readable outline."""
    if n <= 1:
        return mask.astype(float)
    k = n // 2
    a = np.pad(mask.astype(float), ((k + 1, k), (k + 1, k)))
    c = a.cumsum(0).cumsum(1)
    s = c[n:, n:] - c[:-n, n:] - c[n:, :-n] + c[:-n, :-n]
    return s / (n * n)


# ---------------------------------------------------------------- data for one figure
def prepare(key, target, df):
    """Grids and test pixels for one sensor and one target, with the same rules as step 5."""
    t = TARGETS[target]
    abund = read_tif(t["emit_layer"])
    ratio = read_tif(f"{key}_{t['suffix']}")
    emit_valid = read_tif("emit_valid") == 1
    used = (read_tif(f"{key}_ok") == 1) & np.isfinite(ratio) & np.isfinite(abund)
    if target == "fe" and C.EXCLUDE_OTHER_IRON:
        used &= read_tif("fe_other") == 0

    # test pixels and AUC exactly as in step 5
    if df is None:                                   # every usable pixel (paper tables)
        a, r = abund[used], ratio[used]
    else:                                            # the random pixels of step 5
        sel = selection(df, key, t["suffix"])
        a = sel[t["emit_col"]].to_numpy(float)
        r = sel[f"{key}_{t['suffix']}"].to_numpy(float)
    keep = np.isfinite(a) & np.isfinite(r)
    a, r = a[keep], r[keep]
    auc, n1, n0 = auc_quantile(a, r)
    hi, lo = np.quantile(a, 1 - C.QUANTILE), np.quantile(a, C.QUANTILE)
    pos = (a >= hi) & (a > 0)
    neg = (a <= lo) & ~pos

    # background classes: 0 = shown, 1 = excluded from the test, 2 = no EMIT data
    bg = np.where(used, 0, np.where(emit_valid, 1, 2))
    top = majority(used & (abund >= hi) & (abund > 0), OUTLINE_WINDOW) >= 0.5
    vals = abund[used & (abund > 0)]
    vmax = float(np.percentile(vals, 99)) if vals.size else 1.0
    rlo, rhi = (np.percentile(ratio[used], (2, 98)) if used.any() else (0.0, 1.0))
    return dict(abund=np.where(used, abund, np.nan), ratio=np.where(used, ratio, np.nan), bg=bg, top=top,
                vmax=float(f"{vmax:.2g}"), rlim=(float(rlo), float(rhi)), auc=auc, n1=n1, n0=n0,
                r_top=r[pos], r_low=r[neg], hi=hi)


# ---------------------------------------------------------------- drawing helpers
def geo_axes(ax, extent):
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    lat0 = (extent[2] + extent[3]) / 2
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(lat0)))      # true proportions at this latitude
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.2f}°{'W' if v < 0 else 'E'}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.2f}°{'S' if v < 0 else 'N'}"))
    ax.tick_params(labelsize=7, length=2)
    for lab in ax.get_yticklabels():
        lab.set_rotation(90)
        lab.set_va("center")
    scale_bar(ax, extent)


def scale_bar(ax, extent):
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    lat0 = (extent[2] + extent[3]) / 2
    km_per_deg = 111.32 * np.cos(np.deg2rad(lat0))
    width_km = (extent[1] - extent[0]) * km_per_deg
    km = max(k for k in (1, 2, 5, 10, 20, 50, 100) if k <= width_km / 4)
    bar = AnchoredSizeBar(ax.transData, km / km_per_deg, f"{km} km", "lower right", pad=0.3, sep=2,
                          borderpad=0.4, frameon=True, size_vertical=(extent[3] - extent[2]) / 150,
                          color=INK, fontproperties={"size": 6.5})
    bar.patch.set(alpha=0.8, edgecolor="none")
    ax.add_artist(bar)


def background(ax, bg, extent):
    from matplotlib.colors import ListedColormap
    ax.imshow(np.ma.masked_equal(bg, 0), extent=extent, cmap=ListedColormap([GREY_MASKED, GREY_NODATA]),
              vmin=1, vmax=2, interpolation="nearest", zorder=0)


def map_panel(fig, ax, grid, extent, bg, cmap, lim, title, cbar_label, outline=None, outline_color=None):
    background(ax, bg, extent)
    im = ax.imshow(grid, extent=extent, cmap=cmap, vmin=lim[0], vmax=lim[1], interpolation="nearest", zorder=1)
    if outline is not None and outline.any():
        h, w = outline.shape
        xs = np.linspace(extent[0], extent[1], w + 1)
        ys = np.linspace(extent[3], extent[2], h + 1)
        xc, yc = (xs[:-1] + xs[1:]) / 2, (ys[:-1] + ys[1:]) / 2
        ax.contour(xc, yc, outline.astype(float), levels=[0.5], colors=[outline_color],
                   linewidths=0.5, zorder=2)
    geo_axes(ax, extent)
    ax.set_title(title, loc="left", color=INK)
    cax = ax.inset_axes([0.0, -0.2, 1.0, 0.035])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="max")
    from matplotlib.ticker import MaxNLocator
    cb.ax.xaxis.set_major_locator(MaxNLocator(5))
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.set_label(cbar_label, fontsize=7, color=INK)
    cb.outline.set_linewidth(0.4)


def dist_panel(ax, p, key, target):
    t, s = TARGETS[target], SENSORS[key]
    bp = ax.boxplot([p["r_top"], p["r_low"]], widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color=INK, lw=1), whiskerprops=dict(color=MUTED, lw=0.7),
                    capprops=dict(color=MUTED, lw=0.7))
    bp["boxes"][0].set(facecolor=s["color"], edgecolor=s["color"], alpha=0.85)
    bp["boxes"][1].set(facecolor="white", edgecolor=s["color"])
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Top 20 %\nn = {p['n1']}", f"Bottom 20 %\nn = {p['n0']}"], fontsize=7)
    ax.set_xlabel("EMIT abundance class", fontsize=7)
    ax.set_ylabel(f"{s['name']} {t['bands'][key]}", fontsize=7)
    ax.tick_params(labelsize=6.5, length=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(f"(c) Test pixels\nAUC = {p['auc']:.2f}, {verdict(p['auc']).lower()}", loc="left", color=INK)


def legend_row(fig, target, key):
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    t = TARGETS[target]
    hand = [Line2D([], [], color=t["outline"], lw=1, label="Areas in the top 20 % of EMIT abundance"),
            Patch(color=GREY_MASKED, label="Excluded from the test"),
            Patch(color=GREY_NODATA, label="No EMIT data")]
    fig.legend(handles=hand, loc="lower center", ncol=3, frameon=False, fontsize=7,
               bbox_to_anchor=(0.5, -0.05))


# ---------------------------------------------------------------- figures
def paper_figure(key, target, df, extent, title=False):
    import matplotlib.pyplot as plt
    t, s = TARGETS[target], SENSORS[key]
    p = prepare(key, target, df)
    fig = plt.figure(figsize=(7.2, 3.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.5])
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))
    map_panel(fig, ax_a, p["abund"], extent, p["bg"], t["cmap"], (0, p["vmax"]),
              f"(a) {t['emit_title']}", ABUND_UNIT)
    map_panel(fig, ax_b, p["ratio"], extent, p["bg"], t["cmap"], p["rlim"],
              f"(b) {t['ratio_title']}, {s['name']} {t['bands'][key]}", f"{t['bands'][key]} (unitless)",
              outline=p["top"], outline_color=t["outline"])
    dist_panel(ax_c, p, key, target)
    legend_row(fig, target, key)
    if title:
        fig.suptitle(f"{site_name()}, EMIT {C.DATE}, {s['name']} {days_text(key)}", fontsize=9)

    stem = os.path.join(C.DIR_OUT, f"fig_{target}_{s['short']}_{C.SITE}_{C.DATE.replace('-', '')}")
    fig.savefig(stem + ".png")
    fig.savefig(stem + ".pdf")
    plt.close(fig)

    caption = (
        f"{site_name()}: EMIT scene {C.DATE} and {s['name']} scene {days_text(key)}, "
        f"on the EMIT grid (about 60 m). "
        f"(a) {t['emit_title']} ({ABUND_UNIT}); colour scale 0 to {p['vmax']} (99th percentile of non-zero values). "
        f"(b) {t['ratio_title']} ({s['name']} {t['bands'][key]}); colour scale {p['rlim'][0]:.3f} to "
        f"{p['rlim'][1]:.3f} (2nd to 98th percentile). Outline: areas where most pixels of a "
        f"{OUTLINE_WINDOW} × {OUTLINE_WINDOW} window are in the top {int(C.QUANTILE * 100)} % of EMIT abundance "
        f"(≥ {p['hi']:.3g}). "
        f"(c) Ratio of the top-{int(C.QUANTILE * 100)} % (n = {p['n1']}) and bottom-{int(C.QUANTILE * 100)} % "
        f"(n = {p['n0']}) pixels{' of the random test sample' if df is not None else ', all usable pixels in the box'}; boxes show the quartiles, whiskers 1.5 × IQR. "
        f"AUC = {p['auc']:.2f} ({verdict(p['auc']).lower()}). "
        f"Mid grey: pixels excluded from the test (cloud or shadow, NDVI ≥ {C.NDVI_MAX}, invalid EMIT"
        f"{', other Fe absorptions' if target == 'fe' else ''}); light grey: no EMIT data.")
    open(stem + "_caption.txt", "w", encoding="utf-8").write(caption + "\n")
    log(f"wrote {os.path.basename(stem)}.png/.pdf  AUC {p['auc']:.3f} ({verdict(p['auc'])}), "
        f"{s['name']} {days_text(key)}")


def overview(key, df, extent):
    """2 x 2 check maps: carbonate and iron, EMIT abundance next to the ratio."""
    import matplotlib.pyplot as plt
    s = SENSORS[key]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.6), constrained_layout=True)
    for row, target in enumerate(("carbonate", "fe")):
        t = TARGETS[target]
        p = prepare(key, target, df)
        map_panel(fig, axes[row, 0], p["abund"], extent, p["bg"], t["cmap"], (0, p["vmax"]),
                  f"({'ac'[row]}) {t['emit_title']}", ABUND_UNIT)
        map_panel(fig, axes[row, 1], p["ratio"], extent, p["bg"], t["cmap"], p["rlim"],
                  f"({'bd'[row]}) {t['ratio_title']}, {t['bands'][key]}  AUC {p['auc']:.2f}",
                  f"{t['bands'][key]} (unitless)", outline=p["top"], outline_color=t["outline"])
    fig.suptitle(f"{site_name()}, EMIT {C.DATE}, {s['name']} {days_text(key)}", fontsize=9)
    out = os.path.join(C.DIR_OUT, f"maps_{s['short']}.png")
    fig.savefig(out)
    plt.close(fig)
    log(f"wrote {os.path.basename(out)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--target", choices=list(TARGETS) + ["all"], default="carbonate")
    ap.add_argument("--sensor", choices=["s2", "ls", "both"], default="both")
    ap.add_argument("--no-overview", action="store_true", help="skip the 2 x 2 check maps")
    ap.add_argument("--title", action="store_true", help="add a title line to the paper figures")
    ap.add_argument("--sample", action="store_true",
                    help="AUC and panel (c) from the random pixels of step 5 instead of every usable pixel")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    style()
    ensure_dirs()
    tr, w, h = load_grid()
    extent = (tr.c, tr.c + w * tr.a, tr.f + h * tr.e, tr.f)
    df = sample() if args.sample else None            # None = every usable pixel
    keys = ["s2", "ls"] if args.sensor == "both" else [args.sensor]
    targets = list(TARGETS) if args.target == "all" else [args.target]
    for key in keys:
        log(f"{SENSORS[key]['name']}: scene day(s) {days_text(key)} (EMIT {C.DATE})")
        for target in targets:
            paper_figure(key, target, df, extent, title=args.title)
        if not args.no_overview:
            overview(key, df, extent)
    log("STEP 6 DONE.")


if __name__ == "__main__":
    main()
