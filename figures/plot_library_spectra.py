"""
Reference spectra of the minerals tested in the paper, with the Sentinel-2 and Landsat 8/9 bands.

Source: USGS Spectral Library Version 7 (Kokaly et al. 2017, https://doi.org/10.5066/F7RR1WDJ).
One-time download: from that page, get  ASCIIdata_splib07a.zip  (measured spectra, about 0.5 GB).
Leave it zipped; the script reads it directly.

Run:
  python plot_library_spectra.py --lib path/to/ASCIIdata_splib07a.zip
  python plot_library_spectra.py --lib ... --pick calcite=WS272 --pick kaolinite=CM9   (choose samples)
  python plot_library_spectra.py --lib ... --list kaolinite                            (list candidates)

Writes to ./library_spectra/:
  Fig_library_spectra.png (300 dpi) and .pdf
      (a) reflectance spectra 0.4-2.5 um, offset for clarity, with the bands used by the ratios
      (b) detail of the SWIR2 band with the diagnostic absorptions
  library_spectra_bands.csv   band reflectances and ratios per mineral and sensor
  library_spectra_samples.txt the USGS sample names used (quote them in the caption)

Bands are simulated with flat-top responses from centre and width, as in the physics control of the paper.
"""
import argparse
import csv
import io
import os
import re
import sys
import zipfile

import numpy as np

# ------------------------------------------------------------------ what to plot
# key: (label, search words in the USGS file name, preferred sample tokens in order)
MINERALS = [
    ("calcite", "Calcite", ["Calcite"], ["WS272", "HS48", "CO2004"]),
    ("dolomite", "Dolomite", ["Dolomite"], ["HS102", "COD2005", "SR-02"]),
    ("kaolinite", "Kaolinite", ["Kaolinite"], ["CM9", "KGa-1"]),
    ("illite", "Illite", ["Illite"], ["IMt-1", "GDS4"]),
    ("muscovite", "Muscovite", ["Muscovite"], ["GDS107", "HS146"]),
    ("montmorillonite", "Montmorillonite", ["Montmorillonite"], ["SWy-1", "SWy-2", "CM20"]),
    ("chlorite", "Chlorite (clinochlore)", ["Clinochlore", "Chlorite"], ["GDS158", "GDS157", "SMR-13"]),
    ("hematite", "Hematite", ["Hematite"], ["GDS27", "FE2602"]),
    ("goethite", "Goethite", ["Goethite"], ["WS222", "GDS134"]),
]
# instruments covering 0.4-2.5 um, best first (file names carry e.g. _BECKa_ or _ASDFRb_)
INSTRUMENT_ORDER = ["BECK", "ASDFR", "ASDNG", "ASDHR", "AVIRIS"]
INSTRUMENT_RE = re.compile(r"_(BECK|ASDFR|ASDNG|ASDHR|AVIRIS)[a-z]_")

# band centre and full width (nm), from Table 1 of the paper
BANDS = {
    "s2": {"name": "Sentinel-2", "blue": (492, 66), "red": (665, 31), "swir1": (1614, 91), "swir2": (2202, 175),
           "labels": {"blue": "B2", "red": "B4", "swir1": "B11", "swir2": "B12"}, "color": "#2a78d6"},
    "ls": {"name": "Landsat 8/9", "blue": (482, 60), "red": (655, 37), "swir1": (1609, 85), "swir2": (2201, 187),
           "labels": {"blue": "B2", "red": "B4", "swir1": "B6", "swir2": "B7"}, "color": "#eb6834"},
}
# diagnostic absorptions (um) marked on panel (a)
ABSORPTIONS = [(0.48, "Fe³⁺"), (0.87, "Fe³⁺"), (2.165, "kaolinite doublet"), (2.205, "Al-OH"),
               (2.25, "Fe-OH (chlorite)"), (2.32, "dolomite"), (2.34, "calcite / Mg-OH")]
INK, MUTED = "#222222", "#666666"
MISSING = -1.0e30          # USGS flags deleted channels with -1.23e34


# ------------------------------------------------------------------ reading the library
class Library:
    def __init__(self, path):
        self.path = path
        self.zip = zipfile.ZipFile(path) if path.lower().endswith(".zip") else None
        names = self.zip.namelist() if self.zip else [
            os.path.relpath(os.path.join(r, f), path) for r, _, fs in os.walk(path) for f in fs]
        self.names = [n.replace("\\", "/") for n in names if n.lower().endswith(".txt")]
        self.wavefiles = [n for n in self.names if "Wavelengths" in os.path.basename(n)]
        self.minerals = [n for n in self.names if "ChapterM_Minerals" in n]
        if not self.minerals:
            sys.exit("No ChapterM_Minerals spectra found: is this ASCIIdata_splib07a.zip?")
        self._wl = {}

    def read(self, name):
        if self.zip:
            txt = self.zip.read(name).decode("latin-1")
        else:
            txt = open(os.path.join(self.path, name), encoding="latin-1").read()
        lines = txt.strip().splitlines()
        vals = np.array([float(x) for x in lines[1:] if x.strip()])
        return lines[0].strip(), vals

    def wavelengths(self, n_values, instrument):
        """Wavelength file (um) with the same channel count, for the same instrument family."""
        key = (n_values, instrument)
        if key not in self._wl:
            fam = "ASD" if instrument.startswith("ASD") else instrument          # BECK, ASD, AVIRIS
            hit = None
            for wf in self.wavefiles:
                if fam.lower() not in os.path.basename(wf).lower():
                    continue
                _, w = self.read(wf)
                if len(w) == n_values:
                    hit = w
                    break
            if hit is None:
                sys.exit(f"No wavelength file with {n_values} channels for {instrument}")
            self._wl[key] = hit
        return self._wl[key]

    def candidates(self, words):
        out = []
        for n in self.minerals:
            base = os.path.basename(n)
            if "Wavelengths" in base or "Bandpass" in base:
                continue
            if any(re.search(rf"_{w}[_ ]", base, re.I) for w in words):
                m = INSTRUMENT_RE.search(base)
                if m and base.endswith("AREF.txt"):
                    out.append((n, m.group(1)))
        return out


def choose(lib, key, words, prefs, pick):
    cands = lib.candidates(words)
    if not cands:
        sys.exit(f"No spectrum found for {key} ({words})")
    if pick:
        cands = [c for c in cands if pick.lower() in os.path.basename(c[0]).lower()] or sys.exit(
            f"--pick {key}={pick} matches nothing; use --list {key}")

    def score(c):
        base = os.path.basename(c[0])
        pref = next((i for i, p in enumerate(prefs) if p.lower() in base.lower()), len(prefs))
        return (pref, INSTRUMENT_ORDER.index(c[1]), base)
    return sorted(cands, key=score)[0]


def spectrum(lib, name, inst):
    header, refl = lib.read(name)
    wl = lib.wavelengths(len(refl), inst)
    ok = (refl > MISSING) & np.isfinite(refl) & (wl >= 0.35) & (wl <= 2.55)
    return header, wl[ok], refl[ok]


def band_value(wl_um, refl, centre_nm, width_nm):
    lo, hi = (centre_nm - width_nm / 2) / 1000, (centre_nm + width_nm / 2) / 1000
    sel = (wl_um >= lo) & (wl_um <= hi)
    return float(np.mean(refl[sel])) if sel.any() else float(np.interp(centre_nm / 1000, wl_um, refl))


def smooth(y, n=1):
    """Running mean for display only. Off by default: on Beckman sampling (about 0.006 um) even a
    5-channel mean blurs narrow features such as the kaolinite doublet."""
    if n <= 1 or len(y) < n:
        return y
    k = np.ones(n) / n
    return np.convolve(np.pad(y, n // 2, mode="edge"), k, mode="valid")[: len(y)]


# ------------------------------------------------------------------ figure
def make_figure(spectra, table, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": "#999999",
                         "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.labelcolor": INK, "savefig.dpi": 300, "savefig.bbox": "tight"})
    n = len(spectra)
    # offset between spectra: most of the largest reflectance range, so curves rarely overlap
    step = round(0.8 * max(float(np.ptp(smooth(r)[(wl >= 0.4) & (wl <= 2.5)])) for _, _, wl, r in spectra), 2)
    step = max(step, 0.25)
    lo_all = min(float(np.min(smooth(r))) for _, _, _, r in spectra)
    fig = plt.figure(figsize=(6.3, 1.0 + 0.62 * n), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0])
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    top = (n - 1) * step + max(float(np.max(smooth(r))) for _, _, _, r in spectra) + 0.12 * step + 0.25
    head = 0.18 * step + 0.12           # room for the band labels

    def draw(ax, x0, x1, bands, absorptions):
        for k, s in BANDS.items():
            for b in bands:
                c, w = s[b]
                lo, hi = (c - w / 2) / 1000, (c + w / 2) / 1000
                if k == "s2":                       # Sentinel-2: shaded column
                    ax.axvspan(lo, hi, color=s["color"], alpha=0.13, lw=0, zorder=0)
                    ax.text((lo + hi) / 2, top + head * 0.55, s["labels"][b], color=s["color"],
                            ha="center", va="bottom", fontsize=6.5)
                else:                               # Landsat: bar along the top
                    ax.plot([lo, hi], [top + head * 0.4] * 2, color=s["color"], lw=2.4,
                            solid_capstyle="butt")
                    ax.text((lo + hi) / 2, top + head * 0.3, s["labels"][b], color=s["color"], ha="center",
                            va="top", fontsize=6.5)
        if absorptions:
            for j, (x, lab) in enumerate(ABSORPTIONS):
                if x0 < x < x1:
                    ax.axvline(x, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=0.5)
                    ax.text(x, 1.01, lab, rotation=90, ha="center", va="bottom", fontsize=6, color=MUTED,
                            transform=ax.get_xaxis_transform())
        for i, (key, label, wl, r) in enumerate(spectra):
            off = (n - 1 - i) * step
            sel = (wl >= x0 - 0.02) & (wl <= x1 + 0.02)
            ax.plot(wl[sel], smooth(r)[sel] + off, color=INK, lw=0.8)
            for k, s in BANDS.items():
                for b in bands:
                    c, _ = s[b]
                    if x0 <= c / 1000 <= x1:
                        ax.plot(c / 1000, table[(key, k)][b] + off, "o", ms=3.2, mfc=s["color"], mec="white",
                                mew=0.4, zorder=3)
        ax.set_xlim(x0, x1)
        ax.set_ylim(lo_all - 0.05, top + head * 1.2)
        ax.set_yticks([])
        ax.set_xlabel("Wavelength (\u00b5m)")
        ax.tick_params(labelsize=7, length=2)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)

    draw(ax_a, 0.40, 2.50, ("blue", "red", "swir1", "swir2"), absorptions=False)
    ax_a.set_title("(a) Library spectra, 0.4\u20132.5 \u00b5m", loc="left", color=INK, fontsize=8.5, pad=58)
    ax_a.set_ylabel("Reflectance (offset for clarity)")
    for i, (key, label, wl, r) in enumerate(spectra):
        near = (wl >= 0.40) & (wl <= 1.00)        # label sits above the curve over the width of the text
        ax_a.text(0.41, (n - 1 - i) * step + float(np.max(r[near])) + 0.03, label, fontsize=6.8,
                  color=INK, va="bottom")
    draw(ax_b, 2.05, 2.42, ("swir2",), absorptions=True)
    ax_b.set_title("(b) SWIR2 detail", loc="left", color=INK, fontsize=8.5, pad=58)

    fig.savefig(os.path.join(out, "Fig_library_spectra.png"))
    fig.savefig(os.path.join(out, "Fig_library_spectra.pdf"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot USGS library spectra with Sentinel-2 and Landsat bands.")
    ap.add_argument("--lib", required=True, help="ASCIIdata_splib07a.zip or its extracted folder")
    ap.add_argument("--out", default="library_spectra")
    ap.add_argument("--pick", action="append", default=[], help="mineral=text in the file name, e.g. calcite=WS272")
    ap.add_argument("--list", metavar="MINERAL", help="list the candidate spectra for one mineral and stop")
    ap.add_argument("--only", help="comma-separated subset, e.g. calcite,dolomite,kaolinite,chlorite,hematite")
    args = ap.parse_args()

    lib = Library(args.lib)
    if args.list:
        m = next(m for m in MINERALS if m[0] == args.list)
        for n, inst in lib.candidates(m[2]):
            print(inst.ljust(8), os.path.basename(n))
        return
    picks = dict(p.split("=", 1) for p in args.pick)
    wanted = [m for m in MINERALS if not args.only or m[0] in args.only.split(",")]
    os.makedirs(args.out, exist_ok=True)

    spectra, table, used = [], {}, []
    for key, label, words, prefs in wanted:
        name, inst = choose(lib, key, words, prefs, picks.get(key))
        header, wl, r = spectrum(lib, name, inst)
        spectra.append((key, label, wl, r))
        used.append(f"{label}: {os.path.basename(name)}  |  {header}")
        for k, s in BANDS.items():
            table[(key, k)] = {b: band_value(wl, r, *s[b]) for b in ("blue", "red", "swir1", "swir2")}
        print(f"{label:<24} {os.path.basename(name)}")

    with open(os.path.join(args.out, "library_spectra_bands.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mineral", "sensor", "blue", "red", "swir1", "swir2", "iron_red_blue", "clay_swir1_swir2",
                    "carbonate_swir2_swir1"])
        for key, label, _, _ in spectra:
            for k in BANDS:
                v = table[(key, k)]
                w.writerow([label, BANDS[k]["name"]] + [f"{v[b]:.4f}" for b in ("blue", "red", "swir1", "swir2")]
                           + [f"{v['red'] / v['blue']:.3f}", f"{v['swir1'] / v['swir2']:.3f}",
                              f"{v['swir2'] / v['swir1']:.3f}"])
    open(os.path.join(args.out, "library_spectra_samples.txt"), "w", encoding="utf-8").write("\n".join(used) + "\n")
    make_figure(spectra, table, args.out)
    print(f"wrote {args.out}/Fig_library_spectra.png/.pdf, library_spectra_bands.csv, library_spectra_samples.txt")


if __name__ == "__main__":
    main()
