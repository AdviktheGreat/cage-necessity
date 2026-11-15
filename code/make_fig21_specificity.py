"""Is the designed accessibility a peak AT each promoter, or just a more open locus?

This is the test behind correction 22. Ranking designs by DNase height at their weakest
promoter rewards a design for raising accessibility everywhere, so height has to be read
alongside locality. Background is taken on sequence no objective term touched: outside
every CDS and outside the +-1,057 bp windows ProCapNet saw.

Writes data/tables/table16_promoter_specificity.csv and
figures/fig21_peak_vs_specificity.png.
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRED = os.path.join(REPO, "data", "predictions")
Z = np.load(os.path.join(PRED, "borzoi_profiles.npz"))
T = np.load(os.path.join(PRED, "dnase_traces.npz"))
MARGIN, N = int(T["margin"][0]), int(T["design_len"][0])
TRK = [str(t) for t in Z["tracks"]]
iD = TRK.index("DNase")
GENES = {"eGFP": (4462, 5182), "mScarlet": (9644, 10343), "TagBFP": (14805, 15507)}
TSS = {g: a - 128 for g, (a, z) in GENES.items()}
KEY = ["control", "final", "B_cage_rna", "B_cage_rna_floor", "PC_rna_procap",
       "D_cage_rna_procap_floor"]
NICE = {"control": "CDS-only control", "final": "lab's final", "B_cage_rna": "CAGE+RNA",
        "B_cage_rna_floor": "CAGE+RNA+floor", "PC_rna_procap": "RNA+PRO-cap",
        "D_cage_rna_procap_floor": "CAGE+RNA+PRO-cap\n+floors"}
PC_HALF = 1100          # ProCapNet input is 2,114 bp centred on each TSS


def trace(key, model):
    """DNase across the design window at 1 bp, for one model."""
    if model == "borzoi":
        return np.repeat(Z[key][iD], 32)[:N]
    return T[f"{model}|{key}"][MARGIN:MARGIN + N]


def untouched_mask():
    m = np.ones(N, dtype=bool)
    for g, (a, z) in GENES.items():
        m[max(0, TSS[g] - PC_HALF):min(N, TSS[g] + PC_HALF)] = False
        m[a:z] = False
    return m


BG = untouched_mask()
MODELS = [("borzoi", "Borzoi"), ("ag1", "AlphaGenome 1"), ("ag2", "AlphaGenome 2")]
R = {}
for key in KEY:
    for mid, _ in MODELS:
        a = trace(key, mid)
        bg = float(np.median(a[BG]))
        per = [float(a[max(0, TSS[g] - 300):TSS[g] + 300].max()) for g in GENES]
        R[(key, mid)] = (min(per), bg, min(per) / bg if bg > 0 else float("nan"), per)

print(f"background: {BG.sum()} bp of {N} ({100 * BG.sum() / N:.0f}%) untouched by any term")
print(f"{'design':30s} {'model':14s} {'worst':>7s} {'bg':>7s} {'specificity':>12s}")
for key in KEY:
    for mid, nm in MODELS:
        w, bg, sp, _ = R[(key, mid)]
        print(f"{key if mid == 'borzoi' else '':30s} {nm:14s} {w:7.3f} {bg:7.3f} {sp:12.2f}")

with open(os.path.join(REPO, "data", "tables", "table16_promoter_specificity.csv"),
          "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["design", "model", "dnase_egfp", "dnase_mscarlet", "dnase_tagbfp",
                "dnase_worst_promoter", "background_median_untouched", "specificity_ratio"])
    for key in KEY:
        for mid, _ in MODELS:
            worst, bg, sp, per = R[(key, mid)]
            w.writerow([key, {"borzoi": "borzoi", "ag1": "alphagenome1",
                              "ag2": "alphagenome2"}[mid]]
                       + [round(v, 3) for v in per] + [round(worst, 3), round(bg, 4), round(sp, 2)])

fig, ax = plt.subplots(1, 2, figsize=(14.5, 6.2))
for j, (mid, nm) in enumerate([("ag1", "AlphaGenome 1"), ("ag2", "AlphaGenome 2")]):
    A = ax[j]
    for key in KEY:
        h, _, sp, _ = R[(key, mid)]
        c = ("#cc2222" if key == "final"
             else "#228833" if key == "D_cage_rna_procap_floor" else "#4477aa")
        A.scatter(h, sp, s=150, color=c, zorder=3, edgecolor="k", linewidth=.6)
        off = (-18, -26) if key == "D_cage_rna_procap_floor" else (10, 6)
        A.annotate(NICE[key], (h, sp), textcoords="offset points", xytext=off, fontsize=8.5)
    A.set_xscale("log")
    A.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
    A.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    A.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    A.set_xlim(0.09, 12)
    A.set_xlabel("DNase at the WORST promoter  (how strong)", fontsize=10)
    A.set_ylabel("promoter DNase ÷ untouched background  (how LOCAL)", fontsize=10)
    A.set_title(f"{nm} (in no design loss)", fontsize=11, fontweight="bold")
    A.grid(alpha=.25)
    A.axhspan(0, 10, color="grey", alpha=.08)
    A.text(.02, .03, "below here the peak barely stands out\nfrom the rest of the construct",
           transform=A.transAxes, fontsize=7.5, color="dimgrey")
fig.suptitle("\"DNase peaks at promoters\" has two axes, and no design wins both.\n"
             "The lab's own design is the most LOCAL; the dual-oracle arm is the STRONGEST "
             "but partly by opening the whole locus.", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, .90])
out = os.path.join(REPO, "figures", "fig21_peak_vs_specificity.png")
fig.savefig(out, dpi=140)
print("wrote", out)
