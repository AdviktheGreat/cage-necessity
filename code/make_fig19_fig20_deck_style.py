import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

D = "/tmp/dualres/"
REPO = "/Users/advikanand/Downloads/borzoi-cage-necessity/"
Z = np.load(D + "borzoi_profiles.npz")
TRK = list(Z["tracks"]); BINW = 32
iC, iR, iD = TRK.index("CAGE+"), TRK.index("RNA+"), TRK.index("DNase")
tgt_c, tgt_r = Z["target_cage"], Z["target_rna"]
GENES = {"eGFP": (4462, 5182), "mScarlet": (9644, 10343), "TagBFP": (14805, 15507)}

ROWS = [("control", "CDS-only control\n(no edits)", {"CAGE": 0, "RNA": 0, "DNase": 0}),
        ("final", "lab's final\nCAGE + RNA", {"CAGE": 1, "RNA": 1, "DNase": 0}),
        ("B_cage_rna", "CAGE + RNA\n(my reproduction)", {"CAGE": 1, "RNA": 1, "DNase": 0}),
        ("PC_rna_procap", "RNA + PRO-cap\n(no CAGE)", {"CAGE": 0, "RNA": 1, "DNase": 0}),
        ("D_cage_rna_procap_floor", "CAGE + RNA + PRO-cap\n+ per-gene floors", {"CAGE": 1, "RNA": 1, "DNase": 0})]
COLS = [("CAGE", iC, tgt_c, "CAGE (initiation)"),
        ("RNA", iR, tgt_r, "RNA-seq (abundance)"),
        ("DNase", iD, None, "DNase (accessibility)")]

fig, axes = plt.subplots(len(ROWS), 3, figsize=(15.5, 13.2), sharex=True)
x = np.arange(624)
for ri, (key, lab, inloss) in enumerate(ROWS):
    prof = Z[key]
    for ci, (cname, ti, tgt, ctitle) in enumerate(COLS):
        ax = axes[ri, ci]
        used = inloss[cname]
        if tgt is not None:
            ax.fill_between(x, 0, tgt, step="mid", color="#ff9900", alpha=.30,
                            label="requested pattern", zorder=1)
            ax.step(x, tgt, where="mid", color="#ff7700", lw=1.0, zorder=2)
        ax.plot(x, prof[ti], color="#1f4e9c" if used else "#cc2222", lw=.85,
                label="Borzoi prediction", zorder=3)
        for g, (a, _) in GENES.items():
            ax.axvline((a - 128) // BINW, color="grey", ls=":", lw=.6, zorder=0)
        tag = "IN THE LOSS" if used else "HELD OUT"
        ax.text(.985, .93, tag, transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=.3", fc="#1f4e9c" if used else "#cc2222", ec="none"))
        if ri == 0: ax.set_title(ctitle, fontsize=11, fontweight="bold")
        if ci == 0: ax.set_ylabel(lab, fontsize=9)
        if ri == len(ROWS) - 1: ax.set_xlabel("Prediction bin across construct", fontsize=9)
        ax.grid(alpha=.18); ax.tick_params(labelsize=7)
axes[0, 0].legend(fontsize=6.5, loc="upper left")
for ax, g in zip(axes[0], ["", "", ""]): pass
fig.suptitle("Borzoi prediction against the requested pattern, across the construct\n"
             "Blue = this readout was in the design objective.  Red = it was not.  "
             "DNase was never specified in any objective in this project.",
             fontsize=12.5, y=.995)
fig.tight_layout(rect=[0, 0, 1, .965])
fig.savefig(REPO + "figures/fig19_designed_vs_heldout.png", dpi=140)
print("wrote fig19_designed_vs_heldout.png")

# ---------------- fig20: three-model DNase traces, extending deck slide 53 ----------------
A1 = np.load(D + "ag_tracks.npz"); A2 = np.load(D + "ag2_tracks.npz")
MARG = 1000
SHOW = [("control", "CDS-only control"), ("final", "lab's final (CAGE+RNA)"),
        ("B_cage_rna_floor", "CAGE+RNA+floor"),
        ("D_cage_rna_procap_floor", "CAGE+RNA+PRO-cap+floors")]
fig, axes = plt.subplots(3, len(SHOW), figsize=(16.5, 7.6), sharex=True)
bp = np.arange(624) * BINW + BINW / 2
agx = np.arange(19968)
for ci, (key, lab) in enumerate(SHOW):
    axes[0, ci].plot(bp, Z[key][iD], color="#4477aa", lw=.9)
    axes[1, ci].plot(agx, A1[f"{key}|dnase"][MARG:MARG + 19968, 0], color="#ee6677", lw=.5)
    axes[2, ci].plot(agx, A2[f"{key}|dnase"][MARG:MARG + 19968, 0], color="#228833", lw=.5)
    axes[0, ci].set_title(lab, fontsize=10, fontweight="bold")
    for r in range(3):
        for g, (a, _) in GENES.items():
            axes[r, ci].axvline(a - 128, color="grey", ls=":", lw=.7)
        axes[r, ci].grid(alpha=.18); axes[r, ci].tick_params(labelsize=7)
    axes[2, ci].set_xlabel("Position across construct (bp)", fontsize=9)
for r, nm in enumerate(["Borzoi DNase\n(designed against)",
                        "AlphaGenome 1 DNase\n(independent)",
                        "AlphaGenome 2 DNase\n(independent)"]):
    axes[r, 0].set_ylabel(nm, fontsize=9)
fig.suptitle("Accessibility across the construct in three models. Dotted lines are the three designed TSSs.\n"
             "DNase is in no design objective anywhere in this project, so all three are reporting on "
             "something they were never asked to produce — and they agree.", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(REPO + "figures/fig20_dnase_three_models.png", dpi=140)
print("wrote fig20_dnase_three_models.png")
