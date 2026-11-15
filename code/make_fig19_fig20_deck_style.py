"""Figures 19 and 20: the source deck's own plot grammar.

The deck plots signal across the whole construct with the requested pattern overlaid,
split by which readout the objective actually asked for. Every other figure in this
repo reduces a profile to scalars, and that reduction is what hid correction 22.

Inputs are committed: data/predictions/borzoi_profiles.npz (whole-construct Borzoi
profiles, regenerate with `modal run code/borzoi_modal_app.py --step profiles`) and
data/predictions/dnase_traces.npz (AlphaGenome 1 and 2 DNase at 1 bp).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRED = os.path.join(REPO, "data", "predictions")
Z = np.load(os.path.join(PRED, "borzoi_profiles.npz"))
T = np.load(os.path.join(PRED, "dnase_traces.npz"))
MARGIN = int(T["margin"][0])
N = int(T["design_len"][0])

TRK = [str(t) for t in Z["tracks"]]
BINW = 32
iC, iR, iD = TRK.index("CAGE+"), TRK.index("RNA+"), TRK.index("DNase")
GENES = {"eGFP": (4462, 5182), "mScarlet": (9644, 10343), "TagBFP": (14805, 15507)}

ROWS = [("control", "CDS-only control\n(no edits)", {"CAGE": 0, "RNA": 0, "DNase": 0}),
        ("final", "lab's final\nCAGE + RNA", {"CAGE": 1, "RNA": 1, "DNase": 0}),
        ("B_cage_rna", "CAGE + RNA\n(my reproduction)", {"CAGE": 1, "RNA": 1, "DNase": 0}),
        ("PC_rna_procap", "RNA + PRO-cap\n(no CAGE)", {"CAGE": 0, "RNA": 1, "DNase": 0}),
        ("D_cage_rna_procap_floor", "CAGE + RNA + PRO-cap\n+ per-gene floors",
         {"CAGE": 1, "RNA": 1, "DNase": 0})]
COLS = [("CAGE", iC, Z["target_cage"], "CAGE (initiation)"),
        ("RNA", iR, Z["target_rna"], "RNA-seq (abundance)"),
        ("DNase", iD, None, "DNase (accessibility)")]

fig, axes = plt.subplots(len(ROWS), 3, figsize=(15.5, 13.2), sharex=True)
x = np.arange(Z["control"].shape[1])
for ri, (key, lab, inloss) in enumerate(ROWS):
    prof = Z[key]
    for ci, (cname, ti, tgt, ctitle) in enumerate(COLS):
        ax = axes[ri, ci]
        if tgt is not None:
            ax.fill_between(x, 0, tgt, step="mid", color="#ff9900", alpha=.30,
                            label="requested pattern", zorder=1)
            ax.step(x, tgt, where="mid", color="#ff7700", lw=1.0, zorder=2)
        used = inloss[cname]
        ax.plot(x, prof[ti], color="#1f4e9c" if used else "#cc2222", lw=.85,
                label="Borzoi prediction", zorder=3)
        for a, _ in GENES.values():
            ax.axvline((a - 128) // BINW, color="grey", ls=":", lw=.6, zorder=0)
        ax.text(.985, .93, "IN THE LOSS" if used else "HELD OUT", transform=ax.transAxes,
                ha="right", va="top", fontsize=7.5, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=.3", fc="#1f4e9c" if used else "#cc2222", ec="none"))
        if ri == 0:
            ax.set_title(ctitle, fontsize=11, fontweight="bold")
        if ci == 0:
            ax.set_ylabel(lab, fontsize=9)
        if ri == len(ROWS) - 1:
            ax.set_xlabel("Prediction bin across construct", fontsize=9)
        ax.grid(alpha=.18)
        ax.tick_params(labelsize=7)
axes[0, 0].legend(fontsize=6.5, loc="upper left")
fig.suptitle("Borzoi prediction against the requested pattern, across the construct\n"
             "Blue = this readout was in the design objective.  Red = it was not.  "
             "DNase was never specified in any objective in this project.",
             fontsize=12.5, y=.995)
fig.tight_layout(rect=[0, 0, 1, .965])
out = os.path.join(REPO, "figures", "fig19_designed_vs_heldout.png")
fig.savefig(out, dpi=140)
print("wrote", out)

SHOW = [("control", "CDS-only control"), ("final", "lab's final (CAGE+RNA)"),
        ("B_cage_rna_floor", "CAGE+RNA+floor"),
        ("D_cage_rna_procap_floor", "CAGE+RNA+PRO-cap+floors")]
fig, axes = plt.subplots(3, len(SHOW), figsize=(16.5, 7.6), sharex=True)
bp = np.arange(Z["control"].shape[1]) * BINW + BINW / 2
agx = np.arange(N)
for ci, (key, lab) in enumerate(SHOW):
    axes[0, ci].plot(bp, Z[key][iD], color="#4477aa", lw=.9)
    axes[1, ci].plot(agx, T[f"ag1|{key}"][MARGIN:MARGIN + N], color="#ee6677", lw=.5)
    axes[2, ci].plot(agx, T[f"ag2|{key}"][MARGIN:MARGIN + N], color="#228833", lw=.5)
    axes[0, ci].set_title(lab, fontsize=10, fontweight="bold")
    for r in range(3):
        for a, _ in GENES.values():
            axes[r, ci].axvline(a - 128, color="grey", ls=":", lw=.7)
        axes[r, ci].grid(alpha=.18)
        axes[r, ci].tick_params(labelsize=7)
    axes[2, ci].set_xlabel("Position across construct (bp)", fontsize=9)
for r, nm in enumerate(["Borzoi DNase\n(designed against)",
                        "AlphaGenome 1 DNase\n(independent)",
                        "AlphaGenome 2 DNase\n(independent)"]):
    axes[r, 0].set_ylabel(nm, fontsize=9)
fig.suptitle("Accessibility across the construct in three models. Dotted lines are the three "
             "designed TSSs.\nDNase is in no design objective anywhere in this project, so all "
             "three are reporting on something they were never asked to produce — and they agree.",
             fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, .93])
out = os.path.join(REPO, "figures", "fig20_dnase_three_models.png")
fig.savefig(out, dpi=140)
print("wrote", out)
