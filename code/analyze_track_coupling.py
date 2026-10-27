"""Can a chromatin assay stand in for CAGE as a validation readout?

Answers this in the three regimes where the question can be asked, because the
answer differs sharply between them and conflating them was an error in an
earlier version of this analysis.

  regime 1  within one design family: `final` and eight perturbations of it
  regime 2  across designs produced by different objectives
  regime 3  across designs whose loss contained CAGE (the production case)

Writes data/tables/table8_track_coupling.csv
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T2 = os.path.join(HERE, "data/tables/table2_promoter_necessity_sufficiency.csv")
T3 = os.path.join(HERE, "data/tables/table3_objective_ablations.csv")
DST = os.path.join(HERE, "data/tables/table8_track_coupling.csv")
GENES = ["eGFP", "mScarlet", "RFP3"]


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


rows2 = list(csv.DictReader(open(T2)))
rows3 = list(csv.DictReader(open(T3)))
out = []

# regime 1: perturbations of a single design
for gene in GENES:
    sub = [r for r in rows2 if r["gene"] == gene]
    cage = [float(r["cage_peak"]) for r in sub]
    for label, key in [("DNase", "dnase_peak"), ("H3K4me3", "h3k4me3"),
                       ("RNA", "rna_mean_cds"), ("Pol2S5", "pol2s5")]:
        v = [float(r[key]) for r in sub]
        out.append({"regime": "1_within_design_family", "scope": gene, "track": label,
                    "n": len(sub), "pearson_vs_cage": round(pearson(cage, v), 3),
                    "spearman_vs_cage": round(spearman(cage, v), 3),
                    "min_value": round(min(v), 3), "max_value": round(max(v), 3)})

# regimes 2 and 3: across designs from different objectives
for regime, keep in [("2_across_all_objectives", lambda r: True),
                     ("3_cage_containing_only", lambda r: r["cage_in_loss"] == "True")]:
    cage, vals = [], {"DNase": [], "RNA": []}
    for r in rows3:
        if not keep(r):
            continue
        for gene in GENES:
            cage.append(float(r[f"{gene}_cage"]))
            vals["DNase"].append(float(r[f"{gene}_dnase"]))
            vals["RNA"].append(float(r[f"{gene}_rna"]))
    for label, v in vals.items():
        out.append({"regime": regime, "scope": "all_genes", "track": label,
                    "n": len(cage), "pearson_vs_cage": round(pearson(cage, v), 3),
                    "spearman_vs_cage": round(spearman(cage, v), 3),
                    "min_value": round(min(v), 3), "max_value": round(max(v), 3)})

with open(DST, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)

for regime in ["1_within_design_family", "2_across_all_objectives", "3_cage_containing_only"]:
    print(f"\n--- {regime} ---")
    print(f"{'scope':10s} {'track':9s} {'Pearson':>9s} {'Spearman':>9s}   range")
    for r in [x for x in out if x["regime"] == regime]:
        print(f"{r['scope']:10s} {r['track']:9s} {r['pearson_vs_cage']:+9.3f} "
              f"{r['spearman_vs_cage']:+9.3f}   {r['min_value']}-{r['max_value']}")

# An ATAC-only screen, scored against CAGE as ground truth.
print("\n--- would an ATAC-only screen agree with CAGE? (CAGE>5 works, DNase>1.0 works) ---")
miscalls = total = 0
for r in rows3:
    if r["design"] in ("initial", "control"):
        continue
    for gene in GENES:
        c, d = float(r[f"{gene}_cage"]), float(r[f"{gene}_dnase"])
        total += 1
        agree = (c > 5) == (d > 1.0)
        miscalls += not agree
        if not agree:
            print(f"  MISCALL  {r['design']:14s} {gene:9s} "
                  f"CAGE {c:6.2f} {'works' if c > 5 else 'dead ':5s}  "
                  f"DNase {d:6.3f} {'works' if d > 1.0 else 'dead ':5s}")
print(f"\n  {miscalls}/{total} gene-design combinations disagree")
print(f"\nwrote {DST}")
