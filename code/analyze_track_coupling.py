"""Element-level coupling between CAGE and candidate surrogate assays.

Computes both Pearson and Spearman correlation between CAGE and each other readout
across the nine perturbation arms, per gene. Pearson is reported because it is the
intuitive "do these track each other" number, Spearman because Pearson is inflated
by the dead-versus-working dynamic range of these arms.

Reads  data/tables/table2_promoter_necessity_sufficiency.csv
Writes data/tables/table8_track_coupling.csv
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(HERE, "data/tables/table2_promoter_necessity_sufficiency.csv")
DST = os.path.join(HERE, "data/tables/table8_track_coupling.csv")

TRACKS = [("DNase", "dnase_peak"), ("H3K4me3", "h3k4me3"),
          ("RNA", "rna_mean_cds"), ("Pol2S5", "pol2s5")]
GENES = ["eGFP", "mScarlet", "RFP3"]


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
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


rows = list(csv.DictReader(open(SRC)))
out = []
for gene in GENES:
    sub = [r for r in rows if r["gene"] == gene]
    cage = [float(r["cage_peak"]) for r in sub]
    for label, key in TRACKS:
        vals = [float(r[key]) for r in sub]
        out.append({"gene": gene, "track": label, "n_arms": len(sub),
                    "pearson_vs_cage": round(pearson(cage, vals), 3),
                    "spearman_vs_cage": round(spearman(cage, vals), 3),
                    "min_value": round(min(vals), 3), "max_value": round(max(vals), 3)})

with open(DST, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)

print(f"{'gene':10s} {'track':9s} {'Pearson':>9s} {'Spearman':>9s}   range")
for r in out:
    print(f"{r['gene']:10s} {r['track']:9s} {r['pearson_vs_cage']:+9.3f} "
          f"{r['spearman_vs_cage']:+9.3f}   {r['min_value']}-{r['max_value']}")

print("\nmScarlet false negative in the lab's `final` design:")
final = {r["gene"]: r for r in rows if r["design"] == "final"}
for gene in GENES:
    print(f"  {gene:10s} CAGE {float(final[gene]['cage_peak']):6.2f}   "
          f"DNase {float(final[gene]['dnase_peak']):6.2f}")
print(f"\nwrote {DST}")
