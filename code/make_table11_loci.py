#!/usr/bin/env python3
"""Write table11 from the random-genomic-loci scoring run.

Tests Yash's suggested practice -- pad with "random genomic loci" -- against the
two options tested earlier (composition-matched synthetic random, and the real
AAVS1 integration site). Input is /data/out/scores_loci.json produced by
`modal run borzoi_modal_app.py::main --step scoreloci`.
"""
import json, csv, os, sys, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = sys.argv[1] if len(sys.argv) > 1 else "/tmp/scores_loci.json"
OUT  = os.path.join(HERE, "..", "data", "tables", "table11_random_loci.csv")
GENES = ["eGFP", "mScarlet", "RFP3"]

rows = json.load(open(SRC))
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["design", "locus", "context_gc", "context_cpg_oe",
                "cage_eGFP", "cage_mScarlet", "cage_RFP3", "cage_sum",
                "rna_eGFP", "rna_mScarlet", "rna_RFP3", "rna_sum",
                "dnase_max", "h3k4me3_max", "h3k4me1_max", "is_real_site"])
    for r in sorted(rows, key=lambda r: (r["design"], r["context_gc"] if "context_gc" in r else r["gc"])):
        gc = r.get("gc"); cpg = r.get("cpg_oe")
        w.writerow([r["design"], r["locus"], gc, cpg,
            *[round(r[f"CAGE+|{g}"], 3) for g in GENES],
            round(sum(r[f"CAGE+|{g}"] for g in GENES), 3),
            *[round(r[f"RNA+|{g}"], 3) for g in GENES],
            round(sum(r[f"RNA+|{g}"] for g in GENES), 3),
            round(max(r[f"DNase|{g}"] for g in GENES), 3),
            round(max(r[f"H3K4me3|{g}"] for g in GENES), 3),
            round(max(r[f"H3K4me1|{g}"] for g in GENES), 3),
            int(r["locus"].startswith("chr19_548"))])
print("wrote", OUT)

for design in sorted({r["design"] for r in rows}):
    rs = [r for r in rows if r["design"] == design]
    cs = [sum(r[f"CAGE+|{g}"] for g in GENES) for r in rs]
    print(f"{design:14s} n={len(rs):2d}  CAGE sum {min(cs):7.2f}-{max(cs):7.2f}  "
          f"fold {max(cs)/max(min(cs),1e-9):5.2f}x  CV {st.stdev(cs)/st.mean(cs)*100:5.1f}%")
