"""AlphaGenome scoring of the synthetic locus designs.

Runs on Google's servers via the AlphaGenome API, so no GPU is needed here.
Gives CAGE, PRO-cap, DNase and RNA-seq from one model at 1 bp, which lets us
ask whether the Borzoi/ProCapNet disagreement is real or a cross-model artifact.
"""
import json, os, modal

app = modal.App("alphagenome-locus")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("numpy", "pandas", "alphagenome"))
vol = modal.Volume.from_name("borzoi-data")
SEC = modal.Secret.from_name("alphagenome-key")

AG_LEN = 524288          # AlphaGenome's 500KB window, identical to our Borzoi context
DESIGN_START = 252160    # where the 19,968 bp design sits inside that window
DLEN = 19968
GENES = {"eGFP": (4462, 5182), "mScarlet": (9644, 10343), "RFP3": (14805, 15507)}
K562 = "EFO:0002067"


@app.function(image=image, volumes={"/data": vol}, secrets=[SEC], timeout=1800)
def probe():
    import numpy as np
    from alphagenome.models import dna_client
    client = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
    md = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    out = {}
    for name in ("cage", "procap", "dnase", "rna_seq", "atac"):
        df = getattr(md, name, None)
        if df is None:
            out[name] = "absent"
            continue
        tot = len(df)
        k = df[df["ontology_curie"] == K562] if "ontology_curie" in df.columns else df.iloc[:0]
        rows = []
        for i, r in k.iterrows():
            rows.append({"row": int(i), "name": str(r.get("name", "")),
                         "strand": str(r.get("strand", ".")),
                         "biosample": str(r.get("biosample_name", ""))})
        out[name] = {"total_tracks": tot, "k562_tracks": len(k), "rows": rows[:12],
                     "columns": list(df.columns)[:14]}
    out["supported_lengths"] = sorted(int(x) for x in dna_client.SUPPORTED_SEQUENCE_LENGTHS.values())
    for name in ("cage", "procap", "dnase", "atac", "rna_seq"):
        d = out[name]
        if d == "absent":
            print(f"{name:8s} ABSENT"); continue
        print(f"{name:8s} total={d['total_tracks']:5d}  K562={d['k562_tracks']:3d}")
        for r in d["rows"]:
            print(f"           row {r['row']:5d}  strand {r['strand']:2s}  {r['name']}")
    print("lengths:", out["supported_lengths"])
    return out


# ---- real promoters, as the yardstick (hg38 TSS) ----
REAL_PROM = {"GAPDH": ("chr12", 6534517), "RPL13A": ("chr19", 49487608),
             "PGK1": ("chrX", 78104248)}
OUTS = ("cage", "procap", "dnase", "atac", "rna_seq")


def _summarise(arrs, tss, label, gene):
    """Metrics matched to the readout, as in docs/12: sense at the TSS, antisense
    upstream where divergent initiation actually peaks, accessibility nearby."""
    import numpy as np
    s = {}
    for k in ("cage", "procap"):
        a = arrs[k]
        plus = a[:, 0]
        minus = a[:, 1] if a.shape[1] > 1 else np.zeros_like(plus)
        sense = float(plus[tss - 25:tss + 26].sum())
        anti = float(np.abs(minus[tss - 250:tss - 25]).sum())
        s[f"{k}_sense"] = sense
        s[f"{k}_anti"] = anti
        s[f"{k}_ratio"] = anti / sense if sense > 1e-9 else float("inf")
        s[f"{k}_area160"] = float(plus[tss:tss + 160].sum())
        s[f"{k}_peak"] = float(plus[tss - 25:tss + 26].max())
    for k in ("dnase", "atac"):
        s[f"{k}_max"] = float(arrs[k][tss - 300:tss + 300, 0].max())
        s[f"{k}_mean"] = float(arrs[k][tss - 300:tss + 300, 0].mean())
    s["rna_mean"] = float(arrs["rna_seq"][:, 0][tss:tss + 1920].mean())
    s["design"], s["gene"], s["tss"] = label, gene, int(tss)
    return s


@app.function(image=image, volumes={"/data": vol}, secrets=[SEC], timeout=3600)
def score(designs="control,final,B_cage_rna,B_cage_rna_floor,PC_rna_procap,A_rna_cpg1,R_rna3_antienh_floor"):
    import numpy as np
    from alphagenome.models import dna_client
    from alphagenome.data import genome
    client = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
    req = [getattr(dna_client.OutputType, o.upper()) for o in OUTS]

    base = json.load(open("/data/ctx/designs.json"))
    flank = open("/data/ctx/aavs1_flank.txt").read().strip()
    L, R = flank[:DESIGN_START], flank[DESIGN_START + DLEN:]

    rows, keep = [], {}
    # margin so promoter upstream regions are plotted with context
    lo, hi = DESIGN_START - 1000, DESIGN_START + DLEN + 1000

    for name in designs.split(","):
        if name in base:
            dsn = base[name]
        else:
            dsn = json.load(open(f"/data/out/design_{name}.json"))["design"]
        assert len(dsn) == DLEN, (name, len(dsn))
        seq = L + dsn + R
        assert len(seq) == AG_LEN
        out = client.predict_sequence(sequence=seq, requested_outputs=req,
                                      ontology_terms=[K562])
        arrs = {o: np.asarray(getattr(out, o).values, dtype=np.float32) for o in OUTS}
        for o in OUTS:
            keep[f"{name}|{o}"] = arrs[o][lo:hi]
        for gene, (atg, _) in GENES.items():
            tss = DESIGN_START + atg - 128
            rows.append(_summarise(arrs, tss, name, gene))
        print(f"{name:24s} done   shapes " +
              " ".join(f"{o}{arrs[o].shape[1]}" for o in OUTS), flush=True)

    # real promoters through the same pipeline, via genomic intervals
    for gene, (chrom, tss_g) in REAL_PROM.items():
        iv = genome.Interval(chrom, tss_g - AG_LEN // 2, tss_g + AG_LEN // 2)
        out = client.predict_interval(interval=iv, requested_outputs=req,
                                      ontology_terms=[K562])
        arrs = {o: np.asarray(getattr(out, o).values, dtype=np.float32) for o in OUTS}
        c = AG_LEN // 2
        for o in OUTS:
            keep[f"REAL_{gene}|{o}"] = arrs[o][c - 1000:c + 1000]
        rows.append(_summarise(arrs, c, "REAL", gene))
        print(f"REAL {gene:18s} done", flush=True)

    np.savez_compressed("/data/out/alphagenome_tracks.npz", **keep)
    json.dump(rows, open("/data/out/alphagenome_summary.json", "w"), indent=1)
    vol.commit()
    print(f"\nwrote {len(keep)} arrays, {len(rows)} summary rows")
    return rows
