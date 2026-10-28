#!/usr/bin/env python3
"""Regenerate table10: every designed arm scored on both CAGE readings.

The CAGE column is reported twice on purpose. Slide 50 of the source deck
describes the objective as an area over the 5'UTR plus 32 bp (~160 bp, 5 bins),
but this repository originally scored it as a peak height in a 10-bin window.
The two differ by more than 2x for a broad initiation zone, which is what
correction 11 is about, so both are kept side by side.
"""
import json, glob, os, sys, csv

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEQD = os.path.join(REPO, "data", "sequences")
OUT = os.path.join(REPO, "data", "tables", "table10_multitrack_arms.csv")

GENES = {"eGFP": (4462, 5182), "mScarlet": (9644, 10343), "RFP3": (14805, 15507)}
TGT_CAGE = {"eGFP": 28.0, "mScarlet": 20.0, "RFP3": 12.0}
TGT_RNA = {"eGFP": 7.0, "mScarlet": 5.0, "RFP3": 3.0}
G = ["eGFP", "mScarlet", "RFP3"]

OBJ = {
    "initial": "none (Markov backbone, no CDSs spliced)",
    "control": "none (backbone + final CDSs)",
    "final": "CAGE + RNA (the lab's own run)",
    "A_rna": "RNA",
    "C_rna_dnase": "RNA + DNase",
    "F_rna_pol2": "RNA + Pol2S5",
    "E_rna_chrrna": "RNA + chrRNA",
    "H_rna_multichrom": "RNA + 6 chromatin marks (scale-normalised)",
    "K_rna_multichrom_raw": "RNA + 6 chromatin marks (raw equal weights)",
    "I_rna_multichrom_antienh": "RNA + 6 marks + H3K4me1 suppressed",
    "G_cage_rna_dw": "CAGE (x0.25) + RNA",
    "B_cage_rna": "CAGE + RNA",
    "B_cage_rna_area": "CAGE + RNA, CAGE as a 5-bin area",
    "B_cage_rna_i2400": "CAGE + RNA, 2400 iterations",
    "B_cage_rna_l0.0001": "CAGE + RNA, edit penalty 1e-4",
    "B_cage_rna_l1e-05": "CAGE + RNA, edit penalty 1e-5",
    "J_cage_rna_chrrna": "CAGE + RNA + chrRNA",
    "J_cage_rna_chrrna_i2400": "CAGE + RNA + chrRNA, 2400 iterations",
    "J_cage_rna_chrrna_l1e-05": "CAGE + RNA + chrRNA, edit penalty 1e-5",
    "M_cage_rna_antienh": "CAGE + RNA + H3K4me1 suppressed",
    "L_cage_rna_chrrna_antienh": "CAGE + RNA + chrRNA + H3K4me1 suppressed",
    "N_cage_rna_promoteridentity": "CAGE + RNA + H3K4me3 up + H3K4me1 down",
}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def load_sequences():
    """Every arm's sequence, keyed the same way the score file keys them."""
    seqs = {}
    three = os.path.join(SEQD, "three_arms.fasta")
    if os.path.exists(three):
        m = read_fasta(three)
        for k, v in (("markov_initial", "initial"), ("cds_only_control", "control"),
                     ("ledidi_final", "final")):
            if k in m:
                seqs[v] = m[k]
    for fp in glob.glob(os.path.join(SEQD, "designed_arms", "*.fasta")):
        for k, v in read_fasta(fp).items():
            seqs[k] = v
    return seqs


def n_edits(seq, ref):
    """Edits are substitutions against the CDS-only control, which is what every
    designed arm was initialised from."""
    if seq is None or ref is None or len(seq) != len(ref):
        return ""
    return sum(1 for a, b in zip(seq, ref) if a != b)


def main():
    sp = os.path.join(REPO, "data", "predictions", "scores_aavs1.json")
    scored = {r["design"]: r for r in json.load(open(sp))}
    seqs = load_sequences()
    ref = seqs.get("control")

    has_area = any("CAGE+area|eGFP" in r for r in scored.values())
    if not has_area:
        sys.exit("scores_aavs1.json predates the area readout; rerun --step score")

    cols = (["design", "objective", "n_edits",
             "cage_pct_of_target_AREA", "cage_pct_of_target_PEAK", "rna_pct_of_target"]
            + [f"cage_area_{g}" for g in G] + [f"cage_peak_{g}" for g in G]
            + ["initiation_breadth", "genes_lost", "h3k4me3_max", "h3k4me1_max",
               "me3_over_me1", "ep300_max", "dnase_max", "pol2a_max", "tbp_max",
               "strand_ratio"])

    rows = []
    for d, o in OBJ.items():
        if d not in scored:
            continue
        r = scored[d]
        area = [r[f"CAGE+area|{g}"] for g in G]
        peak = [r[f"CAGE+|{g}"] for g in G]
        mx = lambda t: max(r[f"{t}|{g}"] for g in G)
        rp = sum(r[f"RNA+|{g}"] for g in G)
        rm = sum(r[f"RNA-|{g}"] for g in G)
        # Breadth near 1.0 means the whole CAGE signal sits in one 32 bp bin.
        breadth = sum(a / max(p, 1e-9) for a, p in zip(area, peak)) / 3
        lost = sum(1 for g in G if r[f"CAGE+|{g}"] < 0.10 * TGT_CAGE[g])
        rows.append({
            "design": d, "objective": o,
            "n_edits": n_edits(seqs.get(d), ref) if d != "final" else n_edits(seqs.get("final"), ref),
            "cage_pct_of_target_AREA": round(100 * sum(area) / 60.0, 1),
            "cage_pct_of_target_PEAK": round(100 * sum(peak) / 60.0, 1),
            "rna_pct_of_target": round(100 * rp / 15.0, 1),
            **{f"cage_area_{g}": round(a, 2) for g, a in zip(G, area)},
            **{f"cage_peak_{g}": round(p, 2) for g, p in zip(G, peak)},
            "initiation_breadth": round(breadth, 2),
            "genes_lost": lost,
            "h3k4me3_max": round(mx("H3K4me3"), 2),
            "h3k4me1_max": round(mx("H3K4me1"), 2) if "H3K4me1|eGFP" in r else "",
            "me3_over_me1": round(mx("H3K4me3") / max(mx("H3K4me1"), 1e-9), 2) if "H3K4me1|eGFP" in r else "",
            "ep300_max": round(mx("EP300"), 2) if "EP300|eGFP" in r else "",
            "dnase_max": round(mx("DNase"), 2),
            "pol2a_max": round(mx("Pol2A"), 2) if "Pol2A|eGFP" in r else "",
            "tbp_max": round(mx("TBP"), 2),
            "strand_ratio": round(rp / max(rm, 1e-9), 1),
        })

    rows.sort(key=lambda x: x["cage_pct_of_target_AREA"])
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}: {len(rows)} designs")
    miss = [r["design"] for r in rows if r["n_edits"] == ""]
    if miss:
        print("no sequence on disk, edit count blank:", ", ".join(miss))


if __name__ == "__main__":
    main()
