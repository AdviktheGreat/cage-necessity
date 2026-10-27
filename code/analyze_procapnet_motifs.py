#!/usr/bin/env python3
"""Cross-model check: do the Ledidi-designed promoters contain the transcription
initiation motifs that ProCapNet independently discovered from K562 PRO-cap data?

ProCapNet (Cochran et al.) shares no training data, architecture or author set with
Borzoi, so agreement here cannot be a Borzoi-specific artifact. Motifs come from
ENCODE ENCFF804ZPG (TF-MoDISco patterns for ENCSR261KBX), annotated via the
accompanying report ENCFF608IHS.

Detection thresholds are calibrated per motif on the lab's own order-5 Markov
backbone with the promoter windows removed, so the null is composition-matched
rather than an arbitrary score cutoff.

Usage:
  # one-time: fetch and unpack the ProCapNet motifs
  curl -sLO https://www.encodeproject.org/files/ENCFF804ZPG/@@download/ENCFF804ZPG.tar.gz
  tar xzf ENCFF804ZPG.tar.gz -C <MOTIF_DIR>
  python3 analyze_procapnet_motifs.py <MOTIF_DIR>/tfmodisco.motifs.raw/profile/meme_files
"""
import glob, os, sys, csv
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

HERE = os.path.dirname(os.path.abspath(__file__))
SEQ  = os.path.join(HERE, "..", "data", "sequences")
OUT  = os.path.join(HERE, "..", "data", "tables", "table9_procapnet_motifs.csv")
MD   = sys.argv[1] if len(sys.argv) > 1 else "/tmp/pcn/memes/tfmodisco.motifs.raw/profile/meme_files"

GENES = {"eGFP": 4462, "mScarlet": 9644, "RFP3": 14805}   # ATG, 0-based
PROMOTER = 600                                            # bp upstream scanned
# Composition of the lab's design, used as PWM background.
BG  = np.array([0.2973, 0.1977, 0.2091, 0.2959])
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}

# Human annotations from the TF-MoDISco report (ENCFF608IHS), positive profile patterns.
NAMES = {1:"BRE/SP", 2:"CA-Inr", 3:"ETS", 4:"NFY", 5:"NRF1", 6:"ATF1", 7:"TATA",
    8:"THAP11", 9:"YY1", 10:"AP1", 14:"TA-Inr", 16:"CTCF", 18:"NRF1-like",
    20:"ZBTB33", 22:"TCT", 23:"BRE/SP_TE", 24:"TATATA", 25:"CA-Inr_TE",
    28:"TATA_TE", 29:"CA-Inr_dimer", 32:"ZBTB33_TATA_TE", 33:"TE", 34:"GC-rich",
    35:"CA-Inr_TE", 36:"YY1-like", 37:"NFY_TE", 38:"ZBTB33_TE", 39:"ETS_TE",
    40:"ETS_dimer", 41:"NFY_C-rich", 42:"ETS_CA-Inr_dimer", 43:"YY1-like"}
for i in (11, 12, 13, 15, 17, 19, 21, 26, 27, 30, 31):
    NAMES[i] = "DPR_CGG"          # CpG-spacing patterns; 11 of 43 positive motifs

# CAGE attainment per design, from table2 (percent of the 28/20/12 request).
CAGE_PCT = {"markov_initial":0.1, "cds_only_control":0.1, "A_rna":1.3,
    "C_rna_dnase":1.4, "F_rna_pol2":0.4, "E_rna_chrrna":12.1,
    "G_cage_rna_dw":53.7, "B_cage_rna":100.0, "ledidi_final":55.0}


def read_fasta(path):
    seqs, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]; seqs[k] = []
        elif k:
            seqs[k].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def load_designs():
    cand = dict(read_fasta(os.path.join(SEQ, "three_arms.fasta")))
    for p in sorted(glob.glob(os.path.join(SEQ, "designed_arms", "*.fasta"))):
        name = os.path.basename(p).replace(".fasta", "")
        for v in read_fasta(p).values():
            cand[name] = v
    return {k: v for k, v in cand.items() if len(v) == 19968}


def parse_meme(path):
    rows, started = [], False
    for line in open(path):
        if line.startswith("letter-probability"):
            started = True; continue
        if started:
            t = line.split()
            if len(t) == 4:
                try: rows.append([float(x) for x in t])
                except ValueError: break
            elif rows: break
    return np.array(rows)


def info_content(m):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 + np.nansum(np.where(m > 0, m * np.log2(m), 0), axis=1)


def load_motifs():
    """TF-MoDISco emits fixed-width windows; trim to the informative core."""
    mot = {}
    for i in range(1, 44):
        fp = os.path.join(MD, f"ENCSR261KBX.profile.pos.motif_{i}.ppm.fwd.meme")
        if not os.path.exists(fp): continue
        m = parse_meme(fp)
        if not m.size: continue
        keep = np.where(info_content(m) > 0.25)[0]
        if len(keep): m = m[keep[0]:keep[-1] + 1]
        if 6 <= len(m) <= 25:
            mot[i] = m
    return mot


def pwm(m):      return np.log2(np.maximum(m, 1e-4) / BG)
def ohe(s):
    X = np.zeros((len(s), 4), dtype=np.float32)
    for j, b in enumerate(s):
        if b in IDX: X[j, IDX[b]] = 1
    return X
def scan(X, P):
    w = P.shape[0]
    if X.shape[0] < w: return np.array([])
    W = sliding_window_view(X, (w, 4))[:, 0]
    return np.einsum("nwk,wk->n", W, P.astype(np.float32))
def both_strands(X, P):
    return np.concatenate([scan(X, P), scan(X, P[::-1, ::-1])])


def spearman(a, b):
    def pearson(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        xc, yc = x - x.mean(), y - y.mean()
        return float((xc * yc).sum() / np.sqrt((xc ** 2).sum() * (yc ** 2).sum()))
    return pearson(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))


def main():
    cand = load_designs()
    mot  = load_motifs()
    if not mot:
        sys.exit(f"no motifs found under {MD}; see docstring for the download step")

    # Calibrate each motif on the Markov backbone minus the promoter windows.
    null = cand["markov_initial"]
    mask = np.ones(len(null), bool)
    for a in GENES.values():
        mask[a - PROMOTER:a] = False
    Xnull = ohe("".join(c for c, m in zip(null, mask) if m))
    thr = {i: float(np.quantile(both_strands(Xnull, pwm(m)), 0.999)) for i, m in mot.items()}
    print(f"calibrated {len(mot)} motifs at FPR 1e-3 on {Xnull.shape[0]} bp of Markov null")

    hits = {}
    for name, seq in cand.items():
        for g, a in GENES.items():
            X = ohe(seq[a - PROMOTER:a])
            for i, m in mot.items():
                hits[(name, g, i)] = int((both_strands(X, pwm(m)) >= thr[i]).sum())

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design", "cage_pct_of_target", "motif_id", "motif_name",
                    "eGFP", "mScarlet", "RFP3", "total"])
        for name in sorted(cand):
            for i in sorted(mot):
                per = [hits[(name, g, i)] for g in ("eGFP", "mScarlet", "RFP3")]
                w.writerow([name, CAGE_PCT.get(name, ""), i, NAMES.get(i, "unknown"),
                            *per, sum(per)])
    print(f"wrote {OUT}")

    order  = [k for k in CAGE_PCT if k in cand]
    totals = {k: sum(hits[(k, g, i)] for g in GENES for i in mot) for k in order}
    print(f"\n{'design':20s} {'CAGE %tgt':>10s} {'motifs':>7s}")
    for k in sorted(order, key=lambda k: CAGE_PCT[k]):
        print(f"{k:20s} {CAGE_PCT[k]:10.1f} {totals[k]:7d}")
    xs = [CAGE_PCT[k] for k in order]
    ys = [totals[k] for k in order]
    print(f"\nSpearman(CAGE attainment, initiation-motif count) = {spearman(xs, ys):+.3f}  (n={len(xs)})")
    cage_free = [k for k in order if "cage" not in k and k != "ledidi_final"]
    cage_in   = [k for k in order if "cage" in k or k == "ledidi_final"]
    print(f"  CAGE-free arms: {[totals[k] for k in cage_free]}")
    print(f"  CAGE-containing arms: {[totals[k] for k in cage_in]}")


if __name__ == "__main__":
    main()
