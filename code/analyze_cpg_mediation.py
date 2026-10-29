#!/usr/bin/env python3
"""Mechanism behind the CAGE-to-accessibility link: CpG density mediation.

Reads only committed artifacts (the designed FASTAs and scores_aavs1.json) and
needs no GPU. Reproduces the correlations and the mediation result in docs/11.
"""
import json, glob, os, math, re

REPO = "/Users/advikanand/Downloads/borzoi-cage-necessity"
SEQD = os.path.join(REPO, "data", "sequences")
GENES = {"eGFP": (4462, 5182), "mScarlet": (9644, 10343), "RFP3": (14805, 15507)}
G = list(GENES)

def read_fasta(p):
    s, n, b = {}, None, []
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            if n: s[n] = "".join(b)
            n, b = line[1:].split()[0], []
        elif line: b.append(line)
    if n: s[n] = "".join(b)
    return s

seqs = {}
m = read_fasta(os.path.join(SEQD, "three_arms.fasta"))
for k, v in (("markov_initial","initial"),("cds_only_control","control"),("ledidi_final","final")):
    if k in m: seqs[v] = m[k]
for fp in glob.glob(os.path.join(SEQD, "designed_arms", "*.fasta")):
    seqs.update(read_fasta(fp))

scored = {r["design"]: r for r in json.load(open(os.path.join(REPO,"data","predictions","scores_aavs1.json")))}

def cpg_oe(s):
    c, g = s.count("C"), s.count("G")
    cg = s.count("CG")
    if c == 0 or g == 0: return 0.0
    return cg * len(s) / (c * g)

def gc(s):
    return 100.0 * (s.count("G") + s.count("C")) / max(len(s), 1)

# Motifs found independently by the ISM scan and the ProCapNet pattern set.
MOT = {
 "SP1/KLF GC-box": ["GGGCGG","CCGCCC","GGGGCGGGG","CCCCGCCCC"],
 "NRF1":           ["GCGCATGCGC","CGCATGCG","GCGCCTGCGC"],
 "YY1":            ["AAATGGCG","ATGGCGGC","GCCATNTT".replace("N","")],
 "CREB/ATF":       ["TGACGTCA","TGACGCA","GACGTCA"],
 "NFY/CCAAT":      ["CCAAT","ATTGG"],
 "ETS":            ["CCGGAAGT","ACTTCCGG","CGGAAG"],
 "TATA":           ["TATAAA","TTTATA"],
 "CpG dinucleotide": ["CG"],
}
def count(s, pats):
    return sum(len(re.findall("(?=" + p + ")", s)) for p in pats if p)

rows = []
for dn, s in seqs.items():
    if dn not in scored or len(s) < 16000: continue
    r = scored[dn]
    for gn, (a, z) in GENES.items():
        win = s[a-600:a]          # the 600 bp promoter window used throughout
        if f"DNase|{gn}" not in r: continue
        rows.append(dict(design=dn, gene=gn, win=win,
                         dnase=r[f"DNase|{gn}"],
                         cage=r.get(f"CAGE+area|{gn}", r.get(f"CAGE+|{gn}", 0.0)),
                         me3=r[f"H3K4me3|{gn}"],
                         gc=gc(win), cpg=cpg_oe(win),
                         **{k: count(win, v) for k, v in MOT.items()}))

print(f"{len(rows)} promoter windows from {len(set(x['design'] for x in rows))} designs\n")

def pear(a, b):
    n=len(a); ma=sum(a)/n; mb=sum(b)/n
    ca=[x-ma for x in a]; cb=[y-mb for y in b]
    d=math.sqrt(sum(x*x for x in ca)*sum(y*y for y in cb))
    return sum(x*y for x,y in zip(ca,cb))/d if d else 0.0
def spear(a,b):
    def rk(v):
        o=sorted(range(len(v)), key=lambda i: v[i]); p=[0]*len(v)
        for i,j in enumerate(o): p[j]=i
        return p
    return pear(rk(a), rk(b))

Y = [x["dnase"] for x in rows]
print("WHAT PREDICTS A DNase PEAK AT A PROMOTER?")
print(f"{'sequence feature':22s} {'Pearson':>8s} {'Spearman':>9s}")
print("-"*42)
feats = ["cage","me3","gc","cpg"] + list(MOT)
res=[]
for f in feats:
    X=[x[f] for x in rows]
    res.append((abs(spear(X,Y)), f, pear(X,Y), spear(X,Y)))
for _, f, p, s in sorted(res, reverse=True):
    print(f"{f:22s} {p:+8.3f} {s:+9.3f}")
print("\n(cage and me3 are predictions, not sequence; shown for comparison)")

def pear(a,b):
    n=len(a); ma=sum(a)/n; mb=sum(b)/n
    ca=[x-ma for x in a]; cb=[y-mb for y in b]
    d=math.sqrt(sum(x*x for x in ca)*sum(y*y for y in cb))
    return sum(x*y for x,y in zip(ca,cb))/d if d else 0.0
def partial(x,y,z):
    """correlation of x,y with the influence of z removed"""
    rxy,rxz,ryz = pear(x,y), pear(x,z), pear(y,z)
    den = math.sqrt((1-rxz**2)*(1-ryz**2))
    return (rxy-rxz*ryz)/den if den>1e-12 else float('nan')

CAGE_FREE = {"A_rna","C_rna_dnase","E_rna_chrrna","F_rna_pol2","H_rna_multichrom",
             "I_rna_multichrom_antienh","K_rna_multichrom_raw"}
ctrl = {x["gene"]: x for x in rows if x["design"]=="control"}

print("=== STEP 1: does putting CAGE in the objective add CpG to promoters? ===")
print("CpG dinucleotides in the 600 bp upstream window, change vs the unedited control\n")
print(f"{'objective':34s} {'eGFP':>7s} {'mScar':>7s} {'TagBFP':>7s} {'mean':>7s}")
print("-"*68)
def deltas(dn):
    out=[]
    for g in G:
        hit=[x for x in rows if x["design"]==dn and x["gene"]==g]
        out.append(hit[0]["CpG dinucleotide"]-ctrl[g]["CpG dinucleotide"] if hit else float('nan'))
    return out
groups={"CAGE-free":[], "CAGE in objective":[]}
for dn in sorted(set(x["design"] for x in rows)):
    if dn in ("control","initial"): continue
    dl=deltas(dn); mu=sum(dl)/3
    cls="CAGE-free" if dn in CAGE_FREE else "CAGE in objective"
    groups[cls].append(mu)
    if dn in CAGE_FREE or dn in ("final","B_cage_rna","R_rna3_antienh_floor","B_cage_rna_floor"):
        print(f"{'['+('no CAGE' if dn in CAGE_FREE else 'CAGE')+'] '+dn:34s} {dl[0]:+7.0f} {dl[1]:+7.0f} {dl[2]:+7.0f} {mu:+7.1f}")
for k,v in groups.items():
    print(f"\n  {k:20s} n={len(v):2d}  mean CpG gain {sum(v)/len(v):+6.1f}  range {min(v):+.0f} to {max(v):+.0f}")

print("\n=== STEP 2: is CpG the mediator, or does CAGE act independently? ===")
Y=[x["dnase"] for x in rows]; C=[x["cage"] for x in rows]; P=[x["CpG dinucleotide"] for x in rows]
print(f"  r(CAGE, DNase)                      = {pear(C,Y):+.3f}")
print(f"  r(CpG,  DNase)                      = {pear(P,Y):+.3f}")
print(f"  r(CAGE, CpG)                        = {pear(C,P):+.3f}")
print(f"  r(CAGE, DNase | CpG held fixed)     = {partial(C,Y,P):+.3f}   <- CAGE's own contribution")
print(f"  r(CpG,  DNase | CAGE held fixed)    = {partial(P,Y,C):+.3f}   <- CpG's own contribution")

print("\n=== STEP 3: the control case. Did ASKING for DNase produce DNase? ===")
for dn,lbl in (("control","no design at all"),("A_rna","RNA only"),
               ("C_rna_dnase","RNA + DNase in the objective"),
               ("H_rna_multichrom","RNA + 6 marks incl. DNase"),
               ("B_cage_rna","CAGE + RNA"),("final","lab's final")):
    hit=[x for x in rows if x["design"]==dn]
    if not hit: continue
    dnv=[x["dnase"] for x in hit]; cp=[x["CpG dinucleotide"] for x in hit]
    print(f"  {lbl:30s} DNase {min(dnv):5.2f}-{max(dnv):5.2f}   CpG {sum(cp)/3:5.0f}")
