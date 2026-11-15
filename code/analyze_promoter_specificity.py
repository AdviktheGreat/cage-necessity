import numpy as np
D="/tmp/dualres/"
Z=np.load(D+"borzoi_profiles.npz"); A1=np.load(D+"ag_tracks.npz"); A2=np.load(D+"ag2_tracks.npz")
TRK=list(Z["tracks"]); iD=TRK.index("DNase"); MARG=1000; N=19968
GENES={"eGFP":(4462,5182),"mScarlet":(9644,10343),"TagBFP":(14805,15507)}
TSS={g:a-128 for g,(a,z) in GENES.items()}
KEY=["control","final","B_cage_rna","B_cage_rna_floor","PC_rna_procap","D_cage_rna_procap_floor"]
def get(k,src):
    if src=="bz": return np.repeat(Z[k][iD],32)[:N]
    return (A1 if src=="ag1" else A2)[f"{k}|dnase"][MARG:MARG+N,0]

# background = editable, outside CDS, and outside the widest window ANY term acted on
# (ProCapNet input is 2114 bp centred on each TSS => +-1057)
untouched=np.ones(N,bool)
for g,(a,z) in GENES.items():
    untouched[max(0,TSS[g]-1100):min(N,TSS[g]+1100)]=False
    untouched[a:z]=False
print(f"untouched editable background: {untouched.sum()} bp of {N} "
      f"({100*untouched.sum()/N:.0f}%)")
print()
print("Promoter DNase vs background on sequence NO objective term touched")
print(f"{'design':30s} {'model':5s} {'worst_prom':>10s} {'bg_med':>8s} {'specificity':>12s}")
res={}
for k in KEY:
    for src in ("bz","ag1","ag2"):
        a=get(k,src); bg=float(np.median(a[untouched]))
        worst=min(float(a[max(0,TSS[g]-300):TSS[g]+300].max()) for g in GENES)
        res[(k,src)]=(worst,bg,worst/bg if bg>0 else float("nan"))
        print(f"{k if src=='bz' else '':30s} {src:5s} {worst:10.3f} {bg:8.3f} {worst/bg:12.2f}")
print()
print("="*94)
print("Winner by specificity, per model (background = untouched sequence):")
for src in ("bz","ag1","ag2"):
    order=sorted(KEY,key=lambda k:-res[(k,src)][2])
    print(f"  {src:4s}: " + " > ".join(f"{k}({res[(k,src)][2]:.1f})" for k in order[:4]))
print()
print("The two axes do not agree, which is the honest result:")
print(f"{'design':30s} {'AG1 height':>11s} {'AG1 specif':>11s} {'AG2 height':>11s} {'AG2 specif':>11s}")
for k in KEY:
    print(f"{k:30s} {res[(k,'ag1')][0]:11.2f} {res[(k,'ag1')][2]:11.1f} "
          f"{res[(k,'ag2')][0]:11.2f} {res[(k,'ag2')][2]:11.1f}")
np.save("/tmp/dualres/spec_res.npy", np.array([[k,s,*map(str,res[(k,s)])] for k in KEY for s in ("bz","ag1","ag2")],dtype=object), allow_pickle=True)
