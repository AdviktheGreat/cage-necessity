import numpy as np, csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D="/tmp/dualres/"; REPO="/Users/advikanand/Downloads/borzoi-cage-necessity/"
Z=np.load(D+"borzoi_profiles.npz"); A1=np.load(D+"ag_tracks.npz"); A2=np.load(D+"ag2_tracks.npz")
TRK=list(Z["tracks"]); iD=TRK.index("DNase"); MARG=1000; N=19968
GENES={"eGFP":(4462,5182),"mScarlet":(9644,10343),"TagBFP":(14805,15507)}
TSS={g:a-128 for g,(a,z) in GENES.items()}
KEY=["control","final","B_cage_rna","B_cage_rna_floor","PC_rna_procap","D_cage_rna_procap_floor"]
NICE={"control":"CDS-only control","final":"lab's final","B_cage_rna":"CAGE+RNA",
      "B_cage_rna_floor":"CAGE+RNA+floor","PC_rna_procap":"RNA+PRO-cap",
      "D_cage_rna_procap_floor":"CAGE+RNA+PRO-cap\n+floors"}
def get(k,src):
    if src=="bz": return np.repeat(Z[k][iD],32)[:N]
    return (A1 if src=="ag1" else A2)[f"{k}|dnase"][MARG:MARG+N,0]
untouched=np.ones(N,bool)
for g,(a,z) in GENES.items():
    untouched[max(0,TSS[g]-1100):min(N,TSS[g]+1100)]=False; untouched[a:z]=False
R={}
for k in KEY:
    for src in ("bz","ag1","ag2"):
        a=get(k,src); bg=float(np.median(a[untouched]))
        per=[float(a[max(0,TSS[g]-300):TSS[g]+300].max()) for g in GENES]
        R[(k,src)]=(min(per),bg,min(per)/bg,per)

fig,ax=plt.subplots(1,2,figsize=(14.5,6.2))
for j,(src,nm) in enumerate((("ag1","AlphaGenome 1"),("ag2","AlphaGenome 2"))):
    A=ax[j]
    for k in KEY:
        h,bg,sp,_=R[(k,src)]
        c="#cc2222" if k=="final" else ("#228833" if k=="D_cage_rna_procap_floor" else "#4477aa")
        A.scatter(h,sp,s=150,color=c,zorder=3,edgecolor="k",linewidth=.6)
        A.annotate(NICE[k],(h,sp),textcoords="offset points",xytext=(10,6) if k!="D_cage_rna_procap_floor" else (-18,-26),fontsize=8.5)
    A.set_xscale("log")
    A.set_xticks([0.1,0.2,0.5,1,2,5,10])
    A.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    A.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    A.set_xlim(0.09,12)
    A.set_xlabel("DNase at the WORST promoter  (how strong)",fontsize=10)
    A.set_ylabel("promoter DNase ÷ untouched background  (how LOCAL)",fontsize=10)
    A.set_title(f"{nm} (in no design loss)",fontsize=11,fontweight="bold")
    A.grid(alpha=.25)
    A.axhspan(0,10,color="grey",alpha=.08)
    A.text(.02,.03,"below here the peak barely stands out\nfrom the rest of the construct",
           transform=A.transAxes,fontsize=7.5,color="dimgrey")
fig.suptitle("\"DNase peaks at promoters\" has two axes, and no design wins both.\n"
             "The lab's own design is the most LOCAL; the dual-oracle arm is the STRONGEST "
             "but partly by opening the whole locus.",fontsize=12)
fig.tight_layout(rect=[0,0,1,.90])
fig.savefig(REPO+"figures/fig21_peak_vs_specificity.png",dpi=140)
print("wrote fig21_peak_vs_specificity.png")

with open(REPO+"data/tables/table16_promoter_specificity.csv","w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["design","model","dnase_egfp","dnase_mscarlet","dnase_tagbfp",
                "dnase_worst_promoter","background_median_untouched","specificity_ratio"])
    for k in KEY:
        for src,nm in (("bz","borzoi"),("ag1","alphagenome1"),("ag2","alphagenome2")):
            h,bg,sp,per=R[(k,src)]
            w.writerow([k,nm]+[round(v,3) for v in per]+[round(h,3),round(bg,4),round(sp,2)])
print("wrote table16_promoter_specificity.csv")
for k in KEY:
    print(f"  {k:28s} ag1 {R[(k,'ag1')][0]:6.2f} / {R[(k,'ag1')][2]:6.1f}x    "
          f"ag2 {R[(k,'ag2')][0]:6.2f} / {R[(k,'ag2')][2]:6.1f}x")
