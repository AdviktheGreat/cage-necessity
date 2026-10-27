
import json, re, numpy as np, statistics as st
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S=json.load(open('/tmp/ism_aavs1_rep0.json'))
GENES={'eGFP':(4462,5182),'mScarlet':(9644,10343),'RFP3':(14805,15507)}
D=json.load(open('/tmp/vol/ctx/designs.json')); fin=D['final']; ctl=D['control']
per={}
for sh in S:
    g=sh['gene']; per.setdefault(g,{'ref':sh['ref'],'rows':[]}); per[g]['rows']+=sh['rows']
RES={}
for g,(a,z) in GENES.items():
    rows=per[g]['rows']; ref=per[g]['ref']; pos=sorted({r[0] for r in rows})
    d={k:{} for k in ('cage','dnase','rna')}
    for p,rb,al,c,dn,rn in rows:
        d['cage'].setdefault(p,[]).append(c-ref['cage']); d['dnase'].setdefault(p,[]).append(dn-ref['dnase']); d['rna'].setdefault(p,[]).append(rn-ref['rna'])
    RES[g]={'pos':np.array(pos),'ref':ref,**{k:np.array([st.mean(d[k][p]) for p in pos]) for k in d}}

IUPAC={'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','N':'[ACGT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]'}
def rx(c): return ''.join(IUPAC[x] for x in c)
def rc(s): return s[::-1].translate(str.maketrans('ACGTRYSWKMNBDHV','TGCAYRSWMKNVHDB'))
MOTIFS={'SP1/KLF GC-box':'GGGGCGGGG','SP1 core':'GGGCGG','NRF1':'YGCGCAYGCGCR','NRF1 half':'GCGCATGC',
 'YY1':'AANATGGCG','NFY CCAAT':'CCAATCA','CREB/ATF':'TGACGTCA','ETS/GABPA':'ACCGGAAGT','ETS core':'CCGGAAGT',
 'E-box':'CACGTG','TATA':'TATAAA','GATA1':'WGATAR','KLF1 erythroid':'CCMCRCCC','NFE2/AP1':'TGCTGAGTCA',
 'AP1':'TGASTCA','RUNX1':'TGTGGTT','CTCF':'CCRCCASYAGRKGG','ZNF143':'TTCCCAGAATG','THAP11':'ACTAYRNNNCCCR'}
print("=== motifs PRESENT in the designed 600bp promoter window but ABSENT in the control (both strands) ===")
for g,(a,z) in GENES.items():
    lo,hi=a-500,a+100
    f=fin[lo:hi]; c=ctl[lo:hi]
    print(f"\n-- {g}")
    found=False
    for nm,cons in MOTIFS.items():
        for strand,pat in (('+',rx(cons)),('-',rx(rc(cons)))):
            nf=len(re.findall(pat,f)); nc=len(re.findall(pat,c))
            if nf>nc:
                hits=[m.start() for m in re.finditer(pat,f)]
                attrib=[]
                for h in hits:
                    gp=lo+h; m=(RES[g]['pos']>=gp)&(RES[g]['pos']<gp+len(cons))
                    if m.any(): attrib.append(float(np.abs(RES[g]['cage'][m]).mean()))
                locs=', '.join(f"{lo+h-a:+d}" for h in hits)
                print(f"   {nm:<16} {strand} new={nf-nc} (final {nf}, ctrl {nc}) at {locs}   mean|dCAGE| at site = {np.mean(attrib) if attrib else float('nan'):.3f}")
                found=True
    if not found: print("   (none newly created)")
    bg=float(np.abs(RES[g]['cage']).mean()); print(f"   [background mean|dCAGE| over window = {bg:.3f}]")

fig,axes=plt.subplots(3,1,figsize=(13,9))
for gi,(g,(a,z)) in enumerate(GENES.items()):
    ax=axes[gi]; x=RES[g]['pos']-a
    ax.plot(x,RES[g]['cage'],lw=0.9,color='#d62728',label=f"CAGE (ref {RES[g]['ref']['cage']:.1f})")
    ax.plot(x,RES[g]['dnase'],lw=0.9,color='#1f77b4',label=f"DNase (ref {RES[g]['ref']['dnase']:.2f})")
    ax.plot(x,RES[g]['rna'],lw=0.9,color='#2ca02c',label=f"RNA (ref {RES[g]['ref']['rna']:.1f})")
    ax.axvline(-128,color='k',ls='--',lw=0.8); ax.axvline(0,color='gray',ls=':',lw=0.8)
    ax.text(-128,ax.get_ylim()[1]*0.82,' target TSS',fontsize=7)
    ax.text(0,ax.get_ylim()[1]*0.82,' ATG',fontsize=7)
    ax.set_title(f"{g}  -  ISM (mean over 3 alt bases); negative = base is activating",fontsize=10)
    ax.set_ylabel('d readout'); ax.legend(fontsize=7,loc='lower left'); ax.grid(alpha=0.25,lw=0.4)
axes[-1].set_xlabel('position relative to ATG (bp)')
fig.suptitle('Borzoi in-silico mutagenesis of the Ledidi-designed promoters, AAVS1 context',fontsize=12)
fig.tight_layout(); fig.savefig('/tmp/borzoi_out/fig_ism.png',dpi=130)
print("\nwrote fig_ism.png")

