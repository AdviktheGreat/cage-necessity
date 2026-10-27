
import json, math, statistics as st
import numpy as np
S=json.load(open('/tmp/ism_aavs1_rep0.json'))
GENES={'eGFP':(4462,5182),'mScarlet':(9644,10343),'RFP3':(14805,15507)}
import re
designs=json.load(open('/tmp/vol/ctx/designs.json'))
fin=designs['final']; ctl=designs['control']
per={}
for sh in S:
    g=sh['gene']; per.setdefault(g,{'ref':sh['ref'],'rows':[]}); per[g]['rows']+=sh['rows']
print("ISM: per-base, 3 alts, window ATG-500..ATG+100, AAVS1 context, Borzoi rep0\n")
RES={}
for g,(a,z) in GENES.items():
    ref=per[g]['ref']; rows=per[g]['rows']
    pos=sorted({r[0] for r in rows})
    d={k:{} for k in ('cage','dnase','rna')}
    for p,refb,alt,c,dn,rn in rows:
        d['cage'].setdefault(p,[]).append(c-ref['cage'])
        d['dnase'].setdefault(p,[]).append(dn-ref['dnase'])
        d['rna'].setdefault(p,[]).append(rn-ref['rna'])
    prof={k:np.array([st.mean(d[k][p]) for p in pos]) for k in d}
    RES[g]={'pos':np.array(pos),'prof':prof,'ref':ref}
    def gini(x):
        x=np.sort(np.abs(x)); n=len(x); cx=np.cumsum(x)
        return float((n+1-2*np.sum(cx)/cx[-1])/n) if cx[-1]>0 else float('nan')
    def top(x,k=20):
        idx=np.argsort(-np.abs(x))[:k]; return float(np.abs(x[idx]).sum()/np.abs(x).sum())
    print(f"-- {g}   ref: CAGE {ref['cage']:.3f}  DNase {ref['dnase']:.3f}  RNA {ref['rna']:.3f}")
    for k in ('cage','dnase','rna'):
        x=prof[k]
        print(f"     {k:<6} mean|d|={np.abs(x).mean():.4f}  max|d|={np.abs(x).max():.4f}  Gini={gini(x):.3f}  top20/total={100*top(x):.1f}%")
    def pear(x,y):
        if x.std()==0 or y.std()==0: return float('nan')
        return float(np.corrcoef(x,y)[0,1])
    print(f"     r(CAGE,DNase)={pear(prof['cage'],prof['dnase']):+.3f}  r(CAGE,RNA)={pear(prof['cage'],prof['rna']):+.3f}  r(DNase,RNA)={pear(prof['dnase'],prof['rna']):+.3f}")
    # spatial centroid / spread of |attribution|
    for k in ('cage','dnase','rna'):
        w=np.abs(prof[k]); w=w/w.sum()
        cen=float((RES[g]['pos']*w).sum()); sd=float(np.sqrt(((RES[g]['pos']-cen)**2*w).sum()))
        print(f"     {k:<6} centroid={cen-a:+7.1f} bp from ATG   spread(sd)={sd:6.1f} bp")
np.save('/tmp/ism_prof.npy', {g:{'pos':RES[g]['pos'],**{k:RES[g]['prof'][k] for k in RES[g]['prof']}} for g in RES}, allow_pickle=True)

print("\n=== hotspots: top 12 positions by |d CAGE| and the 16bp designed sequence there ===")
MOT={'SP1/GC-box':'GGGGCGGGG','SP1-core':'GGGCGG','NRF1':'GCGCATGCGC','NFY/CCAAT':'CCAAT','CREB':'TGACGTCA',
     'ETS':'GGAAGT','E-box':'CACGTG','TATA':'TATAAA','YY1':'AAGATGGC','NRF1-core':'GCGCA','KLF':'GGGGNGGGG','CTCF':'CCGCGNGGNGGCAG'}
for g,(a,z) in GENES.items():
    x=RES[g]['prof']['cage']; pos=RES[g]['pos']
    idx=np.argsort(-np.abs(x))[:12]
    print(f"\n-- {g}")
    for i in sorted(idx, key=lambda i: pos[i]):
        p=int(pos[i]); s=fin[p-8:p+8]; sc=ctl[p-8:p+8]
        hits=[m for m,c in MOT.items() if re.sub('N','.',c) and re.search(re.sub('N','.',c), s)]
        print(f"   pos {p-a:+6d} from ATG  dCAGE={x[i]:+8.3f}  dDNase={RES[g]['prof']['dnase'][i]:+7.3f}  dRNA={RES[g]['prof']['rna'][i]:+7.3f}  final={s}  ctrl={sc}  {','.join(hits)}")

