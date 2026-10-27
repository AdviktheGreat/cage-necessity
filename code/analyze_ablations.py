
import json, math, sys
import numpy as np
F='/tmp/scores_aavs1.json'
R=json.load(open(F))
G=['eGFP','mScarlet','RFP3']
TC={'eGFP':28.0,'mScarlet':20.0,'RFP3':12.0}; TR_={'eGFP':7.0,'mScarlet':5.0,'RFP3':3.0}
ORDER=['initial','control','final','A_rna','B_cage_rna','G_cage_rna_dw','C_rna_dnase','E_rna_chrrna','F_rna_pol2']
LBL={'initial':'Markov init','control':'CDS-only ctrl','final':'their design',
     'A_rna':'A RNA only','B_cage_rna':'B CAGE+RNA','G_cage_rna_dw':'G CAGE(x0.25)+RNA',
     'C_rna_dnase':'C RNA+DNase','E_rna_chrrna':'E RNA+chromRNA','F_rna_pol2':'F RNA+Ser5Pol2'}
d={r['design']:r for r in R}
have=[k for k in ORDER if k in d]
print("designs scored:", have)
print("\n" + "="*104)
print("CAGE+ promoter peak  (CAGE was in the objective ONLY for arms B and G)")
print("="*104)
print(f"{'design':<20}{'eGFP/28':>12}{'mScarl/20':>12}{'RFP3/12':>12}{'mean %tgt':>12}   objective")
for k in have:
    r=d[k]; v=[r[f'CAGE+|{g}'] for g in G]; pct=np.mean([100*r[f'CAGE+|{g}']/TC[g] for g in G])
    inobj='CAGE in loss' if k in ('B_cage_rna','G_cage_rna_dw') else ('--' if k in ('initial','control') else 'CAGE NOT in loss')
    if k=='final': inobj='CAGE in loss (theirs)'
    print(f"{LBL[k]:<20}{v[0]:12.3f}{v[1]:12.3f}{v[2]:12.3f}{pct:11.1f}%   {inobj}")
print("\n" + "="*104)
print("RNA+ mean over CDS  (RNA in every arm's objective)")
print("="*104)
print(f"{'design':<20}{'eGFP/7':>12}{'mScarl/5':>12}{'RFP3/3':>12}{'mean %tgt':>12}")
for k in have:
    r=d[k]; v=[r[f'RNA+|{g}'] for g in G]; pct=np.mean([100*r[f'RNA+|{g}']/TR_[g] for g in G])
    print(f"{LBL[k]:<20}{v[0]:12.3f}{v[1]:12.3f}{v[2]:12.3f}{pct:11.1f}%")
print("\n" + "="*104)
print("DNase promoter peak  (in the objective ONLY for arm C)")
print("="*104)
print(f"{'design':<20}{'eGFP':>12}{'mScarlet':>12}{'RFP3':>12}{'mean':>12}   objective")
for k in have:
    r=d[k]; v=[r[f'DNase|{g}'] for g in G]
    inobj='DNase in loss' if k=='C_rna_dnase' else '--' if k in ('initial','control') else 'DNase NOT in loss'
    print(f"{LBL[k]:<20}{v[0]:12.3f}{v[1]:12.3f}{v[2]:12.3f}{np.mean(v):12.3f}   {inobj}")
print("\n" + "="*104)
print("HELD-OUT tracks never in ANY objective (promoter-window peak, mean of 3 genes)")
print("="*104)
ho=['CAGE2+','DNase2','DNase3','DNase4','H3K4me3','TBP','Pol2A','Pol2S5']
print(f"{'design':<20}" + "".join(f"{t:>10}" for t in ho))
for k in have:
    r=d[k]; print(f"{LBL[k]:<20}" + "".join(f"{np.mean([r[f'{t}|{g}'] for g in G]):10.3f}" for t in ho))
print("\n" + "="*104)
print("RNA- (antisense) promoter-window peak - strand cleanliness")
print("="*104)
for k in have:
    r=d[k]; print(f"{LBL[k]:<20} RNA- mean over genes = {np.mean([r[f'RNA-|{g}'] for g in G]):8.3f}   RNA+/RNA- = {np.mean([r[f'RNA+|{g}'] for g in G])/max(np.mean([r[f'RNA-|{g}'] for g in G]),1e-9):8.2f}")
def spear(x,y):
    rx=np.argsort(np.argsort(x)); ry=np.argsort(np.argsort(y))
    return float(np.corrcoef(rx,ry)[0,1])
print("\n" + "="*104)
print("CAN RNA-SEQ ALONE RANK DESIGNS THE WAY CAGE DOES?  (across all scored designs x 3 genes)")
print("="*104)
cg=[d[k][f'CAGE+|{g}'] for k in have for g in G]
for t in ['RNA+','DNase','DNase2','H3K4me3','chrRNA+','Pol2S5','TBP']:
    v=[d[k][f'{t}|{g}'] for k in have for g in G]
    print(f"   CAGE+ vs {t:<9} Pearson={np.corrcoef(cg,v)[0,1]:+.3f}   Spearman={spear(cg,v):+.3f}")
