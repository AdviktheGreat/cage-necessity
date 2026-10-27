
import json, statistics as st
d=json.load(open('/tmp/borzoi_out/final/results_0_predictions.json'))
BIN0,NB,BINW=2760,624,32
LOC={'CAGE+':0,'CAGE-':1,'DNase':2,'RNA+':3,'RNA-':4}
genes={'eGFP':(4462,5182),'mScarlet':(9644,10343),'RFP3':(14805,15507)}
tgtC={'eGFP':28,'mScarlet':20,'RFP3':12}; tgtR={'eGFP':7,'mScarlet':5,'RFP3':3}
print("Target model from slides 50-51: TSS = ATG-128 (4 bins); RNA block = TSS -> stop+1056 (33 bins)")
print(f"{'gene':<10}{'TSS bp':>8}{'TSS bin':>9}{'RNA block bins':>17}{'CAGE tgt':>10}{'RNA tgt':>9}")
B={}
for g,(a,z) in genes.items():
    tss=a-128; b0=tss//BINW; b1=(z+1056)//BINW
    B[g]=(b0,b1); print(f"{g:<10}{tss:8d}{b0:9d}{str(b0)+'-'+str(b1):>17}{tgtC[g]:10d}{tgtR[g]:9d}")

print("\n=== per-replicate achievement (is the design tuned to ONE replicate?) ===")
for g,(a,z) in genes.items():
    b0,b1=B[g]
    print(f"\n-- {g}  (CAGE target {tgtC[g]}, RNA target {tgtR[g]})")
    print(f"   {'rep':<6}{'CAGE peak':>11}{'% of tgt':>10}{'at bin':>8}{'RNA mean(CDS)':>15}{'% of tgt':>10}")
    for r in range(4):
        cp=d[r][LOC['CAGE+']][BIN0:BIN0+NB]
        w=cp[max(0,b0-4):b0+6]; pk=max(w); at=cp.index(pk)
        rn=d[r][LOC['RNA+']][BIN0:BIN0+NB][a//BINW:-(-z//BINW)]; rm=st.mean(rn)
        print(f"   {r:<6}{pk:11.3f}{100*pk/tgtC[g]:9.1f}%{at:8d}{rm:15.4f}{100*rm/tgtR[g]:9.1f}%")
    cp=[st.mean(d[r][LOC['CAGE+']][BIN0+i] for r in range(4)) for i in range(NB)]
    w=cp[max(0,b0-4):b0+6]
    print(f"   {'ens':<6}{max(w):11.3f}{100*max(w)/tgtC[g]:9.1f}%{cp.index(max(w)):8d}")

