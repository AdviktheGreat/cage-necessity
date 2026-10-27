
import json
rows=json.load(open('/tmp/flank_test.json'))
G=['eGFP','mScarlet','RFP3']
def get(fl,dn): 
    for r in rows:
        if r['flank']==fl and r['design']==dn: return r
flanks=[('orig',0),('rand',101),('rand',202),('rand',303),('aavs1',0),('N',0)]
names={'orig':'orig random','rand':'random','aavs1':'AAVS1 real','N':'all-N'}
print("=== CAGE+ promoter peak for the FINAL design, across flanking contexts ===")
print(f"{'context':<16}{'eGFP':>9}{'mScarlet':>10}{'RFP3':>9}   (targets 28 / 20 / 12)")
vals={g:[] for g in G}
for fl,sd in flanks:
    r=get(fl,'final'); lab=names[fl]+(f" s{sd}" if fl=='rand' else "")
    print(f"{lab:<16}"+"".join(f"{r['CAGE+|'+g]:>9.3f}" if g=='eGFP' else f"{r['CAGE+|'+g]:>10.3f}" if g=='mScarlet' else f"{r['CAGE+|'+g]:>9.3f}" for g in G))
    for g in G: vals[g].append(r['CAGE+|'+g])
print("\nfold-range across contexts (max/min):")
for g in G: print(f"   {g:<10} {max(vals[g])/max(min(vals[g]),1e-9):8.1f}x   min={min(vals[g]):.3f} max={max(vals[g]):.3f}")

print("\n=== same, RNA+ mean over CDS (targets 7 / 5 / 3) ===")
print(f"{'context':<16}{'eGFP':>9}{'mScarlet':>10}{'RFP3':>9}")
for fl,sd in flanks:
    r=get(fl,'final'); lab=names[fl]+(f" s{sd}" if fl=='rand' else "")
    print(f"{lab:<16}{r['RNA+|eGFP']:>9.3f}{r['RNA+|mScarlet']:>10.3f}{r['RNA+|RFP3']:>9.3f}")

print("\n=== DNase promoter peak (never specified) ===")
print(f"{'context':<16}{'eGFP':>9}{'mScarlet':>10}{'RFP3':>9}")
for fl,sd in flanks:
    r=get(fl,'final'); lab=names[fl]+(f" s{sd}" if fl=='rand' else "")
    print(f"{lab:<16}{r['DNase|eGFP']:>9.3f}{r['DNase|mScarlet']:>10.3f}{r['DNase|RFP3']:>9.3f}")

print("\n=== design vs control vs initial, IN THE AAVS1 REAL LOCUS ===")
print(f"{'readout':<18}{'initial':>10}{'control':>10}{'final':>10}{'final/ctrl':>12}")
for k in ['CAGE+|mScarlet','CAGE+|RFP3','CAGE+|eGFP','DNase|mScarlet','DNase|RFP3','RNA+|mScarlet','RNA+|RFP3','RNA+|eGFP']:
    i,c,f=get('aavs1','initial')[k],get('aavs1','control')[k],get('aavs1','final')[k]
    print(f"{k:<18}{i:10.4f}{c:10.4f}{f:10.4f}{(f/c if c>1e-9 else float('nan')):11.1f}x")

