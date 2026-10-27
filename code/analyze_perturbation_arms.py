
import json, math
A=json.load(open('/tmp/arms_aavs1.json'))
def g(n): return [r for r in A if r['design']==n][0]
G=['eGFP','mScarlet','RFP3']
print("=== AAVS1 real-locus context: SUFFICIENCY (transplant 600bp promoter into CDS-only control) ===")
for tn in ['CAGE+','DNase','RNA+']:
    print(f"\n-- {tn}")
    print(f"{'gene':<10}{'control':>10}{'+its prom':>11}{'+all proms':>12}{'final':>10}{'recovered':>11}")
    for gg in G:
        k=f'{tn}|{gg}'
        c,s1,sa,fi=g('control')[k],g(f'ctrl_prom_{gg}')[k],g('ctrl_prom_all')[k],g('final')[k]
        rec=(s1-c)/(fi-c)*100 if abs(fi-c)>1e-9 else float('nan')
        print(f"{gg:<10}{c:10.4f}{s1:11.4f}{sa:12.4f}{fi:10.4f}{rec:10.0f}%")
print("\n=== AAVS1: NECESSITY (scramble the same windows inside the final design) ===")
for tn in ['CAGE+','DNase','RNA+']:
    print(f"\n-- {tn}")
    print(f"{'gene':<10}{'final':>10}{'scram all':>11}{'scram mSc':>11}{'control':>10}{'lost':>8}")
    for gg in G:
        k=f'{tn}|{gg}'
        fi,sa,sm,c=g('final')[k],g('final_scram_all')[k],g('final_scram_mScarlet')[k],g('control')[k]
        lost=(fi-sa)/(fi-c)*100 if abs(fi-c)>1e-9 else float('nan')
        print(f"{gg:<10}{fi:10.4f}{sa:11.4f}{sm:11.4f}{c:10.4f}{lost:7.0f}%")
def pear(x,y):
    n=len(x); mx,my=sum(x)/n,sum(y)/n
    sx=math.sqrt(sum((a-mx)**2 for a in x)); sy=math.sqrt(sum((b-my)**2 for b in y))
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/(sx*sy) if sx*sy>0 else float('nan')
order=[r['design'] for r in A]
print("\n=== AAVS1: track coupling across the 9 arms (per gene) ===")
for gg in G:
    c=[r[f'CAGE+|{gg}'] for r in A]; d=[r[f'DNase|{gg}'] for r in A]; n=[r[f'RNA+|{gg}'] for r in A]
    p=[r[f'Pol2S5|{gg}'] for r in A]; h=[r[f'H3K4me3|{gg}'] for r in A]
    print(f"  {gg:<10} r(CAGE,DNase)={pear(c,d):+.3f}  r(CAGE,RNA)={pear(c,n):+.3f}  r(CAGE,Pol2S5)={pear(c,p):+.3f}  r(CAGE,H3K4me3)={pear(c,h):+.3f}")
print("\n=== context comparison, FINAL design, achievement vs target ===")
tg={'eGFP':(28,7),'mScarlet':(20,5),'RFP3':(12,3)}
for fl in ['aavs1','N','orig']:
    R=json.load(open(f'/tmp/arms_{fl}.json')); f=[r for r in R if r['design']=='final'][0]
    print(f"-- flank={fl}")
    for gg in G:
        tc,tr=tg[gg]
        print(f"     {gg:<10} CAGE {f[f'CAGE+|{gg}']:7.3f}/{tc:<3}({100*f[f'CAGE+|{gg}']/tc:5.1f}%)   RNA {f[f'RNA+|{gg}']:6.3f}/{tr:<2}({100*f[f'RNA+|{gg}']/tr:6.1f}%)   DNase {f[f'DNase|{gg}']:6.3f}")

