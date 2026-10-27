
import json, re
D=json.load(open('/tmp/vol/ctx/designs.json'))
GENES={'eGFP':4462,'mScarlet':9644,'RFP3':14805}
print("=== TATA check: canonical core-promoter position is about -35..-20 from the TSS ===")
print("    designed TSS = ATG-128, so canonical TATA would sit near ATG-163..-148\n")
for g,a in GENES.items():
    print(f"-- {g}")
    for nm in ('initial','control','final'):
        w=D[nm][a-500:a+100]
        hits=[(m.start()-500) for m in re.finditer('TATAAA',w)]
        hits_rc=[(m.start()-500) for m in re.finditer('TTTATA',w)]
        can=D[nm][a-175:a-140]
        print(f"   {nm:<8} TATAAA(+) at {hits if hits else 'none'}   TTTATA(-) at {hits_rc if hits_rc else 'none'}   seq[-175..-140]={can}")
print("\n=== eGFP strong RNA-attribution region (-390..-350): what did Ledidi build? ===")
a=GENES['eGFP']
for nm in ('control','final'):
    print(f"   {nm:<8} {D[nm][a-390:a-350]}")
diff=[i-390 for i in range(40) if D['control'][a-390+i]!=D['final'][a-390+i]]
print(f"   edited offsets in that 40bp window: {diff}  ({len(diff)}/40 = {100*len(diff)/40:.0f}%)")
print("\n=== GC and CpG o/e of the ISM windows (ATG-500..ATG+100), control vs final ===")
def stats(s):
    gc=(s.count('G')+s.count('C'))/len(s)
    cg=s.count('CG'); exp=s.count('C')*s.count('G')/len(s)
    return gc, (cg/exp if exp>0 else 0)
for g,a in GENES.items():
    sc=stats(D['control'][a-500:a+100]); sf=stats(D['final'][a-500:a+100])
    print(f"   {g:<10} control GC={100*sc[0]:.1f}% CpG o/e={sc[1]:.2f}   ->   final GC={100*sf[0]:.1f}% CpG o/e={sf[1]:.2f}")
print("\n=== edit density in the ISM window vs whole non-CDS design ===")
tot=sum(1 for i in range(19968) if D['control'][i]!=D['final'][i])
print(f"   total control->final diffs across 19968bp: {tot}")
for g,a in GENES.items():
    n=sum(1 for i in range(a-500,a+100) if D['control'][i]!=D['final'][i])
    print(f"   {g:<10} {n}/600 = {100*n/600:.1f}% edited in ISM window")

