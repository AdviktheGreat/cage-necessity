
import re, random, json
def load(fn): return re.search(r"Outputs:\s*\['([ACGTN]+)'\]", open(fn).read()).group(1)
D='/Users/advikanand/Downloads/'
ini=load(D+'output_X_0_initial.txt'); ctl=load(D+'output_X_0_initial_with_cds.txt'); fin=load(D+'output_X_0_final.txt')
genes={'eGFP':4462,'mScarlet':9644,'RFP3':14805}
W=600  # promoter window: ATG-600 .. ATG (contains the designed TSS at ATG-108..-117 and the 128bp 5'UTR)
rng=random.Random(7)
def transplant(base, which):
    s=list(base)
    for g in which:
        a=genes[g]-W; b=genes[g]
        s[a:b]=list(fin[a:b])
    return ''.join(s)
def scramble(base, which):
    s=list(base)
    for g in which:
        a=genes[g]-W; b=genes[g]
        w=list(base[a:b]); rng.shuffle(w); s[a:b]=w
    return ''.join(s)
arms={
 'ctrl_prom_eGFP':    transplant(ctl,['eGFP']),
 'ctrl_prom_mScarlet':transplant(ctl,['mScarlet']),
 'ctrl_prom_RFP3':    transplant(ctl,['RFP3']),
 'ctrl_prom_all':     transplant(ctl,['eGFP','mScarlet','RFP3']),
 'final_scram_all':   scramble(fin,['eGFP','mScarlet','RFP3']),
 'final_scram_mScarlet': scramble(fin,['mScarlet']),
}
# sanity
for n,s in arms.items():
    assert len(s)==19968
    for g,atg in genes.items():
        assert s[atg:atg+3]=='ATG', (n,g)   # CDS starts untouched
print("built arms:", list(arms.keys()))
print("transplant/scramble window: ATG-600..ATG (600bp) per gene; CDSs untouched")
for n,s in arms.items():
    d_ctl=sum(1 for i in range(19968) if s[i]!=ctl[i]); d_fin=sum(1 for i in range(19968) if s[i]!=fin[i])
    print(f"  {n:<22} diffs vs control {d_ctl:5d} | vs final {d_fin:5d}")

# embed in the SAME flanks as before
meta=json.load(open('/tmp/borzoi_in/meta.json')); start=meta['design_start']; CTX=meta['CTX']
ref=open('/tmp/borzoi_in/control_ctx.txt').read()
left=ref[:start]; right=ref[start+19968:]
for n,s in arms.items():
    full=left+s+right
    assert len(full)==CTX
    open(f'/tmp/borzoi_in/{n}_ctx.txt','w').write(full)
print("wrote contexts for", len(arms), "arms (identical flanks to the first run)")

