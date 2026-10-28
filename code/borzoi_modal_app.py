import modal, os, json, time

app = modal.App("borzoi-ledidi")
data = modal.Volume.from_name("borzoi-data", create_if_missing=True)
hfc  = modal.Volume.from_name("borzoi-hf-cache", create_if_missing=True)
image = (modal.Image.debian_slim(python_version="3.11")
    .pip_install("numpy","torch","transformers","huggingface_hub","borzoi-pytorch==0.5.1","ledidi==2.2.0")
    .env({"HF_HOME":"/cache/hf"}))

BIN0, NB, BINW = 2760, 624, 32
DESIGN_START, DLEN = 252160, 19968
TR = {"CAGE+":306,"CAGE-":307,"DNase":1397,"RNA+":6122,"RNA-":6123,
      "chrRNA+":6254,"chrRNA-":6255,"Pol2S5":2450,"H3K4me3":2199,"TBP":2852,
      "DNase2":1398,"DNase3":1399,"CAGE2+":872,"DNase4":1901,"Pol2A":2477,
      "nucRNA+":6250,"H3K9ac":2200,"H3K4me1":2883,"H3K27me3":2196,
      "H3K36me3":2197,"EP300":2834}
GENES = {"eGFP":(4462,5182),"mScarlet":(9644,10343),"RFP3":(14805,15507)}

def ohe_np(seq):
    import numpy as np
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    x = np.zeros((4,len(arr)), dtype="float32")
    for j,c in enumerate(b"ACGT"): x[j] = (arr==c)
    return x

def _flank(kind):
    if kind=="aavs1": s=open("/data/ctx/aavs1_flank.txt").read().strip()
    elif kind=="orig": s=open("/data/ctx/final_ctx.txt").read().strip()
    elif kind=="N":    s="N"*524288
    elif kind.startswith("loc:"): s=open(f"/data/ctx/loci/{kind[4:]}.txt").read().strip()
    else: raise ValueError(kind)
    return s[:DESIGN_START], s[DESIGN_START+DLEN:]

def build(design, kind):
    L,R=_flank(kind); s=L+design+R; assert len(s)==524288; return s

def load(rep=0):
    from borzoi_pytorch import Borzoi
    return Borzoi.from_pretrained(f"johahi/borzoi-replicate-{rep}").cuda().eval()

def readouts(y, tracks=("CAGE+","CAGE-","DNase","DNase2","RNA+","RNA-","chrRNA+","Pol2S5","H3K4me3","TBP")):
    out={}
    for n in tracks:
        p = y[TR[n], BIN0:BIN0+NB].float().cpu().numpy()
        for g,(a,z) in GENES.items():
            tb=(a-128)//BINW; lo,hi=max(0,tb-4),tb+6
            out[f"{n}|{g}"] = float(p[a//BINW:-(-z//BINW)].mean()) if n in ("RNA+","RNA-","chrRNA+") else float(p[lo:hi].max())
            if n in ("CAGE+","CAGE-"):
                # Slide 50 describes the CAGE objective as an AREA over the 5'UTR plus
                # 32 bp (~160 bp, 5 bins), not a peak height at one bin. Reporting both
                # so "percent of target" can be checked against either reading.
                out[f"{n}area|{g}"] = float(p[tb:tb+5].sum())
        out[f"{n}|_mean"]=float(p.mean()); out[f"{n}|_max"]=float(p.max())
    return out

# ---------------------------------------------------------- all arms x flanks
@app.function(image=image, gpu="A100-80GB", volumes={"/data":data,"/cache":hfc}, timeout=3600)
def arms(flank="aavs1"):
    import torch, numpy as np
    model=load(0); D=json.load(open("/data/ctx/designs.json")); rows=[]
    for nm,dsn in D.items():
        X=torch.from_numpy(ohe_np(build(dsn,flank)))[None].cuda()
        with torch.no_grad(), torch.autocast("cuda",dtype=torch.bfloat16):
            y=model(X)[0].float()
        r=readouts(y); r.update({"design":nm,"flank":flank}); rows.append(r)
        del X,y; torch.cuda.empty_cache()
    json.dump(rows, open(f"/data/out/arms_{flank}.json","w")); data.commit()
    return {"n":len(rows),"flank":flank}

# ---------------------------------------------------------- ISM
@app.function(image=image, gpu="A100-80GB", volumes={"/data":data,"/cache":hfc},
              timeout=3600, max_containers=20)
def ism_shard(spec):
    import torch, numpy as np
    gene, p0, p1, flank, rep = spec["gene"], spec["p0"], spec["p1"], spec["flank"], spec["rep"]
    model=load(rep)
    design=json.load(open("/data/ctx/designs.json"))["final"]
    full=build(design, flank)
    a,z=GENES[gene]; tb=(a-128)//BINW; lo,hi=max(0,tb-4),tb+6
    cb0,cb1 = a//BINW, -(-z//BINW)
    base=ohe_np(full)
    BATCH=4
    def run(mats):
        X=torch.from_numpy(np.stack(mats)).cuda()
        with torch.no_grad(), torch.autocast("cuda",dtype=torch.bfloat16):
            y=model(X).float()
        c=y[:,TR["CAGE+"],BIN0+lo:BIN0+hi].max(dim=1).values
        d=y[:,TR["DNase"],BIN0+lo:BIN0+hi].max(dim=1).values
        r=y[:,TR["RNA+"],BIN0+cb0:BIN0+cb1].mean(dim=1)
        o=torch.stack([c,d,r],dim=1).cpu().numpy()
        del X,y; torch.cuda.empty_cache(); return o
    ref=run([base])[0]
    res=[]; pend=[]; meta=[]
    for pos in range(p0,p1):
        gp=DESIGN_START+pos
        refb=int(base[:,gp].argmax())
        for alt in range(4):
            if alt==refb: continue
            m=base.copy(); m[:,gp]=0.0; m[alt,gp]=1.0
            pend.append(m); meta.append((pos,refb,alt))
            if len(pend)==BATCH:
                for o,mt in zip(run(pend),meta): res.append(list(mt)+[float(o[0]),float(o[1]),float(o[2])])
                pend=[]; meta=[]
    if pend:
        for o,mt in zip(run(pend),meta): res.append(list(mt)+[float(o[0]),float(o[1]),float(o[2])])
    return {"gene":gene,"p0":p0,"p1":p1,"rep":rep,"flank":flank,
            "ref":{"cage":float(ref[0]),"dnase":float(ref[1]),"rna":float(ref[2])},"rows":res}

# ---------------------------------------------------------- Ledidi design
OBJ = {
  "A_rna":        {"CAGE+":0.0, "RNA+":1.0},
  "B_cage_rna":   {"CAGE+":1.0, "RNA+":1.0},
  "C_rna_dnase":  {"RNA+":1.0, "DNase":1.0},
  "E_rna_chrrna": {"RNA+":1.0, "chrRNA+":1.0},
  "G_cage_rna_dw":{"CAGE+":0.25, "RNA+":1.0},
  "F_rna_pol2":   {"RNA+":1.0, "Pol2S5":1.0},
  # Does RNA + many chromatin modifications substitute for CAGE?
  "H_rna_multichrom": {"RNA+":2.0, "DNase":1.0, "H3K4me3":1.0, "H3K9ac":1.0,
                       "TBP":1.0, "Pol2S5":1.0, "Pol2A":1.0},
  # Same, plus suppression of H3K4me1, the enhancer mark.
  "I_rna_multichrom_antienh": {"RNA+":2.0, "DNase":1.0, "H3K4me3":1.0, "H3K9ac":1.0,
                               "TBP":1.0, "Pol2S5":1.0, "Pol2A":1.0, "H3K4me1":2.0},
  # Combination recommended from the first round of ablations but never run.
  "J_cage_rna_chrrna": {"CAGE+":1.0, "RNA+":1.0, "chrRNA+":1.0},
  # Same track set as H but with raw equal weights, so the high-amplitude ChIP
  # marks dominate. Directly comparable to arms A/C/F.
  "K_rna_multichrom_raw": {"RNA+":1.0, "DNase":1.0, "H3K4me3":1.0, "H3K9ac":1.0,
                           "TBP":1.0, "Pol2S5":1.0, "Pol2A":1.0},
  # Arm I showed that suppressing H3K4me1 (the enhancer mark) gives RNA+/RNA- of
  # 654 vs the lab's 59.7, but kills initiation. Arm J showed CAGE+RNA+chrRNA
  # gets initiation and chromatin realism but the worst antisense of any arm.
  # These two combine the halves that worked.
  "L_cage_rna_chrrna_antienh": {"CAGE+":1.0, "RNA+":1.0, "chrRNA+":1.0, "H3K4me1":3.0},
  "M_cage_rna_antienh":        {"CAGE+":1.0, "RNA+":1.0, "H3K4me1":3.0},
  # Tests the recommendation that falls out of correction 9: make promoter
  # identity an explicit requirement (H3K4me3 up, H3K4me1 down) rather than
  # something that emerges from a long enough run.
  "N_cage_rna_promoteridentity": {"CAGE+":1.0, "RNA+":1.0, "H3K4me3":2.0, "H3K4me1":2.0},
}
# Multi-track arms normalise each track by its scale so high-amplitude marks
# (Pol2A ~22) cannot swamp RNA (~7) through units alone.
NORM_ARMS = {"H_rna_multichrom", "I_rna_multichrom_antienh"}
TGT_CAGE = {"eGFP":28.0, "mScarlet":20.0, "RFP3":12.0}
TGT_RNA  = {"eGFP":7.0,  "mScarlet":5.0,  "RFP3":3.0}

# Peak amplitude asked of each point-like chromatin track, set to roughly what
# the lab's own optimised design achieves there, so the request is realistic.
# H3K4me1 is 0.0 because the enhancer mark should be suppressed, not produced.
TGT_AMP = {"DNase":3.2, "H3K4me3":17.0, "H3K9ac":10.0, "TBP":16.7,
           "Pol2S5":16.7, "Pol2A":22.3, "H3K4me1":0.0}
# Relative ladder across genes, matching the 28/20/12 CAGE request.
LADDER  = {"eGFP":1.4, "mScarlet":1.0, "RFP3":0.6}
# Typical scale of each track, used only for weight normalisation.
TSCALE  = {"CAGE+":20.0,"RNA+":5.0,"chrRNA+":5.0,"DNase":3.2,"H3K4me3":17.0,
           "H3K9ac":10.0,"TBP":16.7,"Pol2S5":16.7,"Pol2A":22.3,"H3K4me1":5.0}

def target_profiles(track_names):
    import numpy as np
    prof = np.zeros((len(track_names), NB), dtype="float32")
    for ti,tn in enumerate(track_names):
        for g,(a,z) in GENES.items():
            tb=(a-128)//BINW; bz=(z+1056)//BINW
            if tn in ("CAGE+",):
                prof[ti, tb] = TGT_CAGE[g]
            elif tn in ("RNA+","chrRNA+"):
                prof[ti, tb:bz+1] = TGT_RNA[g]
            elif tn in TGT_AMP:
                prof[ti, tb] = TGT_AMP[tn] * LADDER[g]
    return prof

@app.function(image=image, gpu="A100-80GB", volumes={"/data":data,"/cache":hfc}, timeout=10800)
def design(arm="B_cage_rna", flank="aavs1", init="control", max_iter=600,
           l=1e-3, lr=0.5, batch_size=2, rep=0, seed=0):
    import torch, numpy as np
    from ledidi import Ledidi
    weights = OBJ[arm]; track_names = list(weights.keys())
    tids = [TR[t] for t in track_names]
    wv = [weights[t] for t in track_names]
    if arm in NORM_ARMS:
        # Per-track scale normalisation, then rescaled so the summed objective
        # weight matches an unnormalised arm. Without the rescale the whole
        # objective shrinks by ~1/scale^2 and the edit penalty blocks every edit.
        raw = [weights[t] / (TSCALE[t] ** 2) for t in track_names]
        k = sum(weights[t] for t in track_names) / sum(raw)
        wv = [x * k for x in raw]
    w = torch.tensor(wv, dtype=torch.float32).cuda()[None,:,None]

    base = load(rep)
    for p in base.parameters(): p.requires_grad_(False)

    class Head(torch.nn.Module):
        def __init__(self, m): super().__init__(); self.m=m
        def forward(self, X):
            with torch.autocast("cuda", dtype=torch.bfloat16):
                y = self.m(X)
            return (y[:, tids, BIN0:BIN0+NB].float()) * w.sqrt()

    head = Head(base).cuda().eval()
    init_design = json.load(open("/data/ctx/designs.json"))[init]
    X = torch.from_numpy(ohe_np(build(init_design, flank)))[None].cuda()

    mask = np.ones(524288, dtype=bool)
    cds = set()
    for a,z in GENES.values(): cds.update(range(a,z))
    for p in range(DLEN):
        if p not in cds: mask[DESIGN_START+p] = False
    n_mut = int((~mask).sum())
    mask_t = torch.from_numpy(mask).cuda()

    y_bar = torch.from_numpy(target_profiles(track_names))[None].cuda() * w.sqrt()

    torch.manual_seed(seed)
    led = Ledidi(head, shape=(4,524288), output_loss=torch.nn.MSELoss(),
                 input_mask=mask_t, l=l, lr=lr, batch_size=batch_size,
                 max_iter=max_iter, early_stopping_iter=max_iter,
                 report_iter=50, verbose=True, random_state=seed).cuda()
    t0=time.time()
    out = led.fit_transform(X, y_bar)
    dt = time.time()-t0

    # pick the best of the returned batch by output loss
    out = out.detach()
    with torch.no_grad():
        yy = head(out.cuda())
        err = ((yy - y_bar)**2).mean(dim=(1,2))
        best = int(err.argmin())
    seq_oh = out[best].detach().cpu().numpy()
    idx = seq_oh.argmax(axis=0)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[idx].tobytes().decode()
    dsn = seq[DESIGN_START:DESIGN_START+DLEN]
    n_edits = sum(1 for i in range(DLEN) if dsn[i]!=init_design[i])

    res = {"arm":arm,"flank":flank,"init":init,"seconds":round(dt,1),"n_mutable":n_mut,
           "n_edits":n_edits,"edit_frac":round(n_edits/n_mut,4),"max_iter":max_iter,
           "l":l,"lr":lr,"batch_size":batch_size,"tracks":track_names,
           "final_output_mse":float(err[best])}
    # Longer runs or different edit penalties must not overwrite earlier results.
    label = arm
    if max_iter != 600: label += f"_i{max_iter}"
    if abs(l - 1e-3) > 1e-12: label += f"_l{l:g}"
    res["arm"] = label
    with open(f"/data/out/design_{label}.json","w") as f: json.dump({**res,"design":dsn}, f)
    data.commit()
    return res

HELDOUT = ("CAGE+","CAGE2+","DNase","DNase2","DNase3","DNase4","RNA+","RNA-",
           "chrRNA+","nucRNA+","Pol2S5","Pol2A","H3K4me3","TBP",
           "H3K4me1","H3K9ac","H3K27me3","H3K36me3","EP300")

@app.function(image=image, gpu="A100-80GB", volumes={"/data":data,"/cache":hfc}, timeout=5400)
def score_designs(flank="aavs1"):
    import torch, numpy as np, glob
    model=load(0); rows=[]
    base=json.load(open("/data/ctx/designs.json"))
    cand={k:base[k] for k in ("initial","control","final")}
    for fp in sorted(glob.glob("/data/out/design_*.json")):
        d=json.load(open(fp)); cand[d["arm"]]=d["design"]
    for nm,dsn in cand.items():
        X=torch.from_numpy(ohe_np(build(dsn,flank)))[None].cuda()
        with torch.no_grad(), torch.autocast("cuda",dtype=torch.bfloat16):
            y=model(X)[0].float()
        r=readouts(y, tracks=HELDOUT); r["design"]=nm; rows.append(r)
        del X,y; torch.cuda.empty_cache()
    json.dump(rows, open(f"/data/out/scores_{flank}.json","w")); data.commit()
    return {"n":len(rows),"designs":[r["design"] for r in rows]}


# ------------------------------------------------- random genomic loci as context
# Yash's suggestion is to pad with "random genomic loci". That is a third option,
# distinct from both composition-matched synthetic DNA and the real integration
# site: real sequence, but not the site the construct actually lands in. These
# loci are sampled uniformly over autosome length with a fixed seed.
LOCI_SEED = 20261027
CHROM_LEN = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,
"chr5":181538259,"chr6":170805979,"chr7":159345973,"chr8":145138636,
"chr9":138394717,"chr10":133797422,"chr11":135086622,"chr12":133275309,
"chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,
"chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,
"chr21":46709983,"chr22":50818468}

def _sample_loci(n):
    import random
    rng = random.Random(LOCI_SEED)
    names = list(CHROM_LEN); tot = sum(CHROM_LEN.values()); out = []
    while len(out) < n:
        x = rng.randrange(tot); acc = 0
        for c in names:
            if acc + CHROM_LEN[c] > x: ch, off = c, x - acc; break
            acc += CHROM_LEN[c]
        # keep clear of telomeres/centromere-adjacent assembly gaps
        if off < 5_000_000 or off + 524288 > CHROM_LEN[ch] - 5_000_000: continue
        out.append((ch, off))
    return out

@app.function(image=image, volumes={"/data":data}, timeout=5400)
def fetch_loci(n=10):
    import urllib.request, os
    os.makedirs("/data/ctx/loci", exist_ok=True)
    picks = _sample_loci(n) + [("chr19", 54862840)]   # last one = AAVS1, as a positive control
    meta = []
    for ch, s in picks:
        nm = f"{ch}_{s}"; fp = f"/data/ctx/loci/{nm}.txt"
        if os.path.exists(fp) and len(open(fp).read().strip()) == 524288:
            seq = open(fp).read().strip()
        else:
            url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38"
                   f";chrom={ch};start={s};end={s+524288}")
            seq = None
            for _ in range(4):
                try:
                    with urllib.request.urlopen(url, timeout=240) as r:
                        seq = json.loads(r.read().decode())["dna"].upper(); break
                except Exception as e:
                    print("retry", nm, type(e).__name__); time.sleep(5)
            if seq is None or len(seq) != 524288:
                print("FAILED", nm); continue
            open(fp, "w").write(seq)
        nn = seq.count("N")
        gc = (seq.count("G") + seq.count("C")) / max(len(seq) - nn, 1)
        exp = (seq.count("C") * seq.count("G")) / max(len(seq), 1)
        meta.append({"name":nm,"chrom":ch,"start":s,"gc":round(gc,4),
                     "n_frac":round(nn/len(seq),4),
                     "cpg_oe":round(seq.count("CG")/max(exp,1e-9),3)})
        print(f"  {nm:22s} GC={gc*100:5.2f}%  N={nn/len(seq)*100:5.2f}%")
    json.dump(meta, open("/data/ctx/loci/meta.json","w"), indent=1)
    data.commit()
    return {"n":len(meta),"loci":[m["name"] for m in meta]}

@app.function(image=image, gpu="A100-80GB", volumes={"/data":data,"/cache":hfc}, timeout=5400)
def score_loci(designs="final,B_cage_rna,control"):
    import torch, numpy as np, glob
    model = load(0)
    meta = json.load(open("/data/ctx/loci/meta.json"))
    base = json.load(open("/data/ctx/designs.json"))
    cand = {}
    for k in designs.split(","):
        if k in base: cand[k] = base[k]
    for fp in sorted(glob.glob("/data/out/design_*.json")):
        d = json.load(open(fp))
        if d["arm"] in designs.split(","): cand[d["arm"]] = d["design"]
    rows = []
    for m in meta:
        if m["n_frac"] > 0.01:
            print("skip high-N locus", m["name"]); continue
        for nm, dsn in cand.items():
            X = torch.from_numpy(ohe_np(build(dsn, f"loc:{m['name']}")))[None].cuda()
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
                y = model(X)[0].float()
            r = readouts(y, tracks=("CAGE+","RNA+","RNA-","DNase","H3K4me3","H3K4me1","Pol2A"))
            r.update({"design":nm,"locus":m["name"],"gc":m["gc"],"cpg_oe":m["cpg_oe"]})
            rows.append(r)
            del X, y; torch.cuda.empty_cache()
        print("done", m["name"])
    json.dump(rows, open("/data/out/scores_loci.json","w")); data.commit()
    return {"n":len(rows),"loci":len({r["locus"] for r in rows})}

@app.local_entrypoint()
def main(step: str="arms", flank: str="aavs1", rep: int=0, arm: str="B_cage_rna",
         max_iter: int=600, l: float=1e-3, lr: float=0.5):
    if step=="arms":
        for fl in flank.split(","):
            print(json.dumps(arms.remote(fl)))
    elif step=="ism":
        specs=[]
        for g,(a,z) in GENES.items():
            for p0 in range(a-500, a+100, 100):
                specs.append({"gene":g,"p0":p0,"p1":min(p0+100,a+100),"flank":flank,"rep":rep})
        print(f"launching {len(specs)} ISM shards")
        t0=time.time(); out=list(ism_shard.map(specs))
        json.dump(out, open(f"/tmp/ism_{flank}_rep{rep}.json","w"))
        import shutil
        print(f"done in {time.time()-t0:.0f}s; shards={len(out)}; rows={sum(len(o['rows']) for o in out)}")
        with open(f"/tmp/ism_{flank}_rep{rep}.json") as f: blob=f.read()
        vol_path=f"/out/ism_{flank}_rep{rep}.json"
        with data.batch_upload(force=True) as b: b.put_file(f"/tmp/ism_{flank}_rep{rep}.json", vol_path)
        print("uploaded", vol_path)
    elif step=="design":
        handles=[design.spawn(arm=a, flank=flank, max_iter=max_iter, l=l, lr=lr, rep=rep)
                 for a in arm.split(",")]
        for h in handles:
            try: print(json.dumps(h.get()))
            except Exception as e: print("ARM FAILED:", repr(e)[:300])
    elif step=="score":
        print(json.dumps(score_designs.remote(flank), indent=2))
    elif step=="fetchloci":
        print(json.dumps(fetch_loci.remote(10), indent=2))
    elif step=="scoreloci":
        print(json.dumps(score_loci.remote(arm), indent=2))
