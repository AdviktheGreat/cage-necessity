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
      "nucRNA+":6250,"H3K9ac":2199}
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
}
TGT_CAGE = {"eGFP":28.0, "mScarlet":20.0, "RFP3":12.0}
TGT_RNA  = {"eGFP":7.0,  "mScarlet":5.0,  "RFP3":3.0}

def target_profiles(track_names):
    import numpy as np
    prof = np.zeros((len(track_names), NB), dtype="float32")
    for ti,tn in enumerate(track_names):
        for g,(a,z) in GENES.items():
            tb=(a-128)//BINW; bz=(z+1056)//BINW
            if tn in ("CAGE+",):
                prof[ti, tb] = TGT_CAGE[g]
            elif tn in ("DNase","Pol2S5"):
                prof[ti, tb] = TGT_RNA[g]           # modest accessibility target at the TSS
            elif tn in ("RNA+","chrRNA+"):
                prof[ti, tb:bz+1] = TGT_RNA[g]
    return prof

@app.function(image=image, gpu="A100-80GB", volumes={"/data":data,"/cache":hfc}, timeout=10800)
def design(arm="B_cage_rna", flank="aavs1", init="control", max_iter=600,
           l=1e-3, lr=0.5, batch_size=2, rep=0, seed=0):
    import torch, numpy as np
    from ledidi import Ledidi
    weights = OBJ[arm]; track_names = list(weights.keys())
    tids = [TR[t] for t in track_names]
    w = torch.tensor([weights[t] for t in track_names], dtype=torch.float32).cuda()[None,:,None]

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
    with open(f"/data/out/design_{arm}.json","w") as f: json.dump({**res,"design":dsn}, f)
    data.commit()
    return res

HELDOUT = ("CAGE+","CAGE2+","DNase","DNase2","DNase3","DNase4","RNA+","RNA-",
           "chrRNA+","nucRNA+","Pol2S5","Pol2A","H3K4me3","TBP")

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
