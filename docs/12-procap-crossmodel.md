# PRO-cap cross-model validation: the designs initiate divergently
> **Superseded in part by [docs/13](13-alphagenome-arbitration.md).** AlphaGenome predicts
> CAGE, PRO-cap, DNase and RNA from one sequence, and it overturns three claims below:
> §3's verdict on the lab's `final` (corr 20), §5's 3-of-3 result (corr 20), and §6's
> whole premise (corr 19). §1, §2 and §4 stand. Each affected section is flagged inline.
>
> **Followed up in [docs/14](14-dual-oracle.md).** The two-model design loop in §5 was extended
> with CAGE and a per-gene floor on both initiation terms. That fixed its gene-abandonment
> problem and left its overfitting problem exactly where it was.


**The authors' questions.** (1) If we design with CAGE and RNA-seq, do we get
reasonable PRO-cap signal? (2) If we design with RNA-seq and PRO-cap, do we get
good CAGE signal and good DNase peaks at promoters? And the stated primary
criterion: *"the main thing we want to see is DNase peaks at promoters. If the
optimization has converged to a plausible solution, we will see DNase peaks at
promoters."*

Borzoi has **zero** PRO-cap tracks out of 7,611, so answering either question
needs a second model. [ProCapNet](https://www.encodeproject.org/files/ENCFF976FHE/)
is a BPNet trained on ENCODE K562 PRO-cap (ENCSR261KBX): 6.4M parameters, a
2,114 bp input, both strands at **single-base** resolution. It shares no
architecture, training data or authors with Borzoi, and at 2 kb every promoter
window sits entirely inside the 19,968 bp design, so the flank-padding problem that
dominates every Borzoi result here does not arise at all.

Reproduce with `--step procapscore`, `--step procapprof`, `--step procapdesign`.

## 1. The method validates on real promoters

Windows centred on the annotated hg38 TSS of three housekeeping genes:

| | predicted count | reads at TSS (±25 bp) | peak offset | antisense / sense | focus |
|---|---|---|---|---|---|
| GAPDH | 665 | 554 | **+1** | 0.050 | 0.834 |
| RPL13A | 1366 | 1239 | **0** | 0.021 | 0.907 |
| PGK1 | 1205 | 572 | **+7** | 0.160 | 0.474 |
| 8 random genomic windows | 14–138 | **0.0–0.3** | — | 5.6–1003 | 0.000–0.011 |

ProCapNet places the initiation peak within 7 bp of the annotated TSS in all three
cases, and random genomic DNA produces essentially no initiation. **So real
promoters are sharp, sense-stranded, and hold antisense to 2–16% of sense.**

## 2. Q1: a CAGE+RNA design gives weak and inconsistent PRO-cap signal

Reads at the designed TSS, by objective:

| Design / promoter | reads at TSS | antisense / sense | verdict |
|---|---|---|---|
| CDS-only control | 0.3 – 2.3 | 1.7 – 10.4 | background, as expected |
| RNA only | 0.4 – 2.1 | 1.1 – 6.2 | background |
| **CAGE + RNA** (arm B) | **2.9 – 6.0** | 0.5 – 2.9 | barely above background |
| CAGE + RNA + floor | 0.9 – 22.8 | 0.7 | one promoter works |
| **CAGE + RNA + chrRNA** (arm J) | **4.8 – 60.9** | 0.13 – 1.8 | two promoters work |
| **lab's `final`** | **4.1 – 55.8** | 1.7 – 3.7 | one promoter strong |

**Plain CAGE + RNA does not transfer.** Arm B hits 107% of its Borzoi CAGE target
and still produces only 2.9–6.0 reads at the designed TSS, against 0.3–2.3 for a
construct with no promoter at all and 554–1239 for a real promoter. Hitting the
CAGE target is not sufficient for an independent model to see initiation.

The best designs do transfer: arm J reaches 60.9 reads at eGFP and 42.2 at TagBFP,
and the lab's `final` reaches 55.8 at TagBFP. That is **~90× the no-promoter level
and ~9× below a housekeeping promoter**, which for a construct targeting moderate
expression is a plausible place to land.

## 3. The divergence test: 1 of 24 CAGE-designed promoters passes

> **Partly retracted — [correction 20](08-corrections.md).** These verdicts come from
> ProCapNet alone. AlphaGenome passes the lab's `final` at mScarlet (0.108) and TagBFP
> (0.236), so **0 of 3 becomes 2 of 3** for that design. Only the eGFP promoter fails
> under both models. The divergence *method* is sound and the real-promoter range is
> reproduced by both; the per-design verdicts are model-dependent and were sent to the
> authors before a second model had been consulted.


Real promoters hold antisense initiation to 2–16% of sense. Measuring the same
ratio on the designs, with antisense summed over −250 to −25 bp where divergent
initiation actually peaks:

| Design / promoter | sense at TSS | antisense upstream | ratio | |
|---|---|---|---|---|
| arm J / eGFP | 60.1 | **110.1** | **1.83** | divergent |
| lab's `final` / TagBFP | 55.8 | **121.3** | **2.17** | divergent |
| **arm J / TagBFP** | **42.2** | **5.5** | **0.131** | **passes** |
| CAGE+RNA+floor / mScarlet | 21.8 | 16.1 | 0.74 | divergent |
| lab's `final` / mScarlet | 11.1 | 19.4 | 1.74 | divergent |
| lab's `final` / eGFP | 4.1 | 15.1 | 3.66 | divergent |
| real promoter range | — | — | **0.021 – 0.160** | |

**Of 24 promoter-design pairs from CAGE objectives, exactly one has both
real-promoter-like directionality and more than 10 reads of sense initiation: arm
J's TagBFP promoter** (ratio 0.131, focus 0.491, both inside the real range).

Arm J's eGFP promoter has *more antisense than sense initiation*. So does the
lab's `final` at all three of its promoters. These are enhancer signatures,
measured by a model that has never seen a CAGE objective.

**This is the third independent line of evidence for the same conclusion**, and
the first from outside Borzoi: H3K4me1 exceeding H3K4me3
([docs/10](10-multitrack-chromatin.md)), the near-continuous enhancer blanket in
the authors' own NTv3 segmentation on slide 54, and now divergent initiation at
base resolution.
## 4. DNase predicts the independent model, but only inside one family of designs

Over the 24 pairs from CAGE-based objectives:

| Borzoi readout | vs ProCapNet reads at TSS | |
|---|---|---|
| **DNase at the promoter** (never optimised) | Pearson +0.755 | Spearman +0.595 |
| CAGE area (directly optimised) | Pearson +0.365 | Spearman +0.538 |

On that evidence the authors' instinct looks confirmed: the DNase peak is a better
indicator of genuine initiation than the quantity actually being optimised.

**Then the PRO-cap-designed sequence from §5 was added, and it collapses:**

| n = 27, including the PRO-cap-designed sequence | Pearson | Spearman |
|---|---|---|
| DNase vs ProCapNet | **+0.466** (was +0.755) | +0.486 |
| CAGE vs ProCapNet | **+0.040** (was +0.365) | +0.315 |

Three new points did that, because they are extreme: 85.4 PRO-cap reads with DNase
0.18, and 32.4 reads with DNase 0.10.

**So DNase is a useful convergence signal for designs produced by a CAGE
objective, and not a general law about accessibility and initiation.** This is the
same failure mode as corrections 7, 8 and 16 in a new place: a correlation measured
inside one family of sequences, then read as if it held everywhere. The authors'
heuristic is sound for their current pipeline; it should not be trusted once the
objective changes.

## 5. Q2: designing with RNA + PRO-cap gives real promoters, but not CAGE and not DNase

> **Partly retracted — [correction 20](08-corrections.md).** ProCapNet was the model being
> optimised in this loop, so its 3-of-3 verdict is not independent evidence. AlphaGenome,
> which was in no design loop, scores **1 of 3** (eGFP 0.072 passes; mScarlet 0.430 and
> TagBFP 2.170 fail) and drops mScarlet from 84.8 reads to 7.4. **What survives
> independent judging is positional precision**, −4 / +17 / +0 bp from the designed TSSs
> against −172 to +148 bp for CAGE-designed arms ([docs/13 §6](13-alphagenome-arbitration.md)).
> The negative half of this section — no CAGE, no DNase at two of three — is confirmed.


A combined objective with **no CAGE anywhere**: Borzoi supplies RNA abundance over
the 60-bin window per gene, ProCapNet supplies sense initiation at the TSS, an
antisense penalty over -250 to -25 bp, and a magnitude floor. Two models in one
autograd graph, 919 edits, 106 seconds.

| | CAGE % of target | RNA % | DNase min | DNase max | best PRO-cap | passing divergence |
|---|---|---|---|---|---|---|
| no design | 0.3 | 1.3 | 0.07 | 0.07 | 1.5 | 0 / 3 |
| RNA only | 5.0 | 98.3 | 0.11 | 0.15 | 1.3 | 0 / 3 |
| CAGE + RNA | **107.6** | 98.2 | 0.15 | 0.80 | 5.2 | 0 / 3 |
| **RNA + PRO-cap** | **14.2** | 95.9 | 0.10 | **2.43** | **84.8** | **3 / 3** |
| lab's `final` | 112.1 | 131.6 | 0.37 | 5.65 | 55.8 | 0 / 3 |

Per promoter:

| RNA + PRO-cap | DNase | CAGE area (% of target) | PRO-cap reads | antisense / sense |
|---|---|---|---|---|
| eGFP | **2.43** | 7.59 (27%) | **68.9** | **0.049** |
| mScarlet | 0.18 | 0.59 (3%) | **84.8** | **0.070** |
| TagBFP | 0.10 | 0.32 (3%) | **32.2** | **0.051** |

**The answer to Q2 is no on both counts asked, and yes on something better.**

- **CAGE does not follow.** 14.2% of target. And this is not a sharp-peak-versus-
  area artifact: the CAGE *peak* is 2.67 / 0.29 / 0.13 and the *area* is 7.59 /
  0.59 / 0.32, an area/peak ratio of 2.0-2.8, broader than arm B's 1.0-1.1. Borzoi
  simply does not see initiation in this sequence.
- **DNase mostly does not follow.** 2.43 at eGFP is the best any CAGE-free
  objective has produced, three times arm B's 0.80, but the other two promoters sit
  at 0.18 and 0.10.
- **The promoters themselves are the best in this entire project.** 68.9 / 84.8 /
  32.2 reads at all three designed start sites, against 55.8 as the previous best at
  any single promoter, and antisense held to 0.049 / 0.070 / 0.051, all three inside
  the real-promoter range of 0.021-0.160.

**It is the only design where all three promoters pass the divergence test.**

| Design | promoters with >10 sense reads and real-promoter directionality |
|---|---|
| **RNA + PRO-cap** | **3 / 3** |
| CAGE + RNA + chrRNA (arm J) | 1 / 3 |
| lab's `final` | 0 / 3 |
| every other design profiled | 0 / 3 |

## 6. The two models disagree about what a promoter is, and that is the real finding

> **RETRACTED — [correction 19](08-corrections.md).** This section attributes the
> disagreement to the two assays. Hold the model fixed and they agree: AlphaGenome CAGE
> vs its own PRO-cap is Pearson **+0.953**, Spearman +0.918 over the same 21 pairs. The
> disagreement was **Borzoi versus ProCapNet**, a variable I had not controlled. The
> conclusion that external filtering is load-bearing survives and is strengthened; the
> reason given here is wrong. See [docs/13 §3](13-alphagenome-arbitration.md).


ProCapNet sees 84.8 reads of sharp, unidirectional initiation at this design's
mScarlet TSS. Borzoi, looking at the same bases, reports CAGE at 3% of target and
DNase at 0.18. Run it the other way and the disagreement is just as large: arm B
hits 107% of Borzoi's CAGE target and ProCapNet sees 2.9 reads, against 0.3-2.3 for
a construct with no promoter at all.

**CAGE and PRO-cap measure the same physical quantity, capped 5' ends, and
optimising either one yields roughly 100% of its own target and roughly 10% of the
other's.** Ledidi is exploiting whichever model is in the loop.

The practical consequence is the strongest recommendation in this repository:
**the authors' post-design filtering with an external model is not polish, it is
load-bearing.** A design can satisfy one oracle's initiation target almost
perfectly while a second oracle sees nothing. Neither model is obviously the
trustworthy one here; what matters is that agreement between them is the only
evidence that survives.

## Figures

- `figures/fig06_procap_profiles.png` — base-resolution profiles: GAPDH, arm J,
  the lab's `final`, and the no-promoter control, forward above and reverse below
  the axis. The divergence is visible directly.
- `figures/fig07_crossmodel_procap.png` — Borzoi CAGE and Borzoi DNase against
  ProCapNet reads.
- `figures/fig08_dnase_per_promoter.png` — the primary criterion, per promoter,
  across objectives.
- `figures/fig09_procap_vs_real.png` — where the designs sit between background
  and a real promoter.
- `figures/fig10_procap_designed.png` — the PRO-cap-designed locus against the
  lab's CAGE-designed one, all three promoters, with GAPDH for scale.

## Limits

- One PRO-cap model, one cell type, fold 0 of 7.
- The positive controls are housekeeping promoters, among the strongest in the
  genome. A fair magnitude comparison needs real promoters at matched expression,
  which was not done.
- ProCapNet sees a synthetic Markov backbone, which is out of its training
  distribution. The random-genomic negative controls behave sensibly, which is
  reassuring but not conclusive.
- The PRO-cap design is a single run at one weighting. Given how badly weight
  interpolation behaved for the Borzoi objectives
  ([docs/10 §11](10-multitrack-chromatin.md)), it should not be assumed that
  nearby weights give nearby designs.
