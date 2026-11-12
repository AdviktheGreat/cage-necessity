# 13 — AlphaGenome as an arbiter: what survived, what did not

> **Followed up in [docs/14](14-dual-oracle.md).** The oracle overfitting diagnosed in §7 was
> tested directly, by putting a second independent initiation oracle into the design loss. It did
> not help, which turns §7's correlation into a result that reproduces on demand.

[docs/12](12-procap-crossmodel.md) used ProCapNet to judge Borzoi-designed loci and found the
two models disagreed sharply. That left the central question unresolved: was the disagreement a
fact about CAGE and PRO-cap, or an artifact of comparing two unrelated models?

AlphaGenome settles it. It predicts **CAGE, PRO-cap, DNase, ATAC and RNA-seq from one sequence
at 1 bp**, so the assays can be compared with the model held constant. It was used in neither
design loop, so it is an independent judge of every design here.

All numbers below: K562, the same 524,288 bp AAVS1 context used for Borzoi, so the comparison is
sequence-for-sequence. Code in [`code/alphagenome_app.py`](../code/alphagenome_app.py).

| modality | AlphaGenome K562 track | note |
|---|---|---|
| CAGE | rows 150 (+), 423 (−) | |
| PRO-cap | rows 2 (+), 8 (−), **`ENCSR740IPL`** | a *different* experiment from ProCapNet's `ENCSR261KBX` |
| DNase | row 118 | |
| ATAC | row 60 | |
| RNA-seq | rows 119 (+ polyA), 390 (−) | |

## 1. It validates on real promoters

| promoter | PRO-cap sense | antisense/sense | peak offset | mass within ±10 bp |
|---|---|---|---|---|
| GAPDH | 9277.2 | 0.103 | **−1** | 0.847 |
| RPL13A | 5648.2 | 0.012 | **−1** | 0.883 |
| PGK1 | 4556.6 | 0.038 | **−1** | 0.663 |

Sharp, sense-stranded, at the annotated base. Same picture ProCapNet gave, so both PRO-cap models
are calibrated on real sequence. They part company only on designed sequence.

## 2. Resolved: Borzoi's DNase can be trusted

docs/12 §6 ended with a caveat I could not close — maybe Borzoi reported low DNase on the
PRO-cap-designed locus simply because it could not *see* that sequence, the same way it reported
CAGE at 14% of target.

**It is not blind.** Over 21 design-gene pairs:

| | Pearson | Spearman |
|---|---|---|
| **Borzoi DNase vs AlphaGenome DNase** | **+0.872** | **+0.899** |

Two architecturally unrelated models agree closely on accessibility, including on the
PRO-cap-designed locus (Borzoi 2.43 / 0.18 / 0.10; AlphaGenome 4.34 / 0.42 / 0.23 — same
promoter ranked first, same two failures). So the earlier conclusion stands and is now
independently supported: **designing against PRO-cap does not deliver accessible promoters at
all three genes.** See [fig13](../figures/fig13_model_agreement.png), panel 2.

## 3. Retracted: the two assays do not fundamentally disagree

docs/12 §6 claimed CAGE and PRO-cap capture different sequence features, since optimising either
gave ~100% of its own target and ~10% of the other's. **Within AlphaGenome that is false:**

| | Pearson | Spearman |
|---|---|---|
| **AlphaGenome CAGE vs AlphaGenome PRO-cap** | **+0.953** | **+0.918** |

Hold the model fixed and the two assays track each other almost perfectly across the full range,
from the CDS-only control to real housekeeping promoters. The disagreement in docs/12 was
**between Borzoi and ProCapNet**, not between CAGE and PRO-cap. Recorded as
[correction 19](08-corrections.md).

Consistent with this, the weakest agreement in the whole comparison is Borzoi CAGE against
AlphaGenome CAGE — the *same assay*, two models, Spearman only **+0.516** — while the same assay
pair for DNase, which no objective optimised, reaches +0.899. Cross-model agreement degrades
exactly where a model was optimised against.

## 4. Retracted: the PRO-cap design's 3-for-3 was overfitting to its own oracle

The headline of docs/12 §5 was that RNA + PRO-cap produced properly unidirectional promoters at
all three genes. ProCapNet was **in that design loop**. AlphaGenome was not:

| design / promoter | ProCapNet ratio | AlphaGenome ratio | |
|---|---|---|---|
| RNA+PRO-cap / eGFP | 0.049 ✓ | **0.072 ✓** | holds |
| RNA+PRO-cap / mScarlet | 0.070 ✓ | **0.430 ✗** | fails independently |
| RNA+PRO-cap / TagBFP | 0.051 ✓ | **2.170 ✗** | fails independently |

**3 of 3 becomes 1 of 3** once the judge is not the oracle that was optimised. The magnitudes
collapse too: mScarlet 84.8 reads by ProCapNet, 7.4 by AlphaGenome. Recorded as
[correction 20](08-corrections.md).

Across all 21 pairs the two PRO-cap models agree on **17 of 21** divergence verdicts
(Pearson +0.669 on magnitude, +0.546 on the ratio). All four disagreements fall on the two
designs where one model had a stake.

## 5. Softened: the lab's `final` is not divergent at all three promoters

docs/12 §3 reported `final` running 1.7–3.7× more antisense than sense at every promoter, which
I sent to the lab as a headline. AlphaGenome does not reproduce that:

| `final` / promoter | ProCapNet | AlphaGenome | |
|---|---|---|---|
| eGFP | 3.662 ✗ | 1.233 ✗ | divergent in **both** |
| mScarlet | 1.742 ✗ | **0.108 ✓** | disagree |
| TagBFP | 2.173 ✗ | **0.236 ✓** | disagree |

**0 of 3 becomes 2 of 3.** The directionality claim is model-dependent and should not have been
sent as settled. What survives both models is the **eGFP promoter specifically**: divergent by
ProCapNet and by AlphaGenome, and the one promoter where fig12 shows a large antisense peak
270 bp upstream. The three-independent-lines argument that the designs are enhancer-like
(H3K4me1 > H3K4me3, the lab's own NTv3 segmentation, divergent initiation) now rests on two
lines plus a model-dependent third.

## 6. New: designing against PRO-cap places initiation at the intended base

This is the clearest thing PRO-cap buys, and it survives independent judging.
[fig15](../figures/fig15_positional_precision.png):

| | peak offset from designed TSS | mass within ±10 bp |
|---|---|---|
| real promoters | −1 / −1 / −1 | 0.663 – 0.883 |
| **RNA + PRO-cap design** | **−4 / +17 / +0** | **0.152 – 0.430** |
| lab's `final` | +41 / +27 / +46 | 0.032 – 0.272 |
| 15 CAGE-designed promoters | **−172 to +148** | 0.003 – 0.264 |

Every CAGE-designed arm scatters; the PRO-cap arm lands within 17 bp at all three genes and
exactly on the base at TagBFP. Even so, no design reaches real-promoter focus — designed
initiation is spread over 100–200 bp where a real promoter is a single spike
([fig12](../figures/fig12_alphagenome_promoters.png)).

So PRO-cap's value is **position and direction**, not accessibility and not magnitude. That is
also the one job CAGE in Borzoi does worst.

## 7. New, and the most consequential: the expression ladder does not survive an independent model

The deliverable is three genes at 7 / 5 / 3, a 2.33 : 1.67 : 1.00 ratio. Every table in this
repo has reported RNA as a percentage of that target, computed from Borzoi — **the model that
had RNA in its loss.** [fig14](../figures/fig14_rna_ladder_two_models.png):

| design | Borzoi RNA+ (eGFP / mScarlet / TagBFP) | order | AlphaGenome polyA+ | order |
|---|---|---|---|---|
| control | 0.081 / 0.059 / 0.054 | ✓ | 0.029 / 0.004 / 0.002 | ✓ |
| **lab's `final`** | 7.763 / 7.086 / 4.888 | ✓ | **0.050 / 0.270 / 1.820** | **inverted** |
| CAGE+RNA (arm B) | 6.975 / 4.937 / 2.820 | ✓ | 0.757 / 0.357 / 0.057 | ✓ |
| CAGE+RNA+floor | 6.425 / 5.335 / 2.387 | ✓ | 0.312 / 0.885 / 0.346 | ✗ |
| RNA+CpG | 6.921 / 4.899 / 3.158 | ✓ | 3.275 / 0.051 / 0.120 | ✗ |
| arm R | 6.011 / 4.016 / 3.086 | ✓ | 0.016 / 0.007 / 0.083 | ✗ |
| RNA+PRO-cap | 6.763 / 4.809 / 2.812 | ✓ | 0.189 / 0.046 / 0.048 | ✗ |

**Borzoi says 7 of 7 designs have the ladder in the right order. AlphaGenome says 2 of 7.** For
the lab's `final` the order is reversed, with TagBFP at 37× eGFP where it should sit at 0.43×.

Two caveats that keep this from being a verdict:

- The absolute scales differ (Borzoi ~7, AlphaGenome ~0.05–1.8), so only **order and ratio** are
  comparable, not magnitude.
- `final`'s AlphaGenome profile rises monotonically toward the 3′ end and stays high past the
  last gene (pre-eGFP 0.001, eGFP 0.050, gap 0.032, mScarlet 0.270, gap 0.149, TagBFP 1.820,
  **post-TagBFP 0.631**). That is a read-through pattern, not three independent genes. Arm B by
  contrast gives three discrete declining blocks ending at 0.013, which is what an independent
  three-gene locus should look like. The premature `AATAAA` 34 bp after the mScarlet stop
  ([docs/02](02-sequence-analysis.md)) is a candidate cause, and **Borzoi is blind to it**.

Recorded as [correction 21](08-corrections.md). Arm B getting the order right is **not** a claim
that arm B is the better design — its ratio is 13.2 : 6.2 : 1.0, far steeper than the target, and
[correction 9](08-corrections.md) already retracted one "the short arms win" claim.

## 8. What this does to the method

The pattern across §2–§7 is one thing stated four ways: **a readout that was in the loss is not
evidence.** Quantified here:

| readout | in any design loss? | cross-model agreement (Spearman) |
|---|---|---|
| DNase | no | **+0.899** |
| PRO-cap magnitude | only for the PRO-cap arm | +0.648 |
| CAGE | yes, in most arms | **+0.516** |
| RNA ladder order | yes, in every arm | **2 of 7 reproduce** |

Agreement falls monotonically with how hard the objective pushed on that readout. Every headline
in docs/01–12 that rested on an optimised readout needs an independent check before it is
believed, and the two that rested on unoptimised readouts (DNase, chromatin marks) are the two
that survived.

This is the strongest possible argument for the lab's existing external-filtering step, and it
argues for one change: **filter on readouts that were not in the loss**, and treat the loss
readouts as diagnostics only.


## 10. AlphaGenome 2 (released 2026-11-11): the ladder result replicates

AG2 is a conditional model trained on >100 M single-cell profiles across 68 tissues. It is
**not** a drop-in upgrade for this project, for one reason that matters here:

| | AG1 | AG2 |
|---|---|---|
| **PRO-cap** | 12 tracks | **removed** |
| CAGE | 546 | 550 |
| DNase | 305 | **385** |
| ATAC | 167 | 179 |
| scATAC | — | **731 (new)** |
| methylation | — | **new** |
| contact maps | 28 | **73** |
| released weights | no | no |
| input gradients | no | yes, but keyed to an annotated `gene_id` |

So **everything about initiation direction in §3–§5 has to stay on AG1**, since AG2 cannot predict
PRO-cap at all. And the new gradient endpoint does not unlock design: it takes a `gene_id`, not an
arbitrary sequence, so it cannot be pointed at a synthetic construct. AlphaGenome remains a filter,
not a design oracle, in both versions.

What AG2 is good for is the question §7 left open, because expression is exactly where it improved.

**It replicates the inversion, and the two independent models agree on every design**
([fig16](../figures/fig16_ladder_three_models.png)):

| design | Borzoi | AlphaGenome 1 | AlphaGenome 2 |
|---|---|---|---|
| control | ✓ | ✓ | ✓ |
| **lab's `final`** | ✓ | **wrong** | **wrong** |
| CAGE+RNA (arm B) | ✓ | ✓ | ✓ |
| CAGE+RNA+floor | ✓ | wrong | wrong |
| RNA+CpG | ✓ | wrong | wrong |
| arm R | ✓ | wrong | wrong |
| RNA+PRO-cap | ✓ | ✓ | ✓ |
| **order wrong in** | **0 of 7** | **4 of 7** | **4 of 7** |

Two models that share no training pipeline with Borzoi, and only partially with each other, return
**identical verdicts on all seven designs**, and disagree with Borzoi on four. That is much harder
to dismiss than one model's opinion, and it moves [correction 21](08-corrections.md) from
"one model disagrees" to "the only two independent checks available both disagree."

**They differ on severity, which matters for interpretation.** For `final`, AG1 gives
0.081 / 0.340 / 1.849 (a 23× spread, badly inverted) while AG2 gives 0.854 / 1.062 / 1.709 (a 2×
spread, mildly inverted). The target spread is 2.33×. So AG2's reading is that the three genes
land at roughly the right *dynamic range* and in the wrong *order* — a milder failure than AG1
suggested, but a failure of the thing being specified.

One caveat on `PC_rna_procap`: its mScarlet and TagBFP values sit within ~10% of each other, so
whether it counts as correctly ordered flips with the choice of averaging window. Treat it as
"flat", not "correct".

### DNase, now confirmed by a third model

| design | Borzoi | AG1 | AG2 |
|---|---|---|---|
| lab's `final` | 3.70 / 0.37 / 5.65 | 1.30 / 3.56 / 13.38 | 7.64 / 1.40 / 18.05 |
| CAGE+RNA+floor | 1.09 / 4.33 / 1.90 | 2.52 / 8.88 / 4.81 | 7.28 / **43.60** / 8.58 |
| RNA+CpG | 3.05 / 0.13 / 0.12 | 11.50 / 0.96 / 0.45 | **26.17** / 0.43 / 0.55 |
| RNA+PRO-cap | 2.43 / 0.18 / 0.10 | 4.34 / 0.42 / 0.23 | **25.80** / 0.60 / 0.23 |

All three models rank the promoters within each design the same way, including the PRO-cap arm's
one-good-two-dead pattern. Accessibility is now the single most robust readout in this project,
across three unrelated models, and it remains the one no objective ever optimised.

**A data-file bug caught after pushing, recorded here rather than quietly fixed.** The first
version of the AG2 scoring helper selected a track by taking the *first* match on the requested
strand. K562 carries two CAGE tracks per strand in AG2 — `hCAGE` and a low-quality `LQhCAGE` — so
that silently returned `LQhCAGE`, up to 2.5x off (the lab's `final` at eGFP read 21.4 instead of
53.8). The CAGE columns of `data/predictions/alphagenome2_summary.json` were wrong in commit
`9919eb5` and are corrected here; the selector now matches on assay name and raises if the
intended track is absent, rather than falling back.

Nothing above changes. `table14_alphagenome2.csv` was computed separately with explicit column
indices and was always right, and the RNA, DNase, scATAC and methylation values were unaffected
because those modalities have one track per condition. But it is the same species of error as
[corrections 7, 8 and 16](08-corrections.md): a selection step that silently picked the wrong
thing, where a fallback hid the problem instead of failing loudly.

### New: the designed CpG islands are not predicted to be methylated

[docs/11](11-what-cage-is-doing.md) found that CAGE optimisation works by installing CpG islands.
An obvious worry follows: a methylated CpG island is a *silenced* one. AG2 predicts methylation
directly, over the 600 bp upstream of each ATG:

| | eGFP | mScarlet | TagBFP |
|---|---|---|---|
| CDS-only control (7–9 CpGs) | 0.004 | 0.003 | 0.003 |
| lab's `final` (25/20/30 CpGs) | 0.019 | 0.018 | 0.013 |
| RNA+CpG (23/29/27 CpGs) | 0.014 | 0.024 | 0.020 |
| **real GAPDH / RPL13A / PGK1** | **0.006 / 0.009 / 0.015** | | |

Every designed promoter sits in the same low range as real housekeeping CpG-island promoters.
More CpGs does raise predicted methylation slightly (control 0.003 → designs 0.013–0.024), which
is mechanically unsurprising since methylation requires CpGs to occur at, but nothing approaches a
silenced level.

**This is weak evidence and should not be leaned on.** Transgene silencing is driven by
integration site, copy number and the host's silencing machinery — none of which a sequence model
sees. It rules out the narrow worry that the designed sequence is *intrinsically* a methylation
target. It says nothing about what happens to an integrated construct.

## 9. Limits

- One AlphaGenome PRO-cap experiment (`ENCSR740IPL`, 12 tracks across 6 biosamples). Thin, and
  likely overlapping ProCapNet's training data at the assay level even though the accession differs.
- Seven designs, 21 design-gene pairs. The correlations are estimates, not laws — the mistake
  [corrected in docs/12 §4](12-procap-crossmodel.md) was reading one family's correlation as general.
- AlphaGenome is also just a model. Where it and Borzoi disagree, nothing here says which is right;
  the claim is only that **agreement between independent models is weaker evidence than it looked,
  and disagreement is concentrated exactly on the optimised readouts.**
- No wet-lab data. All of this is model-versus-model.
- AlphaGenome 2 cannot predict PRO-cap, so no directionality claim here is checked by more
  than ProCapNet and AG1. AG1 and AG2 also share a lineage, so their agreement on the
  expression ladder is not fully independent evidence.
- The methylation result is a prediction about naked sequence, not about an integrated
  transgene, and should not be read as evidence against silencing.
