# Multi-track objectives, promoter identity, and random genomic context

Fourth round of experiments, prompted by two questions from the lab:

> "curious if designing with multiple tracks such as RNA + many chromatin
> modifications can help?"

and, on flanking context:

> "use some random genomic loci"

Seven new Ledidi arms, one cross-model motif check, and an 11-locus context
sweep. One earlier conclusion is retracted (see
[corrections](08-corrections.md#correction-9)).

Tables: [table9](../data/tables/table9_procapnet_motifs.csv),
[table10](../data/tables/table10_multitrack_arms.csv),
[table11](../data/tables/table11_random_loci.csv).

---

## 1. Chromatin modifications do not substitute for CAGE

Three new arms, all initialised from the CDS-only control with CDSs and flanks
frozen, 600 iterations, AAVS1 context, identical to the earlier ablations.

| Arm | Objective | edits | **CAGE % of target** | RNA % | DNase | H3K4me3 | TBP |
|---|---|---|---|---|---|---|---|
| — | CDS-only control | 0 | 0.1 | 2.4 | 0.07 | 1.16 | 2.02 |
| A | RNA | 555 | 1.3 | 104.9 | 0.15 | 1.74 | 1.66 |
| C | RNA + DNase | 646 | 1.4 | 101.9 | 0.11 | 1.61 | 1.65 |
| **H** | RNA + DNase + H3K4me3 + H3K9ac + TBP + Pol2S5 + Pol2A, scale-normalised | 804 | **1.5** | 101.7 | 0.72 | 3.86 | 4.78 |
| **K** | same, raw equal weights | 1509 | **2.1** | 83.4 | 0.86 | 6.25 | 4.92 |
| B | CAGE + RNA | 689 | **100.0** | 102.8 | 0.80 | 10.67 | 4.30 |

Six chromatin marks move CAGE attainment from 1.3% to 1.5%. Weighting them
harder reaches 2.1% and costs 17 points of RNA accuracy.

The control that makes this interpretable: **the chromatin targets themselves
were met.** H3K4me3 rose 1.74 → 6.25, TBP 1.66 → 4.92, DNase 0.15 → 0.86. This
is not a failure to optimise. Producing the chromatin signature of an active
promoter does not produce initiation, because those marks are consequences of
transcription rather than causes of it.

### Weight normalisation matters, and got this wrong on the first attempt

The first run of arms H and I returned `n_edits: 0` for both. Cause: per-track
weights were divided by each track's squared scale so high-amplitude ChIP tracks
(Pol2A, amplitude ~22) could not swamp RNA (~7), but the result was not rescaled.
Summed objective weight fell from 8.0 to 0.20, a 40× weakening, and the edit
penalty (`l=1e-3`) then blocked every candidate edit. Fixed by rescaling so the
normalised weights sum to the same total as the unnormalised ones:

```python
raw = [weights[t] / (TSCALE[t] ** 2) for t in track_names]
k   = sum(weights[t] for t in track_names) / sum(raw)
wv  = [x * k for x in raw]
```

A zero-edit result from Ledidi should always be read as a weighting bug before it
is read as a biological statement.

## 2. Suppressing H3K4me1 improves strand specificity 3×, and costs a gene

H3K4me1 marks enhancers; H3K4me3 marks promoters. Asking for one down while the
other goes up attacks the promoter-versus-enhancer failure mode directly.

Antisense leakage, measured as summed RNA+ over summed RNA− across the three genes:

| Arm | Objective | RNA+ | RNA− | **RNA+/RNA−** |
|---|---|---|---|---|
| J | CAGE+RNA+chrRNA | 15.66 | 6.85 | **2.3** |
| E | RNA+chrRNA | 16.86 | 6.55 | 2.6 |
| A | RNA | 15.73 | 4.09 | 3.8 |
| B | CAGE+RNA | 15.41 | 3.64 | 4.2 |
| — | **lab's `final`** | 20.00 | 0.335 | **59.7** |
| K | RNA+6 marks, raw | 12.50 | 0.091 | 136.9 |
| **M** | **CAGE+RNA + suppress H3K4me1** | 15.02 | 0.087 | **173.6** |
| **I** | **RNA+6 marks + suppress H3K4me1** | 15.02 | 0.023 | **653.8** |

The term does what it claims. Measured directly against each arm's
unconstrained sibling:

| Comparison | H3K4me1 |
|---|---|
| I vs H | 23.8 → 2.1 (−91%) |
| M vs B | 26.3 → 1.3 (−95%) |
| L vs J | 24.1 → 1.7 (−93%) |

And H3K4me1 is the best single predictor of strand specificity among the marks
scored, over 13 designed arms:

| | r with log₁₀(RNA+/RNA−) |
|---|---|
| **H3K4me1** | **−0.76** |
| H3K4me3 | −0.42 |
| DNase | −0.40 |
| EP300 | −0.38 |

Mechanistically this is the expected result: enhancers initiate divergently and
mRNA promoters unidirectionally, so suppressing the enhancer mark suppresses
antisense transcription.

**The cost is severe.** H3K4me1 and H3K4me3 are correlated in Borzoi's training
data, so the optimiser satisfies "less H3K4me1" partly by flattening H3K4me3 too.
Arm I reached the lowest CAGE of any arm run (0.7%). Arm M kept two genes at
19.0 and 18.7 but drove RFP3 to 0.05.

A divergence penalty applied to a model that measures strand directly
(PRO-cap, base resolution) would not have this entanglement.

**Chromatin-fraction RNA actively causes antisense leakage.** Arms E and J are
the only two with chrRNA in the objective and the only two below 3.0. Chromatin-
associated RNA includes nascent bidirectional transcripts, so requesting it
requests divergent transcription. This is a cost in the recommendation made in
[conclusions](07-conclusions.md) that had not been measured.

## 3. Added constraints make the optimiser abandon a gene

Counting genes whose CAGE falls below 10% of that gene's own target:

| Arm | Objective | genes lost |
|---|---|---|
| B | CAGE+RNA | **0** |
| J | CAGE+RNA+chrRNA | **0** |
| — | lab's `final` | **0** |
| G | CAGE(×0.25)+RNA | 1 (mScarlet) |
| M | CAGE+RNA+anti-H3K4me1 | 1 (RFP3) |
| L | CAGE+RNA+chrRNA+anti-H3K4me1 | 2 |
| E | RNA+chrRNA | 2 |
| A, C, F, H, K, I | no CAGE | 3 |

**Three of twelve objectives deliver all three genes.** The failures are not
under-convergence: arm M made 1,354 edits, nearly twice arm B's 689, and still
zeroed RFP3. The optimiser actively trades one gene for the others.

The cause is that the loss is summed over positions and genes, so abandoning one
gene and overdelivering on the remaining two is a legal way to reduce total MSE.
**A per-gene floor — a penalty that activates when any single gene falls below
some fraction of its own target — would close this.** This is the cheapest
useful change identified in this project.

It also reframes arm G's mScarlet collapse, previously reported as a curiosity in
[ablations](06-objective-ablations.md). It is the generic behaviour of this
objective under added constraint, observed in four independent arms.

## 4. Promoter identity: `final` is the only design that has it

This reverses an earlier claim. H3K4me1 was not in the held-out scoring set until
this round; adding it changes the ranking.

The ratio H3K4me3 / H3K4me1 distinguishes promoter-like from enhancer-like
chromatin character:

| Design | **H3K4me3/H3K4me1** | H3K4me3 | H3K4me1 | CAGE % |
|---|---|---|---|---|
| **lab's `final`** | **3.60** | 27.38 | 7.61 | 55.0 |
| Markov initial (untouched) | 0.91 | 1.10 | 1.21 | 0.1 |
| J CAGE+RNA+chrRNA | 0.89 | 21.40 | 24.11 | 86.7 |
| E RNA+chrRNA | 0.89 | 24.16 | 27.02 | 12.1 |
| M CAGE+RNA+anti-H3K4me1 | 0.78 | 1.02 | 1.31 | 63.0 |
| **B CAGE+RNA** | **0.41** | 10.67 | 26.32 | **100.0** |
| G CAGE(×0.25)+RNA | 0.34 | 10.56 | 31.27 | 53.7 |
| A RNA | 0.10 | 1.74 | 16.69 | 1.3 |

`final` is a 4× outlier. **Every short arm produced here is enhancer-like**, and
the arm with the best CAGE attainment (B, 100% of target) has the second-worst
promoter character of any design.

This is the failure mode this project diagnosed for CAGE-only objectives, now
found in its own CAGE+RNA arms. Scoring designs on the objective they were
optimised for, without checking promoter identity, concealed it.

Consequences:

- **The lab's 2,909 edits are not inefficiency.** That is approximately what
  promoter identity costs. The 689-edit arm looks better on CAGE because it
  skipped it.
- **CAGE + RNA is sufficient for the CAGE number, not for a promoter.** The
  longer, heavier optimisation reached something the short arms did not.
- **Promoter identity is nowhere in the objective.** It emerged. There is
  currently no mechanism guaranteeing it, which suggests adding an explicit
  H3K4me3-up / H3K4me1-down ratio term.
- The anti-H3K4me1 arms suppress H3K4me1 by flattening rather than converting:
  their ratios are 0.78, 0.62 and 0.31, all far below 3.60.

Nothing produced here reaches `final`'s 3.60. Second place in the table is the
untouched Markov backbone.

## 5. Cross-model check: ProCapNet's motifs, independent of Borzoi

Every result above is Borzoi predicting Borzoi, which invites the objection that
Ledidi is exploiting model-specific blind spots rather than designing real
regulatory DNA.

[ProCapNet](https://doi.org/10.1101/2024.05.28.596294) shares no training data,
architecture or authors with Borzoi. ENCODE publishes the initiation motifs it
discovered by TF-MoDISco from K562 PRO-cap (`ENCSR261KBX`): patterns
`ENCFF804ZPG`, annotations `ENCFF608IHS`. Detection thresholds were calibrated
per motif at FPR 1e-3 on 18,168 bp of the lab's own order-5 Markov backbone with
promoter windows removed, so the null is composition-matched. Script:
[`analyze_procapnet_motifs.py`](../code/analyze_procapnet_motifs.py).

| Design | CAGE % of target | initiation-motif instances in 3×600 bp |
|---|---|---|
| CDS-only control | 0.1 | **137** (baseline) |
| F RNA+Pol2S5 | 0.4 | 120 |
| A RNA | 1.3 | 141 |
| C RNA+DNase | 1.4 | 148 |
| E RNA+chrRNA | 12.1 | 189 |
| G CAGE(×0.25)+RNA | 53.7 | 209 |
| lab's `final` | 55.0 | **262** |
| B CAGE+RNA | 100.0 | 183 |

**Spearman +0.850** (n=9). Every CAGE-free arm sits between 120 and 189, against
a 137 baseline; every CAGE-containing arm sits between 183 and 262.

Which motifs, relative to the CDS-only control:

| ProCapNet motif | control | B (CAGE+RNA) | A (RNA only) | `final` |
|---|---|---|---|---|
| DPR_CGG family (11 of 43 patterns, CpG spacing) | 51 | 81 (+30) | 50 (−1) | 113 (+62) |
| NRF1 / NRF1-like | 13 | 16 (+3) | 15 (+2) | 45 (+32) |
| ZBTB33 (Kaiso, methyl-CpG) | 5 | 10 (+5) | 6 (+1) | 15 (+10) |
| YY1 | 3 | 6 (+3) | 1 (−1) | 10 (+7) |
| BRE/SP (SP1 GC-box, ProCapNet's top pattern) | 7 | 9 (+2) | 6 (−1) | 8 (+1) |
| ATF1/CREB | 1 | 1 (0) | 1 (0) | 6 (+5) |

Four independent lines converge:

1. The DPR_CGG family dominates the gain, matching the separately measured rise
   in CpG observed/expected from 0.35–0.55 to 0.82–0.89 with GC% moving only
   2–4 points ([ISM](05-ism-mechanism.md)).
2. BRE/SP gains match the SP1 GC-box found by ISM in RFP3 at −123 and −92,
   installed by a single G→C turning `GGGAGGG` into `GGGCGGG`.
3. YY1 and ATF1/CREB gains match the `ATGGCG` and `GACGTCA` elements found by ISM
   in eGFP's upstream RNA hotspot.
4. NRF1 and ZBTB33 are CpG-island promoter factors, consistent with (1).

**So the CAGE-driven edits install genuine initiation syntax rather than a
Borzoi artifact.** This is the strongest independent evidence in the project and
required no GPU time.

Limits: motif matches are not functional elements; n=9 designs; 24 of 43 patterns
survived the length filter; the separation is not perfect, since arm E at 189
exceeds arm B at 183. It validates that initiation syntax is present, not that it
is present in 28/20/12 proportion.

## 6. Random genomic loci are safe; synthetic random padding is not

[Context robustness](04-context-robustness.md) reported a 437× swing between
composition-matched synthetic random padding and the real AAVS1 locus, and raised
the question of what the lab's pipeline pads with. The answer from the lab is
"random genomic loci" — real sequence, but not the integration site. That is a
third option, and it had not been tested.

Eleven 524,288 bp windows: ten sampled uniformly over autosome length with a
fixed seed (20261027), plus AAVS1 as a positive control. All had 0% N. Context GC
spans 34.42% to 52.77%. Same `final` sequence spliced into each.

| Context | eGFP CAGE | mScarlet CAGE | RFP3 CAGE | CAGE sum |
|---|---|---|---|---|
| composition-matched synthetic random | 0.02 | 4.68 | 0.04 | ~4.7 |
| all-N | 9.66 | 11.94 | 9.45 | 31.1 |
| **random real loci (n=10)** | **median 9.93** | **median 13.88** | **median 10.02** | **26.24–40.72** |
| real AAVS1 integration site | 10.75 | 13.51 | 8.75 | 33.02 |

**Across ten random real loci the CAGE total varies 1.55× (CV 13.1%), against
437× for synthetic padding.** AAVS1 lands at the median (33.02 is exactly the
median of the eleven). Per gene the spread is 1.61× (eGFP), 1.49× (mScarlet),
2.10× (RFP3).

So the 437× effect was never "the integration site versus elsewhere." It was
**real DNA versus synthetic DNA.** Composition-matched random padding is
uniquely pathological: at ~40.7% GC with no genes it reads as a gene desert, and
the model suppresses whatever is placed in it. Literal N, carrying no
information, is far safer. Any real locus is safer still.

**Practical conclusion: the lab's practice is sound, and the earlier warning
should be narrowed.** It applies to synthetic padding, not to random genomic
padding. The remaining caveat is that RNA is more context-sensitive than CAGE
(CV 26.5% vs 13.1%, 2.48× range), so RNA-based comparisons across runs with
different padding deserve care.

The CDS-only control stays at 0.02–0.18 CAGE across all eleven loci, confirming
the signal comes from the design rather than the context.

The GC correlations are not trustworthy at this n. For `final`, r(GC, CAGE) =
+0.391; for arm B it is −0.465, and that sign flip is driven by the two high-GC
loci landing at opposite extremes for arm B (chr10 at 42.50, chr19 at 60.00).
With only two loci above 44% GC, these coefficients should not be interpreted.

## What would resolve the remaining uncertainty

1. **A per-gene floor in the loss.** Cheapest change; prevents gene abandonment.
2. **An explicit H3K4me3/H3K4me1 ratio term.** Makes promoter identity a
   requirement rather than an emergent property of long runs.
3. **Longer runs of arms B and J** to test whether promoter identity is a
   function of optimisation budget or of the objective. If 2,909 edits at 600
   iterations reaches 3.60, the objective is fine and the arms were simply short.
4. **PRO-cap-based design** for strand specificity without the H3K4me1/H3K4me3
   entanglement. See the ProCapNet route in [conclusions](07-conclusions.md).
