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

## 7. Initiation directionality: the discriminator CAGE does constrain

Sections 1–6 used H3K4me3/H3K4me1 as the promoter-versus-enhancer test. That was
a proxy, chosen because it was available. The published discriminator is the
**orientation index** — the fraction of initiation on the stronger strand at a
given site. Mikhaylichenko et al. 2018 (*Genes Dev*) put the boundaries at **≥0.8
unidirectional**, characteristic of mRNA promoters, and **≤0.6 divergent**,
characteristic of enhancers, and note that bidirectional promoters often act as
strong enhancers while unidirectional ones generally cannot.

This repository had never measured it, because the reverse-strand CAGE track was
missing from the scored set. Adding it (track 307, K562 CAGE−) gives:

| Objective | eGFP | mScarlet | RFP3 | mean OI | verdict | CAGE % |
|---|---|---|---|---|---|---|
| starting point | 0.671 | 0.518 | 0.544 | **0.578** | divergent | 0.3 |
| F RNA+Pol2S5 | 0.650 | 0.527 | 0.624 | 0.600 | divergent | 1.1 |
| A RNA only | 0.718 | 0.580 | 0.599 | 0.632 | intermediate | 5.0 |
| E RNA+chrRNA | 0.753 | 0.586 | 0.658 | 0.666 | intermediate | 42.3 |
| C RNA+DNase | 0.766 | 0.634 | 0.650 | 0.683 | intermediate | 5.3 |
| K RNA+6 marks raw | 0.798 | 0.575 | 0.725 | 0.699 | intermediate | 4.5 |
| H RNA+6 marks | 0.820 | 0.640 | 0.768 | 0.743 | intermediate | 3.9 |
| I RNA+6 marks+me1↓ | 0.849 | 0.662 | 0.769 | **0.760** | intermediate | 2.4 |
| B CAGE+RNA | 0.963 | **0.566** | 0.984 | **0.838** | unidirectional | 107.6 |
| J CAGE+RNA+chrRNA | 0.900 | 0.754 | 0.959 | 0.871 | unidirectional | 110.9 |
| N +me3↑/me1↓ | 0.922 | 0.911 | 0.828 | 0.887 | unidirectional | 75.5 |
| L J+me1↓ | 0.988 | 0.805 | 0.878 | 0.890 | unidirectional | 46.6 |
| G CAGE(×0.25)+RNA | 0.962 | 0.801 | 0.986 | 0.916 | unidirectional | 61.7 |
| M CAGE+RNA+me1↓ | 0.996 | 0.989 | 0.815 | 0.933 | unidirectional | 70.1 |
| J, `l`=1e-5 | 0.972 | 0.926 | 0.954 | 0.951 | unidirectional | 113.4 |
| B, `l`=1e-5 | 0.983 | 0.963 | 0.986 | 0.977 | unidirectional | 104.0 |
| **lab's `final`** | **0.994** | **0.993** | **0.991** | **0.993** | unidirectional | 112.1 |

**The split is clean and has no overlap.** Every objective without CAGE tops out at
0.760. Every objective with CAGE starts at 0.838. Eight CAGE-free objectives, nine
CAGE-containing, no exceptions in either direction.

So CAGE does constrain a genuine promoter property, and it is this one. Without
it the designed sites initiate in both directions, which is enhancer behaviour.

### Three consequences

**1. Arm B has a defect nothing else detected.** Its mScarlet promoter scores
**0.566**, below the divergent threshold. Per-strand: forward 20.77, reverse
**15.91**. The CAGE target was met at 104% of its own value, the arm holds the
intended ladder in 11 of 11 genomic contexts, and the middle gene is still firing
almost equally in both directions. In a cell that is an antisense transcript at
comparable strength to the sense one. A single summed-over-genes CAGE term does not
notice.

**2. Directionality and chromatin character are separate axes.** Over 17 designs
the orientation index and H3K4me3/H3K4me1 correlate at only Spearman **+0.505**.
Arm B at `l`=1e-5 scores 0.977 on directionality, the second best in the set, and
0.21 on chromatin, the worst in the set. Scoring one does not tell you the other.

**3. The lab's `final` is the only design good on both.** 0.993 and 3.60. Its
antisense CAGE is 0.11–0.19 against 14.71–31.57 sense, so forward:reverse ratios of
112 to 178. On directionality the short arms are closer to it than the chromatin
proxy suggested, but nothing matches it on both axes at once.

### What this corrects

[Correction 9](08-corrections.md) concluded, from H3K4me3/H3K4me1 alone, that every
arm here is enhancer-like. On the orientation index that is **too harsh**: the
CAGE-containing arms are unidirectional by the published criterion. The accurate
statement is that they are promoter-like in initiation directionality and
enhancer-like in chromatin character, and the two properties come apart. `final`
remains the only design that satisfies both, and the reason it does is still
unexplained.

## 8. What the objective actually specifies, window by window

Prompted by a direct question from the authors: does this reimplementation enforce
5′ and 3′ UTR lengths, since the original designs set UTR coordinates in order to
specify RNA level across 5′UTR + CDS + 3′UTR?

Audit of `target_profiles()`:

| Component | This work | Source deck | Match? |
|---|---|---|---|
| TSS position | ATG − 128 | slide 51: 128 bp 5′UTR | **yes** |
| RNA window | TSS → CDS end + 1056, 60 bins ≈ 1,920 bp | slide 50: ~1,900 bp, ~60 bins; slide 51: 1,056 bp 3′UTR | **yes** |
| **CAGE window** | **1 bin (32 bp) at the TSS** | **slide 50: area over 5′UTR + 32 bp, 5 bins** | **no** |

So the RNA half of the objective does enforce the UTR geometry and matches the
deck. The CAGE half does not: it asks for a single-bin spike where the deck
describes an area. That is the mechanical cause of the breadth gap in
[§7](#7-initiation-directionality-the-discriminator-cage-does-constrain) — the short
arms here score 1.04–1.07 on area ÷ peak, meaning all initiation in one 32 bp bin,
against 2.14 for the lab's `final`.

A reporting inconsistency was also found and fixed: RNA was **optimised** over the
full 60-bin window but **reported** as a mean over the CDS only, about a third of
it. Reporting both windows shows this does not matter — the two agree within 7% on
every design (ratio 1.01–1.07) — so no RNA conclusion changes. Both are now
recorded.

### Testing the area target directly

Arm B rerun with CAGE spread over 5 bins at one fifth the amplitude, preserving the
requested area, everything else identical:

| | edits | CAGE % | **breadth** | **me3/me1** | orientation | genes lost |
|---|---|---|---|---|---|---|
| B, CAGE as a 1-bin spike | 689 | 107.6 | **1.07** | 0.41 | 0.838 | 0 |
| **B, CAGE as a 5-bin area** | 648 | 52.6 | **2.98** | **0.49** | 0.844 | **1** |
| lab's `final` | 2,909 | 112.1 | 2.14 | **3.60** | 0.993 | 0 |

**The shape hypothesis is confirmed and the promoter-identity hypothesis is not.**
Broadening the target does broaden the result, from 1.07 to 2.98, overshooting the
lab's 2.14. But promoter character moves only 0.41 → 0.49, against 3.60. Initiation
breadth is now ruled out **causally**, not just correlationally: breadth was
manipulated directly and chromatin character did not follow.

Two side effects worth recording:

**It fixed the divergent-initiation defect.** Arm B's mScarlet promoter scores 0.566
on the orientation index with a spike target (forward 20.77, reverse 15.91). With
the area target it scores **0.866** (forward 8.82, reverse 1.36). Asking for a zone
rather than a point removed the bidirectionality, which is consistent with a spike
being satisfiable by an enhancer.

**It cost a gene and most of the magnitude.** TagBFP collapsed to 0.73 against a
target of 12, and total CAGE fell to 52.6%. A flat 5-bin target at one fifth
amplitude is a harder request than a single bin at full amplitude, and the optimiser
paid for it by abandoning the third gene — the same failure mode as
[§3](#3-gene-abandonment-is-the-default-not-an-accident).

So the correct CAGE window is worth asking the authors about, and it is **not** the
explanation for promoter identity. That now survives five tests: budget, edit
penalty, explicit chromatin terms, initiation breadth by correlation, and
initiation breadth by direct manipulation.

## 9. A per-gene floor, and why five earlier negative results were confounded

Every term added to the objective in §1–§8 cost a gene: H3K4me1 suppression,
H3K4me3↑/H3K4me1↓, the 5-bin area target, and reverse-strand suppression all
dropped one gene below 10% of its own CAGE target. [§3](#3-gene-abandonment-is-the-default-not-an-accident)
identified the cause — the loss sums error over genes, so abandoning one and
overdelivering on two is legal — and proposed a per-gene floor as the fix. It was
never tested. It is now.

The floor is a one-sided hinge: a gene contributes nothing while it is above 60% of
its own CAGE target, and a quadratic penalty below it.

```python
got = y_pred[:, cage_row, tb:tb + 5].sum(-1)
pen = pen + torch.relu(floor_g - got).pow(2).mean()
```

### It rescues both terms that had been written off

| Objective | floor | CAGE % | RNA % | **me3/me1** | min OI | genes lost |
|---|---|---|---|---|---|---|
| CAGE+RNA | no | 107.6 | 98.2 | 0.41 | 0.566 | 0 |
| CAGE+RNA | **yes** | 77.4 | 94.3 | **0.90** | 0.569 | 0 |
| + reverse-strand suppressed | no | 72.9 | 97.9 | 0.34 | 0.708 | **1** |
| + reverse-strand suppressed | **yes** | **114.1** | 84.3 | 0.83 | **0.915** | **0** |
| + H3K4me1 suppressed | no | 70.1 | 97.7 | 0.78 | 0.815 | **1** |
| + H3K4me1 suppressed | **yes** | 99.1 | 72.9 | **10.07** | 0.601 | **0** |
| lab's `final` | — | 112.1 | 131.6 | 3.60 | 0.991 | 0 |

**Both previously-failing terms work once the optimiser cannot pay for them with a
gene.** Reverse-strand suppression goes from losing TagBFP at 72.9% CAGE to keeping
all three genes at **114.1%** with the best directionality produced here (0.915).
H3K4me1 suppression goes from 0.78 promoter character to **10.07**, nearly three
times the lab's 3.60 and the first design in this work that is strongly
promoter-like rather than enhancer-like.

The floor alone, with no other change, more than doubles promoter character
(0.41 → 0.90).

### What this corrects

[Correction 12](08-corrections.md) retracted the recommendation to constrain
H3K4me3/H3K4me1 explicitly, on the evidence that arm N made things worse. That
retraction was itself confounded: arm N lost a gene, and the loss of the gene is
what destroyed the chromatin signal. With a floor the chromatin term does what was
originally claimed for it. See [correction 15](08-corrections.md).

Five negative results in this document — the area target, reverse-strand
suppression, H3K4me1 suppression, the H3K4me3↑/H3K4me1↓ pair, and the breadth
manipulation — were all measured on designs that had abandoned a gene. **Gene
abandonment was the confound, not the terms.**

### The terms still do not compose

| Objective | CAGE % | RNA % | me3/me1 | min OI | criteria met |
|---|---|---|---|---|---|
| + reverse-strand, floor | 114.1 | 84.3 | 0.83 | 0.915 | 3 / 5 |
| + H3K4me1, floor | 99.1 | 72.9 | **10.07** | 0.601 | 3 / 5 |
| + H3K4me3↑/H3K4me1↓, floor | 89.2 | 54.1 | 0.68 | 0.640 | 2 / 5 |
| **+ reverse-strand AND H3K4me1, floor** | 75.5 | **54.3** | 0.80 | 0.877 | **2 / 5** |
| **lab's `final`** | 112.1 | 131.6 | 3.60 | 0.991 | **5 / 5** |

Criteria: CAGE 85–125% of target, RNA 85–135%, all three genes retained, min
orientation index ≥ 0.8, H3K4me3/H3K4me1 > 1.5.

Combining the two terms that each work individually produces something worse than
either. And the cost is now visible in one column: **RNA falls monotonically as
terms are added** — 98.2, 94.3, 84.3, 72.9, 54.1, 54.3. With the floor blocking
gene abandonment, the optimiser pays for each new requirement out of transcript
level instead.

**Across 27 designs and 19 objective variants, the lab's `final` is the only one
that satisfies all five criteria at once.** This work can now match or beat it on
any single axis — 114.1% versus 112.1% on CAGE, 10.07 versus 3.60 on promoter
character — and cannot get them together. `final` sits at RNA 131.6%,
*overdelivering* while every loaded objective here undershoots, which suggests it
found a sequence that satisfies everything rather than a point on a trade-off
curve. How it did that is still the open question.
