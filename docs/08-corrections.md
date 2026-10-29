# 08. Corrections log

Claims made during this work that later proved wrong, and why. Kept because the reason for each
error is methodologically instructive. Corrections 1-5 all trace to one bad choice, the synthetic
flanking padding. Corrections 7-9 share a different root: measuring only what was being optimized,
and declaring success before checking anything outside the objective.

## Retracted

### 1. "eGFP effectively failed"
**Claimed:** the eGFP gene in the lab's design produced essentially no CAGE and had RNA *lower* than
the CDS-only control, so it failed.

**Actually:** in real AAVS1 context eGFP reaches CAGE 10.75 and RNA 8.42, which is 120% of its RNA
target of 7. It works.

**Cause:** random padding suppressed it 437×, more than any other gene.

### 2. "The height ladder is inverted"
**Claimed:** achieved RNA was ascending (0.067 / 0.265 / 0.527) against a descending 7/5/3 target,
reproducing the deck's slide 33–34 anti-correlation.

**Actually:** in AAVS1, RNA is 8.42 / 7.14 / 4.44 against targets 7 / 5 / 3. The **ordering is
correct**. The real defects are milder: everything overshoots 1.2–1.5× and the dynamic range is
compressed from 2.33× requested to 1.90× delivered.

**Cause:** same padding artifact.

### 3. "Borzoi reads RNA coverage off CDS-like composition, and that is the null space RNA-only optimization exploits"
**Claimed:** the CDS-only control produced RNA signal (0.162 at eGFP) with no promoter, showing that
the RNA objective has a promoter-free solution available.

**Actually:** in AAVS1 the control gives RNA 0.154 against the design's 8.42, i.e. 1.8%. Coding
sequence produces essentially no RNA signal in realistic context.

**Cause:** padding artifact. This was the most consequential error because it was a *mechanistic*
claim, and it was wrong.

**Replaced by:** the correct mechanism, established by ablation in
[06](06-objective-ablations.md): the RNA objective is trivially satisfiable (555 edits, 3.1% of
mutable positions, 106% of target) and therefore does not constrain the solution. Not "RNA signal
comes from CDS" but "RNA targets are cheap to hit by many routes".

### 4. "CAGE and DNase share a causal substrate" (as originally stated)
**Claimed:** from random-flank data, $r(\text{CAGE},\text{DNase})$ = +0.91/+0.93/+0.67 across
perturbation arms, so the two are driven by the same bases and accessibility comes free with CAGE.

**Actually:** the element-level claim is correct and stronger in AAVS1 (+0.942 to +0.998). But the
implied base-level claim is wrong: per-base ISM correlation is only −0.04 to +0.35.

**Refined to:** same element, different nucleotides. This distinction turned out to be the key to
the whole project, because it explains why DNase works as a validation surrogate but fails as a
design objective. The strength of the element-level claim is itself qualified in item 7 below.

### 5. Promoter sufficiency numbers
**Claimed:** transplanting the 600 bp promoter into the control recovered 46–86% of CAGE.

**Actually:** in AAVS1 it recovers 0–49% of CAGE and ~1% of RNA. The promoters are necessary but
clearly **not** sufficient.

### 9. "The short CAGE+RNA arms beat the lab's design" — retracted

**Claim retracted:** that the short CAGE+RNA arms produced here outperform the
lab's `final` design. Specifically, that arm J "beats the lab's design on CAGE
(86.7% vs 55.0%), on RNA accuracy, and on parsimony (1,072 edits vs 2,909)," and
that arm B at 100% of the CAGE target was the strongest design in the set.

**What was wrong:** H3K4me1 was not in the held-out scoring set, so promoter
versus enhancer character was never measured. Adding it reverses the ranking.

| Design | H3K4me3/H3K4me1 | CAGE % of target |
|---|---|---|
| lab's `final` | **3.60** | 55.0 |
| Markov initial (untouched) | 0.91 | 0.1 |
| J CAGE+RNA+chrRNA | 0.89 | 86.7 |
| B CAGE+RNA | **0.41** | **100.0** |

`final` is the only design in the set with promoter-like chromatin character, by
a factor of four over anything else. Every arm produced here is enhancer-like,
and the arm with the best CAGE attainment has nearly the worst promoter
character. Second place in the table is the untouched starting sequence.

**Why it happened:** the arms were scored on the objective they were optimised
for. This project's central diagnosis is that a CAGE-only objective is satisfied
by enhancers because enhancers are cheaper than promoters. The same mechanism
operates, less extremely, under CAGE+RNA — and it went unnoticed for exactly the
reason it was predicted to be easy to miss: nothing in the scoring set
distinguished the two.

**Revised conclusion:** CAGE + RNA is sufficient to hit the CAGE target, not to
produce a promoter. The lab's longer run reached promoter identity that none of
the 600-iteration arms did, which means its larger edit count is a cost of
building the right object rather than evidence of inefficiency. Promoter identity
does not appear anywhere in the objective, so it is currently emergent rather
than guaranteed; an explicit H3K4me3/H3K4me1 term would make it a requirement.

Detail: [multi-track objectives](10-multitrack-chromatin.md#4-promoter-identity-final-is-the-only-design-that-has-it).

## Corrected for imprecision, not retracted

### 6. "No TATA boxes anywhere"
**Stated too broadly.** Correct version: no `TATAAA` sits at the canonical core-promoter position in
any gene in any of the three sequences. However mScarlet's design did acquire a `TATAAA` 17 bp
upstream of its designed TSS, and RFP3 has one at −258 on the minus strand. Both carry
below-background ISM attribution, so the functional conclusion (TATA-less, CpG-enriched promoters)
stands.

### 7. "ATAC/DNase is a near-perfect CAGE surrogate"
**Stated too strongly.** The quoted +0.942 to +0.998 are Pearson coefficients over arms spanning
0.067 to 5.6, so the on/off contrast carries them. Rank correlation over the same data is +0.733 to
+0.833.

**Refined to:** accessibility is a near-perfect *detector* of whether a promoter fired and a moderate
*ruler* of how strongly. The correction matters because mScarlet in the lab's own design is a true
false negative: highest CAGE of the three genes (13.51), background DNase (0.37). An ATAC-only
screen would have scored the strongest promoter as dead. The recommendation is unchanged in
direction but now explicitly includes CAGE calibration on a few constructs rather than dropping it
entirely.
### 8. "ATAC/DNase can replace CAGE in wet-lab validation" — retracted
**Claimed:** across nine perturbation arms, $r(\text{CAGE}, \text{DNase})$ = +0.942 to +0.998, so
accessibility can stand in for CAGE and the lab can drop CAGE from routine validation. Item 6 above
was a first, insufficient softening of this to "a detector, not a ruler".

**Actually:** the claim fails even as a detector. Scored as a screen (CAGE > 5 works, DNase > 1.0
works), accessibility disagrees with CAGE on **7 of 21** gene-design combinations, including **all
three genes of arm B**, the best design in the study. And among designs whose loss contained CAGE,
which is the production case, the correlation is Pearson **−0.212**, Spearman **−0.033**.

**Cause: I generalized from the wrong population.** The +0.99 figures came from nine arms that are
all `final`, the CDS-only control, or constructions from them. Inside one sequence family, breaking
a promoter kills CAGE and accessibility together, so they co-vary almost perfectly. I read that as a
general property of the two assays. It is not. Across designs produced by *different objectives*,
the coupling drops to +0.218, and among CAGE-containing designs it vanishes.

Arm E is the clean counterexample and it was in the data the whole time: no CAGE term anywhere in
its loss, 7× the accessibility of arm B, and 10× worse initiation. A design can have excellent
chromatin and no initiation. That single fact is simultaneously why chromatin RNA is the most
interesting CAGE-free objective *and* why accessibility cannot report on CAGE.

**Replaced by:** [07](07-conclusions.md) recommendation 2, which drops the surrogate framing
entirely and asks which specification is being verified instead.

**What survives:** the regime-1 coupling is still a real finding about promoter integrity within one
design family, and it remains the reason DNase *gradients* fail as a design objective. The error was
in the scope of the claim, not the measurement.
### 10. "Any result computed with random padding is suspect" — narrowed

**Claim narrowed:** that "any result computed with random padding is suspect"
and that flanking context causes a 437× swing.

**What was imprecise:** "random padding" conflated two very different things.
The 437× swing is specific to *composition-matched synthetic* random DNA. Testing
ten *real* genomic loci sampled uniformly over the autosomes gives a CAGE spread
of only **1.55×** (CV 13.1%), with AAVS1 at the median.

The mechanism is real DNA versus synthetic DNA, not integration site versus
elsewhere. Synthetic padding at ~40.7% GC with no genes reads as a gene desert;
literal N, carrying no information, is much safer; any real locus is safer still.

**Revised conclusion:** padding with random genomic loci, which is the lab's
practice, is sound. The warning applies to synthetic padding only. The residual
caveat is that RNA is roughly twice as context-sensitive as CAGE (CV 26.5% vs
13.1%), so RNA comparisons across differently padded runs still need care.

Detail: [random genomic loci](10-multitrack-chromatin.md#6-random-genomic-loci-are-safe-synthetic-random-padding-is-not).

### 11. "The lab's design reaches only 55–60% of its CAGE target" — retracted

**Claim retracted:** that the lab's `final` design attains 55.0% (AAVS1, this
work) or 59.6% (earlier runs) of its 28/20/12 CAGE request, and the repeated
framing that the short arms here "beat the lab's design on CAGE."

**What was wrong:** the readout. This work measured CAGE as a **peak height** in a
10-bin window around the TSS. Slide 50 of the source deck describes the
objective as an **area** over the 5′UTR plus 32 bp, roughly 160 bp, which is 5
bins. Those are different quantities, and for a broad initiation zone they differ
by a factor of two or more.

Re-scored as a 5-bin area:

| Design | CAGE as peak | CAGE as area |
|---|---|---|
| **lab's `final`** | **55.0%** | **112.1%** |
| B CAGE+RNA | 100.0% | 107.6% |
| J CAGE+RNA+chrRNA | 86.7% | 110.9% |

**The lab's design is on target in aggregate.** It only looked like a 45% shortfall
because a peak was compared against an area.

**What survives:** the per-gene ladder is still not attained, but differently from
what was reported. As an area, `final` gives eGFP 112.8%, mScarlet 73.6% and
RFP3 **174.7%** of their individual targets, so RFP3 overshoots by 1.75× and
mScarlet undershoots. The aggregate 112.1% averages those out.

Arm B as an area is 109.5% / 103.9% / 109.6%, so it does attain the ladder on
both readings. But see [correction 9](#9-the-short-cagerna-arms-beat-the-labs-design--retracted):
it attains it by building single-bin spikes.

**The lesson, again:** a percent-of-target number is only meaningful once the
target's functional form is pinned down. Every "% of CAGE target" in this
repository before this correction is a peak reading.

### 12. "An explicit H3K4me3-up / H3K4me1-down term would make promoter identity a requirement" — retracted

**Recommendation retracted.** It was proposed in
[correction 9](#9-the-short-cagerna-arms-beat-the-labs-design--retracted) and in
[docs/10](10-multitrack-chromatin.md), and it fails when run.

Arm N (`CAGE + RNA + H3K4me3↑ + H3K4me1↓`, weights 1/1/2/2, 600 iterations):

| | H3K4me3 | H3K4me1 | **H3K4me3/H3K4me1** | genes lost |
|---|---|---|---|---|
| B CAGE+RNA (no such term) | 10.67 | 26.32 | **0.41** | 0 |
| **N, with the term** | **0.69** | 2.15 | **0.32** | 1 |
| lab's `final` | 27.38 | 7.61 | **3.60** | 0 |

Asking for H3K4me3 up and H3K4me1 down **drives both down.** H3K4me3 fell from
10.67 to 0.69 despite being explicitly requested at a target near 17. The
optimiser satisfies the easy half of the instruction (flatten H3K4me1) and
abandons the hard half, because the two marks are too correlated in Borzoi's
training data to be moved in opposite directions. Promoter character ends up
*worse* than not asking at all, and RFP3 is lost.

This strengthens rather than weakens the case for PRO-cap-based design: the
promoter-versus-enhancer distinction appears to be unreachable through Borzoi's
H3K4 tracks, whatever weighting is used.

### 13. Two further hypotheses for the promoter-identity gap, both ruled out

Neither is a retracted claim, but both were stated as the likely explanation and
both are wrong.

**Optimisation budget.** Arm B at 2,400 iterations is **byte-identical** to arm B
at 600 iterations: 0 differing bases out of 19,968, same 689 edits, same MSE
0.0633. Arm J differs by 494 bases with promoter character essentially unchanged
(0.89 → 0.87). Ledidi converges at or before 600 iterations here, so the lab's
longer runtime cannot be buying promoter identity through iterations alone.

**Edit penalty.** Loosening `l` increases edits but makes promoter character
*worse*, the opposite of the prediction:

| Arm | `l` | edits | H3K4me3/H3K4me1 |
|---|---|---|---|
| B | 1e-3 | 689 | 0.41 |
| B | 1e-4 | 1,677 | 0.28 |
| B | 1e-5 | 2,352 | 0.21 |
| J | 1e-5 | 2,625 | 0.60 |
| lab's `final` | ? | 2,909 | **3.60** |

So the lab's 2,909 edits are not what produces its 3.60. This work reached 2,625
edits and got 0.60.

**What is left.** Initiation-zone breadth (5-bin area ÷ peak height) correlates
with promoter character at Spearman **+0.679** across 14 designs, and the single-bin
spikes produced here are the most enhancer-like. But it is not sufficient: arm L
has the same breadth as `final` (2.13 vs 2.14) and only 0.62 promoter character.
`final` remains an outlier that none of budget, edit penalty, chromatin terms or
initiation breadth explains. **The cause is still unidentified**, and the most
likely remaining candidate is a difference in the objective's functional form
beyond peak-versus-area — which is a question for the authors, not something that
can be settled from the outputs alone.

### 14. The third gene identified: it is TagBFP, and "RFP3" is a misnomer

**Resolved, not retracted** — this was listed as an open question since the first
commit. The third protein is **TagBFP**, a blue fluorescent protein.

The 702 bp CDS at 14805–15507 translates to 233 amino acids, and that sequence is
an **exact, character-for-character match** to TagBFP in
[FPbase](https://www.fpbase.org/protein/tagbfp/): 233 of 233 residues, zero
substitutions, aligned against all 1,042 proteins in the database. The runner-up,
mTagBFP2, is 98.3%; the best non-TagBFP hit is 51%.

| | TagBFP |
|---|---|
| length | 233 aa (matches exactly) |
| excitation / emission | **402 nm / 457 nm** (blue) |
| oligomeric state | monomeric |
| brightness | 32.76 |
| PDB | 3M24 |

**Why it looked like a red protein.** TagBFP is engineered *from* TagRFP and is
94.3% identical to it (13 substitutions: N18D, M61L, R65K, P125T, N141F, M144T,
S156N, G166S, C170A, F172I, H195Y, D204N, K205N). Same scaffold, shifted colour.
Reading the sequence by eye put it in the eqFP578/TagRFP lineage, which was correct
but one step short, and the provisional label `RFP3` followed from that. It is blue,
not red.

**Two things this confirms.** The lab's files label this gene HaloTag, which it
is not — HaloTag is 894 bp / 297 aa with no peptide-level similarity. And slide 45's
schematic shows **GFP / BFP / RFP**, which was right all along; the schematic should
have been trusted over the file name.

**The construct is therefore:**

| Position | Gene | Colour | Excite / emit | RNA target |
|---|---|---|---|---|
| 4462–5182 | eGFP | green | 488 / 507 nm | 7 |
| 9644–10343 | mScarlet | red | 569 / 594 nm | 5 |
| 14805–15507 | **TagBFP** | **blue** | **402 / 457 nm** | 3 |

**Practical consequence.** Emission peaks at 457, 507 and 594 nm are three
well-separated channels, so **flow cytometry can read all three levels
simultaneously in single cells** with standard filter sets. That strengthens
[recommendation 2](../README.md): the 7/5/3 expression ladder does not need any
sequencing assay to verify, and the three-colour readout is cleaner than this
work assumed when the third protein was thought to be a second red.

**Note on naming.** The identifier `RFP3` is kept throughout the code, data files
and tables for continuity with earlier results. It refers to TagBFP.

### 15. "An explicit chromatin term cannot work in Borzoi" — retracted, it was gene abandonment

[Correction 12](#12-an-explicit-h3k4me3-up--h3k4me1-down-term-would-make-promoter-identity-a-requirement--retracted)
retracted the recommendation to constrain H3K4me3/H3K4me1 explicitly, concluding
that "the promoter-versus-enhancer distinction appears to be unreachable through
Borzoi's H3K4 tracks, whatever weighting is used." **That conclusion was wrong,
and the reason is instructive.**

Arm N lost TagBFP. So did arm M, the 5-bin area arm, and the reverse-strand arm.
Every chromatin and directionality term tested dropped one gene below 10% of its
target, because the loss sums over genes and sacrificing one is a cheap way to
fund a new requirement. The chromatin readouts were then measured on a construct
with a dead gene.

Adding a per-gene floor — a hinge penalty that activates only below 60% of each
gene's own CAGE target — changes the result completely:

| Objective | floor | genes lost | **H3K4me3/H3K4me1** |
|---|---|---|---|
| CAGE+RNA | no | 0 | 0.41 |
| CAGE+RNA + H3K4me1 suppressed | no | **1** | 0.78 |
| **CAGE+RNA + H3K4me1 suppressed** | **yes** | **0** | **10.07** |
| lab's `final` | — | 0 | 3.60 |

**10.07, against the lab's 3.60.** The chromatin term does exactly what was
originally claimed for it. It had been judged on a run where the optimiser spent
the gene instead.

**What this means for the earlier corrections.** Five negative results in
[docs/10](10-multitrack-chromatin.md) were measured on gene-abandoning designs and
are therefore unsafe: the area target, reverse-strand suppression, H3K4me1
suppression, the H3K4me3↑/H3K4me1↓ pair, and the breadth manipulation. The
per-gene floor should have been the *first* thing tested after
[§3](10-multitrack-chromatin.md) identified abandonment, not the last.

**What survives.** The terms still do not compose — combining reverse-strand and
H3K4me1 suppression gives 2 of 5 criteria, worse than either alone — and RNA falls
monotonically as terms are added (98.2 → 54.3%), because with abandonment blocked
the optimiser pays out of transcript level. The lab's `final` remains the only
design of 27 that satisfies all five criteria simultaneously, and that is still
unexplained.

### 16. "CAGE and DNase correlate at −0.21, so accessibility is not evidence of a promoter" — retracted

**Claim retracted.** This work reported that among CAGE-containing designs the
CAGE-DNase Pearson correlation is −0.212, and used it to argue that accessibility
is a consequence rather than evidence, and that asking for DNase peaks is the
weakest available check.

**What was wrong:** the number. It was computed over **9 designs**, using the
**maximum across the three genes** of each track, and the **peak** CAGE readout
that [correction 11](#11-the-labs-design-reaches-only-5560-of-its-cage-target--retracted)
superseded. Recomputed per gene-design pair across 20 CAGE-containing designs with
the area readout, it is **+0.418**, not −0.21. Collapsing three genes to one
maximum destroyed most of the signal and inverted the sign.

**What is actually true:**

| Objective class | n | worst-promoter DNase |
|---|---|---|
| no CAGE in the objective | 8 | 0.04 – 0.22 |
| CAGE in the objective | 20 | 0.04 – 1.09 |

CAGE-free objectives never get their weakest promoter above 0.22, against a control
baseline of 0.07. CAGE-containing objectives reach 6.92 at individual promoters. So
including CAGE in the objective does produce DNase peaks at promoters, and nothing
else tested does. The authors' reasoning was right.

**What survives, in weaker form.** DNase remains unreliable *per promoter*: only 1
of 28 designs clears 1.0 at all three, and the lab's own `final` has 0.37 at
mScarlet while its neighbours are at 3.70 and 5.65. So accessibility is sound as an
aggregate signal that initiation was designed in, and unsound as a per-promoter
pass/fail test. The per-gene floor is what makes it pass at all three
([docs/10 §10](10-multitrack-chromatin.md)).

**The lesson, which is now the third instance of it:** corrections 7, 8 and 16 all
come from collapsing per-gene data to a single aggregate before correlating. The
three genes behave differently enough that any statistic taking a maximum or a mean
across them should be treated as suspect.

### 17. "CpG density can replace CAGE in the objective" — refuted by direct test, same day it was proposed

[docs/11](11-what-cage-is-doing.md) showed that CAGE's effect on accessibility is
almost entirely mediated by CpG dinucleotide density, and concluded that the
objective could therefore be RNA abundance plus a direct CpG constraint, with CAGE
dropped. That was named as untested. It was then tested, and **the substitution
does not work.**

Two arms added a differentiable CpG hinge (one-sided, target 25 per 600 bp
promoter window, the count the lab's own design reached) to an RNA-only objective,
with no CAGE:

| | eGFP CpG / DNase | mScarlet CpG / DNase | TagBFP CpG / DNase | CAGE % |
|---|---|---|---|---|
| control | 7 / 0.07 | 8 / 0.07 | 9 / 0.07 | 0.3 |
| RNA only | 7 / 0.15 | 9 / 0.11 | 8 / 0.11 | 5.0 |
| **RNA + CpG (w=0.25)** | 24 / **3.16** | 22 / 0.17 | 21 / 0.14 | 14.8 |
| **RNA + CpG (w=1.0)** | 23 / **3.05** | **29** / **0.13** | **27** / **0.12** | 10.7 |
| CAGE + RNA | 22 / 0.80 | 8 / 0.15 | 12 / 0.27 | 107.6 |
| lab's `final` | 25 / 3.70 | 20 / 0.37 | 30 / 5.65 | 112.1 |

**The constraint was satisfied and accessibility still did not follow.** At w=1.0
mScarlet reached 29 CpGs, *more* than the 20 in the lab's `final`, and its DNase
stayed at 0.13 while `final` reaches 0.37. TagBFP reached 27 against `final`'s 30
and stayed at 0.12 while `final` reaches 5.65. This is not gene abandonment — the
hinge is per window and all three were filled. **CpG count is simply not
sufficient.**

**What does survive, and it is not nothing.** eGFP went from 0.07 to **3.16** with
no CAGE anywhere in the loss, against 3.70 for the lab's design at the same
promoter. The previous best worst-case for any CAGE-free objective was 0.86, and
across eight CAGE-free arms none had ever exceeded that. So **an accessible
promoter does not require CAGE in the objective**, which is the first direct
evidence on the authors' actual question, and it took one sequence constraint
rather than a new model.

**But accessibility arrived without initiation.** CAGE stayed at 10.7–14.8% of
target in both CpG arms. The result is an open, CpG-rich region with no transcription
start site, which is the enhancer-like failure mode seen throughout
[docs/10](10-multitrack-chromatin.md). So the mediation result in docs/11 was read
too strongly: CpG density is what carries CAGE's effect *on accessibility*, but
CAGE is also specifying where and in which direction initiation happens, and that
part has no CpG substitute. The two have to be specified separately.

**Weak lead, not a finding.** Over these 12 windows, CpGs falling within 150 bp of
the designed TSS track DNase at +0.764 while total CpGs in the window do not. That
is 12 points with obvious confounds and should not be built on without a proper
test, but it suggests *position relative to the TSS* matters more than density.

**The process lesson.** This correction exists because a published implication was
tested within the hour rather than being left to stand as a plausible suggestion.
Correction 15 took fifteen experiments to find; this took one.

## The lesson

Every substantive error above has one cause: padding a 19,968 bp design into Borzoi's required
524,288 bp input with composition-matched random DNA. That felt like the neutral, fair choice. It is
in fact the most suppressive of the three options tested, worse than literally filling the flanks
with N.

Anyone running short designs through Borzoi, Enformer or AlphaGenome should treat the padding as an
experimental variable and report it, not as boilerplate.
