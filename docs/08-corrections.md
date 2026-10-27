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

## The lesson

Every substantive error above has one cause: padding a 19,968 bp design into Borzoi's required
524,288 bp input with composition-matched random DNA. That felt like the neutral, fair choice. It is
in fact the most suppressive of the three options tested, worse than literally filling the flanks
with N.

Anyone running short designs through Borzoi, Enformer or AlphaGenome should treat the padding as an
experimental variable and report it, not as boilerplate.
