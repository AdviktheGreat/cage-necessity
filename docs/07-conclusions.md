# 07. Conclusions and recommendations

## The question was two questions

The project's framing was "CAGE is expensive, can we do this with RNA-seq only?" That conflates two
things with opposite answers.

**Design time.** Borzoi predicts all 7,611 tracks in a single forward pass. A CAGE term in the
Ledidi loss is just an index into a tensor that has already been computed. It costs nothing.

**Validation time.** Running CAGE on a cell line costs real money.

Separating them resolves the project:

| | Is CAGE needed? | Evidence |
|---|---|---|
| In the design loss | **Yes, and irreplaceably** | CAGE+RNA reaches 100.8% of target; RNA-only 1.3%; best CAGE-free substitute 9.7% |
| In wet-lab validation | **Largely no** | element-level $r(\text{CAGE},\text{DNase})$ = +0.942 to +0.998 Pearson, +0.73 to +0.83 rank |

## Why CAGE is irreplaceable in the loss

Two independent experiments converge.

**The RNA objective is nearly free to satisfy.** Every one of six arms reached 99–112% of the RNA
target, and the RNA-only arm did it with 555 edits out of 17,847 mutable positions. An objective that
can be satisfied from almost any direction cannot constrain the solution, so gradient descent is
never forced to build a transcription start site. CAGE is the only term that forces it.

**No neighbouring readout can stand in, because they share the element but not the nucleotides.**
Scrambling a designed promoter destroys CAGE and DNase together ($r \approx 0.99$ at element level),
but per-base ISM correlation between them is only −0.04 to +0.35, and their attribution centroids
are 60–150 bp apart. Gradient descent follows nucleotides, so a DNase or Pol II term pushes the
sequence somewhere else and never assembles a core promoter. The ablations confirm this directly:
DNase 1.4%, Ser5-P Pol II 0.4%, chromatin RNA 9.7%.

This is an **identifiability** property of the objective, not a capacity limit of the model. Borzoi
can represent the right answer under every objective. The question is whether the objective picks it
out uniquely, and only CAGE does.

## Recommendations

### 1. Keep CAGE in the design loss
Free, irreplaceable, and removing it costs the entire phenotype.

### 2. Validate with RNA-seq + ATAC, calibrated against CAGE
Element-level Pearson coupling is 0.942–0.998 for DNase and 0.974–0.999 for H3K4me3, but those
coefficients are inflated by the dead-versus-working dynamic range. Rank correlation over the same
arms is +0.73 to +0.83. Accessibility is therefore a **near-perfect detector of whether a promoter
fired and a moderate predictor of how strongly**.

Two reasons to keep a little CAGE rather than none:

- **mScarlet is a real false negative.** In the lab's `final` design it has the highest CAGE of the
  three genes (13.5) and background DNase (0.37). An ATAC-only screen would have called it dead.
- **Across-design rank correlation is only 0.74–0.81**, and the discordance comes entirely from the
  deliberately pathological arms (RNA-only, RNA+DNase, RNA+Pol II) that hit RNA targets with no
  initiation. A pipeline that keeps CAGE in the loss should not generate those.

So: RNA-seq + ATAC for everything, with CAGE on two or three constructs to fix the relationship.

One note for K562 specifically: **Borzoi has no ATAC track for K562** (232 ATAC tracks exist, none
in K562; K562 has 4 CAGE, 7 DNase, 424 ChIP, 94 RNA). In-silico accessibility for K562 means DNase.

### 3. Try CAGE + RNA + chromatin-fraction RNA
The clearest gap in the current objective. Arm B (CAGE+RNA) overshoots CAGE 44.0 against the lab's
26.9 while undershooting every chromatin mark by 3–6×. Arm E (RNA + chromatin RNA), with no CAGE at
all, matches the lab's design on DNase, H3K4me3 and Pol II. Combining all three is untested and
should give both initiation precision and chromatin realism. Track index 6254/6255,
`RNA:K562 chromatin fraction` (`ENCFF124AUH`).

### 4. Add an explicit antisense penalty
Every arm here produced RNA+/RNA− of 2.6–4.9 against the lab design's 59.7. Specifying only RNA+
lets the optimizer build bidirectional promoters. Free to add, and the designs are confirmed
(+) strand so the constraint is well-posed.

### 5. Design in the intended integration locus
A 437× swing from padding alone, and composition-matched random padding is the *worst* of the three
options tested. If the production pipeline pads with random sequence, switching to AAVS1 or the
actual target locus is potentially a large and free improvement.

### 6. Do not add DNase to the objective
It makes accessibility worse, not better. The arm that targeted DNase achieved the least DNase of
any designed arm.

### 7. Report replicate spread
Borzoi's four replicates disagree 2.7× on the same designed CAGE peak (7.59, 7.33, 19.86, 11.42).
Replicate 2 lands suspiciously close to the exact target. Any single-replicate magnitude should
carry an uncertainty band, and it is worth checking which replicate the production pipeline uses.

## What this does not establish

This is a model auditing a model. Everything here is Borzoi explaining Borzoi, which is the correct
level of analysis for "why does this objective behave this way" and tells you nothing about whether
cells agree. The deck's own experimental results are the relevant check, and they already report
that positions transfer while heights do not.

The two design liabilities found by direct sequence inspection ([02](02-sequence-analysis.md)) are
also invisible to this analysis: Borzoi does not penalize the premature polyadenylation signal 34 bp
after the mScarlet stop (RNA+ 0.269 before versus 0.250 after), so that risk survives optimization
and must be caught by hand.
