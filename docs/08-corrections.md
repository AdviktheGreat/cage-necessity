# 08. Corrections log

Claims made during this work that later proved wrong, and why. Kept because the reason for each
error is methodologically instructive: all but one trace to a single bad choice, the random
flanking padding.

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
design objective.

### 5. Promoter sufficiency numbers
**Claimed:** transplanting the 600 bp promoter into the control recovered 46–86% of CAGE.

**Actually:** in AAVS1 it recovers 0–49% of CAGE and ~1% of RNA. The promoters are necessary but
clearly **not** sufficient.

## Corrected for imprecision, not retracted

### 6. "No TATA boxes anywhere"
**Stated too broadly.** Correct version: no `TATAAA` sits at the canonical core-promoter position in
any gene in any of the three sequences. However mScarlet's design did acquire a `TATAAA` 17 bp
upstream of its designed TSS, and RFP3 has one at −258 on the minus strand. Both carry
below-background ISM attribution, so the functional conclusion (TATA-less, CpG-enriched promoters)
stands.

## The lesson

Every substantive error above has one cause: padding a 19,968 bp design into Borzoi's required
524,288 bp input with composition-matched random DNA. That felt like the neutral, fair choice. It is
in fact the most suppressive of the three options tested, worse than literally filling the flanks
with N.

Anyone running short designs through Borzoi, Enformer or AlphaGenome should treat the padding as an
experimental variable and report it, not as boilerplate.
