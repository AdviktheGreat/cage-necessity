# 06. Objective ablations: six Ledidi runs

## Design of the experiment

Everything held constant except the objective:

| Parameter | Value |
|---|---|
| Initialization | CDS-only control (Markov backbone + three CDSs) |
| Context | real AAVS1 locus |
| Mutable positions | 17,847 (all non-CDS positions inside the 19,968 bp design) |
| Frozen | all three CDSs, all 504 kb of flanking sequence |
| Oracle | Borzoi replicate 0, bf16 autocast |
| Iterations | 600 |
| Batch size | 2 |
| Learning rate | 0.5 |
| Input-loss weight $l$ | $10^{-3}$ |
| Targets | CAGE 28/20/12 at the TSS bin; RNA 7/5/3 over TSS → stop + 1,056 bp |
| Loss | MSE over all 624 design bins per included track, including the zeros |

Including the zero bins means off-target signal is penalized, matching the deck's "desired pattern"
plots.

## Results

![Ablations](../figures/fig_ablations.png)

| Arm | Objective | CAGE % of target | RNA % of target | DNase peak | Edits | Output MSE | Runtime |
|---|---|---|---|---|---|---|---|
| A | RNA only | **1.3%** | 106.2% | 0.121 | 555 (3.1%) | 0.0691 | 212 s |
| C | RNA + DNase | **1.4%** | 100.5% | 0.100 | 646 (3.6%) | 0.1336 | 213 s |
| F | RNA + Ser5-P Pol II | **0.4%** | 99.3% | 0.055 | 972 (5.4%) | 0.3014 | 202 s |
| E | RNA + chromatin RNA | **9.7%** | 112.2% | **2.774** | 1004 (5.6%) | 0.1179 | 215 s |
| G | CAGE(×0.25) + RNA | 55.8% | 102.1% | 0.446 | 566 (3.2%) | 0.1543 | 202 s |
| **B** | **CAGE + RNA** | **100.8%** | 102.7% | 0.408 | 689 (3.9%) | 0.0633 | 200 s |
| — | lab's design (reference) | 59.6% | 137.0% | 3.241 | 2,909 | — | ~84 min |
| — | CDS-only control | 0.1% | 2.4% | 0.068 | — | — | — |
| — | Markov initial | 0.1% | 0.2% | 0.083 | — | — | — |

Per-gene detail in `data/tables/table3_objective_ablations.csv`.

## Findings

### 1. The central claim reproduces exactly

Arm A against arm B is a clean reproduction of slides 47 and 49. CAGE goes from 1.3% to 100.8% on
the strength of the CAGE term alone. RNA reaches target either way.

### 2. The RNA objective is nearly free to satisfy

**Every** arm reached 99–112% of the RNA target, and arm A did it with 555 edits (3.1% of mutable
positions). An objective satisfiable from almost any direction cannot constrain the solution, so
gradient descent is never forced to build initiation machinery. This is the mechanism behind slide 47.

### 3. Nothing substitutes for CAGE

DNase (1.4%), Ser5-P Pol II (0.4%) and chromatin RNA (9.7%) all fail. Ser5-P Pol II is *worse* than
RNA alone despite being the most initiation-specific track available, and it has by far the worst
output MSE (0.3014), meaning Ledidi could not satisfy that objective either.

The ISM result in [05](05-ism-mechanism.md) explains it: these tracks share the designed *element*
with CAGE but not the *nucleotides*, so their gradients point elsewhere.

### 4. Undertraining is not the explanation

Arm B reached 100.8% in the same 600 iterations and 200 seconds. The budget is demonstrably
sufficient when the objective is right. The CAGE-free arms are not simply starved of compute.

### 5. DNase is essentially not optimizable

Arm C put DNase directly in the loss and achieved the **lowest DNase of any designed arm** (0.100),
below arm B which never mentioned it (0.408) and far below arm E (2.774). Its output MSE was the
second worst.

So "add DNase to the objective" is not an available fix, and the accessibility seen in these designs
is a byproduct of the CAGE term rather than something that can be requested.

### 6. Loss weighting matters, and failure is gene-specific

Downweighting CAGE 4× (arm G) dropped CAGE achievement from 100.8% to 55.8%. The collapse was not
uniform: mScarlet fell from 20.4 to 0.236 while eGFP (27.2 → 21.1) and RFP3 (12.4 → 10.9) largely
survived. So part of the joint objective's success is relative gradient magnitude, and genes fail
one at a time rather than together. This rules out the simplest "CAGE just dominates the loss"
objection only partially: the effect is real but insufficient to explain a 77× difference.

### 7. Chromatin-fraction RNA buys chromatin realism, not initiation precision

Held-out tracks, never in any objective:

| Design | CAGE2+ | DNase2 | DNase3 | H3K4me3 | TBP | Pol2A | Ser5-P Pol2 |
|---|---|---|---|---|---|---|---|
| lab's design | 26.9 | 2.44 | 2.59 | 17.1 | 16.7 | 22.3 | 16.7 |
| A RNA only | 0.74 | 0.10 | 0.12 | 1.54 | 1.58 | 5.39 | 4.47 |
| **B CAGE+RNA** | **44.0** | 0.32 | 0.45 | 5.02 | 2.75 | 6.94 | 5.65 |
| C RNA+DNase | 0.78 | 0.09 | 0.10 | 1.41 | 1.63 | 6.11 | 5.23 |
| **E RNA+chromRNA** | 4.05 | **2.34** | **2.67** | **15.6** | **7.12** | **19.5** | **18.5** |
| F RNA+Pol2 | 0.30 | 0.07 | 0.06 | 1.12 | 1.87 | 0.72 | 0.72 |

Arm B overshoots CAGE (44.0 against the lab's 26.9) while badly undershooting every chromatin mark.
Arm E, with **no CAGE anywhere in its loss**, essentially matches the lab's design on accessibility,
H3K4me3 and Pol II occupancy.

The lab's design is the only one that gets both, which is probably what its much longer optimization
bought. **CAGE + RNA + chromatin RNA is the untested combination most worth running next.**

### 8. All my arms have an antisense defect

| Design | RNA− peak | RNA+ / RNA− |
|---|---|---|
| lab's design | 0.112 | **59.7** |
| A RNA only | 1.364 | 3.8 |
| B CAGE+RNA | 1.212 | 4.2 |
| C RNA+DNase | 1.040 | 4.9 |
| E RNA+chromRNA | 2.183 | 2.6 |
| F RNA+Pol2 | 0.025 | 196.2 |

Specifying only RNA+ lets the optimizer build bidirectional promoters. The lab's design is clean, so
something in their setup already handles this. If it is not deliberate, an explicit RNA− penalty is
free and worth adding.

## Can RNA-seq alone rank designs the way CAGE does?

Across all nine scored designs × 3 genes:

| Comparison | Pearson | Spearman |
|---|---|---|
| CAGE+ versus H3K4me3 | +0.329 | **+0.812** |
| CAGE+ versus DNase | +0.218 | **+0.791** |
| CAGE+ versus DNase2 | +0.202 | +0.765 |
| CAGE+ versus RNA+ | +0.441 | +0.736 |
| CAGE+ versus chromatin RNA | +0.157 | +0.661 |
| CAGE+ versus Ser5-P Pol2 | +0.242 | +0.661 |
| CAGE+ versus TBP | +0.245 | +0.429 |

Rank correlation is decent but not high, and Pearson is poor. **The discordance comes entirely from
arms A, C and F**, which hit RNA targets with no initiation. Those are exactly the designs a
CAGE-free validation panel would wave through.

Since a production pipeline keeps CAGE in the loss, such designs should not arise. But this is the
reason to calibrate against CAGE once rather than never.

Data: `data/tables/table3_objective_ablations.csv`, `data/tables/table4_heldout_tracks.csv`,
`data/sequences/designed_arms/*.fasta`.
