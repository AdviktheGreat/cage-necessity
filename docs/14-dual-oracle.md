# 14. The dual-oracle arm: two initiation models in one loss

Run on 2026-11-12 at the lab's request. This is the last design arm of the project and the
only one built specifically to attack a failure the project had already diagnosed, rather than
to explore the objective space.

## Why this arm exists

Two findings pointed at the same fix.

[Correction 20](08-corrections.md) found that the RNA+PRO-cap arm **overfits the oracle inside
its own loss**. ProCapNet, which was in the loss, called all three promoters cleanly divergent.
AlphaGenome 1, which was in no loss, called one of three. The design had learned ProCapNet, not
initiation.

The [per-gene floor result](10-multitrack-chromatin.md) found that every CAGE-free arm builds
one good promoter and abandons two, because a loss that sums error over genes makes that trade
legal.

So: put **two initiation oracles that share no training data** in the same loss, and put a
per-gene floor on each. If the two models have to agree at the same position for the same gene,
agreement should be harder to fake. That was the hypothesis.

## The objective

`design_dual()` in [`../code/borzoi_modal_app.py`](../code/borzoi_modal_app.py). One autograd
graph holds both models: Borzoi sees the whole 524,288 bp input, and ProCapNet sees three
2,114 bp promoter windows sliced out of the same one-hot tensor.

| term | model | form |
|---|---|---|
| CAGE | Borzoi, row 306 | MSE to a 5-bin area target, 28/20/12 ([correction 11](08-corrections.md)) |
| CAGE floor | Borzoi | one-sided hinge per gene at 60% of that gene's own target, strength 10 |
| RNA | Borzoi, row 6122 | MSE to 7/5/3 over the transcript |
| PRO-cap sense | ProCapNet | one-sided hinge per gene, 50% of mass in the central 51 bp |
| PRO-cap antisense | ProCapNet | one-sided hinge per gene, antisense held under 5% of sense |
| PRO-cap count | ProCapNet | one-sided hinge per gene, log-count 5.0 / 4.65 / 4.14 |

Nine one-sided hinges, three per gene, plus two MSE terms. 600 iterations, 2,025 edits of
17,847 mutable positions, 115 s on an A100. Output loss fell from 1.59e3 to 1.075.

## Result 1: the best design in the project on the lab's primary target

Yash's stated priority is *"the main thing we want to see is dnase peaks at promoters."* The
honest statistic for that is **DNase at the worst of the three promoters**, because a design that
abandons one gene is not usable. DNase appears in no design loss anywhere in this project, so all
three models below are reporting on something they were never asked to produce.

| design | Borzoi worst | **AG1 worst** | AG2 worst |
|---|---|---|---|
| CDS-only control | 0.07 | 0.15 | 0.15 |
| lab's `final` | 0.37 | 1.30 | 1.40 |
| CAGE+RNA | 0.15 | 0.50 | 0.72 |
| CAGE+RNA+floor | 1.09 | 2.52 | **7.28** |
| RNA+PRO-cap | 0.10 | 0.23 | 0.23 |
| **CAGE+RNA+PRO-cap+floors** | **1.86** | **4.59** | 6.25 |

By AlphaGenome 1 the margin is 1.8x over the previous best and 3.5x over the lab's `final`. AG2
puts it a close second to the CAGE+floor arm.

**It wins on evenness, not on peak height** (figure 18, panel C). AlphaGenome 1 DNase per
promoter:

| design | eGFP | mScarlet | TagBFP |
|---|---|---|---|
| CAGE+RNA+floor | 2.52 | 8.88 | 4.81 |
| RNA+PRO-cap | 4.34 | 0.42 | 0.23 |
| **CAGE+RNA+PRO-cap+floors** | **4.59** | **6.00** | **5.59** |

Every other design in the project has at least one weak promoter. This one does not. The same
pattern holds in AlphaGenome 1's independent CAGE readout, where the arm puts 664 / 1468 / 664
at the three promoters against the CAGE+floor arm's 38.5 / 1519.7 / 41.9 (one real promoter and
two failures) and the PRO-cap arm's 267 / 40 / 50.

ProCapNet reads also even out, which is the clearest sign the floors did their job: 46.3 / 82.7 /
99.3, against the PRO-cap-only arm's 68.9 / 84.8 / **32.2**.

The cost is RNA, which falls to 83.2% of target from 95.9%. That is the monotonic decline already
documented as terms accumulate, and it is the expected price of nine extra hinges.

## Result 2: directionality still does not survive an independent model

The hypothesis was wrong.

| design | ProCapNet (in the loss) | AlphaGenome 1 (in no loss) |
|---|---|---|
| lab's `final` (never saw PRO-cap) | 0 of 3 | **2 of 3** |
| RNA+PRO-cap | **3 of 3** | 1 of 3 |
| **CAGE+RNA+PRO-cap+floors** | 2 of 3 | **0 of 3** |

Antisense/sense ratios for the new arm are 0.308 / 0.049 / 0.149 by ProCapNet and
3.295 / 1.215 / 0.475 by AlphaGenome 1. Two independent oracles in the same loss did not make
the result transfer. They were fooled together.

Worse, the ranking inverts completely. The design that looks most divergent to the independent
model is the lab's `final`, which never saw a PRO-cap track in its life, and the two arms that
were explicitly optimised for divergence come last. **Optimising a readout destroys that
readout's value as evidence, and adding a second oracle to the loss destroys both.**

This is the project's [unifying result](../README.md) confirmed prospectively rather than
retrospectively. Up to now the evidence was a correlation: cross-model agreement fell
monotonically with how hard the objective pushed on a readout. This arm was a direct test, and
it reproduced the effect on demand.

## What this changes for the pipeline

Nothing about the recommendation, which gets stronger. **Filter on readouts that are absent from
the loss.** The fix for oracle overfitting is not a better loss or more oracles in it; it is
holding a model back. If the lab wants PRO-cap to be a *check*, PRO-cap must stay out of the
design objective.

The practical recommendation does change in one way. If the goal is three usable promoters,
**this arm is the one to build**: it is the only design in the project where no promoter is
abandoned according to two models that saw none of its objective. Its directionality should be
treated as unmeasured, not as good.

## Caveats

- One run, one seed, one initialisation. The weight space is not smooth, and this project has
  already shown that terms do not compose.
- AG1 and AG2 share a lineage, so the two "independent" checks are not fully independent.
- The ladder order is split: Borzoi says yes (5.90 / 3.86 / 2.71), AG2 says yes
  (8.14 / 0.798 / 0.619), AG1 says no (0.806 / 0.856 / 0.609, and those values are nearly flat,
  a 1.4x spread against a 2.33x target). Treat the ladder as unresolved, exactly as
  [correction 21](08-corrections.md) requires.
- DNase is a prediction, not a measurement. Flow cytometry still comes first.

## Reproduce

```
modal run code/borzoi_modal_app.py --step dualdesign --arm D_cage_rna_procap_floor
modal run code/borzoi_modal_app.py --step score
modal run code/borzoi_modal_app.py --step procapprof
modal run code/borzoi_modal_app.py --step procapscore
modal run code/alphagenome_app.py::score  --designs "...,D_cage_rna_procap_floor"
modal run code/alphagenome_app.py::score2 --designs "...,D_cage_rna_procap_floor"
```

Data: [`../data/tables/table15_dual_oracle.csv`](../data/tables/table15_dual_oracle.csv),
[`../data/sequences/designed_arms/D_cage_rna_procap_floor.fasta`](../data/sequences/designed_arms/D_cage_rna_procap_floor.fasta),
[`../figures/fig18_dual_oracle.png`](../figures/fig18_dual_oracle.png).
