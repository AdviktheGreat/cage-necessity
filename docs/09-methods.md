# 09. Methods

## Compute

| Item | Detail |
|---|---|
| Platform | Modal, workspace `advikanand7`, A100-80GB |
| Total compute | ~1 GPU-hour, about $3 |
| Model | `borzoi-pytorch` 0.5.1, weights `johahi/borzoi-replicate-0` |
| Editor | `ledidi` 2.2.0 |
| Precision | bf16 autocast for forward and backward; fp32 checked and equivalent |

### Benchmarks (524,288 bp input, A100-80GB)

| Operation | Time | Peak memory |
|---|---|---|
| Forward, fp32 | 0.87 s | 3.5 GB |
| Forward, bf16 | 0.25 s | 3.1 GB |
| Forward + backward, batch 1 | 0.60 s | 15.4 GB |
| Forward + backward, batch 2 | 0.47 s | 30.3 GB |

A full 600-iteration Ledidi run at batch 2 takes about 200 s, not the ~84 min reported in the deck.

## Model fidelity

`borzoi-pytorch` replicate 0 was validated against the independently served Borzoi used for the
initial runs, on the same 524,288 bp input:

| Metric | Value |
|---|---|
| Pearson over the 624 design bins | **0.9988** |
| argmax bin | 298 in both |
| Peak amplitude | 5.943 versus 7.589 (ratio 1.28) |

Reverse-complement averaging was tested as an explanation for the amplitude gap and made the match
*worse* (plain forward $r$ = 0.9988; rc with strand swap 0.9978; rc without swap 0.9898), so plain
forward passes are used throughout. The residual 1.28× offset is small next to Borzoi's own
replicate spread on the same sequence (7.59, 7.33, 19.86, 11.42), and all comparisons here are
internal to one model.

## Sequence geometry

| Quantity | Value |
|---|---|
| Design length | 19,968 bp |
| Borzoi input length | 524,288 bp |
| Design offset in input | 252,160 |
| Borzoi crop | 163,840 |
| Output bins | 6,144 at 32 bp |
| Design output bins | **2,760 to 3,384** (624 bins) |

CDS coordinates (0-based, half-open): eGFP 4462–5182, mScarlet 9644–10343, RFP3 14805–15507.

## Flanking contexts

| Name | Source |
|---|---|
`orig` | composition-matched random DNA, fixed seed, A 0.2973 / C 0.1977 / G 0.2091 / T 0.2959 |
`rand` seeds 101/202/303 | same composition, different seeds |
`N` | all-N padding |
**`aavs1`** | **hg38 `chr19:54,862,840-55,387,128`, retrieved via UCSC REST API in 9 chunks; 524,288 bp, 0 N, 52.8% GC** |

All substantive results use `aavs1`. `table1_context_robustness.csv` is the exception noted below.

## Borzoi track indices used

| Name | Index | Identifier |
|---|---|---|
| CAGE+ | 306 | `CNhs11250` K562 |
| CAGE− | 307 | `CNhs11250` K562 |
| CAGE2+ (held out) | 872 | `CNhs12336` K562 ENCODE |
| DNase | 1397 | `ENCFF413AHU` K562 |
| DNase2/3/4 (held out) | 1398 / 1399 / 1901 | `ENCFF868NHV` / `ENCFF565YDB` / `ENCFF971AHO` |
| RNA+ | 6122 | `ENCFF980ZHM` K562 |
| RNA− | 6123 | `ENCFF980ZHM` K562 |
| chromatin RNA+ | 6254 | `ENCFF124AUH` K562 chromatin fraction |
| nuclear RNA+ | 6250 | `ENCFF669RVI` K562 nuclear fraction |
| Ser5-P Pol II | 2450 | `ENCFF811SWR` `POLR2AphosphoS5` K562 |
| Pol II total | 2477 | `ENCFF547IOW` `POLR2A` K562 |
| H3K4me3 | 2199 | `ENCFF814IYI` K562 |
| TBP | 2852 | `ENCFF727KHF` K562 |

Track survey result worth recording: Borzoi's human head has 7,611 tracks (ChIP 3,886, RNA 1,543,
CAGE 1,276, DNase 674, ATAC 232). **No ATAC track exists for K562.** There are no PRO-seq, NET-seq,
GRO-seq or TT-seq tracks anywhere in the panel.

## Readout definitions

- **CAGE / DNase / ChIP peak:** maximum over bins $[t-4, t+6)$ where $t$ is the gene's TSS bin,
  $t = \lfloor(\text{ATG} - 128)/32\rfloor$.
- **RNA mean:** mean over the gene's CDS bins.
- **ISM value:** $\Delta = \text{mutant} - \text{reference}$, averaged over the three alternative
  bases at that position. Negative means activating.

## Ledidi configuration

```python
Ledidi(head, shape=(4, 524288), output_loss=MSELoss(),
       input_mask=mask,          # True = frozen; CDSs and all flanks frozen
       l=1e-3, lr=0.5, batch_size=2,
       max_iter=600, early_stopping_iter=600,
       random_state=seed)
```

The oracle is wrapped so that `forward` returns only the selected tracks over the 624 design bins,
scaled by $\sqrt{w}$ per track so that per-track weights enter the MSE linearly. Targets are full
624-bin profiles that are zero everywhere except the specified peaks and blocks, so off-target
signal is penalized.

Mutable positions: 17,847 (all non-CDS positions inside the design). The best member of the returned
batch is selected by output MSE.

## Target construction

From slides 50–51: TSS at ATG − 128 (4 bins), 3′UTR 1,056 bp (33 bins).

- CAGE target: single bin at the TSS, height 28 / 20 / 12 for eGFP / mScarlet / RFP3.
- RNA target: block from the TSS bin to $\lfloor(\text{stop} + 1056)/32\rfloor$, height 7 / 5 / 3.
- DNase and Pol II targets (arms C and F): single bin at the TSS, height 7 / 5 / 3.

## Known inconsistency

`table1_context_robustness.csv` and the context-robustness table in
[04](04-context-robustness.md) were computed with strand-swapped reverse-complement averaging
(the default argument at the time). Everything else uses plain forward passes. Within-table
comparisons are valid; absolute values should not be compared across that boundary. For reference,
the AAVS1/final CAGE peaks are 10.021 / 15.846 / 9.887 under rc-swap and 10.752 / 13.512 / 8.753
under plain forward.

## Reproducibility

All randomness seeded: flank generation (`numpy.random.RandomState`), scrambles
(`random.Random(7)`), Ledidi (`random_state=0`, `torch.manual_seed(0)`).

See the [README](../README.md#reproducing) for the command sequence.
