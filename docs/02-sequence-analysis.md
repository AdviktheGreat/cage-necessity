# 02. Sequence analysis of the lab's designs

Three files were provided:

| File | Contents | Runtime reported |
|---|---|---|
`output_X_0_initial.txt` | order-5 Markov initialization, 19,968 bp | 0.02 s |
`output_X_0_final.txt` | "Ledidi with Borzoi" output, 19,968 bp | 5,048 s (84 min) |
`gfp_scarlet_halo.fasta` | eGFP 720 bp, mScarlet 699 bp, HaloTag7 894 bp | — |

All sequences are (+) strand (confirmed by the lab). All three CDSs are on the + strand and in
frame, consistent with CAGE+ peaks measuring 7× to 293× stronger than CAGE− in this audit.

## Basic comparison

- Both sequences are exactly **19,968 bp**.
- **All changes are substitutions**, so coordinates never shift.
- **4,525 bases differ (22.7%)**.
- Ts/Tv = 0.51, i.e. random-like, not mutation-spectrum-like.
- GC gains 1,108 versus losses 864.

## Inserted CDSs (0-based, half-open)

| Gene | Coordinates | Length | Match to FASTA |
|---|---|---|---|
| eGFP | 4462–5182 | 720 bp | 100% exact |
| mScarlet | 9644–10343 | 699 bp | 100% exact |
| **third gene** | **14805–15507** | **702 bp / 233 aa** | **no match to anything provided** |

### The third gene is not HaloTag

The FASTA's HaloTag7 is 894 bp / 297 aa. The third insert is 702 bp / 233 aa and gives **zero
peptide hits in any of six reading frames**. It encodes an RFP-family protein:

```
MSELIKENMHMKLYMEG...CDLPSKLGHKLN
```

Unverified guess: the TagRFP/mKate family.

Gene spacing is near-perfect (gaps of 4462, 4462, 4462, 4461 bp). A genuine 894 bp HaloTag would
have produced ~4414 bp gaps, so the run was configured with a ~700 bp RFP from the start rather than
substituted later. **This remains unresolved and should be checked with whoever generated the run.**
Referred to as `RFP3` throughout this repository.

## Where the edits concentrate

Edits are enriched exactly where regulatory elements belong.

| Region | Edit rate |
|---|---|
| Outside CDSs, overall | 16.3% |
| 200 bp immediately upstream of an ATG | **22.8% and 26.7%** |
| random 200 bp non-CDS bins (null) | 16.1% ± 4.7 |
| +1,000 to +1,400 bp after each stop | ~23% |

The post-stop enrichment matches the 1,056 bp (33-bin) 3′UTR target from slide 51.

## CpG and promoter character

| Region | initial | final |
|---|---|---|
| Outside CDSs, CpG o/e | 0.23 | **0.50** |
| 400 bp upstream of gene 1 | 0.37 | **0.83** |
| 400 bp upstream of gene 2 | 0.33 | **0.82** |
| 400 bp upstream of gene 3 | 0.29 | **1.07** |

GC content in those windows reaches only 46–48%, just short of the usual CpG-island cutoff. GC-boxes
appear upstream of two genes.

**TATA:** no TATAAA sits at the canonical core-promoter position in any gene, in any of the three
sequences. mScarlet's design did acquire a TATAAA 17 bp upstream of its designed TSS, but ISM shows
it carries below-average attribution ([05](05-ism-mechanism.md)), so these are functionally
TATA-less, CpG-enriched promoters.

## Target model, recovered and verified

From slides 50–51: TSS = ATG − 128 (4 bins), RNA block from TSS to stop + 1,056 bp (33 bins).

| Gene | Predicted TSS bin | Measured CAGE peak bin | CAGE target | RNA target |
|---|---|---|---|---|
| eGFP | 135 | 137–140 | 28 | 7 |
| mScarlet | 297 | **297–298** | 20 | 5 |
| RFP3 | 458 | **458–460** | 12 | 3 |

mScarlet and RFP3 land within one bin. eGFP is off by 2–5 bins (64–160 bp).

## Liabilities found in the design

1. **Premature polyadenylation signal.** `AATAAA` sits **34 bp after the mScarlet stop**, unchanged
   from the Markov backbone, which would truncate the intended 1,056 bp 3′UTR. Borzoi is blind to
   it: RNA+ is 0.269 before versus 0.250 after the signal, so there is no truncation penalty in the
   oracle and the liability survives optimization. **Open experimental risk.**
2. **Upstream in-frame ATG at 14607.** Ledidi removed the two in-frame stop codons that the Markov
   and control sequences retain, creating a potential 66-aa N-terminal extension on RFP3.
3. **Kozak contexts** around the eGFP and mScarlet ATGs are strong but identical in the initial and
   final sequences, so they are incidental rather than designed.

## Borzoi replicate disagreement

On the same sequence, the four Borzoi replicates disagree nearly 3-fold on the designed mScarlet
CAGE peak (random-flank context):

| Replicate | CAGE peak | % of target (20) |
|---|---|---|
| 0 | 7.589 | 37.9% |
| 1 | 7.328 | 36.6% |
| 2 | **19.859** | **99.3%** |
| 3 | 11.419 | 57.1% |
| ensemble mean | 11.422 | 57.1% |

Replicate 2 is consistently the most favourable and lands suspiciously close to the exact target.
Worth checking whether the design was optimized against a single replicate. See
`data/tables/table7_borzoi_replicate_spread.csv`.
