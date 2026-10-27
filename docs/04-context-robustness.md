# 04. Context robustness and promoter necessity/sufficiency

## The problem

Borzoi requires a 524,288 bp input. The designs are 19,968 bp. So ~504 kb of flanking sequence must
be invented, and that choice turns out to matter more than almost anything else measured here.

Padding scheme used throughout: design centred at offset **252,160**, crop 163,840, design occupying
output bins **2,760–3,384** (624 bins of 32 bp).

## Result: a 437× swing from padding alone

Same `final` design, three padding strategies.

| Context | eGFP CAGE | mScarlet CAGE | RFP3 CAGE | eGFP RNA | mScarlet RNA | RFP3 RNA | eGFP DNase | RFP3 DNase |
|---|---|---|---|---|---|---|---|---|
| random, seed 101 | 0.023 | 4.682 | 0.038 | 0.063 | 0.065 | 0.060 | 0.035 | 0.064 |
| random, seed 202 | 0.023 | 4.704 | 0.038 | 0.063 | 0.065 | 0.060 | 0.035 | 0.064 |
| random, seed 303 | 0.023 | 4.700 | 0.038 | 0.063 | 0.065 | 0.060 | 0.035 | 0.064 |
| all-N | 9.663 | 11.936 | 9.454 | 4.594 | 5.359 | 3.299 | 3.736 | 5.991 |
| **real AAVS1** | **10.021** | **15.846** | **9.887** | **6.496** | **7.182** | **7.690** | 3.529 | 6.865 |

Fold-range across contexts: **eGFP 437×**, RFP3 262×, mScarlet 3.4×.

Three independent random seeds agree to within about 2% on every one of 30 measured quantities, so
this is a systematic property of AT-rich random DNA, not sampling noise.

### Interpretation

The Markov backbone is 41% GC. Half a megabase of composition-matched random DNA apparently reads to
Borzoi as inactive, gene-poor, heterochromatin-like territory and suppresses the construct inside
it. All-N padding is uninformative and therefore permissive. The real AAVS1 locus (52.8% GC,
gene-rich, within the actively transcribed `PPP1R12C`) is the most permissive of the three.

**Consequence: composition-matched random padding is the worst of the three options**, which is
counterintuitive, because it is the one that looks most "fair".

## Why AAVS1

AAVS1 (`PPP1R12C` intron 1) is the standard human safe-harbour integration site, so it is where a
construct like this would actually end up. Region used: hg38 `chr19:54,862,840-55,387,128`,
524,288 bp, zero N, 52.8% GC, retrieved via the UCSC REST API in nine chunks.

Stored at `data/context/aavs1_flank.txt`.

## Promoter necessity and sufficiency

Two complementary perturbations on the 600 bp window upstream of each ATG (which contains the
designed TSS at ATG − 128 and the 128 bp 5′UTR). About 120–140 edited bases per gene are involved.

- **Transplant (sufficiency):** copy the designed window into the CDS-only control.
- **Scramble (necessity):** seeded shuffle of the same window inside the design, preserving
  mononucleotide composition.

All in AAVS1 context. CDSs untouched in every arm.

### Necessity: near-total

| Track | eGFP | mScarlet | RFP3 |
|---|---|---|---|
| CAGE lost on scramble | 96% | **100%** | 99% |
| DNase lost on scramble | 87% | 99% | **100%** |
| RNA lost on scramble | 92% | 84% | 87% |

### Sufficiency: poor

| Track | eGFP | mScarlet | RFP3 |
|---|---|---|---|
| CAGE recovered by transplant | 0% | 49% | 5% |
| DNase recovered by transplant | 0% | 61% | 36% |
| RNA recovered by transplant | ~0% | 1% | 1% |

### Interpretation

The designed promoters are **necessary but not sufficient**. Removing them collapses everything;
installing them in an otherwise undesigned backbone recovers little, and almost nothing for RNA.
Ledidi's other ~2,800 edits across the construct are doing required work, and the design behaves as
a distributed, cooperative system rather than three isolated promoters.

This pushes back on the deck's attribution-map reading that the oracles are "essentially local".

## Track coupling across the nine arms

Correlating readouts across all nine perturbation arms (which span dead to fully working):

| Gene | r(CAGE, DNase) | r(CAGE, H3K4me3) | r(CAGE, RNA) | r(CAGE, Ser5-P Pol2) |
|---|---|---|---|---|
| eGFP | **+0.998** | +0.999 | +0.997 | +0.971 |
| mScarlet | **+0.994** | +0.984 | +0.633 | +0.586 |
| RFP3 | **+0.942** | +0.974 | +0.993 | +0.971 |

**This is the practical basis for dropping CAGE from validation.** At the level of whole elements,
accessibility and H3K4me3 track CAGE almost perfectly.

Note that mScarlet is an outlier: highest CAGE of the three (13.5) but by far the lowest DNase
(0.37, against 3.7 and 5.6). A real dissociation worth understanding.

## Caveats

- Scrambling preserves mononucleotide but not dinucleotide composition, so it destroys CpG
  dinucleotides, and these promoters are CpG-enriched. Some measured necessity may reflect CpG
  content rather than specific motifs. A dinucleotide-preserving shuffle is the better control.
- The 600 bp window is a choice; elements further upstream were not individually tested.
- eGFP necessity figures for DNase and RNA are unstable because the transplant arms recovered
  nothing, making the denominator small.

Data: `data/tables/table1_context_robustness.csv`,
`data/tables/table2_promoter_necessity_sufficiency.csv`.
