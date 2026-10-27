# 03. The CDS-only control

## Purpose

The lab's `final` sequence differs from `initial` in two ways at once: three CDSs were inserted and
~2,900 regulatory edits were made. Any comparison between them therefore confounds "adding coding
sequence" with "adding designed regulation".

The control isolates the second effect: the Markov backbone with the three final CDSs spliced in at
their exact coordinates, and **nothing else changed**.

## Construction

Files produced: `data/sequences/initial_with_cds.fasta` and
`data/sequences/output_X_0_initial_with_cds.txt` (identical payload, two formats for convenience).

Because every Ledidi change is a substitution, CDS coordinates in `final` apply unchanged to
`initial`, so the splice is a direct window copy at 4462–5182, 9644–10343 and 14805–15507.

## Verification (all checks passed)

| Check | Result |
|---|---|
| Length | 19,968 bp |
| Three CDS windows byte-identical to `final` | yes |
| eGFP, mScarlet match provided FASTA | 100% |
| Diffs versus `initial` **outside** the windows | **0** |
| Diffs versus `initial` **inside** the windows | 1,616 |
| Diffs versus `final` **inside** the windows | **0** |
| Diffs versus `final` **outside** the windows | 2,909 (the Ledidi regulatory edits) |
| All ORFs ATG-start, in frame, terminal TAA, no internal stops | yes |
| Alphabet | A/C/G/T only |
| GC (non-CDS) | 40.8%, matching `initial`; `final` is 42.2% |
| CpG o/e (non-CDS) | 0.23, matching `initial`; `final` is 0.50 |

So the control is a clean two-way decomposition: `initial` → `control` adds coding sequence only,
`control` → `final` adds regulation only.

## Known asymmetries to remember when interpreting

1. The control retains an uninterrupted **upstream ATG at 4393**, giving a 23-aa N-terminal
   extension on eGFP. Ledidi removed this in `final`.
2. Conversely, the control retains the two in-frame stop codons before RFP3 that `final` lost.
3. The mScarlet premature poly(A) signal at +34 is present in **all three** sequences, so it does
   not confound any comparison.

## Headline result from the control

In real AAVS1 context, coding sequence alone buys essentially nothing regulatory:

| Readout | initial | control | final | final / control |
|---|---|---|---|---|
| CAGE+ mScarlet | 0.0238 | 0.0189 | 15.846 | **839×** |
| CAGE+ RFP3 | 0.0101 | 0.0139 | 9.887 | **713×** |
| CAGE+ eGFP | 0.0113 | 0.0389 | 10.021 | **258×** |
| DNase RFP3 | 0.0933 | 0.0716 | 6.865 | **96×** |
| DNase mScarlet | 0.1438 | 0.0635 | 0.796 | 12.5× |
| RNA+ RFP3 | 0.0256 | 0.0610 | 7.691 | **126×** |
| RNA+ mScarlet | 0.0492 | 0.1073 | 7.182 | **67×** |
| RNA+ eGFP | 0.0124 | 0.1174 | 6.496 | **55×** |

The CDS-only control is indistinguishable from the bare Markov backbone on CAGE and DNase. All of
the regulatory signal comes from the Ledidi edits, and DNase was never in the objective.
